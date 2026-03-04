"""
Extract acquisition parameters from DICOM files
Includes Siemens CSA header parser for MRSI data
"""
import pydicom
import os
import struct
import re


def parse_csa_header(csa_header_bytes):
    """
    Parse Siemens CSA header and extract parameters
    CSA headers start with 'SV10' and contain name-value pairs
    """
    if not csa_header_bytes or len(csa_header_bytes) < 8:
        return {}
    
    # Check for SV10 magic signature
    if csa_header_bytes[:4] != b'SV10':
        return {}
    
    params = {}
    
    try:
        # Skip first 8 bytes (SV10 + unused)
        pos = 8
        
        # Read number of tags
        num_tags = struct.unpack('<I', csa_header_bytes[pos:pos+4])[0]
        pos += 4
        
        # Skip unused bytes
        pos += 4
        
        # Parse each tag
        for i in range(num_tags):
            if pos + 84 > len(csa_header_bytes):
                break
            
            # Tag name (64 bytes, null-terminated)
            name_bytes = csa_header_bytes[pos:pos+64]
            name = name_bytes.split(b'\x00')[0].decode('latin-1', errors='ignore')
            pos += 64
            
            # VM (value multiplicity)
            vm = struct.unpack('<I', csa_header_bytes[pos:pos+4])[0]
            pos += 4
            
            # VR (value representation) - 4 bytes
            vr = csa_header_bytes[pos:pos+4].decode('latin-1', errors='ignore').strip('\x00')
            pos += 4
            
            # Syngo DT (4 bytes)
            pos += 4
            
            # Number of items
            nitems = struct.unpack('<I', csa_header_bytes[pos:pos+4])[0]
            pos += 4
            
            # Skip reserved bytes
            pos += 4
            
            # Read values
            values = []
            for j in range(nitems):
                if pos + 8 > len(csa_header_bytes):
                    break
                
                # Length of this item (4 bytes, each)
                item_length = struct.unpack('<4I', csa_header_bytes[pos:pos+16])
                pos += 16
                
                if item_length[1] > 0 and pos + item_length[1] <= len(csa_header_bytes):
                    value_bytes = csa_header_bytes[pos:pos+item_length[1]]
                    value = value_bytes.decode('latin-1', errors='ignore').strip('\x00')
                    values.append(value)
                    pos += (item_length[1] + 3) // 4 * 4  # Align to 4-byte boundary
            
            if values:
                params[name] = values if len(values) > 1 else values[0]
        
    except Exception as e:
        pass  # Silently fail on parsing errors
    
    return params


def extract_siemens_params(ds):
    """Extract acquisition parameters from Siemens DICOM with CSA headers"""
    params = {}
    
    # Try to get CSA Series Header Info (0029, 1120)
    csa_series_tag = (0x0029, 0x1120)
    if csa_series_tag in ds:
        csa_data = ds[csa_series_tag].value
        csa_params = parse_csa_header(csa_data)
        
        # Extract specific parameters
        if 'UsedPatientWeight' in csa_params:
            params['patient_weight'] = csa_params['UsedPatientWeight']
        if 'MrPhoenixProtocol' in csa_params:
            # Parse Phoenix protocol (contains TR, TE, etc.)
            phoenix = csa_params['MrPhoenixProtocol']
            
            # Extract TR
            tr_match = re.search(r'lRepetitionTime\s*=\s*(\d+)', phoenix)
            if not tr_match:
                tr_match = re.search(r'alTR\[0\]\s*=\s*(\d+)', phoenix)
            if tr_match:
                params['tr'] = float(tr_match.group(1)) / 1000.0  # Convert µs to ms
            
            # Extract TE
            te_match = re.search(r'alTE\[0\]\s*=\s*(\d+)', phoenix)
            if not te_match:
                te_match = re.search(r'lEchoTime\s*=\s*(\d+)', phoenix)
            if te_match:
                params['te'] = float(te_match.group(1)) / 1000.0  # Convert µs to ms
            
            # Extract TI
            ti_match = re.search(r'alTI\[0\]\s*=\s*(\d+)', phoenix)
            if ti_match:
                params['ti'] = float(ti_match.group(1)) / 1000.0  # Convert µs to ms
            
            # Extract flip angle
            fa_match = re.search(r'adFlipAngleDegree\[0\]\s*=\s*([\d.]+)', phoenix)
            if fa_match:
                params['flip_angle'] = float(fa_match.group(1))
            
            # Extract bandwidth
            bw_match = re.search(r'alDwellTime\[0\]\s*=\s*(\d+)', phoenix)
            if bw_match:
                dwell_time_ns = float(bw_match.group(1))
                params['dwell_time_us'] = dwell_time_ns / 1000.0
                params['bandwidth_hz'] = 1000000.0 / dwell_time_ns if dwell_time_ns > 0 else None
            
            # Extract spectral points
            spec_pts_match = re.search(r'lVectorSize\s*=\s*(\d+)', phoenix)
            if spec_pts_match:
                params['spectral_points'] = int(spec_pts_match.group(1))
            
            # Extract voxel size
            fov_match = re.search(r'sSliceArray\.asSlice\[0\]\.dReadoutFOV\s*=\s*([\d.]+)', phoenix)
            phase_fov_match = re.search(r'sSliceArray\.asSlice\[0\]\.dPhaseFOV\s*=\s*([\d.]+)', phoenix)
            thickness_match = re.search(r'sSliceArray\.asSlice\[0\]\.dThickness\s*=\s*([\d.]+)', phoenix)
            cols_match = re.search(r'sKSpace\.lBaseResolution\s*=\s*(\d+)', phoenix)
            
            if fov_match:
                params['fov_readout'] = float(fov_match.group(1))
            if phase_fov_match:
                params['fov_phase'] = float(phase_fov_match.group(1))
            if thickness_match:
                params['slice_thickness'] = float(thickness_match.group(1))
            if cols_match:
                params['base_resolution'] = int(cols_match.group(1))
    
    return params


def extract_dicom_params(filepath, label, show_debug=False):
    """Extract and print acquisition parameters from a DICOM file"""
    # Check if file exists first
    if not os.path.exists(filepath):
        print(f"\n{'='*80}")
        print(f"FILE NOT FOUND: {label}")
        print(f"{'='*80}")
        return None
    
    # Try to read the DICOM file
    try:
        ds = pydicom.dcmread(filepath, force=True)
    except Exception as e:
        print(f"\n{'='*80}")
        print(f"⚠️ ERROR reading {label}")
        print(f"{'='*80}")
        print(f"Error: {str(e)}")
        return None
    
    # Debug: Show all available tags (optional)
    if show_debug:
        print(f"\n{'='*80}")
        print(f"DEBUG: All available DICOM tags")
        print(f"{'='*80}")
        print(f"Total tags found: {len(ds.keys())}")
        for key in ds.keys():  # Show ALL tags
            try:
                tag_name = ds[key].name
                tag_value = str(ds[key].value)[:150]  # Truncate long values
                print(f"{key} - {tag_name}: {tag_value}")
            except:
                pass
        print(f"\n")
    
    params = {}
    params['label'] = label
    params['protocol'] = getattr(ds, 'ProtocolName', 'N/A')
    params['series'] = getattr(ds, 'SeriesDescription', 'N/A')
    params['tr'] = getattr(ds, 'RepetitionTime', None)
    params['te'] = getattr(ds, 'EchoTime', None)
    params['ti'] = getattr(ds, 'InversionTime', None)
    params['flip_angle'] = getattr(ds, 'FlipAngle', None)
    params['matrix'] = f"{ds.Rows} x {ds.Columns}" if hasattr(ds, 'Rows') and hasattr(ds, 'Columns') else None
    params['pixel_spacing'] = ds.PixelSpacing if hasattr(ds, 'PixelSpacing') else None
    params['slice_thickness'] = getattr(ds, 'SliceThickness', None)
    params['pixel_bandwidth'] = getattr(ds, 'PixelBandwidth', None)
    params['field_strength'] = getattr(ds, 'MagneticFieldStrength', None)
    params['manufacturer'] = getattr(ds, 'Manufacturer', 'N/A')
    params['model'] = getattr(ds, 'ManufacturerModelName', 'N/A')
    
    # If Siemens, try to extract CSA header parameters
    if params['manufacturer'] and 'SIEMENS' in params['manufacturer'].upper():
        siemens_params = extract_siemens_params(ds)
        # Override/add Siemens-specific parameters
        for key, value in siemens_params.items():
            if value is not None:
                params[key] = value
    
    # Calculate derived parameters
    if params['pixel_spacing'] and params['slice_thickness']:
        voxel_vol = params['pixel_spacing'][0] * params['pixel_spacing'][1] * params['slice_thickness']
        params['voxel_size'] = f"{params['pixel_spacing'][0]:.2f} x {params['pixel_spacing'][1]:.2f} x {params['slice_thickness']:.2f} mm³"
        params['voxel_volume'] = voxel_vol
    else:
        params['voxel_size'] = None
        params['voxel_volume'] = None
    
    if params['pixel_spacing'] and hasattr(ds, 'Rows') and hasattr(ds, 'Columns'):
        fov_x = ds.Columns * params['pixel_spacing'][0]
        fov_y = ds.Rows * params['pixel_spacing'][1]
        params['fov'] = f"{fov_x:.1f} x {fov_y:.1f} mm"
    else:
        params['fov'] = None
    
    return params



def print_params(params):
    """Pretty print the parameters"""
    if params is None:
        return
    
    print(f"\n{'='*80}")
    print(f"{params['label']}")
    print(f"{'='*80}")
    print(f"Protocol: {params['protocol']}")
    print(f"Series: {params['series']}")
    
    print(f"\nTIMING PARAMETERS:")
    print(f"  TR: {params['tr']} ms" if params['tr'] else "  TR: N/A")
    print(f"  TE: {params['te']} ms" if params['te'] else "  TE: N/A")
    if params['ti']:
        print(f"  TI: {params['ti']} ms")
    if params['flip_angle']:
        print(f"  Flip Angle: {params['flip_angle']}°")
    
    print(f"\nSPATIAL PARAMETERS:")
    if params['matrix']:
        print(f"  Matrix: {params['matrix']}")
    if 'base_resolution' in params:
        print(f"  Base Resolution: {params['base_resolution']}")
    if params['pixel_spacing']:
        print(f"  In-plane resolution: {params['pixel_spacing'][0]:.3f} x {params['pixel_spacing'][1]:.3f} mm")
    if params['slice_thickness']:
        print(f"  Slice thickness: {params['slice_thickness']} mm")
    if params['voxel_size']:
        print(f"  Voxel size: {params['voxel_size']} (Volume: {params['voxel_volume']:.2f} mm³)")
    if params['fov']:
        print(f"  FOV: {params['fov']}")
    elif 'fov_readout' in params and 'fov_phase' in params:
        print(f"  FOV: {params['fov_readout']:.1f} x {params['fov_phase']:.1f} mm")
    
    if params['pixel_bandwidth']:
        print(f"\nBANDWIDTH:")
        print(f"  {params['pixel_bandwidth']} Hz/pixel")
    
    # Spectroscopy parameters (if available)
    if 'spectral_points' in params or 'bandwidth_hz' in params:
        print(f"\nSPECTROSCOPY PARAMETERS:")
        if 'spectral_points' in params:
            print(f"  Spectral Points: {params['spectral_points']}")
        if 'dwell_time_us' in params:
            print(f"  Dwell Time: {params['dwell_time_us']:.3f} µs")
        if 'bandwidth_hz' in params:
            print(f"  Bandwidth: {params['bandwidth_hz']:.1f} Hz")
    
    print(f"\nSCANNER:")
    print(f"  Field strength: {params['field_strength']} T" if params['field_strength'] else "  Field strength: N/A")
    print(f"  Manufacturer: {params['manufacturer']}")
    print(f"  Model: {params['model']}")
    
    # Additional info
    if 'patient_weight' in params:
        print(f"\nOTHER:")
        print(f"  Patient Weight: {params['patient_weight']} kg")


if __name__ == "__main__":
    import glob
    
    # Directory containing DICOM files (can be .IMA or .dcm)
    dicom_dir = '/Users/bbbartel/ASU Dropbox/BioProteanLab/MRSI/MRS_Clinical_Data/mayo_mrs_data/FLAIR data/RMS01/3D SAG T2 SPC FLAIR C9C1HN013 - 4'
    
    # Find all DICOM files in the directory (.IMA or .dcm)
    ima_files = sorted(glob.glob(os.path.join(dicom_dir, '*.IMA')))
    dcm_files = sorted(glob.glob(os.path.join(dicom_dir, '*.dcm')))
    dicom_files = ima_files + dcm_files
    
    if not dicom_files:
        print(f"No DICOM files found in: {dicom_dir}")
    else:
        print(f"Found {len(dicom_files)} DICOM files")
        
        # Extract parameters from the first file as a representative sample
        print(f"\nAnalyzing first file as representative sample...")
        params = extract_dicom_params(dicom_files[0], f'FLAIR - Slice 1/{len(dicom_files)}', show_debug=False)
        
        if params:
            print_params(params)
            
            # Optionally analyze all files (warning: may produce a lot of output)
            analyze_all = False  # Set to True to see all files
            if analyze_all and len(dicom_files) > 1:
                print(f"\n{'='*80}")
                print(f"Analyzing remaining {len(dicom_files)-1} files...")
                print(f"{'='*80}\n")
                
                for i, filepath in enumerate(dicom_files[1:], start=2):
                    params = extract_dicom_params(filepath, f'FLAIR - Slice {i}/{len(dicom_files)}')
                    if params:
                        print_params(params)
    
    print(f"\n{'='*80}\n")
    print("Summary complete!")
    print("Summary complete!")
