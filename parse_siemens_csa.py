"""
Parse Siemens CSA headers from DICOM files
Based on the CSA header format specification
"""
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
        print("Warning: CSA header doesn't start with 'SV10'")
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
        print(f"Error parsing CSA header: {e}")
    
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
            params['phoenix_protocol'] = phoenix
            
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
    
    # Try CSA Image Header Info (0029, 1110)
    csa_image_tag = (0x0029, 0x1110)
    if csa_image_tag in ds:
        csa_data = ds[csa_image_tag].value
        csa_params = parse_csa_header(csa_data)
        
        if 'ImaAbsTablePosition' in csa_params:
            params['table_position'] = csa_params['ImaAbsTablePosition']
    
    return params


if __name__ == "__main__":
    import pydicom
    import os
    
    # Test with first .IMA file
    dicom_dir = '/Users/bbbartel/ASU Dropbox/BioProteanLab/MRSI/MRS_Clinical_Data/mayo_mrs_data/Original MS Data/RMS01_Anonymized/WIP_3D_MRS_SI_0018'
    ima_files = sorted([f for f in os.listdir(dicom_dir) if f.endswith('.IMA')])
    
    if ima_files:
        filepath = os.path.join(dicom_dir, ima_files[0])
        print(f"Reading: {ima_files[0]}\n")
        
        ds = pydicom.dcmread(filepath, force=True)
        params = extract_siemens_params(ds)
        
        print("="*80)
        print("EXTRACTED SIEMENS MRSI PARAMETERS")
        print("="*80)
        
        print("\nTIMING PARAMETERS:")
        if 'tr' in params:
            print(f"  TR: {params['tr']:.2f} ms")
        if 'te' in params:
            print(f"  TE: {params['te']:.2f} ms")
        if 'ti' in params:
            print(f"  TI: {params['ti']:.2f} ms")
        if 'flip_angle' in params:
            print(f"  Flip Angle: {params['flip_angle']:.1f}°")
        
        print("\nSPATIAL PARAMETERS:")
        if 'fov_readout' in params:
            print(f"  FOV Readout: {params['fov_readout']:.1f} mm")
        if 'fov_phase' in params:
            print(f"  FOV Phase: {params['fov_phase']:.1f} mm")
        if 'slice_thickness' in params:
            print(f"  Slice/Slab Thickness: {params['slice_thickness']:.1f} mm")
        if 'base_resolution' in params:
            print(f"  Base Resolution: {params['base_resolution']}")
        
        print("\nSPECTROSCOPY PARAMETERS:")
        if 'spectral_points' in params:
            print(f"  Spectral Points: {params['spectral_points']}")
        if 'dwell_time_us' in params:
            print(f"  Dwell Time: {params['dwell_time_us']:.3f} µs")
        if 'bandwidth_hz' in params:
            print(f"  Bandwidth: {params['bandwidth_hz']:.1f} Hz")
        
        print("\nOTHER:")
        if 'patient_weight' in params:
            print(f"  Patient Weight: {params['patient_weight']} kg")
    else:
        print("No .IMA files found")
