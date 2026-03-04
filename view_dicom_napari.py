"""
View DICOM images with Napari
"""
import napari
import pydicom
import numpy as np
import glob
import os


def load_dicom_series(directory, pattern='*.dcm'):
    """Load a series of DICOM files as a 3D volume"""
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    
    if not files:
        print(f"No files matching {pattern} found in {directory}")
        return None, None
    
    print(f"Loading {len(files)} DICOM files...")
    
    # Read first file to get dimensions
    ds0 = pydicom.dcmread(files[0], force=True)
    
    # Try to get pixel array
    try:
        img_shape = ds0.pixel_array.shape
        dtype = ds0.pixel_array.dtype
    except:
        print("Warning: Could not read pixel data from DICOM")
        return None, None
    
    # Preallocate volume
    volume = np.zeros((len(files), img_shape[0], img_shape[1]), dtype=dtype)
    
    # Load all slices
    for i, filepath in enumerate(files):
        try:
            ds = pydicom.dcmread(filepath, force=True)
            volume[i] = ds.pixel_array
        except Exception as e:
            print(f"Error reading {os.path.basename(filepath)}: {e}")
            continue
    
    # Get voxel spacing if available
    spacing = [1.0, 1.0, 1.0]  # Default
    if hasattr(ds0, 'SliceThickness'):
        spacing[0] = float(ds0.SliceThickness)
    if hasattr(ds0, 'PixelSpacing'):
        spacing[1] = float(ds0.PixelSpacing[0])
        spacing[2] = float(ds0.PixelSpacing[1])
    
    print(f"Volume shape: {volume.shape}")
    print(f"Spacing: {spacing} mm")
    print(f"Data range: [{volume.min()}, {volume.max()}]")
    
    return volume, spacing


def view_dicom_volume(directory, pattern='*.dcm'):
    """Load and view DICOM volume in Napari"""
    volume, spacing = load_dicom_series(directory, pattern)
    
    if volume is None:
        return
    
    # Create viewer
    viewer = napari.Viewer()
    
    # Add volume as image layer
    viewer.add_image(
        volume,
        name=os.path.basename(directory),
        scale=spacing,
        colormap='gray',
        contrast_limits=[volume.min(), volume.max()]
    )
    
    print("\nNapari viewer opened!")
    print("Controls:")
    print("  - Scroll: Navigate slices")
    print("  - Click/drag: Adjust contrast")
    print("  - Right panel: Layer controls")
    print("  - Close window when done")
    
    napari.run()


if __name__ == "__main__":
    # Example: View FLAIR data
    dicom_dir = '/Users/bbbartel/ASU Dropbox/BioProteanLab/MRSI/MRS_Clinical_Data/mayo_mrs_data/FLAIR data/RMS01/3D SAG T2 SPC FLAIR C9C1HN013 - 4'
    
    # Check if directory exists
    if os.path.exists(dicom_dir):
        view_dicom_volume(dicom_dir, pattern='*.dcm')
    else:
        print(f"Directory not found: {dicom_dir}")
        print("\nPlease update the 'dicom_dir' variable with your data path")
