import numpy as np
import rasterio
from rasterio import features
from scipy.ndimage import label
from osgeo import gdal
import os
import multiprocessing
from joblib import Parallel, delayed

def mpa_by_pixel(i, j, land_cover_block, kernel_radius):
    """
    Computes the Mean Patch Area (MPA) for a given pixel by analyzing a 7x7 neighborhood.
        Parameters:
    - i, j: Coordinates of the pixel within the block.
    - land_cover_block: 2D NumPy array representing the land cover raster block.
    - kernel_radius: Radius of the neighborhood kernel (3 pixels for a 7x7 window).
    Returns:
    - Mean patch area (m²) of the identified patches within the 7x7 window.
    - NaN if the window contains only nodata values.
    """
    # Extract the 7x7 neighborhood around the current pixel
    kernel = land_cover_block[i - kernel_radius:i + kernel_radius + 1,
                              j - kernel_radius:j + kernel_radius + 1]
    # Skip computation if the entire kernel contains only nodata values (-1)
    if np.all(kernel == -1):
        return np.nan
    # Identify unique land cover classes within the neighborhood (excluding nodata)
    unique_classes = np.unique(kernel[kernel >= 0])  
    patch_areas = []  # Store the area of patches found in the kernel
    for cls in unique_classes:
        # Create a binary mask for the specific land cover class
        binary_kernel = (kernel == cls).astype(int)
        # Label connected patches of the same class
        labeled_patches, num_features = label(binary_kernel)
        # Compute the area of each patch (assuming 10x10 m per pixel)
        for patch_id in range(1, num_features + 1):
            patch_size = np.sum(labeled_patches == patch_id)
            patch_areas.append(patch_size * 100)  # Convert to square meters
    # Return the mean patch area for the pixel, or 0 if no patches were found
    return np.mean(patch_areas) if patch_areas else 0

def mpa_by_block(input_raster_path, output_folder_path, block_size=5000, start_block=0):
    """
    Computes the Mean Patch Area (MPA) for an input raster by processing it in blocks.
        Parameters:
    - input_raster_path: Path to the input raster file.
    - output_folder_path: Directory where processed blocks will be saved.
    - block_size: Size of each processing block in pixels (default: 5000x5000), define the block size depend on the input raster size.
    - start_block: Index of the first block to process (useful for resuming interrupted runs).
    Output:
    - Saves each processed block as a separate raster file.
    """
    # Create the output directory if it does not exist
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    # Open the input raster file
    with rasterio.open(input_raster_path) as src:
        profile = src.profile  # Store raster metadata
        nodata_value = src.nodata  # Retrieve the nodata value of the raster
        block_id = 0  # Initialize block counter
        # Loop through raster blocks row-wise
        for block_index in range(0, src.height, block_size):
            # Loop through raster blocks column-wise
            for block_col_index in range(0, src.width, block_size):
                # Skip blocks that have already been processed
                if block_id < start_block:
                    block_id += 1
                    continue
                # Define the actual block size (accounting for raster edges)
                block_height = min(block_size + 6, src.height - block_index)
                block_width = min(block_size + 6, src.width - block_col_index)
                # Read the current block from the raster (including a 3-pixel buffer)
                window = rasterio.windows.Window(
                    max(0, block_col_index - 3), 
                    max(0, block_index - 3), 
                    block_width, 
                    block_height
                )
                land_cover_block = src.read(1, window=window)
                # Convert nodata values to -1 for easier processing
                land_cover_block = np.where(land_cover_block == nodata_value, -1, land_cover_block)
                # Define the 7x7 neighborhood moving window size/kernel size
                kernel_size = 7
                kernel_radius = kernel_size // 2  # 3 pixels (centered window)
                # Initialize an array to store the computed MPA values for the block
                block_result = np.full((block_height - 6, block_width - 6), fill_value=np.nan, dtype=float)
                # Set up parallel processing to compute MPA for each pixel
                num_cores = multiprocessing.cpu_count()
                num_cores = 4  # number of cores for optimal performance, change according the available resources 

                results = Parallel(n_jobs=num_cores)(
                    delayed(mpa_by_pixel)(i, j, land_cover_block, kernel_radius)
                    for i in range(kernel_radius, block_height - kernel_radius)
                    for j in range(kernel_radius, block_width - kernel_radius)
                )
                # Store results in the corresponding block result array
                index = 0
                for i in range(kernel_radius, block_height - kernel_radius):
                    for j in range(kernel_radius, block_width - kernel_radius):
                        block_result[i - kernel_radius, j - kernel_radius] = results[index]
                        index += 1
                # Update raster metadata for the output block
                block_profile = profile.copy()
                block_profile.update(
                    height=block_height - 6,
                    width=block_width - 6,
                    transform=rasterio.windows.transform(window, src.transform),
                    dtype=rasterio.float32,
                    nodata=-9999  # Define a new nodata value for output
                )
                # Define the output path for the processed block
                block_output_path = f"{output_folder_path}/MPA_block_{block_id}.tif"
                # Write the computed MPA values to a new raster file
                with rasterio.open(block_output_path, 'w', **block_profile) as dst:
                    dst.write(block_result, 1)
                print(f"Block {block_id} calculated and saved to: {block_output_path}")
                block_id += 1

# Define input raster path and output folder
input_raster = "/content/drive/My Drive/dw_2022/folder.tif"  
output_folder = "/content/drive/My Drive/dw_2022/folder"
# Start processing from block = number
start_block = 1 #This is especially useful when suddenly the calculation stop. 
# Run the block processing function
mpa_by_block(input_raster, output_folder, start_block=start_block)
