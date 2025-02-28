import numpy as np
import rasterio
from rasterio.windows import Window
import os
import multiprocessing
from joblib import Parallel, delayed

# Function to calculate the Contagion Index for a given 7x7 window. In parameters of the function define the number of land use and land cover class of the input raster
def calculate_contagion(values, n_classes=9):
    """
    Computes the Contagion Index for a given window of land cover values.
    Parameters:
    - values: 1D NumPy array representing a flattened 7x7 pixel neighborhood.
    - n_classes: Number of land cover classes (default is 9).
    Returns:
    - Contagion Index value (float) if valid, otherwise NaN.
    """
    # If the entire window is NaN (no data), return NaN
    if np.all(np.isnan(values)):
        return np.nan  
    # Reshape the flattened array into a 7x7 matrix
    size = int(np.sqrt(len(values)))
    matrix = values.reshape((size, size))
    # Initialize adjacency matrix to count transitions between classes
    adjacency_matrix = np.zeros((n_classes, n_classes))
    # Iterate through each pixel in the 7x7 window
    for i in range(size):
        for j in range(size):
            current_class = matrix[i, j]
            # Ensure valid class values (non-NaN and within the class range)
            if not np.isnan(current_class) and 0 <= current_class < n_classes and cur-rent_class.is_integer():
                neighbors = []
                # Check 8-directional adjacency (including diagonals)
                if i > 0: neighbors.append(matrix[i - 1, j])  # Up
                if i < size - 1: neighbors.append(matrix[i + 1, j])  # Down
                if j > 0: neighbors.append(matrix[i, j - 1])  # Left
                if j < size - 1: neighbors.append(matrix[i, j + 1])  # Right
                if i > 0 and j > 0: neighbors.append(matrix[i - 1, j - 1])  # Up-Left
                if i > 0 and j < size - 1: neighbors.append(matrix[i - 1, j + 1])  # Up-Right
                if i < size - 1 and j > 0: neighbors.append(matrix[i + 1, j - 1])  # Down-Left
                if i < size - 1 and j < size - 1: neighbors.append(matrix[i + 1, j + 1])  # Down-Right

                # Update adjacency counts for neighboring land cover classes
                for neighbor in neighbors:
                    if not np.isnan(neighbor) and 0 <= neighbor < n_classes and neigh-bor.is_integer():
                        adjacency_matrix[int(current_class), int(neighbor)] += 1
    # Compute the total number of adjacency transitions
    total_adjacencies = np.sum(adjacency_matrix)
    # If there are no valid adjacencies, return NaN
    if total_adjacencies == 0:
        return np.nan
    # Calculate the probability of each adjacency transition
    p_ij = adjacency_matrix / total_adjacencies
    # Compute the Contagion Index using entropy formula
    p_ij_log = np.where(p_ij > 0, p_ij * np.log(p_ij), 0)
    contagion_raw = np.sum(p_ij_log)
    contagion_index = 1 + (contagion_raw / (2 * np.log(n_classes)))
    return contagion_index
# Function to process a single pixel in a raster block
def process_pixel(i, j, block, buffer_size, n_classes=9):
    """
    Extracts a 7x7 window around a pixel and computes its Contagion Index.
    Parameters:
    - i, j: Pixel coordinates within the raster block.
    - block: NumPy array representing the land cover raster block.
    - buffer_size: Half-size of the moving window (7x7 → 3-pixel buffer).
    - n_classes: Number of land cover classes.
    Returns:
    - Tuple (i, j, contagion_value) where contagion_value is the computed Contagion Index.
    """
    
    # Extract 7x7 neighborhood window
    window_values = block[i - buffer_size:i + buffer_size + 1, j - buffer_size:j + buffer_size + 1]
    # Compute Contagion Index for the extracted window
    return i, j, calculate_contagion(window_values.flatten(), n_classes=n_classes)

# Block processing parameters
window_size = 7  # size moving window
half_window = window_size // 2  # 3-pixel buffer to account for edge effects
buffer_size = half_window  
block_size = 5000  # Process raster in blocks of 5000x5000 pixels, define depend on the input raster size.

# Specify the starting block (in terms of row and column)
start_row = 0  # Row index to start processing from
start_col = 0  # Column index to start processing from

# Input raster file and output directory
input_file = "/content/drive/My Drive/folder/lulc.tif"
output_dir = "/content/drive/My Drive/folder"
os.makedirs(output_dir, exist_ok=True)  # Create output directory if not exists

# Open the input raster file
with rasterio.open(input_file) as src:
    profile = src.profile  # Copy raster metadata
    profile.update(dtype=rasterio.float32, count=1, nodata=np.nan)  # Update to float32 format

    # Iterate through raster in blocks (starting from the specified row/column)
    for row_start in range(start_row, src.height, block_size - 2 * buffer_size):
        for col_start in range(start_col, src.width, block_size - 2 * buffer_size):
            # Define the block size, considering the image boundaries
            row_end = min(row_start + block_size, src.height)
            col_end = min(col_start + block_size, src.width)

            # Define the block window (including buffer)
            window = Window(col_start, row_start, col_end - col_start, row_end - row_start)
            block = src.read(1, window=window)
            # Skip empty blocks (all NaN)
            if np.all(np.isnan(block)):
                print(f"Skipping empty block at rows {row_start}-{row_end} and cols {col_start}-{col_end}")
                continue

            # Create an output array for storing Contagion index results
            output_block = np.full_like(block, np.nan, dtype=np.float32)
            # Get block dimensions
            block_height, block_width = block.shape
            # Set up parallel processing (use up to 16 CPU cores)
            num_cores = min(16, multiprocessing.cpu_count())  
            results = Parallel(n_jobs=num_cores)(
                delayed(process_pixel)(i, j, block, buffer_size, n_classes=9)
                for i in range(buffer_size, block_height - buffer_size)
                for j in range(buffer_size, block_width - buffer_size)
            )

            # Assign the computed Contagion Index values to the output block
            for i, j, contagion_value in results:
                output_block[i, j] = contagion_value
            # Remove the buffer before saving the block
            core_output_block = output_block[
                buffer_size:block.shape[0] - buffer_size,
                buffer_size:block.shape[1] - buffer_size
            ]

            # Skip saving if there is no valid data
            if np.all(np.isnan(core_output_block)):
                print(f"No valid data in core block at rows {row_start}-{row_end} and cols {col_start}-{col_end}")
                continue

            # Define output transform for the cropped block
            core_transform = rasterio.windows.transform(
                Window(
                    col_start + buffer_size,
                    row_start + buffer_size,
                    core_output_block.shape[1],
                    core_output_block.shape[0]
                ),
                src.transform
            )

            # Update raster metadata for output block
            block_profile = profile.copy()
            block_profile.update(
                height=core_output_block.shape[0],
                width=core_output_block.shape[1],
                transform=core_transform
            )

            # Define the output file path for the block
            output_path = os.path.join(output_dir, f"block_{row_start}_{col_start}.tif")
            # Save the processed block as a new raster file
            with rasterio.open(output_path, "w", **block_profile) as dst:
                dst.write(core_output_block, 1)
            print(f"Saved block: {output_path}")
