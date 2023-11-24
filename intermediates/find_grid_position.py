import os
import time
import argparse
import csv
import pandas as pd
import numpy as np
from spectral import *
from PIL import Image
import matplotlib.pyplot as plt

# Function Definitions
def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def plot_on_axes(ax, data, style, title, legend_label=None, x_data=None):
    if x_data is None:
        x_data = range(len(data))
    ax.plot(x_data, data, style, label=legend_label)
    ax.set_title(title)
    if legend_label:
        ax.legend()

def plot_derivatives(row_avg, row_avg_smooth,
                     col_avg, col_avg_smooth,
                     row_derivative1, col_derivative1,
                     row_derivative2, col_derivative2,
                     window_size, mode, data_path):

    # Plotting
    fig, axs = plt.subplots(3, 2, figsize=(18, 12))
    
    # Average values
    plot_on_axes(axs[0, 0], row_avg, 'o-', mode + ' Row Averages', 'Original')
    plot_on_axes(axs[0, 0], row_avg_smooth, 'x-', mode + ' Row Averages', 'Smoothed', range(window_size-1, len(row_avg)))

    plot_on_axes(axs[0, 1], col_avg, 'o-', mode + ' Column Averages', 'Original')
    plot_on_axes(axs[0, 1], col_avg_smooth, 'x-', mode + ' Column Averages', 'Smoothed', range(window_size-1, len(col_avg)))

    # First derivative
    plot_on_axes(axs[1, 0], row_derivative1, 'x-', mode + ' First Derivative of Smoothed Row Averages')
    plot_on_axes(axs[1, 1], col_derivative1, 'x-', mode + ' First Derivative of Smoothed Column Averages')

    # Second derivative
    plot_on_axes(axs[2, 0], row_derivative2, 's-', mode + ' Second Derivative of Smoothed Row Averages')
    plot_on_axes(axs[2, 1], col_derivative2, 's-', mode + ' Second Derivative of Smoothed Column Averages')

    plt.tight_layout()
    plt.savefig(os.path.join(data_path, f"{mode}_derivatives.png"))
    plt.show()

def plot_heatmap(mode, CLS_matrix, y1, y2, x1, x2, data_path):
    fig, axs = plt.subplots(1, 1)
    im = axs.imshow(CLS_matrix, cmap='viridis', aspect='auto')
    axs.set_title(f'{mode} signal intensity at each position')
    axs.set_xlabel('X position')
    axs.set_ylabel('Y position')

    # Adding red lines
    axs.axhline(y=y1, color='r', linestyle='--')
    axs.axhline(y=y2, color='r', linestyle='--')
    axs.axvline(x=x1, color='r', linestyle='--')
    axs.axvline(x=x2, color='r', linestyle='--')

    # Adding a colorbar for intensity
    cbar = plt.colorbar(im, ax=axs)
    cbar.set_label('Intensity')
    
    plt.savefig(os.path.join(data_path, f"{mode}_heatmap.png"))
    plt.show()

# Main Function
def main(args):
    data = args.data
    cap = args.cap
    index = args.index
    mode = args.mode
    plot = args.plot
    rotation = args.rotation
    flip_horizontal = args.flip_horizontal
    
    coordinates = {}
    
    start_time = time.time()  # Start time

    # List all files in the directory
    all_files = os.listdir(data)

    for mode in ["RGB", "hyperspectral"] if mode == "both" else [mode]:

        if mode == "hyperspectral" :

            selected_file = [f for f in all_files if "Broadband" in f and "hroma" in f and "hdr" in f]
            
            # Print or return the filtered file
            print("Selected file:", selected_file)
            
            # Assert that the length of the filtered files is 1
            assert len(selected_file) == 1, f"Expected 1 file, found {len(selected_file)}"

            hdr_path=os.path.join(data, selected_file[0])

            image1 = envi.open(hdr_path)

            print("Extracting channel at index " + str(index))

            CLS_matrix = np.squeeze(image1[:,:,index],
                                    axis=2)

            CLS_matrix = np.interp(CLS_matrix,
                                   (0,
                                    cap),
                                   (0,
                                    255)).astype(int)
            
            debug_time = time.time()
            print(f"Hyperspectral data loaded, elapsed time: {debug_time - start_time:.2f} seconds")

        if mode == "RGB":

            # Load the image
            selected_file = [f for f in all_files if "hromagrid" in f and "jpg" in f]
            jpg_path=os.path.join(data, selected_file[0])
            image = Image.open(jpg_path)

            # Convert to grayscale
            image_gray = image.convert("L")

            # Create CLS_matrix from the grayscale image
            CLS_matrix = np.array(image_gray)
            
            debug_time = time.time()
            print(f"RGB data loaded, elapsed time: {debug_time - start_time:.2f} seconds")

        if rotation:
            CLS_matrix = np.rot90(CLS_matrix, rotation // 90)  # Rotate by 90/180/270 degrees
        if flip_horizontal:
            CLS_matrix = np.fliplr(CLS_matrix)  # Flip horizontally
            
        # Calculate averages for each row and each column
        row_avg = np.mean(CLS_matrix, axis=1)
        col_avg = np.mean(CLS_matrix, axis=0)

        # Apply a moving average with a window size of 3
        window_size = 4
        row_avg_smooth = moving_average(row_avg, window_size)
        col_avg_smooth = moving_average(col_avg, window_size)

        # Compute the first and second derivatives of the smoothed data
        row_derivative1 = np.diff(row_avg_smooth)
        row_derivative2 = np.diff(row_derivative1)

        col_derivative1 = np.diff(col_avg_smooth)
        col_derivative2 = np.diff(col_derivative1)

        y1 = np.argmax(row_derivative1)
        y2 = np.argmin(row_derivative1)
        x1 = np.argmax(col_derivative1)
        x2 = np.argmin(col_derivative1)
        
        debug_time = time.time()
        print(f"Calculations done, elapsed time: {debug_time - start_time:.2f} seconds")

        
        coordinates[mode] = {"x1": x1, "x2": x2, "y1": y1, "y2": y2}

        csv_path = os.path.join(data, "coordinates.csv")
        df = pd.DataFrame.from_dict(coordinates, orient='index', columns=['x1', 'x2', 'y1', 'y2'])
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'mode'}, inplace=True)

        # Save as CSV
        df.to_csv(csv_path, index=False)
        
        # csv_path = os.path.join(data, "coordinates.csv")
        # with open(csv_path, "w", newline='\n') as csvfile:  # set newline to '\n'
        #     fieldnames = ["mode", "x1", "x2", "y1", "y2"]
        #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        #     writer.writeheader()
        #     for m, coords in coordinates.items():
        #         writer.writerow({"mode": m, **coords})
        
        debug_time = time.time()
        print(f"Coordinates saved, elapsed time: {debug_time - start_time:.2f} seconds")

        # Plotting only if --plot is True
        if plot:
            plot_derivatives(row_avg, row_avg_smooth, col_avg, col_avg_smooth,
                             row_derivative1, col_derivative1, row_derivative2, col_derivative2,
                             window_size, mode, data)

            plot_heatmap(mode, CLS_matrix, y1, y2, x1, x2, data)
            
            debug_time = time.time()
            print(f"Plotting done, elapsed time: {debug_time - start_time:.2f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for plotting and saving data.")
    parser.add_argument("--mode", type=str, choices=['RGB', 'hyperspectral', 'both'], help="Mode: RGB, hyperspectral, or both")
    parser.add_argument("--cap", type=int, help="Cap value for intensity")
    parser.add_argument("--index", type=int, help="Index for hyperspectral channel")
    parser.add_argument("--data", type=str, help="Path to the data directory")
    parser.add_argument("--plot", action='store_true', help="Whether to plot the graphs")
    parser.add_argument("--rotation", type=int, choices=[0, 90, 180, 270], help="Rotate the CLS matrix by the specified angle before computing stats.")
    parser.add_argument("--flip_horizontal", action='store_true', help="Flip the CLS matrix horizontally before computing stats.")


    args = parser.parse_args()
    main(args)