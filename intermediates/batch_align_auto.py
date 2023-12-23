import numpy as np
import argparse
from skimage import transform
import PIL
import PIL.Image as Image
import PIL.ImageOps as ImageOps
from PIL import ImageEnhance
import glob
import os
#from spectral import *

# Function to calculate the moving average used in Version 1
def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

# Gaussian smoothing as used in Version 1
def gaussian_smoothing(data, sigma):
    return np.array(Image.fromarray(data).filter(ImageFilter.GaussianBlur(radius=sigma)))

def compute_homography_matrix(src_borders, dst_borders):
#     print("Source borders before computing matrix: ")
#     print(src_borders)

#     print("Destination borders before computing matrix: ")
#     print(dst_borders)
    
    src_coords = np.array([
        [src_borders[0], src_borders[1]],
        [src_borders[2], src_borders[1]],
        [src_borders[0], src_borders[3]],
        [src_borders[2], src_borders[3]]
    ])

    dst_coords = np.array([
        [dst_borders[0], dst_borders[1]],
        [dst_borders[2], dst_borders[1]],
        [dst_borders[0], dst_borders[3]],
        [dst_borders[2], dst_borders[3]]
    ])

    tform = transform.estimate_transform('projective', src_coords, dst_coords)
    return tform

def process_image(image_path, tform, flip, rotation):
    rgb = Image.open(image_path)
    
    if rotation:
        rgb = rgb.transpose(Image.ROTATE_90)

    if flip:
        rgb = rgb.transpose(Image.FLIP_TOP_BOTTOM)
    
    rgb_np = np.array(rgb)

#     print("Will apply matrix:")
#     print(tform)
    
    tf_img = transform.warp(rgb_np, np.linalg.inv(tform), order=0)
    tf_img_pil = Image.fromarray((tf_img * 255).astype(np.uint8))
    return tf_img_pil

def color_invert_and_transform(image, matrix):
    inverted_image = ImageOps.invert(image)
    transformed_image = inverted_image.convert("RGB", matrix)
    return transformed_image

def create_overlay(rgb_colored, hyp_colored, alpha=0.3):
    return Image.blend(rgb_colored, hyp_colored, alpha)

def main():
    print("pillow version:")
    print(Image.__version__)
    
    parser = argparse.ArgumentParser(description="Script for batch transforming images using homography matrix.")
    parser.add_argument('--left_source', type=int, required=True)
    parser.add_argument('--top_source', type=int, required=True)
    parser.add_argument('--right_source', type=int, required=True)
    parser.add_argument('--bottom_source', type=int, required=True)
    parser.add_argument('--left_target', type=int, required=True)
    parser.add_argument('--top_target', type=int, required=True)
    parser.add_argument('--right_target', type=int, required=True)
    parser.add_argument('--bottom_target', type=int, required=True)
    parser.add_argument('--img_dir', type=str, required=True)
    parser.add_argument('--hyp_path', type=str, required=True)
    parser.add_argument('--channel_index', type=int, default=130)
    parser.add_argument('--overlay_dir', type=str, required=True)

    args = parser.parse_args()

    src_borders = [args.left_source, args.top_source, args.right_source, args.bottom_source]
    dst_borders = [args.left_target, args.top_target, args.right_target, args.bottom_target]

    tform = compute_homography_matrix(src_borders, dst_borders)

    # Combine the lists of file paths for both file types
    image_paths = glob.glob(args.img_dir + "*_rgb.jpg") + glob.glob(args.img_dir + "*_segment_uncropped.png")

    # Loop through the combined list of file paths
    for image_path in image_paths:
        print("Performing alignment for image:")
        print(image_path)
        rgb = Image.open(image_path)
        

        # Process RGB image
        #rgb = rgb.transpose(Image.FLIP_TOP_BOTTOM).transpose(Image.ROTATE_90)
        rgb = rgb.transpose(Image.FLIP_TOP_BOTTOM)

        # Apply transformation using skimage (since PIL doesn't support projective transforms)
        rgb_np = np.array(rgb)
        print("Will apply matrix:")
        print(tform)
        tf_img = transform.warp(rgb_np, np.linalg.inv(tform), order=0)
        tf_img_pil = Image.fromarray((tf_img * 255).astype(np.uint8))
        rgb = tf_img_pil

        hyp_path = os.path.join(args.hyp_path, os.path.basename(image_path).replace('_rgb.jpg', '_Fluorescence.png'))
        
        hyp = Image.open(hyp_path)
        hyp = hyp.transpose(Image.ROTATE_270)

        width, height = hyp.size

        rgb_cropped = rgb.crop((0, 0, width, height))
        # Apply color transformations as in Version 1
        
        rgb_cropped = PIL.ImageOps.invert(rgb_cropped)

        # A 12-value tuple which is a transform matrix for dropping 
        # green channel (in this case)
        matrix = ( 1, 0, 0, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0)

        rgb_colored = rgb_cropped.convert("RGB", matrix)
        
        hyp = PIL.ImageOps.invert(hyp)

        matrix = ( 0, 0, 0, 0,
                   0, 0, 1, 0,
                   0, 0, 1, 0)

        hyp_colored = hyp.convert("RGB", matrix)
        hyp_colored = PIL.ImageEnhance.Contrast(hyp_colored)
        hyp_colored = hyp_colored.enhance(1.5)

        overlay = PIL.Image.blend(rgb_colored, hyp_colored, alpha = 0.3)

        rgb_save_path = image_path.replace('_rgb.jpg', '_rgb_processed.png')
        overlay_save_path = os.path.join(args.overlay_dir, os.path.basename(image_path).replace('_rgb.jpg', '_overlay.jpg'))

        rgb_cropped.save(rgb_save_path)
        overlay.save(overlay_save_path)

if __name__ == "__main__":
    main()
