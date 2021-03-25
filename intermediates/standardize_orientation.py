from PIL import Image
import piexif
import sys
import os

wd = sys.argv[1]

def rotate_jpeg(filename):
    img = Image.open(filename)
    if "exif" in img.info:
        exif_dict = piexif.load(img.info["exif"])

        if piexif.ImageIFD.Orientation in exif_dict["0th"]:
            orientation = exif_dict["0th"].pop(piexif.ImageIFD.Orientation)
            exif_bytes = piexif.dump(exif_dict)

            if orientation == 2:
                img = img.transpose(Image.FLIP_LEFT_RIGHT)
            elif orientation == 3:
                img = img.rotate(180)
            elif orientation == 4:
                img = img.rotate(180).transpose(Image.FLIP_LEFT_RIGHT)
            elif orientation == 5:
                img = img.rotate(-90, expand=True).transpose(Image.FLIP_LEFT_RIGHT)
            elif orientation == 6:
                img = img.rotate(-90, expand=True)
            elif orientation == 7:
                img = img.rotate(90, expand=True).transpose(Image.FLIP_LEFT_RIGHT)
            elif orientation == 8:
                img = img.rotate(90, expand=True)

            img.save(filename, exif=exif_bytes, quality=91)
            
def croploop(wd):
    file_list = os.listdir(wd)
    os.chdir(wd)
    file_list = [x for x in file_list if '.png' in x] + [x for x in file_list if '.jpg' in x]
    
    # If cropped files already overexist, we will overwrite instead of cropping the cropped versions
    file_list = [x for x in file_list if 'cropped' not in x]
    
    # Don't want to crop blend and composite images if they (already) exist
    file_list = [x for x in file_list if 'blend' not in x]
    file_list = [x for x in file_list if 'composite' not in x]
    
    for file in file_list:
        rotate_jpeg(file)
        
def main(wd):
    croploop(wd)

if __name__ == "__main__":
    main(wd = wd)