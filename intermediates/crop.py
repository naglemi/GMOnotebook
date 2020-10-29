import os
import sys
from PIL import Image

wd = sys.argv[1] #"/scratch2/NSF_GWAS/deeplab/deeplab/dataset/JPEG/EJ/"

def cropshrinkIMG(file_path):
    
    try:
        file_in = Image.open(file_path).convert("RGB")
        file_in = file_in.resize((900,900))
        file_in = file_in.rotate(270)
        file_in = file_in.crop((0, 150, 900, 750))

        file_in.save(file_path.replace('rgb', 'rgb_cropped'))
    except: 
        print('Corrupt ' + str(file_path))

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
        cropshrinkIMG(file)

def main(wd):
    croploop(wd)

if __name__ == "__main__":
    main(wd = wd)
