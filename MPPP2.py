"""
Example usage:  python MPPP2.py --sol 123* --input_dir ../data/zcam --output_dir ../images
"""

import glob
import argparse
from src2.Image import Image, Frame
from src2.ImProcessor import ImProcessor


def main(args):
    processor = ImProcessor(args.config_fpath)
    try:
        IMG_paths = []
        suf = 'refs_'+str(args.sol)
        
        IMG_paths     += sorted( glob.glob(args.input_dir + f"/FLF_{args.sol}_*.IMG"))
        IMG_paths     += sorted( glob.glob(args.input_dir + f"/FRF_{args.sol}_*.IMG"))
        IMG_paths     += sorted( glob.glob(args.input_dir + f"/NLF_{args.sol}_*.IMG"))
        # IMG_paths     += sorted( glob.glob(directory_input + f"/NLMV_{sol}*_*.IMG"))
        IMG_paths     += sorted( glob.glob(args.input_dir + f"/NRF_{args.sol}_*.IMG"))
        # IMG_paths     += sorted( glob.glob(directory_input + f"/NRMV_{sol}*_*.IMG"))
        IMG_paths     += sorted( glob.glob(args.input_dir + f"/*/ZL0_{args.sol}_*.IMG"))
        IMG_paths     += sorted( glob.glob(args.input_dir + f"/*/ZR0_{args.sol}_*.IMG"))
    
        IMG_paths = IMG_paths[:]
        print( len(IMG_paths), 'images\n')
        if len(IMG_paths):
            processor.process_images(img_paths = IMG_paths,
                                output_dir = args.output_dir,
                                suf = suf.replace('*',''),
                                find_offset_mode=False,
                                frame = "site",
                                angles = 'opk',
                                save_im = True)
    except:
        print(f"FAILED TO PROCESS! sols:{args.sol}* in {args.input_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sol', type=str, required=True, help="Regex pattern of the sol to be analyse")
    parser.add_argument('--input_dir', type=str, required=True, help="Path to your input directory")
    parser.add_argument('--output_dir', type=str, required=True, help="Path to your output directory")
    parser.add_argument('--config_fpath', type=str, default='./params/image_processing_config.json', 
                        help='filepath to image processor config file')

    args = parser.parse_args()
    main(args)
