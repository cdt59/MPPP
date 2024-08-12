from src2.Image import Image, Frame
from src2.ImProcessor import ImProcessor

# create a processor object
processor = ImProcessor("./params/image_processing_config.json")

processor.process_images(img_paths=["./data/SOL01151\FLF_1151_0769122851_275RAD_N0522638FHAZ00206_0A0295J01.IMG"],
                         output_dir="./test_results",
                         suf='refs_709',
                         find_offset_mode=True,
                         frame=Frame.SITE,
                         angles='opk',
                         save_im=True)
