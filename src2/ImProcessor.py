from src2.Image import Image, Camera
import colour_demosaicing
import numpy as np
from typing import Tuple

"""
Function to process images

Look into image_process()
"""


class ImProcessor():

    def process_image(self, image: Image):
        """
        Function to process image. This funciton modifies the image object in place.
        :param image: np.array, image to process
        """
        # Don't understand what this is
        if self.filename.split('_N')[0][-3:] == 'RZS':
            ftau = np.float32(
                image.label['DERIVED_IMAGE_PARMS']['RAD_ZENITH_SCALING_FACTOR'])
            image *= ftau

        # if the image has one color band, either demosaic or stack the image to
        # make it a color image
        if len(image.shape) == 2:
            if image.cam_type == Camera.NAVCAM_VCE:
                image = np.stack([image, image, image], axis=-1)
            else:
                image = colour_demosaicing.demosaicing_CFA_Bayer_Malvar2004(
                    image, 'RGGB')

        # solar elevation angle
        d = 57.296
        try:
            mu = np.sin(
                image.label['SITE_DERIVED_GEOMETRY_PARMS']['SOLAR_ELEVATION'][0]/d)
        except:
            mu = 1

        # find ftau value
        ftau = self.find_ftau(image, mu)

        down_sample = image.filename.split('_')[-1][3]

        # default padding values
        pad_left, pad_right, pad_top, pad_bottom = 0, 0, 0, 0

    def find_ftau(self, image: Image, mu: float) -> float:
        """
        Find ftau value using the camera type and the solar elevation angle
        :param image: np.array, image to process
        :param mu: float, solar elevation angle
        """
        match image.cam_type:
            case Camera.HSF:
                ftau = 1.0
            case Camera.HNM:
                ftau = 1.0
            case Camera.SHERLOC:
                ftau = 1.0
            case _:
                tau = 0.5 if image.sol <= 700 else 0.8
                tau_ref = 0.3
                tau_min = 0.2
                ftau = np.maximum(mu * np.exp(-(tau-tau_ref)/6/mu), tau_min)

        return ftau

    def calculate_padding(self, image: Image) -> Tuple[int, int, int, int]:
        """
        Caculates the amount of padding for a given image
        :param image: np.array, image to process
        :return: left, right, top, bottom padding values
        """
        pad_left, pad_right, pad_top, pad_bottom = 0, 0, 0, 0

        # padding for Mastcam-Z images with non-standard sizes
        if image.cam_type == Camera.ZCAM_LEFT or \
                image.cam_type == Camera.ZCAM_RIGHT:

            full_height,  full_width = [1200, 1648]

            if image.image.shape[0] != full_height or image.image.shape[1] != full_width:
                pad_left = image.label['MINI_HEADER']['FIRST_LINE_SAMPLE'] - 1
                pad_right = full_width - image.label['MINI_HEADER']['LINE_SAMPLES'] -\
                    image.label['MINI_HEADER']['FIRST_LINE_SAMPLE'] + 1
                pad_top = image.label['MINI_HEADER']['FIRST_LINE'] - 1
                pad_bottom = full_height - image.label['MINI_HEADER']['LINES'] -\
                    image.label['MINI_HEADER']['FIRST_LINE'] + 1

        # padding for HAZCAM front, Hazcam rear and Navcam images
        elif image.cam_type == Camera.HAZCAM_FRONT or\
                image.cam_type == Camera.HAZCAM_REAR or\
                image.cam_type == Camera.NAVCAM:

            # parse the downsample value and calc expected image size
            downsample_char = image.filename.split('_')[-1][3]
            downsample = int(
                downsample_char) if downsample_char.isdigit() else None
            full_h, full_w, c = 3840/2**downsample, 5120/2**downsample, 3

            # if the image is not the expected size, calculate padding
            if (downsample is not None) and (image.image.shape != (full_h, full_w, c)):
                tile_first_line_sample = image.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE_SAMPLE']
                tile_first_line = image.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE']

                pad_left = np.min(tile_first_line_sample) - 1
                pad_right = np.max(
                    full_w - np.max(tile_first_line_sample) - 1280 + 1, 0)
                pad_top = np.min(tile_first_line) - 1
                pad_bottom = np.max(
                    full_h - np.max(tile_first_line) - 960 + 1, 0)

        return pad_left, pad_right, pad_top, pad_bottom
    
    
