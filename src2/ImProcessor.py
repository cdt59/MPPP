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
            image.image *= ftau

        # if the image has one color band, either demosaic or stack the image to
        # make it a color image
        self.single_to_triple_channels(image)

        # solar elevation angle
        mu = self.calculate_solar_elevation_angle(image)

        # find ftau value
        ftau = self.find_ftau(image, mu)

        # padding values
        paddings = self.calculate_padding(image)
        if any([pad != 0 for pad in paddings]) and image.pad_im:
            image.image = self.pad_image(image.image, paddings)

        # create image mask
        mask_image = self.make_image_mask(image, paddings)

    def single_to_triple_channels(self, image: Image) -> None:
        """
        Checks if the image has one color band, either demosaic or stack the
        image to make it a color image
        :param image: np.array, image to process
        """
        if len(image.image.shape) == 2:
            if image.cam_type == Camera.NAVCAM_VCE:
                image.image = np.stack([image.image, image.image, image.image],
                                       axis=-1)
            else:
                image.image = colour_demosaicing.demosaicing_CFA_Bayer_Malvar2004(
                    image.image, 'RGGB')

    def calculate_solar_elevation_angle(self, image: Image) -> float:
        """
        Calculate the solar elevation angle
        :param image: np.array, image to process
        :return: float, solar elevation angle
        """
        d = 57.296
        try:
            mu = np.sin(
                image.label['SITE_DERIVED_GEOMETRY_PARMS']['SOLAR_ELEVATION'][0]/d)
        except:
            mu = 1
        return mu

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
        demosaic or stack the image to make it a color image
        :param image: np.array, image to process
        :return: left, right, top, bottom padding values
        """
        pad_left, pad_right, pad_top, pad_bottom = 0, 0, 0, 0

        # padding for Mastcam-Z images with non-standard sizes
        if image.cam_type == Camera.ZCAM_LEFT or \
                image.cam_type == Camera.ZCAM_RIGHT:

            full_height,  full_width = [1200, 1648]

            if image.image.shape[0] != full_height or image.image.shape[1] != full_width:
                first_line_sample = image.label['MINI_HEADER']['FIRST_LINE_SAMPLE']
                line_samples = image.label['MINI_HEADER']['LINE_SAMPLES']
                first_line = image.label['MINI_HEADER']['FIRST_LINE']
                lines = image.label['MINI_HEADER']['LINES']

                pad_left = first_line_sample - 1
                pad_right = full_width - line_samples - first_line_sample + 1
                pad_top = first_line - 1
                pad_bottom = full_height - lines - first_line + 1

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

    def pad_image(self, im: np.array, paddings: Tuple) -> np.array:
        """
        Function to pad an image
        :param image: np.array, image to process
        :param paddings: padding values (left, right, top, bottom)
        :return: np.array, padded image
        """
        pad_left, pad_right, pad_top, pad_bottom = paddings
        if len(im.shape) == 3:
            padded_image = np.hstack([np.zeros((im.shape[0], pad_left, 3)), im,
                                      np.zeros((im.shape[0], pad_right, 3)), ])
            padded_image = np.vstack([np.zeros((pad_top, padded_image.shape[1], 3)),
                                      padded_image,
                                      np.zeros((pad_bottom, padded_image.shape[1], 3)), ])
        else:
            padded_image = np.hstack([np.zeros((im.shape[0], pad_left)), im,
                                      np.zeros((im.shape[0], pad_right)), ])
            padded_image = np.vstack([np.zeros((pad_top, padded_image.shape[1])),
                                      padded_image,
                                      np.zeros((pad_bottom, padded_image.shape[1])), ])
        return padded_image

    def make_image_mask(self, image: Image, paddings: Tuple) -> np.array:
        """
        Function to create an image mask
        :param image: np.array, image to process.
        :param paddings: padding values (left, right, top, bottom)
        :return: np.array, 8-bit image mask
        """
        mask_image = image.mask_image.copy()
        pad_left, pad_right, pad_top, pad_bottom = paddings

        # if there are non-zero padding values, pad the mask image
        if any([pad != 0 for pad in paddings]) and image.pad_im:
            mask_image = self.pad_image(mask_image, paddings)

            if pad_bottom == 0 and pad_right == 0:
                mask_image[pad_top:, pad_left:][image.image[:, :, 1] == 0] = 0
            elif pad_bottom == 0:
                mask_image[pad_top:, pad_left:-
                           pad_right][image.image[:, :, 1] == 0] = 0
            elif pad_right == 0:
                mask_image[pad_top:-pad_bottom,
                           pad_left:][image.image[:, :, 1] == 0] = 0
            else:
                mask_image[pad_top:-pad_bottom, pad_left:-
                           pad_right][image.image[:, :, 1] == 0] = 0

        else:
            mask_image[image.image[:, :, 1] == 0] = 0

        return mask_image
