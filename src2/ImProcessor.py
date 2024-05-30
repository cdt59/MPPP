from src2.Image import Image, Camera
import os
import cv2
import colour_demosaicing
import numpy as np
from typing import Tuple

"""
Function to process images

TODO: look into image and im in original MPPP, make sure they are referenced properly here
TODO: change function to take in array rather than image object for clearer funciton calls with less assumptions
"""

clip_low = 0.01
scale_red = 1.0
scale_blue = 1.0
scale = 12
gamma = 2


class ImProcessor():

    def process_image(self, image: Image):
        """
        Function to process image. This funciton modifies the image object in place.
        :param image: np.array, image to process
        """
        # Don't understand what this is
        if image.filename.split('_N')[0][-3:] == 'RZS':
            ftau = np.float32(
                image.label['DERIVED_IMAGE_PARMS']['RAD_ZENITH_SCALING_FACTOR'])
            image.image *= ftau

        # if the image has one color band, either demosaic or stack the image to
        # make it a color image
        self.single_to_triple_channels(image)
        image.proc_image = image.image.copy()

        # solar elevation angle
        mu = self.calculate_solar_elevation_angle(image)

        # find ftau value
        ftau = self.find_ftau(image, mu)

        # padding values
        if image.pad_im:
            paddings = self.calculate_padding(image)
            if any([pad != 0 for pad in paddings]):
                image.proc_image = self.pad_image(image.proc_image, paddings)

        # create image mask
        mask_image = self.make_image_mask(image, paddings)
        image.mask_image = mask_image

        # apply color and brightness correction
        self.color_brightness_correction(
            image, ftau, scale, scale_red, scale_blue)

        if image.filename[:3] == 'HNM':
            im_max = np.percentile(image.proc_image, 99.9) * 1.01
            image.proc_image /= im_max

        # clip the image
        self.clip_image(image, clip_low)

        # gamma correction
        self.gamma_correction(image, gamma)

        # rescale image to 8 unsigned bits
        image.im8 = np.clip(255*image.proc_image, 0, 255).astype('uint8')

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
            mu = 1.0
        return mu

    def find_ftau(self, image: Image, mu: float) -> float:
        """
        Find ftau value using the camera type and the solar elevation angle
        :param image: np.array, image to process
        :param mu: float, solar elevation angle
        """
        if (image.filename.startswith('Z') or image.filename.startswith('S')) \
                and 'IOF_N' in image.filename:
            return 1.0

        match image.cam_type:
            case Camera.HSF:
                return 1.0
            case Camera.HNM:
                return 1.0
            case Camera.SHERLOC:
                return 1.0
            case _:
                tau = 0.5 if image.sol <= 700 else 0.8
                tau_ref = 0.3
                ftau_min = 0.2
                return np.maximum(mu * np.exp(-(tau-tau_ref)/6/mu), ftau_min)

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
                image.cam_type == Camera.NAVCAM or\
                image.cam_type == Camera.NAVCAM_VCE:

            # parse the downsample value and calc expected image size
            downsample_char = image.filename.split('_')[-1][3]
            downsample = int(
                downsample_char) if downsample_char.isdigit() else None
            if (downsample is not None) and (downsample in [0, 1, 2]):
                if downsample == 0:
                    full_h, full_w, c = 3840, 5120, 3
                elif downsample == 1:
                    full_h, full_w, c = 1920, 2560, 3
                elif downsample == 2:
                    full_h, full_w, c = 960, 1280, 3

                # if the image is not the expected size, calculate padding
                if image.image.shape != (full_h, full_w, c):

                    tile_first_line_sample = image.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE_SAMPLE']
                    tile_first_line = image.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE']

                    pad_left = np.min(tile_first_line_sample) - 1
                    pad_right = np.max(
                        full_w - np.max(tile_first_line_sample) - 1280 + 1, 0)
                    pad_top = np.min(tile_first_line) - 1
                    pad_bottom = np.max(
                        full_h - np.max(tile_first_line) - 960 + 1, 0)

        return pad_left, pad_right, pad_top, pad_bottom

    def pad_image(self, im_: np.array, paddings: Tuple) -> np.array:
        """
        Function to pad an image
        :param image: np.array, image to process
        :param paddings: padding values (left, right, top, bottom)
        :return: np.array, padded image
        """
        im = im_.copy()
        if len(im.shape) == 3:
            im = np.hstack([np.zeros((im.shape[0], paddings[0], 3)),
                           im, np.zeros((im.shape[0],  paddings[1], 3)), ])
            im = np.vstack([np.zeros((paddings[2], im.shape[1],  3)),
                           im, np.zeros((paddings[3],  im.shape[1], 3)), ])
        else:
            im = np.hstack([np.zeros((im.shape[0],   paddings[0])),
                           im, np.zeros((im.shape[0],  paddings[1])), ])
            im = np.vstack([np.zeros((paddings[2], im.shape[1],)),
                           im, np.zeros((paddings[3],  im.shape[1])), ])
        return im

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
                mask_image[pad_top:,
                           pad_left:][image.image[:, :, 1] == 0] = 0
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

        pro_mask = self._process_mask(mask_image, image)

        return pro_mask

    def _process_mask(self, mask: np.array, image: Image) -> np.array:
        """
        Function to process the mask image
        :param mask: np.array, mask image
        :param image: np.array, image that is tied to this mask. Assumes image was already padded
        :return: np.array, processed mask image
        """
        pro_mask = mask.copy()

        # Mars2020 Mastcam-Z mask processing
        if image.filename[0] in ['Z', 'S']:

            pro_mask[:4, :] = 0
            pro_mask[-1:, :] = 0
            pro_mask[:, :24] = 0
            pro_mask[:, -17:] = 0

            # use pre-saved mask
            parent_path = os.getcwd()
            if image.filename[:2] == 'ZL':
                mask_path = os.path.join(parent_path, 'params/ZL.jpg')
            if image.filename[:2] == 'ZR':
                mask_path = os.path.join(parent_path, 'params/ZL.jpg')
            else:
                mask_path = os.path.join(parent_path, 'params/S.jpg')
            mask = cv2.imread(mask_path)
            pro_mask[mask[:, :, 0] < 100] = 0

        # Mars2020 SuperCam RMI mask processing
        elif image.filename[0] == 'L':

            pro_mask[image.image == 0] = 0
            pro_mask[1800:, :, :] = 0
            pro_mask = cv2.blur(pro_mask, (20, 20))
            pro_mask[pro_mask < 255] = 0

        # HNM
        elif image.filename[:3] == 'HNM':

            parent_path = os.getcwd()
            mask_path = os.path.join(parent_path, 'params/HNM.jpg')
            mask = cv2.imread(mask_path)
            pro_mask[mask[:, :, 0] < 100] = 0

        # Mars2020 Ecam mask processing
        else:
            pro_mask[:2, :] = 0
            pro_mask[-2:, :] = 0
            pro_mask[:, :3] = 0
            pro_mask[:, -3:] = 0

            # use pre-saved mask
            if image.filename[0] == 'F':

                try:
                    downsample = int(image.filename.split('_')[-1][3])
                except:
                    print(
                        f"Error in downsample informaiton for filename:{image.filename}", error=True)
                    exit()

                parent_path = os.getcwd()
                if image.filename[:2] == 'FL':
                    mask_path = os.path.join(
                        parent_path, 'params/FL{}.jpg'.format(downsample))
                else:
                    mask_path = os.path.join(
                        parent_path, 'params/FR{}.jpg'.format(downsample))

                mask = cv2.imread(mask_path)
                pro_mask[mask[:, :, 0] < 100] = 0

            if 'MV' in image.filename or 'M_' in image.filename:

                parent_path = os.getcwd()
                if image.filename[:2] == 'NL':
                    mask_path = os.path.join(parent_path, 'params/NL2_vce.jpg')
                else:
                    mask_path = os.path.join(parent_path, 'params/NR2_vce.jpg')

                mask = cv2.imread(mask_path)
                pro_mask[mask[:, :, 0] < 100] = 0

        return pro_mask

    def color_brightness_correction(self, image: Image, ftau: float, scale: float, scale_red: float, scale_blue: float) -> None:
        """
        Function to apply color and brightness correction to an image
        :param image: np.array, image to process
        :param ftau: float, ftau value
        :param scale: float, scale value
        :param scale_red: float, red scale value
        :param scale_blue: float, blue scale value
        """
        image.proc_image[:, :, 0] *= scale / ftau * scale_red
        image.proc_image[:, :, 1] *= scale / ftau * 1
        image.proc_image[:, :, 2] *= scale / ftau * scale_blue

    def clip_image(self, image: Image, clip_low: float) -> np.array:
        """
        Function to clip the image
        :param image: np.array, image to process
        """
        image.proc_image = (image.proc_image - clip_low)/(1 - clip_low)
        image.proc_image = np.clip(image.proc_image, 0, 1)

    def gamma_correction(self, image: Image, gamma: float):
        """
        Function to apply gamma correction to an image
        :param image: np.array, image to process
        :param gamma: float, gamma value
        """
        if gamma != 1.0:
            image.proc_image = image.proc_image**(1 / gamma)
