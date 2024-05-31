from src2.Image import Image, Camera
import os
import cv2
from colour_demosaicing import demosaicing_CFA_Bayer_Malvar2004
import numpy as np
from typing import Tuple, Dict

"""
Function to process images

TODO: change function to take in array rather than image object for clearer function calls with less assumptions
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
        image.image = self._single_to_triple_channels(
            image.image, image.cam_type)
        image.proc_image = image.image.copy()

        # solar elevation angle
        mu = self._calculate_solar_elevation_angle(image.label)

        # find ftau value
        ftau = self._find_ftau(mu, image.label, image.cam_type, image.filename)

        # padding values
        if image.pad_im:
            paddings = self._calculate_padding(image.image, image.label,
                                               image.cam_type, image.filename)
            if any([pad != 0 for pad in paddings]):
                image.proc_image = self._pad_image(image.proc_image, paddings)

        # create image mask
        mask_image = self._make_image_mask(
            image.image, paddings, image.pad_im, image.filename)
        image.mask_image = mask_image

        # apply color and brightness correction
        image.proc_image = self._color_brightness_correction(
            image.proc_image, ftau, scale, scale_red, scale_blue)

        if image.filename[:3] == 'HNM':
            im_max = np.percentile(image.proc_image, 99.9) * 1.01
            image.proc_image /= im_max

        # clip the image
        image.proc_image = self._clip_image(image.proc_image, clip_low)

        # gamma correction
        image.proc_image = self._gamma_correction(image.proc_image, gamma)

        # rescale image to 8 unsigned bits
        image.im8 = np.clip(255*image.proc_image, 0, 255).astype('uint8')

    def _single_to_triple_channels(self, im_arr: np.array, cam_type: Camera) -> np.array:
        """
        Checks if the image has one color band, either demosaic or stack the
        image to make it a color image. Assumes im_arr is either a 2D or a 3D array.
        :param image: np.array, image to process
        :param cam_type: Camera, camera type
        """
        if len(im_arr.shape) == 2:
            if cam_type == Camera.NAVCAM_VCE:
                return np.stack([im_arr, im_arr, im_arr], axis=-1)
            else:
                return demosaicing_CFA_Bayer_Malvar2004(im_arr, 'RGGB')

        # if it is a three channel image, return the image
        return im_arr

    def _calculate_solar_elevation_angle(self, im_label: Dict) -> float:
        """
        Calculate the solar elevation angle
        :param image_label: meta data associate with the image
        :return: float, solar elevation angle
        """
        d = 57.296
        try:
            mu = np.sin(
                im_label['SITE_DERIVED_GEOMETRY_PARMS']['SOLAR_ELEVATION'][0]/d)
        except:
            mu = 1.0
        return mu

    def _find_ftau(self, mu: float, im_label: Dict, cam_type: Camera, fname: str) -> float:
        """
        Find ftau value using the camera type and the solar elevation angle
        :param mu: float, solar elevation angle
        :param im_label: meta data associate with the image
        :param cam_type: Camera, camera type
        :param fname: str, filename of the image
        """
        if (fname.startswith('Z') or fname.startswith('S')) \
                and 'IOF_N' in fname:
            return 1.0

        match cam_type:
            case Camera.HSF:
                return 1.0
            case Camera.HNM:
                return 1.0
            case Camera.SHERLOC:
                return 1.0
            case _:
                if im_label['LOCAL_TRUE_SOLAR_TIME_SOL'].isdigit() and\
                        int(im_label['LOCAL_TRUE_SOLAR_TIME_SOL']) <= 700:
                    tau = 0.5
                else:
                    tau = 0.8
                tau_ref = 0.3
                ftau_min = 0.2
                return np.maximum(mu * np.exp(-(tau-tau_ref)/6/mu), ftau_min)

    def _calculate_padding(self, im_arr: np.array, im_label: Dict,
                           cam_type: Camera, fname: str) -> Tuple[int, int, int, int]:
        """
        demosaic or stack the image to make it a color image
        :param im_arr: np.array, image to process
        :param im_label: meta data associate with the image
        :param cam_type: Camera, camera type
        :param fname: str, filename of the image
        :return: Tuple[int, int, int, int], padding values (left, right, top, bottom)
        """
        pad_left, pad_right, pad_top, pad_bottom = 0, 0, 0, 0

        # padding for Mastcam-Z images with non-standard sizes
        if cam_type == Camera.ZCAM_LEFT or \
                cam_type == Camera.ZCAM_RIGHT:

            full_height,  full_width = [1200, 1648]

            if im_arr.shape[0] != full_height or im_arr.shape[1] != full_width:
                first_line_sample = im_label['MINI_HEADER']['FIRST_LINE_SAMPLE']
                line_samples = im_label['MINI_HEADER']['LINE_SAMPLES']
                first_line = im_label['MINI_HEADER']['FIRST_LINE']
                lines = im_label['MINI_HEADER']['LINES']

                pad_left = first_line_sample - 1
                pad_right = full_width - line_samples - first_line_sample + 1
                pad_top = first_line - 1
                pad_bottom = full_height - lines - first_line + 1

        # padding for HAZCAM front, Hazcam rear and Navcam images
        elif cam_type in [Camera.HAZCAM_FRONT, Camera.HAZCAM_REAR, Camera.NAVCAM, Camera.NAVCAM_VCE]:

            # parse the downsample value and calc expected image size
            downsample_char = fname.split('_')[-1][3]
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
                if im_arr.shape != (full_h, full_w, c):

                    tile_first_line_sample = im_label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE_SAMPLE']
                    tile_first_line = im_label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE']

                    pad_left = np.min(tile_first_line_sample) - 1
                    pad_right = np.max(
                        full_w - np.max(tile_first_line_sample) - 1280 + 1, 0)
                    pad_top = np.min(tile_first_line) - 1
                    pad_bottom = np.max(
                        full_h - np.max(tile_first_line) - 960 + 1, 0)

        return pad_left, pad_right, pad_top, pad_bottom

    def _pad_image(self, im_: np.array, paddings: Tuple) -> np.array:
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

    def _make_image_mask(self, OG_image: np.array, paddings: Tuple,
                         pad_im: bool, im_fname: str) -> np.array:
        """
        Function to create an image mask
        :param OG_image: np.array, original image
        :param paddings: Tuple, padding values (left, right, top, bottom)
        :param pad_im: bool, if True, pad the mask image
        :param im_fname: str, filename of the image
        :return: np.array, mask image
        """
        mask_image = np.ones(OG_image.shape[:2])*255
        pad_left, pad_right, pad_top, pad_bottom = paddings

        # if there are non-zero padding values, pad the mask image
        if any([pad != 0 for pad in paddings]) and pad_im:

            mask_image = self._pad_image(mask_image, paddings)

            if pad_bottom == 0 and pad_right == 0:
                mask_image[pad_top:,
                           pad_left:][OG_image[:, :, 1] == 0] = 0
            elif pad_bottom == 0:
                mask_image[pad_top:, pad_left:-
                           pad_right][OG_image[:, :, 1] == 0] = 0
            elif pad_right == 0:
                mask_image[pad_top:-pad_bottom,
                           pad_left:][OG_image[:, :, 1] == 0] = 0
            else:
                mask_image[pad_top:-pad_bottom, pad_left:-
                           pad_right][OG_image[:, :, 1] == 0] = 0

        else:
            mask_image[OG_image[:, :, 1] == 0] = 0

        pro_mask = self._process_mask(mask_image, OG_image, im_fname)

        return pro_mask

    def _process_mask(self, mask: np.array, OG_image: np.array, im_fname: str) -> np.array:
        """
        Function to process the mask image
        :param mask: np.array, mask image
        :param OG_image: np.array, original image
        :param im_fname: str, filename of the image
        :return: np.array, processed mask image
        """
        pro_mask = mask.copy()

        # Mars2020 Mastcam-Z mask processing
        if im_fname[0] in ['Z', 'S']:
            pro_mask[:4, :] = 0
            pro_mask[-1:, :] = 0
            pro_mask[:, :24] = 0
            pro_mask[:, -17:] = 0

            # use pre-saved mask
            parent_path = os.getcwd()
            if im_fname[:2] == 'ZL':
                mask_path = os.path.join(parent_path, 'params/ZL.jpg')
            if im_fname[:2] == 'ZR':
                mask_path = os.path.join(parent_path, 'params/ZR.jpg')
            else:
                mask_path = os.path.join(parent_path, 'params/S.jpg')
            stored_mask = cv2.imread(mask_path)
            pro_mask[stored_mask[:, :, 0] < 100] = 0

        # Mars2020 SuperCam RMI mask processing
        elif im_fname[0] == 'L':

            pro_mask[OG_image == 0] = 0
            pro_mask[1800:, :, :] = 0
            pro_mask = cv2.blur(pro_mask, (20, 20))
            pro_mask[pro_mask < 255] = 0

        # HNM
        elif im_fname[:3] == 'HNM':

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
            if im_fname[0] == 'F':

                try:
                    downsample = int(im_fname.split('_')[-1][3])
                except:
                    print(
                        f"Error in downsample informaiton for filename:{im_fname}")
                    exit()

                parent_path = os.getcwd()
                if im_fname[:2] == 'FL':
                    mask_path = os.path.join(
                        parent_path, 'params/FL{}.jpg'.format(downsample))
                else:
                    mask_path = os.path.join(
                        parent_path, 'params/FR{}.jpg'.format(downsample))

                mask = cv2.imread(mask_path)
                pro_mask[mask[:, :, 0] < 100] = 0

            if 'MV' in im_fname or 'M_' in im_fname:

                parent_path = os.getcwd()
                if im_fname[:2] == 'NL':
                    mask_path = os.path.join(parent_path, 'params/NL2_vce.jpg')
                else:
                    mask_path = os.path.join(parent_path, 'params/NR2_vce.jpg')

                mask = cv2.imread(mask_path)
                pro_mask[mask[:, :, 0] < 100] = 0

        return pro_mask

    def _color_brightness_correction(self, proc_image: np.array, ftau: float, scale: float, scale_red: float, scale_blue: float) -> np:
        """
        Function to apply color and brightness correction to an image
        :param proc_image: np.array, image to process
        :param ftau: float, ftau value
        :param scale: float, scale value
        :param scale_red: float, red scale value
        :param scale_blue: float, blue scale value
        """
        proc_image_ = proc_image.copy()
        proc_image_[:, :, 0] *= scale / ftau * scale_red
        proc_image_[:, :, 1] *= scale / ftau * 1
        proc_image_[:, :, 2] *= scale / ftau * scale_blue
        return proc_image_

    def _clip_image(self, proc_image: np.array, clip_low: float) -> np.array:
        """
        Function to clip the image
        :param proc_image: np.array, image to process
        :param clip_low: float, clip low value
        """
        proc_image_ = proc_image.copy()
        proc_image_ = (proc_image_ - clip_low)/(1 - clip_low)
        proc_image_ = np.clip(proc_image_, 0, 1)

        return proc_image_

    def _gamma_correction(self, proc_image: np.array, gamma: float):
        """
        Function to apply gamma correction to an image
        :param proc_image: np.array, image to process
        :param gamma: float, gamma value
        """
        proc_image_ = proc_image.copy()
        if gamma != 1.0:
            proc_image_ = proc_image_**(1 / gamma)

        return proc_image_
