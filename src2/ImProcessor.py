from src2.Image import Image, Camera, Frame
import os
import cv2
from colour_demosaicing import demosaicing_CFA_Bayer_Malvar2004
import numpy as np
import json
import time
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List


class ImProcessor():
    """
    Class that handles image processing

    Example usage:
        processor = ImProcessor()
        img = Image("./path/to/image")
        processed_img = processor.process_image(img)
    """

    def __init__(self, params_path: str):
        """
        :param params_path: str, path to json file where image processing parameters are stored
        """
        with open(params_path, 'r') as f:
            self.params = json.load(f)

    def process_images(self, img_paths: List[str], output_dir: str, suf: str, find_offset_mode: bool = False, frame: Frame = Frame.SITE3, angles: str = 'opk', save_im: bool = True):
        """
        Processes a list of images
        img_paths: list of paths to the images to be processed
        output_dir: str, path to the output directory
        suf: str, suffix used in the saved CSV file
        find_offset_mode: indicate if offsets should be found. If so, save the offsets in a txt file
        frame: str, frame of the image
        angles: one of 'opk' or 'ypr' to indicate angle type to be returned
        save_im: Indicates if the processed image should be saved
        """
        pos_lines = []
        error_lines = []
        veh_XYZs = []
        im_XYZs = []
        veh_azs = []
        im_azs = []
        im_els = []
        sols = []
        rmcs = []

        print(f"Number of images: {len(img_paths)}")

        for i, img_path in enumerate(img_paths):
            try:
                img = Image(img_path, frame=frame)
                print(f"{i} opened {os.path.basename(img_path)}")
                IMG_loaded = True

            except:
                print(f"{os.path.basename(img_path)} failed to process!", end='\n\n')
                error_lines.append(os.path.basename(img_path)+'\n')
                IMG_loaded = False

            if IMG_loaded:
                # color processing parameters
                scale_scale = self.params['scales']['scale_scale']
                img.clip_low = self.params['clip_low']
                img.gamma = self.params['gamma']
                img.pad_im = self.params['pad_im']
                img.save_im = save_im
                img.save_mask = self.params['save_mask']
                img.find_offset_mode = find_offset_mode

                match img.cam_type:
                    case Camera.ZCAM_LEFT | Camera.ZCAM_RIGHT:
                        scales = self.params['scales']['zcam']
                        img.clip_low = self.params['clip_low_z']
                        img.pad_im = self.params['pad_im_z']

                    case Camera.SHERLOC:
                        scales = self.params['scales']['sherloc']
                        img.clip_low = 0.0

                    case Camera.RMI:
                        scales = self.params['scales']['rmi']

                    case Camera.NAVCAM:
                        scales = self.params['scales']['navcam']

                    case Camera.NAVCAM_VCE:
                        scales = self.params['scales']['navcam_vce']
                        img.clip_low = 0.0
                        img.gamma = 1.0

                    case Camera.HAZCAM_FRONT:
                        scales = self.params['scales']['hazcam_front']
                        img.clip_low = img.clip_low/2

                    case Camera.HAZCAM_REAR:
                        scales = self.params['scales']['hazcam_rear']
                        img.clip_low = img.clip_low/2

                    case Camera.HSF:
                        scales = self.params['scales']['hsf']

                    case Camera.HNM:
                        scales = self.params['scales']['hnm']
                        img.clip_low = 0.3
                        img.gamma = 1.0

                img.scale = scales['scale_factor'] * scale_scale
                img.scale_red = scales['scale_red']
                img.scale_blue = scales['scale_blue']

                file_extension = self.params['file_extension']

                img.focus_mc = -1
                img.zoom_mc = -1

                # create save directory
                img.save_path_full = self.make_save_path(
                    img.IMG_path, output_dir, fullpath=True, file_extension=file_extension)
                img.save_path = self.make_save_path(
                    img.IMG_path, output_dir, fullpath=False)
                img.save_name = img.save_path_full.split('/')[-1]
                csv_save_path = img.save_path_full

                # process and save image
                if img.save_im:
                    img = self.process_image(img, img.scale, img.scale_red,
                                             img.scale_blue, img.clip_low, img.gamma)
                    if img.save_mask:
                        img.im8a = cv2.cvtColor(img.im8, cv2.COLOR_BGR2RGBA)
                        img.im8a[:, :, 3] = img.mask_image
                        cv2.imwrite(img.save_path_full, img.im8a)
                    else:
                        cv2.imwrite(img.save_path_full, img.im8[:, :, ::-1])

                # TODO: ADD CMOD Implementation

        csv_save_path = os.path.dirname(
            csv_save_path)+'/positions_'+suf+'_'+str(frame)+'_' +\
            time.strftime("%Y%m%d-%H%M%S") + '.txt'
        with open(csv_save_path, 'w') as file:
            for pos_line in pos_lines:
                file.write(pos_line)

        print(f"saved {csv_save_path}")

        self.plot_image_locations(img_paths, im_XYZs, veh_XYZs,
                                  veh_azs, im_azs, im_els)

        if find_offset_mode:
            sites = [rmcs[i][0] for i in range(len(rmcs))[::-1]]
            drives = [rmcs[i][1] for i in range(len(rmcs))[::-1]]
            Xs = [veh_XYZs[i][0] for i in range(len(veh_XYZs))[::-1]]
            Ys = [veh_XYZs[i][1] for i in range(len(veh_XYZs))[::-1]]
            Zs = [veh_XYZs[i][2] for i in range(len(veh_XYZs))[::-1]]

            table = np.stack([sols[::-1], sites, drives, Xs, Ys, Zs], axis=1)
            np.round(table, 4)

            np.savetxt(output_dir+"/offsets_" +
                       suf+".csv", table, delimiter="\t")

    def process_image(self, image: Image, scale: float, scale_red: float, scale_blue: float, clip_low: float, gamma: float) -> Image:
        """
        Function to process a single image
        :param image: Image, image object to be processed
        :param scale: float, scale value
        :param scale_red: float, red scale value
        :param scale_blue: float, blue scale value
        :param clip_low: float, lower bound value for use in clipping
        :param gamma: float, gamma value for color correction
        :return: Image, processed image object
        """
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
        else:
            paddings = [0, 0, 0, 0]

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

        return image

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

        pad_left = max(0, pad_left)
        pad_right = max(0, pad_right)
        pad_top = max(0, pad_top)
        pad_bottom = max(0, pad_bottom)

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

    def make_save_path(self, IMG_path, directory_output, fullpath=True, file_extension='.png'):
        """
        make_save_path sorts the images into an output directory organized by camera type and each 100 sols of the mission
        """

        filename = os.path.basename(IMG_path)
        sol = int(filename[4:8])
        camera = filename[0]
        mission = 'Mars2020'  # mission name is hardcoded for now

        if camera in ['F', 'N', 'R']:
            camera_type = 'eng'
        elif camera in ['H']:
            camera_type = 'heli'
        elif camera in ['Z', 'L', 'S']:
            camera_type = 'sci'

        sol_floor_100 = int(np.floor(sol/100) * 100)
        sol_range_100 = str(sol_floor_100).zfill(4) + '-' + \
            str(sol_floor_100).zfill(4)[:2] + '99'

        save_path = directory_output + '/sols_' + sol_range_100 + '_' + camera_type

        if not os.path.exists(save_path):
            # Create a new directory because it does not exist
            os.makedirs(save_path)
            print("The new directory is created: ", save_path)

        if fullpath:
            return save_path + '/' + filename.split('.')[0] + file_extension
        else:
            return save_path

    def plot_image_locations(self, IMG_paths, im_xyzs, rover_xyzs, rover_rots, im_azs, im_els):
        '''
        plot_image_locations displays the Northing vs Easting locations of each image and rover position

        future work: replace the input arrays with a single pandas dataframe
        '''

        plt.figure(figsize=[12, 8])

        scale = np.round(np.std(np.array(rover_xyzs), axis=0).max()/4+1)

        for i in range(len(im_xyzs)):

            filename = os.path.basename(IMG_paths[i])

            marker = '*k'
            if filename[:2] in ['FL', 'RL']:
                marker = 'ob'
            if filename[:2] in ['FR', 'RR']:
                marker = 'or'
            if filename[:2] == 'NL':
                marker = 'sb'
            if filename[:2] == 'NR':
                marker = 'sr'
            if filename[:2] == 'ZL':
                marker = '^b'
            if filename[:2] == 'ZR':
                marker = '^r'

            if 'MV' not in IMG_paths[i]:

                plt.plot(rover_xyzs[i][0], rover_xyzs[i][1], color='k',    marker=(
                    4, 0, 45 + rover_rots[i]), ms=30, )
                plt.plot(rover_xyzs[i][0], rover_xyzs[i][1], color='gray', marker=(
                    3, 0, 120 + rover_rots[i]), ms=20, )

                sol = os.path.basename(IMG_paths[i])[4:8]
                if i > 1:
                    if sol != os.path.basename(IMG_paths[i-1])[4:8]:
                        plt.text(rover_xyzs[i][0]+scale/4, rover_xyzs[i][1]+scale/4, 'Sol ' + sol,
                                 bbox=dict(facecolor='w', alpha=0.5, edgecolor='w'), size='large')

                if i > 1 and os.path.basename(IMG_paths[i])[:2] == 'NLF_':
                    plt.plot([rover_xyzs[i][0], rover_xyzs[i][1]], [
                        rover_xyzs[i][0], rover_xyzs[i][1]], '--', color='gray')

                cos_az = np.cos(im_azs[i]/57.3)
                sin_az = np.sin(im_azs[i]/57.3)
                cos_el = np.cos(im_els[i]/57.3)

                if os.path.basename(IMG_paths[i])[0] == 'Z':
                    plt.arrow(im_xyzs[i][0], im_xyzs[i][1], scale*cos_el*sin_az, scale*cos_el*cos_az,
                              color=marker[1], lw=int(scale/32), linestyle='dashed')
                else:
                    plt.arrow(im_xyzs[i][0], im_xyzs[i][1], scale*cos_el*sin_az, scale*cos_el*cos_az,
                              color=marker[1], lw=int(scale/32))

            plt.plot(im_xyzs[i][0], im_xyzs[i][1], marker)

        plt.axis('equal')
        plt.xlim([np.round(plt.gca().get_xlim()[0])-3,
                  np.round(plt.gca().get_xlim()[1])+3])
        plt.ylim([np.round(plt.gca().get_ylim()[0])-3,
                  np.round(plt.gca().get_ylim()[1])+3])

        plt.xlabel('Easting Site Frame')
        plt.ylabel('Northing Site Frame')
