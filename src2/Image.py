from enum import Enum
from planetaryimage import PDS3Image
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import colour_demosaicing
import os
import cv2
import time
import glob
import pandas as pd


class Frame(Enum):
    SITE3 = 'site3'
    RNAV = 'rnav'


class Camera(Enum):
    ZCAM_LEFT = "ZCAM_LEFT"  # Mars2020 Mastcam-Z left camera
    ZCAM_RIGHT = "ZCAM_RIGHT"  # Mars2020 Mastcam-Z right camera
    SHERLOC = "SHERLOC"  # SHERLOC Watson camera
    RMI = "RMI"  # Mars2020 supercam RMI
    NAVCAM = "NAVCAM"  # Mars2020 Navcam
    NAVCAM_VCE = "NAVCAM_VCE"  # Mars2020 Navcam VCE
    HAZCAM_FRONT = "HAZCAM_FRONT"  # Mars2020 Hazcam front
    HAZCAM_REAR = "HAZCAM_REAR"  # Mars2020 Hazcam rear
    HSF = "HELI_RTE"  # Heli Ingenuity RTE
    HNM = "HELI_NAV"  # Heli Ingenuity Navcam
    NONE = "NONE"  # Failure to parse camera type


"""
Class that describes an IMG image. Includes both the image data and the metadata.
"""


class Image:

    def __init__(self, IMG_path: str, just_label: bool = False, frame: Frame = Frame.SITE3):
        """
        Constructor for Image class.

        :param IMG_path: str, path to the IMG file
        :param just_label: bool, if True, the image data will not be loaded
        :param frame: str, frame of the image
        """

        # parse the camera, camera type and sol
        self.IMG_path = IMG_path
        self.filename = os.path.basename(IMG_path)
        self.name = self.filename[:self.filename.find('.')]
        self.cam = self.filename[:2] + \
            self.filename[45:48]  # TODO: remove redundancy
        self.sol = int(self.filename[4:8])
        self.cam_type = self._parse_cam_type(self.filename)

        # PDS image header metadata
        self.label = PDS3Image.open(IMG_path).label

        # If there is an image array
        if not just_label:
            self.image = np.float32(PDS3Image.open(IMG_path).image)
            self.proc_image = None
            self.mask_image = np.ones(
                self.image.shape[:2])*255  # TODO: remove this!!
            self.scale = self.label['DERIVED_IMAGE_PARMS']['RADIANCE_SCALING_FACTOR'][0]
            self.image *= self.scale

            self.pad_im = False

        self.find_offset_mode = None
        self.frame = frame

        # Parse the rover motion counters and the LMST
        try:
            self.site = int(self.label['ROVER_MOTION_COUNTER'][0])
            self.drive = int(self.label['ROVER_MOTION_COUNTER'][1])
            self.LMST = self.label['LOCAL_MEAN_SOLAR_TIME'].split('M')[1]
        except:
            self.site = -1
            self.drive = -1
            self.LMST = -1

        if self.cam_type == Camera.ZCAM_LEFT or self.cam_type == Camera.ZCAM_RIGHT:
            art_dev_pos = self.label['MINI_HEADER']['ARTICULATION_DEV_POSITION']
            self.focus_mc = art_dev_pos[0]
            self.zoom_mc = art_dev_pos[1]
            self.filter_mc = art_dev_pos[2]

    def _parse_cam_type(self, fname: str) -> Camera:
        """
        Parses the camera type from the fname of the camera
        :param fname: str, filename of the image
        """

        if fname[0:2] == 'ZL':
            return Camera.ZCAM_LEFT

        if fname[0:2] == 'ZR':
            return Camera.ZCAM_RIGHT

        if fname[0] == 'S':
            return Camera.SHERLOC

        if fname[0] == 'L':
            return Camera.RMI

        if fname[0] == 'N':
            return Camera.NAVCAM

        if 'MV' in fname or 'M_' in fname:
            return Camera.NAVCAM_VCE

        if fname[0] == 'F':
            return Camera.HAZCAM_FRONT

        if fname[0] == 'R':
            return Camera.HAZCAM_REAR

        if fname[0:3] == 'HSF':
            return Camera.HSF

        if fname[0:3] == 'HNM':
            return Camera.HNM

        return Camera.NONE
