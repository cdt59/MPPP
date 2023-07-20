# MPPP functions

'''

'''

from planetaryimage import PDS3Image
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interp1d
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import colour_demosaicing
import urllib.request, json 
import os
import cv2
import time
import glob

from numpy.linalg import inv, norm, det


class image:
    
    '''
    The image class holds all parameters specific to the IMG file
    '''

    def __init__(self, IMG_path, just_label=False, frame='site3' ):
        
        self.IMG_path    = IMG_path
        self.filename    = os.path.basename( IMG_path )
        self.label       = PDS3Image.open( IMG_path ).label                  # PDS image header metadata
        
        self.cam         = self.filename[:2]
        self.sol         = int( self.filename[4:8] )
        
        if not just_label:
            self.image       = np.float32( PDS3Image.open( IMG_path ).image )    # three band float-32 image array
            self.mask_image  = np.ones( self.image.shape[:2] )*255               # one band boolian image array
            # int to float scaling factor
            self.scale       = self.label['DERIVED_IMAGE_PARMS']['RADIANCE_SCALING_FACTOR'][0]
            self.image      *= self.scale
        
        self.find_offsets_mode = None
        self.frame = frame
        
        try:
            self.site  = int( self.label['ROVER_MOTION_COUNTER'][0] )
            self.drive = int( self.label['ROVER_MOTION_COUNTER'][1] )
            self.LMST  = self.label['LOCAL_MEAN_SOLAR_TIME'].split('M')[1]
        except:
            self.site  = 0
            self.drive = 0
            self.LMST  = 0
            
            
        if self.filename[0] == 'Z':
            self.focus_mc = self.label['MINI_HEADER']['ARTICULATION_DEV_POSITION'][0]
            self.zoom_mc  = self.label['MINI_HEADER']['ARTICULATION_DEV_POSITION'][1]
            self.filt_mc  = self.label['MINI_HEADER']['ARTICULATION_DEV_POSITION'][2]
            
            # if self.focus_mc > 2**10: self.focus_mc -= 2**11
            

    def image_process( self ):


        if self.filename.split('_N')[0][-3:]=='RZS': 
            self.ftau   = np.float32( self.label['DERIVED_IMAGE_PARMS']['RAD_ZENITH_SCALING_FACTOR'] )
            self.image *= self.ftau

        # if the image has one color band, either demosaic or stack the image to make it a color image
        if len(self.image.shape)==2:
            if 'MV0' in self.IMG_path:
                self.image = np.stack( [self.image,self.image,self.image], axis=-1)
            else:
                # self.image = colour_demosaicing.demosaicing_CFA_Bayer_bilinear  ( self.image, 'RGGB' )
                self.image = colour_demosaicing.demosaicing_CFA_Bayer_Malvar2004( self.image, 'RGGB' )
                # self.image = colour_demosaicing.demosaicing_CFA_Bayer_Menon2007 ( im_.image, 'RGGB' )


        d  = 57.296
        try:
            self.mu = np.sin( self.label['SITE_DERIVED_GEOMETRY_PARMS']['SOLAR_ELEVATION'][0]/d )
        except:
            self.mu = 1.0

        self.find_tau()

        self.tau_ref  = 0.3
        self.ftau     = self.mu * np.exp( - ( self.tau - self.tau_ref ) / 6 / self.mu )
        self.ftau_min = 0.2
        if self.ftau  < self.ftau_min: self.ftau = self.ftau_min


        self.im = self.image.copy()
        self.down_sample = self.filename.split('_')[-1][3]

        self.pad_left,self.pad_right,self.pad_top,self.pad_bottom = [0,0,0,0]


        '''
        future work: move these photometric adjustments to a separate function, photo_adjust(...)
        '''
        # Mars2020 Mastcam-Z color processing
        if self.filename[0] == 'Z':

            
            # pad Mastcam-Z images for th non-standard sizes

            self.pad_left,self.pad_right,self.pad_top,self.pad_bottom = [0,0,0,0]

            self.down_sample == '0'
            self.full_height, self.full_width = [ 1200, 1648 ]

            if self.pad_im:
                if self.im.shape != ( self.full_height, self.full_width, 3):                    

                    self.pad_left   =                    self.label['MINI_HEADER']['FIRST_LINE_SAMPLE'] - 1
                    self.pad_right  = self.full_width  - self.label['MINI_HEADER']['LINE_SAMPLES']      - self.label['MINI_HEADER']['FIRST_LINE_SAMPLE']  + 1
                    self.pad_top    =                    self.label['MINI_HEADER']['FIRST_LINE']        - 1
                    self.pad_bottom = self.full_height - self.label['MINI_HEADER']['LINES']             - self.label['MINI_HEADER']['FIRST_LINE']         + 1

                    if self.pad_top!=0 or self.pad_bottom!=0 or self.pad_left!=0 or self.pad_right!=0:
                        
                        print( 'resizing image size {} by padding = [ left, right, top, bottom ] = [ {}, {}, {}, {} ]'.format( \
                            self.image.shape, self.pad_left, self.pad_right, self.pad_top, self.pad_bottom ))
                        
                        self.im = pad_image( self.image, pad = [ self.pad_left, self.pad_right, self.pad_top, self.pad_bottom ] )

        
        # Mars2020 SuperCam RMI color processing
        # elif self.filename[0]=='L':

            # im /= flat

            # w = 400
            # high_scale = np.percentile( im[w:-w,w:-w,:], 99.8 )
            # im /= high_scale
            # clip_low = np.percentile( im[w:-w,w:-w,:], .05 )
            # clip_low = 0.3
            # high_cut = np.percentile( im[(w+300):-w,w:-w,:], 99.5 )
            # print( 'scale',high_scale, 'cut', clip_low)


        # Ingenuity Return-to-Earch (RTE) color processing
        elif self.filename[0:3] == 'HSF':
            self.ftau = 1.0

            # im /= np.load( 'C:/Users/cornell/Mastcam-Z/ws/HSF/HSF_flat_v1.npy' )
            # im /= np.percentile( im[400:-10,100:-100,:], 99.9 )*1.0
            # w = 100
            # clip_low = 0.2  #np.percentile( im[w:-w,w:-w,:], .5 )


        # Ingenuity Navcam color processing
        elif self.filename[0:3] == 'HNM':
            self.ftau = 1.0
#             w = 50
#             self.im /= np.percentile( im[w:-w,w:-w,:], 99.95 )*1.0
#             clip_low = np.percentile( im[w:-w,w:-w,:], 0.01 )*1.0


        # Mars2020 SHERLOC WATSON color processing
        elif self.filename[0]=='S':
            self.ftau = 1.0
       

        if self.filename[0] in [ 'F', 'N', 'R']:
            # Monochromatic VCE Navcam images
            # if 'MV0' in self.IMG_path:
            #     self.clip_low = 0.25

            # Pad to the image's standard dimensions [ full_height, full_width, 3 ]
            if self.pad_im:
                if   ( self.down_sample == '0' and self.im.shape!=(3840, 5120, 3) ) or \
                     ( self.down_sample == '1' and self.im.shape!=(1920, 2560, 3) ) or \
                     ( self.down_sample == '2' and self.im.shape!=( 960, 1280, 3) ):
                    if self.down_sample == '0': self.full_height, self.full_width = [ 3840, 5120 ]
                    if self.down_sample == '1': self.full_height, self.full_width = [ 1920, 2560 ]
                    if self.down_sample == '2': self.full_height, self.full_width = [  960, 1280 ]

                    self.pad_left   =                    np.min(self.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE_SAMPLE']) - 1
                    self.pad_right  = self.full_width  - np.max(self.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE_SAMPLE']) - 1280 + 1
                    self.pad_top    =                    np.min(self.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE'])        - 1
                    self.pad_bottom = self.full_height - np.max(self.label['INSTRUMENT_STATE_PARMS']['TILE_FIRST_LINE'])        - 960  + 1
                    
                    if self.pad_right < 0: self.pad_right = 0
                    if self.pad_bottom < 0: self.pad_bottom = 0
                    
                    
                    
                    if self.pad_top!=0 or self.pad_bottom!=0 or self.pad_left!=0 or self.pad_right!=0:
                        
                        print( 'resizing image size {} by padding = [ left, right, top, bottom ] = [ {}, {}, {}, {} ]'.format( \
                            self.image.shape, self.pad_left, self.pad_right, self.pad_top, self.pad_bottom ))
                        
                        self.im = pad_image( self.im, pad = [self.pad_left,self.pad_right,self.pad_top,self.pad_bottom] )

        # make image mask
        self.mask_im = self.mask_image.copy() 

        if ( self.pad_top!=0 or self.pad_bottom!=0 or self.pad_left!=0 or self.pad_right!=0 ) and self.pad_im:

            self.mask_im = pad_image( self.mask_im, pad = [ self.pad_left, self.pad_right, self.pad_top, self.pad_bottom ] )

            if  self.pad_bottom==0 and self.pad_right==0:
                self.mask_im[ self.pad_top:,                 self.pad_left:                ][ self.image[:,:,1] == 0 ] = 0
            elif self.pad_bottom==0:
                self.mask_im[ self.pad_top:,                 self.pad_left:-self.pad_right ][ self.image[:,:,1] == 0 ] = 0
            elif self.pad_right ==0:
                self.mask_im[ self.pad_top:-self.pad_bottom, self.pad_left:                ][ self.image[:,:,1] == 0 ] = 0
            else:
                self.mask_im[ self.pad_top:-self.pad_bottom, self.pad_left:-self.pad_right ][ self.image[:,:,1] == 0 ] = 0
        else:
            self.mask_im[ self.image[:,:,1] ==0 ] = 0

        # Mars2020 Mastcam-Z mask processing    
        if self.filename[0] in [ 'Z', 'S']:
            
            self.mask_im[ :4,  :] = 0
            self.mask_im[ -1:, :] = 0
            self.mask_im[ : ,:24] = 0
            self.mask_im[ :,-17:] = 0
            
            if 'IOF_N' in self.filename:
                self.ftau = 1.0
                
            # use pre-saved mask
            parent_path  = os.path.split( os.getcwd() )[0]
            if self.filename[:2] == 'ZL':
                mask_path = os.path.join( parent_path, 'params/ZL.jpg' )
            if self.filename[:2] == 'ZR':
                mask_path = os.path.join( parent_path, 'params/ZL.jpg' )
            else:
                mask_path = os.path.join( parent_path, 'params/S.jpg' )
            mask = cv2.imread( mask_path )
            self.mask_im[ mask[:,:,0] < 100 ] = 0
               

        # Mars2020 SuperCam RMI mask processing        
        if self.filename[0] == 'L':

            self.mask_im[ self.image==0 ] = 0
            self.mask_im[1800:,:,:] = 0
            self.mask_im = cv2.blur( self.mask_im, (20,20))
            self.mask_im[ self.mask_im<255]=0 
            
            
        if self.filename[:3] == 'HNM':

            parent_path  = os.path.split( os.getcwd() )[0]
            mask_path = os.path.join( parent_path, 'params/HNM.jpg' )
            mask = cv2.imread( mask_path )
            self.mask_im[ mask[:,:,0] < 100 ] = 0
            

        # Mars2020 Ecam mask processing    
        else:             
            self.mask_im[  :1, :] = 0
            self.mask_im[ -1:, :] = 0
            self.mask_im[ :, :2 ] = 0
            self.mask_im[ :,-2: ] = 0
            
            
            # use pre-saved mask
            if self.filename[0] == 'F':
            
                parent_path  = os.path.split( os.getcwd() )[0]
                if self.filename[:2] == 'FL':
                    mask_path = os.path.join( parent_path, 'params/FL{}.jpg'.format(self.down_sample) )
                else:
                    mask_path = os.path.join( parent_path, 'params/FR{}.jpg'.format(self.down_sample) )
                
                mask = cv2.imread( mask_path )
                self.mask_im[ mask[:,:,0] < 100 ] = 0
                
            if 'MV0' in self.filename or 'M_0' in self.filename:
            
                parent_path  = os.path.split( os.getcwd() )[0]
                if self.filename[:2] == 'NL':
                    mask_path = os.path.join( parent_path, 'params/NL2_vce.jpg' )
                else:
                    mask_path = os.path.join( parent_path, 'params/NR2_vce.jpg' )
                
                mask = cv2.imread( mask_path )
                self.mask_im[ mask[:,:,0] < 100 ] = 0       
            
            
        if self.filename[:3] == 'HNM' and 1:
           
            im_max = np.max( self.im )
            for i in range(3):
                self.im[:,:,i] = cv2.equalizeHist( np.uint8( self.im[:,:,i] * 255 / im_max )  ).astype('float') / 255 * im_max / self.scale_red
            self.clip_low = 0
            
            
        # apply color and brightnesss corrections
        self.im[:,:,0] *= self.scale / self.ftau * self.scale_red
        self.im[:,:,1] *= self.scale / self.ftau * 1
        self.im[:,:,2] *= self.scale / self.ftau * self.scale_blue
                
        # apply clipping
        self.im = ( self.im - self.clip_low )/( 1 - self.clip_low )
        self.im = np.clip( self.im, 0, 1 )
        
        # apply gamma corection
        if self.gamma != 1.0: 
            self.im = self.im**( 1/self.gamma ) 

        # rescale image to 8 unsigned bits
        self.im8 = np.clip( 255*self.im, 0, 255 ).astype('uint8')
        
                
        
    def find_tau( self ):
        
        '''
        return the tau, the opacity of Mars sky at the time of the image
        
        future work:
        * access the M2020 tau record
        * estimate the tue for the specific sol
        '''
        
        self.tau = 0.8
        
#         if self.sol >=700:
#             self.tau = 0.5


    # Cmod
    
    def image_reference( self ):
        
        
        GEOMETRIC_CAMERA_MODEL = self.label['GEOMETRIC_CAMERA_MODEL']
        self.cmod_from_cahvor( GEOMETRIC_CAMERA_MODEL )       
        
        if self.filename[0] == 'H':
            self.find_Rt_veh2site_inginuity( )
        else:
            self.find_Rt_veh2site_perseverance( )
        
        self.R_cam2ned  = R.from_matrix( [[0,-1,0],[1,0,0],[0,0,1]] ) 
        self.R_site2cam = self.R_cam2ned * self.R_veh2cam * self.R_veh2site.inv()
        self.ypr = find_ypr_from_R( self.R_site2cam  )    
        self.opk = find_opk_from_R( self.R_site2cam  )    
        
        self.yaw, self.pitch, self.roll = self.ypr
        self.omega, self.phi, self.kappa = self.opk      
        self.az, self.el = find_azel_from_ypr( self.ypr  )

        
    
        
        
    def cmod_from_cahvor( self, GEOMETRIC_CAMERA_MODEL ):
        
        self.projection = 'frame'
        self.C  = np.array(  GEOMETRIC_CAMERA_MODEL['MODEL_COMPONENT_1'], dtype=np.float64 )
        self.A  = np.array(  GEOMETRIC_CAMERA_MODEL['MODEL_COMPONENT_2'], dtype=np.float64 )
        self.H  = np.array(  GEOMETRIC_CAMERA_MODEL['MODEL_COMPONENT_3'], dtype=np.float64 )
        self.V  = np.array(  GEOMETRIC_CAMERA_MODEL['MODEL_COMPONENT_4'], dtype=np.float64 )
        self.O  = np.array(  GEOMETRIC_CAMERA_MODEL['MODEL_COMPONENT_5'], dtype=np.float64 )
        self.R  = np.array(  GEOMETRIC_CAMERA_MODEL['MODEL_COMPONENT_6'], dtype=np.float64 )

        self.hs = norm( np.cross( self.H, self.A ) )
        self.vs = norm( np.cross( self.V, self.A ) )
        self.hc = np.dot( self.H, self.A ) 
        self.vc = np.dot( self.V, self.A ) 

        self.hp = ( self.H - self.hc* self.A ) / self.hs
        self.vp = ( self.V - self.vc* self.A ) / self.vs

        self.theta = np.arcsin( ( - norm( np.cross( self.vp, self.hp ) )
                                  / norm( self.vp )
                                  / norm( self.hp ) ) )
        self.theta_degrees = np.rad2deg( self.theta )
        
        self.K_cam = np.array([
                    [ -self.hs*np.sin(self.theta), self.hs*np.cos(self.theta), self.hc ],
                    [                           0,                    self.vs, self.vc ],
                    [                           0,                          0,       1 ], ])

        self.rot_cam = np.matmul( inv( self.K_cam ), 
                                  np.array( [ self.H, self.V, self.A ] ) )
        
        self.R_cam = R.from_matrix( self.rot_cam )       
        self.R_veh2cam = self.R_cam
        
        self.R_ned2enu = R.from_matrix( [[0,1,0],[1,0,0],[0,0,-1]] )
        
        self.w, self.h = [ self.label['IMAGE']['LINE_SAMPLES'], self.label['IMAGE']['LINES'] ]
        if self.w==1600: self.w = 1648

        self.k1 = self.R[1]
        self.k2 = self.R[2]
        self.k3 = 0
        self.k4 = 0

        self.p1 = 0
        self.p2 = 0        
        
        self.cxp = self.hc - self.w/2
        self.cyp = self.vc - self.h/2
        
        self.f  =  self.vs
        self.b1 = -self.hs * np.sin( self.theta ) - self.vs
        self.b2 =  self.hs * np.cos( self.theta )
        
        
        
    
    def find_Rt_veh2site_inginuity( self ):
        
        self.pan      = 0
        self.tilt     = 0
        self.az       = 0
        self.az_veh   = 0
        self.C = np.array( self.label['GEOMETRIC_CAMERA_MODEL']['MODEL_COMPONENT_1'] ) - np.array( self.label['GEOMETRIC_CAMERA_MODEL']['MODEL_TRANSFORM_VECTOR'] )

        self.q_HELI_M = q_wxyz2xyzw( self.label['HELI_M_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'] )
        self.R_HELI_M = R.from_quat( self.q_HELI_M )
        self.C_HELI_M = self.R_HELI_M.apply( self.C, inverse=0 ) \
                      + np.array( self.label['HELI_M_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'] )
        
                        
        self.q_HELI_G = q_wxyz2xyzw( self.label['HELI_G_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'] )
        self.R_HELI_G = R.from_quat( self.q_HELI_G )
        self.C_HELI_G = self.R_HELI_G.apply( self.C_HELI_M, inverse=0 ) \
                      + np.array( self.label['HELI_G_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'] )
        
        R_HELI_takeoff = self.R_HELI_G.inv() * self.R_HELI_M
        C_HELI_takeoff = R_HELI_takeoff.apply( self.label['HELI_G_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'] ) 
        
#         print( "C_HELI_M", np.round( self.C_HELI_M,2), "C_HELI_takeoff", np.array( self.label['HELI_G_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'] ) )
        
        
                 
        self.xyz      = self.C_HELI_G.copy()
        self.xyz_veh  = self.C_HELI_G.copy()  # approximation


        self.X_shift, self.Y_shift, self.Z_shift = [0,0,0]
        
#         if self.sol == 174:
#             sclk       = int( self.filename[9:19])
#             sclk_start = 682390500
#             sclk_end   = 682390670
#             sclk_inter = ( sclk - sclk_start ) / ( sclk_end - sclk_start )
#             self.X_shift, self.Y_shift, self.Z_shift = sclk_inter * np.array( [5.449, -11.234, -0.408] )
            
            
        self.X, self.Y, self.Z = xyz_ned2enu( self.xyz )
        self.X_offset, self.Y_offset, self.Z_offset = xyz_ned2enu( self.xyz_veh )
        
        self.X        += self.X_shift
        self.Y        += self.Y_shift
        self.Z        += self.Z_shift
        self.X_offset += self.X_shift
        self.Y_offset += self.Y_shift
        self.Z_offset += self.Z_shift
        
        self.R_veh2site =  self.R_HELI_G * self.R_HELI_M 
        self.t_veh2site = -self.xyz_veh
        
        
    def find_Rt_veh2site_perseverance( self ):
        
        try:
            self.pan  = np.rad2deg( self.label['RSM_ARTICULATION_STATE']['ARTICULATION_DEVICE_ANGLE'][0] )
            self.tilt = np.rad2deg( self.label['RSM_ARTICULATION_STATE']['ARTICULATION_DEVICE_ANGLE'][1] )
        except:
            self.pan  = np.rad2deg( self.label['RSM_ARTICULATION_STATE']['ARTICULATION_DEVICE_ANGLE'][0][0] )
            self.tilt = np.rad2deg( self.label['RSM_ARTICULATION_STATE']['ARTICULATION_DEVICE_ANGLE'][1][0] )
        
        if self.frame == 'site3':
            
            self.az     =   self.label['SITE_DERIVED_GEOMETRY_PARMS' ]['INSTRUMENT_AZIMUTH'][0]
            self.az_veh = ( self.label['ROVER_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'][0] - 
                            self.label['SITE_DERIVED_GEOMETRY_PARMS' ]['INSTRUMENT_AZIMUTH'][0])%360

            self.q  = q_wxyz2xyzw( self.label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'] )
 
            self.R_veh = R.from_quat( self.q )
            self.Cr =  self.R_veh.apply( self.C, inverse=0 )

            self.xyz_veh = np.array( self.label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'] )        
            self.xyz = self.Cr + np.array( self.label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'] )
            
            
        if self.frame == 'rnav':
        
            self.az     = self.label['ROVER_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'][0]
            self.az_veh = 0

            self.q  = [0,0,0,1]

            self.R_veh = R.from_quat( self.q )
            self.Cr =  self.R_veh.apply( self.C, inverse=0 )

            self.xyz_veh = np.array( [0,0,0] )  
            self.xyz = self.Cr + self.xyz_veh
            
            
        self.X, self.Y, self.Z = xyz_ned2enu( self.xyz )
        self.X_offset, self.Y_offset, self.Z_offset = xyz_ned2enu( self.xyz_veh )        

        if not self.find_offsets_mode:
            
            if self.frame == 'site3':
                self.X_shift, self.Y_shift, self.Z_shift = XYZ_shift_offsets( self.site, self.drive )
            else:
                self.X_shift, self.Y_shift, self.Z_shift = [ 0, 0, 0 ]
            
            self.X        += self.X_shift
            self.Y        += self.Y_shift
            self.Z        += self.Z_shift
            self.X_offset += self.X_shift
            self.Y_offset += self.Y_shift
            self.Z_offset += self.Z_shift
        
        self.R_veh2site = self.R_veh
        self.t_veh2site = self.xyz
        
                
    



def pad_image( im, pad = [0,0,0,0] ):
    
    if len( im.shape ) == 3:
        im     = np.hstack( [ np.zeros( (im.shape[0],   pad[0] , 3)), im, np.zeros( ( im.shape[0],  pad[1], 3)), ] )
        im     = np.vstack( [ np.zeros( (pad[2]  , im.shape[1],  3)), im, np.zeros( ( pad[3],  im.shape[1], 3)), ] )
    else:
        im     = np.hstack( [ np.zeros( (im.shape[0],   pad[0]    )), im, np.zeros( ( im.shape[0],  pad[1]   )), ] )
        im     = np.vstack( [ np.zeros( (pad[2]  , im.shape[1],   )), im, np.zeros( ( pad[3],  im.shape[1]   )), ] )
    return im

                    

def make_save_path( IMG_path, directory_output, fullpath = True, file_extension = '.png' ):
    
    '''
    make_save_path sorts the images into an output directory organized by camera type and each 100 sols of the mission
    '''
    
    filename = os.path.basename( IMG_path )
    sol      = int( filename[4:8] )
    camera   = filename[0]
    mission  = 'Mars2020' # mission name is hardcoded for now
    
    if camera in ['F','N','R']:
        camera_type = 'eng'
    elif camera in ['H']:
        camera_type = 'heli'
    elif camera in ['Z','L','S']:
        camera_type = 'sci'

    sol_floor_100 = int(np.floor( sol/100 ) * 100)
    sol_range_100 = str(sol_floor_100).zfill(4) + '-' + str(sol_floor_100).zfill(4)[:2] + '99'

    save_path = directory_output + '/sols_' + sol_range_100 + '_' + camera_type 
    
    if not os.path.exists(save_path):
        # Create a new directory because it does not exist
        os.makedirs(save_path)
        print("The new directory is created: ", save_path )

    if fullpath:
        return save_path + '/' + filename.split('.')[0] + file_extension
    else:
        return save_path 



def plot_image_locations( IMG_paths, im_xyzs, rover_xyzs, rover_rots, im_azs, im_els ):
    
    '''
    plot_image_locations displays the Northing vs Easting locations of each image and rover position
    
    future work: replace the input arrays with a single pandas dataframe
    '''
    
    
    plt.figure( figsize=[12,8])
    
    scale = np.round( np.std( np.array(rover_xyzs), axis=0 ).max()/4+1 )

    for i in range(len(im_xyzs)):
        
        filename = os.path.basename( IMG_paths[i] )
        
        marker = '*k'
        if filename[:2] in ['FL','RL']: marker = 'ob'
        if filename[:2] in ['FR','RR']: marker = 'or'
        if filename[:2] ==  'NL':       marker = 'sb'
        if filename[:2] ==  'NR':       marker = 'sr'
        if filename[:2] ==  'ZL':       marker = '^b'
        if filename[:2] ==  'ZR':       marker = '^r'

            
        if 'MV' not in IMG_paths[i]:
            
            plt.plot( rover_xyzs[i][0], rover_xyzs[i][1], color='k',    marker=(4, 0, 45  + rover_rots[i]), ms=30, )
            plt.plot( rover_xyzs[i][0], rover_xyzs[i][1], color='gray', marker=(3, 0, 120 + rover_rots[i]), ms=20, )
        
              
            
            plt.text( rover_xyzs[i][0]+scale/4, rover_xyzs[i][1]+scale/4, 'Sol '+ os.path.basename( IMG_paths[i] )[4:8], \
                      bbox=dict(facecolor='w', alpha=0.5, edgecolor='w'), size='large' )

            if i>1 and os.path.basename( IMG_paths[i] )[:2]=='NLF_':
                plt.plot( [rover_xyzs[i][0], rover_xyzs[i][1] ], [rover_xyzs[i][0], rover_xyzs[i][1] ], '--', color='gray' )
            
            cos_az = np.cos(im_azs[i]/57.3)
            sin_az = np.sin(im_azs[i]/57.3)
            cos_el = np.cos(im_els[i]/57.3)


            if os.path.basename( IMG_paths[i] )[0]=='Z':
                plt.arrow( im_xyzs[i][0], im_xyzs[i][1], scale*cos_el*sin_az, scale*cos_el*cos_az,
                           color=marker[1], lw = int(scale/32), linestyle='dashed' )
            else:
                plt.arrow( im_xyzs[i][0], im_xyzs[i][1], scale*cos_el*sin_az, scale*cos_el*cos_az,
                           color=marker[1], lw = int(scale/32) )
                
        plt.plot( im_xyzs[i][0], im_xyzs[i][1], marker )
        

    plt.axis('equal')
    plt.xlim( [ np.round(plt.gca().get_xlim()[0])-3, np.round(plt.gca().get_xlim()[1])+3 ] )
    plt.ylim( [ np.round(plt.gca().get_ylim()[0])-3, np.round(plt.gca().get_ylim()[1])+3 ] )

    plt.xlabel( 'Easting Site Frame')
    plt.ylabel( 'Northing Site Frame')
#     plt.tight_layout()
#     plt.savefig( path + '/positions'+suf+'.jpg', dpi=300  )




def XYZ_shift_offsets( site, drive ):
    
    '''
    XYZ_shift_offsets finds most accurate Site-Nav offset for each site index and drive
    
    '''

    # print( site, drive )

    parent_path  = os.path.split( os.getcwd() )[0]
    waypoint_shift_path = os.path.join( parent_path, 'params/Mars2020_waypoint_shifts.csv' )

    shift_params = np.loadtxt( waypoint_shift_path, delimiter=',', skiprows=1 )

    site_shifts  = shift_params[ np.where( shift_params[:,1]==site)[0] ]
    site_drives  = site_shifts[:,2]

    if drive in site_drives:
        drive_site_shift = site_shifts[ np.where( site_shifts[:,2]==drive)[0] ][0,:]

    elif drive > site_drives.min() and drive < site_drives.max():
        drive_site_shift = interp1d( site_shifts[:,2], site_shifts, axis=0)(drive)

    elif drive >= site_drives.max():
        drive_site_shift = site_shifts[-1,:]

    else:
        drive_site_shift = np.zeros(12)

    # print( drive_site_shift )
    X_shift, Y_shift, Z_shift = drive_site_shift[9:]

    # X_shift, Y_shift, Z_shift = [ 0,0,0 ]
 
    return X_shift, Y_shift, Z_shift


def remove_duplicate_IMGs( IMG_paths ):
    
    names = [ IMG_paths[i][:-5] for i in range(len(IMG_paths)) ]
    duplicates = list( set(  [ names[i] for i, x in enumerate(names) if i != names.index(x)] ))
    print( len( duplicates ))
    for i in range( len( duplicates) ):
        all_i_paths = sorted( glob.glob( duplicates[i] +'*.IMG'))[::-1]
        duplicates_i_paths = all_i_paths[1:]
        print( '\nkeeping  ', os.path.basename( all_i_paths[0] ) )
        for j in range(len( duplicates_i_paths )):
            print( 'removing ', os.path.basename( duplicates_i_paths[j] ) )
            os.remove( duplicates_i_paths[j] )
            
            
def image_list_process( IMG_paths, directory_output, suf, find_offsets_mode = 0, frame='site3' ):
    
    
    # File parameters    

    '''
    future work: save these calibration prameters as a text files, which we load for each camera
    '''
    
    file_extension = ''
    
    # save images and thereby overwrite existing images
    save_im    = 1

    # add an alpha channel to the output images
    save_mask  = 1

    # add transparrent pixels to restore the image's full, standard size
    pad_im     = 1
    pad_im_z   = 1

    # turn on when finding the waypoint offsets
#     find_offsets_mode = 1

    # set the color values
    gamma      = 2.2      # gamma value
    gamma      = 2        # gamma value

    # fraction of the dynamic range to clip off the lower values of image 
    clip_low_z = 0.02  # for the Mastcam-Z cameras
    clip_low   = 0.05  # for everything else


    # scale all the scale parameters below bsy the same number
    scale_scale = 18

    # color balance parameters for the Mars 2020 science cameras
    scale_z,  scale_red_z,  scale_blue_z  = [ 1.1*scale_scale, 0.7 , 1.5  ] # Mastcam-Z 
    scale_l,  scale_red_l,  scale_blue_l  = [ 1.0*scale_scale, 0.75, 1.40 ] # SuperCam RMI
    scale_s,  scale_red_s,  scale_blue_s  = [ 1.0*scale_scale, 0.85, 1.40 ] # SHERLOC WATSON 

    # color balance parameters for the Mars 2020 engineering cameras
    scale_n,  scale_red_n,  scale_blue_n  = [ 1.0*scale_scale, 0.75, 1.2  ] # Navcam
    scale_v,  scale_red_v,  scale_blue_v  = [ 1.3*scale_scale, 1.10, 0.93 ] # Grayscale VCE Navcam
    scale_f,  scale_red_f,  scale_blue_f  = [ 1.1*scale_scale, 0.78, 1.25 ] # Front Hazcam
    scale_r,  scale_red_r,  scale_blue_r  = [ 1.1*scale_scale, 0.78, 1.25 ] # Rear Hazcam
    scale_hr, scale_red_hr, scale_blue_hr = [ 1.0*scale_scale, 0.75, 1.43 ] # Inginuity RTE
    scale_hn, scale_red_hn, scale_blue_hn = [ 1.0*scale_scale, 1.08 , 0.92 ] # Inginuity Navcam
    
    pos_lines    = []
    error_lines  = []
    veh_XYZs     = []
    im_XYZs      = []
    veh_azs     = []
    im_azs       = []
    im_els       = []
    sols         = []
    rmcs         = []
    ims          = []
    im_save_path = ''

    print( len(IMG_paths), 'images\n')


    for i in range(len(IMG_paths))[::][:]:

        
        ####################################################
        ################# *** debugging *** ################
        ####################################################
        try:    # catch all the images that fail to process

            # open image
            im = image( IMG_paths[i], frame=frame )
            print( i, im.filename )
            IMG_loaded = True
            
        except:
            print( os.path.basename( IMG_paths[i]), 'failed to process! \n' )
            error_lines.append( os.path.basename( IMG_paths[i]) +'\n' )
            IMG_loaded = False
            
            
            
        if IMG_loaded:

            # Set color processing parameters
            im.scale       = scale_scale
            im.scale_red   = 1
            im.scale_blue  = 1
            im.clip_low    = clip_low
            im.gamma       = gamma
            im.pad_im      = pad_im
            im.save_im     = save_im
            im.save_mask   = save_mask
            im.find_offsets_mode = find_offsets_mode

            # Mars 2020 Mastcam-Z
            if im.cam[0] == 'Z':
                im.scale       = scale_z
                im.scale_red   = scale_red_z
                im.scale_blue  = scale_blue_z
                im.clip_low    = clip_low_z
                im.pad_im      = pad_im_z

    #             if 'IOF_N' in im.IMG_path:
    #                 im.scale       = scale_n*1.4
    #                 im.scale_red   = 0.65
    #                 im.scale_blue  = 1.3

            # Mars 2020 SHERLOC WATSON
            if im.cam[0] == 'S':
                im.scale       = scale_s
                im.scale_red   = scale_red_s
                im.scale_blue  = scale_blue_s
                im.clip_low    = 0.0

            # Mars 2020 SuperCam RMI
            if im.cam[0] == 'L':
                im.scale       = scale_l
                im.scale_red   = scale_red_l
                im.scale_blue  = scale_blue_l

            # Mars 2020 Navcam
            if im.cam[0] == 'N':
                im.scale       = scale_n
                im.scale_red   = scale_red_n
                im.scale_blue  = scale_blue_n

                
                
            # Mars 2020 Navcam VCE images
            if 'MV0' in im.filename or 'M_0' in im.filename:
                im.scale       = scale_v
                im.scale_red   = scale_red_v
                im.scale_blue  = scale_blue_v
                im.clip_low    = 0.1


            # Mars 2020 Front Hazcam
            if im.cam[0] == 'F':
                im.scale       = scale_f
                im.scale_red   = scale_red_f
                im.scale_blue  = scale_blue_f
                im.clip_low    = clip_low/2

            # Mars 2020 Rear Hazcam
            if im.cam[0] == 'R':
                im.scale       = scale_r
                im.scale_red   = scale_red_r
                im.scale_blue  = scale_blue_r
                im.clip_low    = clip_low/2

            # Heli Ingenuity RTE 
            if im.filename[0:3] == 'HSF':
                im.scale       = scale_hr
                im.scale_red   = scale_red_hr
                im.scale_blue  = scale_blue_hr

            # Heli Ingenuity Navcam  
            if im.filename[0:3] == 'HNM':
                im.scale       = scale_hn
                im.scale_red   = scale_red_hn
                im.scale_blue  = scale_blue_hn
                im.clip_low    = 0.3
                
                
            file_extension = '.png'
            
            im.focus_mc = -1
            im.zoom_mc  = -1
            
#             if im.cam[0] == 'Z':
                
#                 im.focus_mc = im.label['MINI_HEADER']['ARTICULATION_DEV_POSITION'][0]
#                 im.zoom_mc  = im.label['MINI_HEADER']['ARTICULATION_DEV_POSITION'][1]
                
#                 D_near, D_med = [ 3.0, 6.0 ] 
                
#                 if im.filename[45:48] == '034':
#                     if im.cam[0:2] == 'ZL':
#                         B,A = [ 1216.2, 2194 ]
#                     if im.cam[0:2] == 'ZR':
#                         B,A = [ 1276.6, 2264 ]
                        
                    
#                     D_mc = A / ( B - im.focus_mc )
#                     if D_mc < 0 : D_mc = 100
                    
#                     if D_mc < D_near:
#                         file_extension = '_near.png'
#                         print( im.focus_mc,D_mc , file_extension )
                        
#                     if D_mc >= D_near and D_mc < D_med:
#                         file_extension = '_med.png'
#                         print( im.focus_mc, D_mc , file_extension )


                    

            # create save directory
            im.save_path_full = make_save_path( im.IMG_path, directory_output, fullpath=True, file_extension = file_extension  ) 
            im.save_path      = make_save_path( im.IMG_path, directory_output, fullpath=False ) 
            im.save_name      = im.save_path_full.split('/')[-1]
            csv_save_path     = im.save_path_full

            # process and save image
            if im.save_im:

                im.image_process( )

                if im.save_mask:
                    im.im8a = cv2.cvtColor( im.im8, cv2.COLOR_BGR2RGBA )
                    im.im8a[:,:,3] = im.mask_im
                    cv2.imwrite( im.save_path_full, im.im8a )                
                else:
                    cv2.imwrite( im_save_path_full, im.im8[:,:,::-1] )  


            # find image position and rotation parameters
            im.image_reference( )

            # save reference data for plotting        
            '''
            future work: replace these lists with pandas dataframes
            '''
            im_XYZs .append( [ im.X, im.Y, im.Z ] )
            veh_XYZs.append( [ im.X_offset, im.Y_offset, im.Z_offset ] )
            veh_azs .append( im.az_veh )
            im_azs  .append( im.az )
            im_els  .append( im.el )
            rmcs    .append( im.label['ROVER_MOTION_COUNTER'])
            sols    .append( int(im.label['LOCAL_TRUE_SOLAR_TIME_SOL']) )

            # create a line for the reference file
#             # Label	 X/East	Y/North	Z/Altitude	Yaw	Pitch	Roll
#             pos_line =  im.save_name+'\t'\
#                  +str( np.round( im.X,4))+'\t'\
#                  +str( np.round( im.Y,4))+'\t'\
#                  +str( np.round( im.Z,4))+'\t'\
#                  +str( np.round( im.yaw,  4))+'\t'\
#                  +str( np.round( im.pitch,4))+'\t'\
#                  +str( np.round( im.roll, 4))+'\n'

            # Label	 X/East	Y/North	Z/Altitude	Omega Phi Kappa  pan tilt
            pos_line =  im.save_name\
                     + '\t'\
                     + str( np.round( im.X    , 5))+'\t'\
                     + str( np.round( im.Y    , 5))+'\t'\
                     + str( np.round( im.Z    , 5))+'\t'\
                     + str( np.round( im.omega, 5))+'\t'\
                     + str( np.round( im.phi  , 5))+'\t'\
                     + str( np.round( im.kappa, 5))+'\t'\
                     + str( np.round( im.pan  , 5))+'\t'\
                     + str( np.round( im.tilt , 5))+'\t'\
                     + str( np.round( im.focus_mc, 5))+'\t'\
                       '\n'

            pos_lines.append( pos_line )

            if 1:
                try:
                    print( 'sol {} site {} drive {}  zenith angle {:0.0f} scale {:0.2f}'.
                                format( im.sol, im.site, im.drive, im.el*57.3, im.ftau ) )
                except:
                    print( 'sol {} site {} drive {}'.
                                format( im.sol, im.site, im.drive, ) )
                print( 'XYZ_ENU = [{:0.2f}, {:0.2f}, {:0.2f}] YPR = [{:0.0f}, {:0.0f}, {:0.0f}]  OPK = [{:0.0f}, {:0.0f}, {:0.0f}]'.format(im.X,im.Y,im.Z,im.yaw,im.pitch,im.roll,im.omega,im.phi,im.kappa) )
                print( )




    current_time = time.strftime("%Y%m%d-%H%M%S")


#     #save failed images list as TXT
#     if len(error_lines) > 0:
#         csv_save_path = os.path.dirname( csv_save_path)+'/failed_'+suf+'_'+current_time+'.txt'
#         with open(csv_save_path,'w') as file:
#             for error_line in error_lines:
#                 file.write(error_line)
#     print( 'saved', csv_save_path )

    #save image positions as CSV file
    csv_save_path = os.path.dirname( csv_save_path)+'/positions_'+suf+'_'+frame+'_'+current_time+ '.txt'
    with open(csv_save_path,'w') as file:
        for pos_line in pos_lines:
            file.write(pos_line)
    print( 'saved', csv_save_path )

    len( pos_lines )
    
    plot_image_locations( IMG_paths, im_XYZs, veh_XYZs, veh_azs, im_azs, im_els )
    
    if find_offsets_mode:
        sites  = [ rmcs[i][0] for i in range(len(rmcs))[::-1] ]
        drives = [ rmcs[i][1] for i in range(len(rmcs))[::-1] ]
        Xs     = [ veh_XYZs[i][0] for i in range(len(veh_XYZs))[::-1] ]
        Ys     = [ veh_XYZs[i][1] for i in range(len(veh_XYZs))[::-1] ]
        Zs     = [ veh_XYZs[i][2] for i in range(len(veh_XYZs))[::-1] ]

        table = np.stack( [sols[::-1], sites, drives, Xs, Ys, Zs], axis=1)
        np.round( table, 4 )

        np.savetxt( directory_output+"/offsets_"+suf+".csv", table, delimiter="\t")


        
def xyz_ned2enu( xyz ):
    return np.array( [ xyz[1], xyz[0], -xyz[2] ] )
        
def q_wxyz2xyzw( q_wxyz ):
     return np.array([ q_wxyz[1], q_wxyz[2], q_wxyz[3], q_wxyz[0]] )
        
# def find_ypr_from_R( R_site2cam ):    
#     R_cam2ned = R.from_matrix( [[0,-1,0],[1,0,0],[0,0,1]] )    
#     angles = ( R_site2cam * R_cam2ned ).as_euler( 'zxy',degrees=1 )
#     return np.array([ ( - angles[0] ) % 360, - angles[1], angles[2] ])

def find_ypr_from_R( R_ ):
    # finds yaw, pitch, roll from rotation matrix R_cam2site
    ypr = -(R_).as_euler('zyx',degrees=1)
    if np.abs( ypr[2] ) > 90:
        ypr = [ (ypr[0]+180)%360, ypr[1], (ypr[2]+180)%360 ]
    return np.array( ypr )

def find_R_from_ypr( ypr ):
    # finds matrix R_cam2site from yaw, pitch, roll angles    
    R_ = R.from_euler( 'zyx', -np.array(ypr) ,degrees=1) 
    return R_

def find_azel_from_ypr( ypr ):   
    return np.array([ ypr[0] % 360, ypr[1] - 90 ])


def find_opk_from_R( R_ ):
    # finds omega, phi, kappa from rotation matrix R_cam2site
    R_cam2ned  = R.from_matrix( [[0,-1,0],[1,0,0] ,[0,0,1] ] )
    angles = ( R_cam2ned.inv()*R_.inv() ).as_euler('XYZ',degrees=1) 
    opk    = [ angles[0], -angles[1], -angles[2]-90 ]
    return np.array( opk )

def find_R_from_opk( opk ):
    # finds matrix R_cam2site from omega, phi, kappa angles    
    R_cam2ned  = R.from_matrix( [[0,-1,0],[1,0,0] ,[0,0,1] ] )   
    R_angles = R.from_euler( 'XYZ', np.array(opk) ,degrees=1) 
    R_ = ( R_cam2ned * R_angles).inv()
    return R_



    
# def find_R_cam_from_ypr_veh( ypr, R_veh ):    
#     R_cam2ned = R.from_matrix( [[0,-1,0],[1,0,0],[0,0,1]] ) 
#     angles = [ -ypr[0], -ypr[1], ypr[2] ]
#     return R.from_euler( 'zxy', angles, degrees=1 ) * R_cam2ned.inv() * R_veh

# def find_R_cam_from_ypr( ypr ):    
#     R_cam2ned = R.from_matrix( [[0,-1,0],[1,0,0],[0,0,1]] ) 
#     angles = [ -ypr[0], -ypr[1], ypr[2] ]
#     return R.from_euler( 'zxy', angles, degrees=1 ) * R_cam2ned.inv() 