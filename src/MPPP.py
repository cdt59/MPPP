# MPPP funcitons

'''
future work: define an image class
'''

import numpy as np
from planetaryimage import PDS3Image
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import colour_demosaicing
from PIL import Image
import os


class image:
    
    '''
    The image class holds all parameters specific to a IMG file
    '''

    def __init__(self, IMG_path ):
        
        self.IMG_path    = IMG_path
        self.filename    = os.path.basename( IMG_path )
        self.label       = PDS3Image.open( IMG_path ).label                  # PDS image header metadata
        self.image       = np.float32( PDS3Image.open( IMG_path ).image )    # three band float-32 image array
        self.mask_image  = np.ones( self.image.shape[:2] )*255               # one band boolian image array
        self.cam         = self.filename[:2]
        
        # int to float scaling factor
        self.scale       = self.label['DERIVED_IMAGE_PARMS']['RADIANCE_SCALING_FACTOR'][0]
        self.image      *= self.scale


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
        self.mu = np.sin( self.label['SITE_DERIVED_GEOMETRY_PARMS']['SOLAR_ELEVATION'][0]/d )

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
            # w = 40
            # im /= np.percentile( im[w:-w,w:-w,:], 99.95 )*1.0
            # clip_low = np.percentile( im[w:-w,w:-w,:], 0.01 )*1.0


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

                    if self.pad_top!=0 or self.pad_bottom!=0 or self.pad_left!=0 or self.pad_right!=0:
                        self.im = pad_image( self.im, pad = [self.pad_left,self.pad_right,self.pad_top,self.pad_bottom] )


        # make image mask

        self.mask_im = self.mask_image.copy() 

        if ( self.pad_top!=0 or self.pad_bottom!=0 or self.pad_left!=0 or self.pad_right!=0 ) and pad_im:

            print( 'resizing image size {} by padding = [ left, right, top, bottom ] = [ {}, {}, {}, {} ]'.format( \
                    self.image.shape, self.pad_left, self.pad_right, self.pad_top, self.pad_bottom ))    

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
            self.mask_im[ : ,:25] = 0
            self.mask_im[ :,-18:] = 0

        # Mars2020 SuperCam RMI mask processing        
        if self.filename[0] == 'L':

            self.mask_im[ self.image==0 ] = 0
            self.mask_im[1800:,:,:] = 0
            self.mask_im = cv2.blur( self.mask_im, (20,20))
            self.mask_im[ self.mask_im<255]=0 

        # Mars2020 Ecam mask processing    
        else:             
            self.mask_im[  :1, :] = 0
            self.mask_im[ -1:, :] = 0
            self.mask_im[ :, :2 ] = 0
            self.mask_im[ :,-2: ] = 0

            
            
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
        
        # print( 'processed image', self.filename)



    def image_reference( self ):


            '''
            future work: move exterior camera calculations to separate function
            '''


            self.az    = self.label['SITE_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'  ][0]
            self.el    = self.label['SITE_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_ELEVATION'][0]
            self.xyz   = np.array( self.label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR'].copy() )
            self.C     = self.label['GEOMETRIC_CAMERA_MODEL']['MODEL_COMPONENT_1'].copy()
            self.el   += 90
            self.rl    = 0
            if self.filename[:2]=='FL': self.rl = + 10
            if self.filename[:2]=='FR': self.rl = - 10


            try: 
                self.rot       = 57.3*np.float32( self.label['RSM_ARTICULATION_STATE']['ARTICULATION_DEVICE_ANGLE'][0] )
                self.rot_rover = (self.label['ROVER_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'][0] - self.label['SITE_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'][0])%360
            except: 
                self.rot       = 57.3*np.float32( self.label['RSM_ARTICULATION_STATE']['ARTICULATION_DEVICE_ANGLE'][0][0] )
                self.rot_rover = (self.label['ROVER_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'][0] - self.label['SITE_DERIVED_GEOMETRY_PARMS']['INSTRUMENT_AZIMUTH'][0])%360

            self.q = self.label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'].copy()
            self.q = [ self.q[1], self.q[2], self.q[3], self.q[0]]
            self.Cr = R.from_quat( self.q ).apply( self.C, inverse=0 )

            self.xyz_rover = self.xyz.copy()


            self.xyz[0] += self.Cr[0]
            self.xyz[1] += self.Cr[1]
            self.xyz[2] += self.Cr[2]

            self.X =  self.xyz[1]
            self.Y =  self.xyz[0]
            self.Z = -self.xyz[2]

            self.X_offset =  self.xyz_rover[1]
            self.Y_offset =  self.xyz_rover[0]
            self.Z_offset = -self.xyz_rover[2]

            self.sol     = int( self.label['LOCAL_TRUE_SOLAR_TIME_SOL'] )
            self.site    = int( self.label['ROVER_MOTION_COUNTER'][0] )
            self.drive   = int( self.label['ROVER_MOTION_COUNTER'][1] )

            self.LMST = self.label['LOCAL_MEAN_SOLAR_TIME'].split('M')[1]

            self.x_shift, self.y_shift, self.z_shift = xyz_shift_offsets( self.site, self.drive )

            if 1 and self.find_offsets_mode==0:
                self.X += self.x_shift
                self.Y += self.y_shift
                self.Z += self.z_shift
                self.X_offset += self.x_shift
                self.Y_offset += self.y_shift
                self.Z_offset += self.z_shift
        
    def find_tau( self ):
        self.tau = 0.6

    

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



def plot_image_locations( RAD_paths, im_xyzs, rover_xyzs, rover_rots, im_azs, im_els ):
    
    '''
    plot_image_locations displays the Northing vs Easting locations of each image and rover position
    
    future work: replace the input arrays with a single pandas dataframe
    '''
    
    
    plt.figure( figsize=[12,8])
    
    scale = np.round( np.std( np.array(rover_xyzs), axis=0 ).max()/4+1 )

    for i in range(len(im_xyzs)):
        
        filename = os.path.basename( RAD_paths[i] )
        
        marker = '*k'
        if filename[:2] in ['FL','RL']: marker = 'ob'
        if filename[:2] in ['FR','RR']: marker = 'or'
        if filename[:2] ==  'NL':       marker = 'sb'
        if filename[:2] ==  'NR':       marker = 'sr'
        if filename[:2] ==  'ZL':       marker = '^b'
        if filename[:2] ==  'ZR':       marker = '^r'

        plt.plot( rover_xyzs[i][0], rover_xyzs[i][1], color='k',    marker=(4, 0, 45  + rover_rots[i]), ms=30, )
        plt.plot( rover_xyzs[i][0], rover_xyzs[i][1], color='gray', marker=(3, 0, 120 + rover_rots[i]), ms=20, )
        plt.text( rover_xyzs[i][0]+scale/5, rover_xyzs[i][1]+scale/5, 'Sol '+ RAD_paths[i].split('\\')[-1][4:8], \
                  bbox=dict(facecolor='w', alpha=0.5, edgecolor='w'), size='large' )

        if i>1 and RAD_paths[i].split('\\')[-1][:2]=='NLF_':
            plt.plot( [rover_xyzs[i][0], rover_xyzs[i][1] ], [rover_xyzs[i][0], rover_xyzs[i][1] ], '--', color='gray' )

        plt.plot( im_xyzs[i][0], im_xyzs[i][1], marker )
        
        cos_az = np.cos(im_azs[i]/57.3)
        sin_az = np.sin(im_azs[i]/57.3)
        cos_el = np.cos(im_els[i]/57.3)

        if RAD_paths[i].split('\\')[-1][0]=='Z':
            plt.arrow( im_xyzs[i][0], im_xyzs[i][1], scale*cos_el*sin_az, scale*cos_el*cos_az,
                       color=marker[1], lw = int(scale/10), linestyle='dashed' )
        else:
            plt.arrow( im_xyzs[i][0], im_xyzs[i][1], scale*cos_el*sin_az, scale*cos_el*cos_az,
                       color=marker[1], lw = int(scale/10) )

    plt.axis('equal')
    plt.xlim( [ np.round(plt.gca().get_xlim()[0])-3, np.round(plt.gca().get_xlim()[1])+3 ] )
    plt.ylim( [ np.round(plt.gca().get_ylim()[0])-3, np.round(plt.gca().get_ylim()[1])+3 ] )

    plt.xlabel( 'Easting Site Frame')
    plt.ylabel( 'Northing Site Frame')
#     plt.tight_layout()
#     plt.savefig( path + '/positions'+suf+'.jpg', dpi=300  )




def xyz_shift_offsets( site, drive ):
    
    '''
    xyz_shift_offsets finds most accurate Site-Nav offset for each site index and drive
    
    future work: 
        * load parameters from an spreadsheet
        * interpolate between drives for more acurate localization for mid-drive images
    '''

    if 0:
        x_shift, y_shift, z_shift = [ 0,0,0 ]
        
    elif site== 33 and drive>= 262: x_shift, y_shift, z_shift = [ -2365.65193, 756.0839, 30.94624 ]
    elif site== 33 and drive>= 0: x_shift, y_shift, z_shift = [ -2365.93, 756.577, 31.424 ]
    elif site== 32 and drive>= 2340: x_shift, y_shift, z_shift = [ -2456.1649, 463.127000000001, 24.9395 ]
    elif site== 32 and drive>= 1394: x_shift, y_shift, z_shift = [ -2456.61, 465.855000000001, 25.042036 ]
    elif site== 32 and drive>= 1184: x_shift, y_shift, z_shift = [ -2457.2022, 466.5075, 24.7525469 ]
    elif site== 32 and drive>= 1174: x_shift, y_shift, z_shift = [ -2457.2047, 466.5164, 24.7574308 ]
    elif site== 32 and drive>= 1038: x_shift, y_shift, z_shift = [ -2457.2206, 466.9109, 24.886537 ]
    elif site== 32 and drive>= 896: x_shift, y_shift, z_shift = [ -2457.2561, 466.9587, 25.002247 ]
    elif site== 32 and drive>= 774: x_shift, y_shift, z_shift = [ -2457.4309, 467.4988, 24.805481 ]
    elif site== 32 and drive>= 672: x_shift, y_shift, z_shift = [ -2456.3776, 466.8088, 25.053988 ]
    elif site== 32 and drive>= 580: x_shift, y_shift, z_shift = [ -2457.2906, 467.257, 24.97048 ]
    elif site== 32 and drive>= 482: x_shift, y_shift, z_shift = [ -2456.9152, 467.1303, 24.987918 ]
    elif site== 32 and drive>= 378: x_shift, y_shift, z_shift = [ -2456.9546, 467.49592, 24.936184 ]
    elif site== 32 and drive>= 250: x_shift, y_shift, z_shift = [ -2457.1925, 467.918052, 24.848803 ]
    elif site== 32 and drive>= 148: x_shift, y_shift, z_shift = [ -2456.996, 467.6681, 24.908544 ]
    elif site== 32 and drive>= 0: x_shift, y_shift, z_shift = [ -2456.636, 467.802, 24.748 ]
    elif site== 31 and drive>= 856: x_shift, y_shift, z_shift = [ -2792.2897, 397.7234, 45.06612 ]
    elif site== 31 and drive>= 500: x_shift, y_shift, z_shift = [ -2791.6302, 398.9164, 44.80022 ]
    elif site== 31 and drive>= 0: x_shift, y_shift, z_shift   = [ -2790.25, 397.382, 44.625 ]
    elif site== 30 and drive>= 2188: x_shift, y_shift, z_shift = [ -2790.2754, 397.0656, 44.28703 ]
    elif site== 30 and drive>= 1776: x_shift, y_shift, z_shift = [ -2788.6738, 397.4148, 44.35993 ]
    elif site== 30 and drive>= 1524: x_shift, y_shift, z_shift = [ -2788.1609, 397.4227, 44.31772 ]
    elif site== 30 and drive>= 1344: x_shift, y_shift, z_shift = [ -2787.764, 397.3139, 44.09667 ]
    elif site== 30 and drive>= 1172: x_shift, y_shift, z_shift = [ -2787.676, 396.852000000001, 43.87697 ]
    elif site== 30 and drive>= 1096: x_shift, y_shift, z_shift = [ -2787.57, 397.021000000001, 43.91903 ]    
    elif site== 30 and drive>= 0: x_shift, y_shift, z_shift = [ -2790.208, 397.307, 44.625 ]
    elif site== 29 and drive>= 144: x_shift, y_shift, z_shift = [ -2779.67751, 360.9082, 43.36229 ]
    elif site== 29 and drive>= 58: x_shift, y_shift, z_shift = [ -2779.57856, 360.9551, 43.389303 ]
    elif site== 28 and drive>= 0: x_shift, y_shift, z_shift = [ -2779.803, 362.544, 43.391 ]
    elif site== 27 and drive>= 230: x_shift, y_shift, z_shift = [ -2759.3294, 337.0496, 42.477692 ]
    elif site== 29 and drive>= 0: x_shift, y_shift, z_shift = [ -2779.953, 361.526, 43.385 ]
    elif site== 27 and drive>= 0: x_shift, y_shift, z_shift = [ -2759.746, 337.659, 42.545 ]
    elif site== 26 and drive>= 5226: x_shift, y_shift, z_shift = [ -2507.6057, 768.672999999999, 35.10683 ]
    elif site== 26 and drive>= 5214: x_shift, y_shift, z_shift = [ -2507.5848, 768.698, 35.07032 ]
    elif site== 26 and drive>= 5150: x_shift, y_shift, z_shift = [ -2507.6548, 769.068, 34.77809 ]
    elif site== 26 and drive>= 4200: x_shift, y_shift, z_shift = [ -2507.8953, 767.539999999999, 34.98978 ]
    elif site== 26 and drive>= 3482: x_shift, y_shift, z_shift = [ -2508.8886, 765.945, 34.82269 ]
    elif site== 26 and drive>= 3008: x_shift, y_shift, z_shift = [ -2509.2189, 764.168999999999, 35.05277 ]
    elif site== 26 and drive>= 2548: x_shift, y_shift, z_shift = [ -2509.657, 762.9411, 35.72451 ]
    elif site== 26 and drive>= 2018: x_shift, y_shift, z_shift = [ -2508.192, 762.0176, 35.63733 ]
    elif site== 26 and drive>= 1652: x_shift, y_shift, z_shift = [ -2507.5925, 762.3143, 35.63282 ]
    elif site== 26 and drive>= 1222: x_shift, y_shift, z_shift = [ -2506.0705, 762.6273, 35.1771000000001 ]
    elif site== 26 and drive>= 1154: x_shift, y_shift, z_shift = [ -2505.8738, 762.4009, 35.0390000000001 ]
    elif site== 26 and drive>= 1004: x_shift, y_shift, z_shift = [ -2505.74932, 762.1176, 35.5139 ]
    elif site== 26 and drive>= 850: x_shift, y_shift, z_shift = [ -2505.9002, 762.1751, 35.0638 ]
    elif site== 26 and drive>= 756: x_shift, y_shift, z_shift = [ -2506.1233, 762.2654, 35.0623000000001 ]
    elif site== 26 and drive>= 630: x_shift, y_shift, z_shift = [ -2506.27595, 762.5299, 35.1129 ]
    elif site== 26 and drive>= 470: x_shift, y_shift, z_shift = [ -2506.1379, 762.73, 35.1706000000001 ]
    elif site== 26 and drive>= 218: x_shift, y_shift, z_shift = [ -2506.8696, 762.903, 35.25684 ]    
    elif site== 26 and drive>= 96: x_shift, y_shift, z_shift = [ -2506.57645, 763.1905, 35.02604 ]
    elif site== 26 and drive>= 0: x_shift, y_shift, z_shift = [ -2506.646, 763.608, 35.132 ]
    elif site== 25 and drive>= 446: x_shift, y_shift, z_shift = [ -2377.102, 773.1197, 34.057027 ]
    elif site== 25 and drive>= 214: x_shift, y_shift, z_shift = [ -2863.463, 353.292, 47.6723 ]
    elif site== 25 and drive>= 94: x_shift, y_shift, z_shift = [ -2378.3308, 772.55288, 34.159531 ]
    elif site== 25 and drive>= 0: x_shift, y_shift, z_shift = [ -2378.559, 772.6, 34.28 ]
    elif site== 24 and drive>= 4470: x_shift, y_shift, z_shift = [ -2784.097, 344.957, 44.7962 ]
    elif site== 24 and drive>= 4412: x_shift, y_shift, z_shift = [ -2784.589, 344.875000000001, 44.6932 ]
    elif site== 24 and drive>= 4358: x_shift, y_shift, z_shift = [ -2784.123, 345.111, 44.9133 ]
    elif site== 24 and drive>= 4232: x_shift, y_shift, z_shift = [ -2784.009, 345.824, 45.1816 ]
    elif site== 24 and drive>= 3900: x_shift, y_shift, z_shift = [ -2784.822, 345.428000000001, 44.9822 ]
    elif site== 24 and drive>= 3290: x_shift, y_shift, z_shift = [ -2785.105, 347.497000000001, 44.6974 ]
    elif site== 24 and drive>= 3170: x_shift, y_shift, z_shift = [ -2784.907, 347.98, 44.5468999999999 ]
    elif site== 24 and drive>= 2770: x_shift, y_shift, z_shift = [ -2784.878, 348.168000000001, 44.494 ]
    elif site== 24 and drive>= 1992: x_shift, y_shift, z_shift = [ -2782.724, 350.576, 44.8594 ]
    elif site== 24 and drive>= 1970: x_shift, y_shift, z_shift = [ -2782.817, 350.504, 44.8395 ]
    elif site== 24 and drive>= 1328: x_shift, y_shift, z_shift = [ -2783.57156458304, 352.3443, 44.4426 ]
    elif site== 24 and drive>= 902: x_shift, y_shift, z_shift = [ -2782.27490275483, 351.9806, 43.99387 ]
    elif site== 24 and drive>= 0: x_shift, y_shift, z_shift = [ -2780.97097343423, 352.043, 43.445 ]
    elif site== 23 and drive>= 3142: x_shift, y_shift, z_shift = [ -2517.86510288296, 436.8147, 25.0537000000001 ]
    elif site== 23 and drive>= 2872: x_shift, y_shift, z_shift = [ -2518.33886184144, 437.198999999999, 25.2623 ]
    elif site== 23 and drive>= 2532: x_shift, y_shift, z_shift = [ -2519.0071656846, 438.075999999999, 24.8130000000001 ]
    elif site== 23 and drive>= 2096: x_shift, y_shift, z_shift = [ -2518.80717808731, 438.161999999999, 24.6291 ]
    elif site== 23 and drive>= 1180: x_shift, y_shift, z_shift = [ -2520.90882184935, 436.184, 25.2165000000001 ]
    elif site== 23 and drive>= 824: x_shift, y_shift, z_shift = [ -2520.79226220751, 435.622, 25.26953 ]
    elif site== 23 and drive>= 0: x_shift, y_shift, z_shift = [ -2520.80223501524, 433.615, 25.12 ]
    elif site== 22 and drive>= 788: x_shift, y_shift, z_shift = [ -2273.7408528018, 553.3168, 23.57345 ]
    elif site== 22 and drive>= 634: x_shift, y_shift, z_shift = [ -2274.2041555917, 553.2385, 23.63658 ]
    elif site== 22 and drive>= 532: x_shift, y_shift, z_shift = [ -2274.30828523706, 552.95073, 23.93847 ]
    elif site== 22 and drive>= 412: x_shift, y_shift, z_shift = [ -2274.02412229753, 552.5224, 24.07615 ]
    elif site== 22 and drive>= 0: x_shift, y_shift, z_shift = [ -2275.14664997495, 552.38, 24.309 ]
    elif site== 21 and drive>= 2746: x_shift, y_shift, z_shift = [ -2100.36943232086, 547.29, 21.314 ]
    elif site== 21 and drive>= 2560: x_shift, y_shift, z_shift = [ -1633.16766411043, 741.719, 20.171649 ]
    elif site== 21 and drive>= 1126: x_shift, y_shift, z_shift = [ -1634.86543382004, 739.830999999999, 20.3523 ]
    elif site== 21 and drive>= 0: x_shift, y_shift, z_shift = [ -1636.2605840328, 738.027, 20.042 ]
    elif site== 20 and drive>= 1422: x_shift, y_shift, z_shift = [ -1173.22753966411, 949.798, 17.30686 ]
    elif site== 20 and drive>= 0: x_shift, y_shift, z_shift = [ -1174.41393238602, 947.725, 17.851 ]
    elif site== 19 and drive>= 2066: x_shift, y_shift, z_shift = [ -661.715855069015, 1154.713, 12.71619 ]
    elif site== 19 and drive>= 1126: x_shift, y_shift, z_shift = [ -663.061026936431, 1153.289, 12.707 ]
    elif site== 19 and drive>= 0: x_shift, y_shift, z_shift = [ -664.104104559188, 1151.764, 12.919 ]
    elif site== 18 and drive>= 1622: x_shift, y_shift, z_shift = [ -276.619254179765, 1350.702, 12.85237 ]
    elif site== 18 and drive>= 320: x_shift, y_shift, z_shift = [ -276.91806500552, 1348.3362, 12.1960758 ]
    elif site== 18 and drive>= 0: x_shift, y_shift, z_shift = [ -277.37471203441, 1348.107, 12.167 ]
    elif site== 17 and drive>= 1064: x_shift, y_shift, z_shift = [ 73.1858874593579, 1397.4184, 7.24982 ]
    elif site== 17 and drive>= 0: x_shift, y_shift, z_shift = [ 71.5034355566901, 1397.758, 7.793 ]
    elif site== 16 and drive>= 0: x_shift, y_shift, z_shift = [ 199.90989036014, 1403.036, 6.72 ]
    elif site== 15 and drive>= 1334: x_shift, y_shift, z_shift = [ 603.820390305231, 1176.952, 1.73934000000001 ]
    elif site== 15 and drive>= 0: x_shift, y_shift, z_shift = [ 601.229265583504, 1176.533, 2.169 ]
    elif site== 14 and drive>= 1176: x_shift, y_shift, z_shift = [ 622.620908260725, 747.927, 1.273707 ]
    elif site== 14 and drive>= 0: x_shift, y_shift, z_shift = [ 621.342422723249, 750.545, 1.089 ]
    elif site== 13 and drive>= 2448: x_shift, y_shift, z_shift = [ 453.568712555918, 384.848, -0.59128999999999 ]
    elif site== 13 and drive>= 1280: x_shift, y_shift, z_shift = [ 452.374842550113, 386.591, -0.535944 ]
    elif site== 13 and drive>= 0: x_shift, y_shift, z_shift = [ 452.860466749729, 389.596, -1.292 ]
    elif site== 12 and drive>= 0: x_shift, y_shift, z_shift = [ 193.349994355496, 258.971, 0.154 ]
    elif site== 11 and drive>= 2472: x_shift, y_shift, z_shift = [ 89.2529943554968, -47.689, 1.41034 ]
    elif site== 11 and drive>= 1110: x_shift, y_shift, z_shift = [ 89.3270381853386, -44.649, 1.1154259 ]
    elif site== 11 and drive>= 108: x_shift, y_shift, z_shift = [ 87.0732292114366, -43.005, 1.188 ]
    elif site== 11 and drive>= 0: x_shift, y_shift, z_shift = [ 87.2619827257403, -43.457, 1.188 ]
    elif site== 10 and drive>= 1114: x_shift, y_shift, z_shift = [ 57.885182391218, -235.108, 0.173082000000001 ]
    elif site== 10 and drive>= 0: x_shift, y_shift, z_shift = [ 56.6108781291436, -233.38, 0.191 ]
    elif site== 9 and drive>= 9314: x_shift, y_shift, z_shift = [ 30.7753594225956, -269.049, 0.15 ]
    elif site== 9 and drive>= 7894: x_shift, y_shift, z_shift = [ 47.4748286327977, -515.725, -0.903 ]
    elif site== 9 and drive>= 6884: x_shift, y_shift, z_shift = [ 65.4253826940879, -697.926, -5.091 ]
    elif site== 9 and drive>= 5676: x_shift, y_shift, z_shift = [ 75.7954438593235, -854.339, -12.072 ]
    elif site== 9 and drive>= 2982: x_shift, y_shift, z_shift = [ -528.570601157329, -709.806000000001, -3.78408 ]
    elif site== 9 and drive>= 2858: x_shift, y_shift, z_shift = [ -539.76716477192, -709.688, -3.97515 ]
    elif site== 9 and drive>= 2818: x_shift, y_shift, z_shift = [ -539.548982780045, -709.758, -3.90320000000001 ]
    elif site== 9 and drive>= 1554: x_shift, y_shift, z_shift = [ -537.558718905617, -709.156000000001, -4.408639 ]
    elif site== 9 and drive>= 364: x_shift, y_shift, z_shift = [ -534.984542891075, -709.13755, -4.727058 ]
    elif site== 9 and drive>= 276: x_shift, y_shift, z_shift = [ -534.984542891075, -709.13755, -4.727058 ]
    elif site== 9 and drive>= 178: x_shift, y_shift, z_shift = [ -535.043353278405, -709.14076, -4.776216 ]
    elif site== 9 and drive>= 82: x_shift, y_shift, z_shift = [ -534.960427769811, -709.04462, -4.604324 ]
    elif site== 9 and drive>= 0: x_shift, y_shift, z_shift = [ -534.82976299399, -709.314, -4.528 ]
    elif site== 8 and drive>= 734: x_shift, y_shift, z_shift = [ -441.686695930272, -630.7351, 0.73767 ]
    elif site== 8 and drive>= 458: x_shift, y_shift, z_shift = [ -441.632947637872, -631.2016, 0.89597 ]
    elif site== 8 and drive>= 410: x_shift, y_shift, z_shift = [ -442.461622935661, -630.911, 0.86209 ]
    elif site== 8 and drive>= 256: x_shift, y_shift, z_shift = [ -442.035662724996, -631.39685, 0.81640999999999 ]
    elif site== 8 and drive>= 64: x_shift, y_shift, z_shift = [ -441.974337550378, -632.05725, 0.884683 ]
    elif site== 8 and drive>= 0: x_shift, y_shift, z_shift = [ -442.276042296893, -632.051, 0.756 ]
    elif site== 7 and drive>= 2378: x_shift, y_shift, z_shift = [ -363.682592710141, -824.91, -4.816 ]
    elif site== 7 and drive>= 2280: x_shift, y_shift, z_shift = [ -363.650158355931, -824.677, -4.714 ]
    elif site== 7 and drive>= 2216: x_shift, y_shift, z_shift = [ -363.765833525645, -824.447, -4.942 ]
    elif site== 7 and drive>= 2050: x_shift, y_shift, z_shift = [ -363.279463749648, -824.218, -5.138 ]
    elif site== 7 and drive>= 1836: x_shift, y_shift, z_shift = [ -363.550000256358, -823.837, -4.906 ]
    elif site== 7 and drive>= 1716: x_shift, y_shift, z_shift = [ -363.482388489528, -824.019, -4.856 ]
    elif site== 7 and drive>= 1556: x_shift, y_shift, z_shift = [ -363.281536978566, -823.758, -5.013 ]
    elif site== 7 and drive>= 1358: x_shift, y_shift, z_shift = [ -363.047164616848, -823.173, -5.116 ]
    elif site== 7 and drive>= 1172: x_shift, y_shift, z_shift = [ -362.788897312806, -822.954, -5.23255 ]
    elif site== 7 and drive>= 346: x_shift, y_shift, z_shift = [ -365.491464390023, -823.204, -4.467 ]
    elif site== 7 and drive>= 0: x_shift, y_shift, z_shift = [ -366.570706899171, -823.192, -4.614 ]
    elif site== 6 and drive>= 2666: x_shift, y_shift, z_shift = [ 49.2552015925663, -1002.496, -15.183 ]
    elif site== 6 and drive>= 2250: x_shift, y_shift, z_shift = [ 48.2690550955953, -1003.131, -14.713 ]
    elif site== 6 and drive>= 1752: x_shift, y_shift, z_shift = [ 47.2862899018857, -1002.662, -14.71 ]
    elif site== 6 and drive>= 1648: x_shift, y_shift, z_shift = [ 47.1228902439419, -1002.741, -15.084 ]
    elif site== 6 and drive>= 1450: x_shift, y_shift, z_shift = [ 45.8303964633962, -1002.484, -14.88 ]
    elif site== 6 and drive>= 892: x_shift, y_shift, z_shift = [ 44.8359770584786, -1002.884, -14.868 ]
    elif site== 6 and drive>= 410: x_shift, y_shift, z_shift = [ 43.4471135109867, -1002.711, -15.069 ]
    elif site== 6 and drive>= 170: x_shift, y_shift, z_shift = [ 42.6807860286917, -1002.702, -14.934 ]
    elif site== 6 and drive>= 0: x_shift, y_shift, z_shift = [ 42.9589719413418, -1003.669, -15.012 ]
    elif site== 5 and drive>= 2512: x_shift, y_shift, z_shift = [ 42.1951371772719, -519.968, -0.946999999999999 ]
    elif site== 5 and drive>= 2296: x_shift, y_shift, z_shift = [ 42.1400458754569, -520.256, -1.036 ]
    elif site== 5 and drive>= 1812: x_shift, y_shift, z_shift = [ 42.7783919369836, -521.464, -1.343 ]
    elif site== 5 and drive>= 1572: x_shift, y_shift, z_shift = [ 43.2635505157394, -522.221, -1.354 ]
    elif site== 5 and drive>= 1388: x_shift, y_shift, z_shift = [ 43.6779013539158, -522.673, -1.187 ]
    elif site== 5 and drive>= 894: x_shift, y_shift, z_shift = [ 44.4290071167894, -523.826, -0.89 ]
    elif site== 5 and drive>= 500: x_shift, y_shift, z_shift = [ 45.7907972365013, -524.803, -0.974 ]
    elif site== 5 and drive>= 0: x_shift, y_shift, z_shift = [ 47.1997606973803, -525.883, -0.724 ]
    elif site== 4 and drive>= 2222: x_shift, y_shift, z_shift = [ 76.5417490013376, -17.429, 0.501 ]
    elif site== 4 and drive>= 1878: x_shift, y_shift, z_shift = [ 77.1254205139582, -17.393, 0.442 ]
    elif site== 4 and drive>= 1860: x_shift, y_shift, z_shift = [ 77.0841330319459, -17.351, 0.413 ]
    elif site== 4 and drive>= 1776: x_shift, y_shift, z_shift = [ 77.0841330319459, -17.351, 0.413 ]
    elif site== 4 and drive>= 1712: x_shift, y_shift, z_shift = [ 77.9003835423752, -18, 0.354 ]
    elif site== 4 and drive>= 1644: x_shift, y_shift, z_shift = [ 78.3710991162154, -18.659, 0.24 ]
    elif site== 4 and drive>= 1422: x_shift, y_shift, z_shift = [ 78.3710991162154, -18.659, 0.24 ]
    elif site== 4 and drive>= 1250: x_shift, y_shift, z_shift = [ 79.1452258194926, -18.453, 0.213 ]
    elif site== 4 and drive>= 1062: x_shift, y_shift, z_shift = [ 79.4493718734854, -18.966, 0.05 ]
    elif site== 4 and drive>= 922: x_shift, y_shift, z_shift = [ 80.0837877886232, -19.831, 0.189 ]
    elif site== 4 and drive>= 822: x_shift, y_shift, z_shift = [ 80.1387312493609, -20.767, 0.192 ]
    elif site== 4 and drive>= 738: x_shift, y_shift, z_shift = [ 79.6605612081629, -20.992, 0.41 ]
    elif site== 4 and drive>= 644: x_shift, y_shift, z_shift = [ 79.3242794897827, -21.108, 0.325 ]
    elif site== 4 and drive>= 592: x_shift, y_shift, z_shift = [ 79.3652015309346, -21.675, 0.495 ]
    elif site== 4 and drive>= 510: x_shift, y_shift, z_shift = [ 79.097083205629, -22.15, 0.447 ]
    elif site== 4 and drive>= 430: x_shift, y_shift, z_shift = [ 79.169992684349, -22.5365000000001, 0.451844 ]
    elif site== 4 and drive>= 372: x_shift, y_shift, z_shift = [ 79.03968547493, -22.919, 0.648 ]
    elif site== 4 and drive>= 218: x_shift, y_shift, z_shift = [ 78.7591294758789, -23.672, 0.806 ]
    elif site== 4 and drive>= 136: x_shift, y_shift, z_shift = [ 78.8072612752607, -23.626, 0.787 ]
    elif site== 4 and drive>= 48: x_shift, y_shift, z_shift = [ 78.3718579316058, -24.027, 0.73 ]
    elif site== 4 and drive>= 0: x_shift, y_shift, z_shift = [ 78.1942259581876, -24.177, 0.61 ]
    elif site== 3 and drive>= 2430: x_shift, y_shift, z_shift = [ -4.29537206183093, 0.84, 1.889 ]
    elif site== 3 and drive>= 2346: x_shift, y_shift, z_shift = [ -3.93510098768898, 0.651999999999998, 1.719 ]
    elif site== 3 and drive>= 2208: x_shift, y_shift, z_shift = [ -3.2421363891512, 0.612, 1.681 ]
    elif site== 3 and drive>= 2150: x_shift, y_shift, z_shift = [ -2.8742231207639, 0.159, 1.748 ]
    elif site== 3 and drive>= 2046: x_shift, y_shift, z_shift = [ -2.71955156142316, -0.316, 1.623 ]
    elif site== 3 and drive>= 1950: x_shift, y_shift, z_shift = [ -2.26360489415657, -0.317, 1.645 ]
    elif site== 3 and drive>= 1850: x_shift, y_shift, z_shift = [ -2.22216259443773, -0.599, 1.509 ]
    elif site== 3 and drive>= 1708: x_shift, y_shift, z_shift = [ -2.12295671931542, -0.543999999999999, 1.336 ]
    elif site== 3 and drive>= 1416: x_shift, y_shift, z_shift = [ -0.730249654556161, -0.305, 0.983 ]
    elif site== 3 and drive>= 1398: x_shift, y_shift, z_shift = [ -0.497864674634571, -0.460999999999999, 0.983 ]
    elif site== 3 and drive>= 1392: x_shift, y_shift, z_shift = [ -0.675545006739909, -0.0589999999999993, 0.973 ]
    elif site== 3 and drive>= 1374: x_shift, y_shift, z_shift = [ -0.531401265733731, -0.112, 0.881 ]
    elif site== 3 and drive>= 1266: x_shift, y_shift, z_shift = [ -1.20790695928578, 0.00499999999999989, 0.842 ]
    elif site== 3 and drive>= 1044: x_shift, y_shift, z_shift = [ -1.13244386078872, -0.613, 0.688 ]
    elif site== 3 and drive>= 828: x_shift, y_shift, z_shift = [ -1.54665842428645, -0.754999999999996, 0.829 ]
    elif site== 3 and drive>= 792: x_shift, y_shift, z_shift = [ -1.83157067707323, -0.774000000000001, 0.817 ]
    elif site== 3 and drive>= 770: x_shift, y_shift, z_shift = [ -1.36344032575406, -0.928000000000001, 0.857 ]
    elif site== 3 and drive>= 578: x_shift, y_shift, z_shift = [ 0.21159490692979, -0.628999999999998, 0.616 ]
    elif site== 3 and drive>= 386: x_shift, y_shift, z_shift = [ 0.210352483635505, -1.12, 0.369 ]
    elif site== 3 and drive>= 38:  x_shift, y_shift, z_shift = [ 2.77076111174037, -4.351, 0.04 ]
    elif site== 3 and drive>= 0:   x_shift, y_shift, z_shift = [ 0, 0, 0 ]
    # else:                          x_shift, y_shift, z_shift = [ 0, 0, 0 ]
        
    return x_shift, y_shift, z_shift

