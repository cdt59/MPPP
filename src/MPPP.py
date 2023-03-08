# MPPP functions

'''

'''

import numpy as np
from planetaryimage import PDS3Image
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interp1d
import colour_demosaicing
from PIL import Image
import os


class image:
    
    '''
    The image class holds all parameters specific to the IMG file
    '''

    def __init__(self, IMG_path ):
        
        self.IMG_path    = IMG_path
        self.filename    = os.path.basename( IMG_path )
        self.label       = PDS3Image.open( IMG_path ).label                  # PDS image header metadata
        self.image       = np.float32( PDS3Image.open( IMG_path ).image )    # three band float-32 image array
        self.mask_image  = np.ones( self.image.shape[:2] )*255               # one band boolian image array
        self.cam         = self.filename[:2]
        self.sol         = int( self.filename[4:8] )
        
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
        
        '''
        return the tau, the opacity of Mars sky at the time of the image
        
        future work:
        * access the M2020 tau record
        * estimate the tue for the specific sol
        '''
        
        self.tau = 0.8
        
        if self.sol >=700:
            self.tau = 0.5

    

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




def xyz_shift_offsets( site, drive ):
    
    '''
    xyz_shift_offsets finds most accurate Site-Nav offset for each site index and drive
    
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
    x_shift, y_shift, z_shift = drive_site_shift[9:]

    # x_shift, y_shift, z_shift = [ 0,0,0 ]
 
    return x_shift, y_shift, z_shift

