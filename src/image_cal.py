file_extension = ''

'''
future work: save these calibration prameters as a text files, which we load for each camera
'''

# save images and thereby overwrite existing images
save_im    = 1

# add an alpha channel to the output images
save_mask  = 1

# add transparrent pixels to restore the image's full, standard size
pad_im     = 1
pad_im_z   = 1

# turn on when finding the waypoint offsets
find_offsets_mode = 0

# set the color values
gamma      = 2.2      # gamma value
gamma      = 2        # gamma value

# fraction of the dynamic range to clip off the lower values of image 
clip_low_z = 0.02  # for the Mastcam-Z cameras
clip_low   = 0.05  # for everything else


# scale all the scale parameters below bsy the same number
scale_scale = 20

# color balance parameters for the Mars 2020 science cameras
scale_z,  scale_red_z,  scale_blue_z  = [ 1.0*scale_scale, 0.7 , 1.5  ] # Mastcam-Z 
scale_l,  scale_red_l,  scale_blue_l  = [ 1.0*scale_scale, 0.75, 1.40 ] # SuperCam RMI
scale_s,  scale_red_s,  scale_blue_s  = [ 1.0*scale_scale, 0.85, 1.40 ] # SHERLOC WATSON 

# color balance parameters for the Mars 2020 engineering cameras
scale_n,  scale_red_n,  scale_blue_n  = [ 1.0*scale_scale, 0.75, 1.2  ] # Navcam
scale_f,  scale_red_f,  scale_blue_f  = [ 1.1*scale_scale, 0.78, 1.25 ] # Front Hazcam
scale_r,  scale_red_r,  scale_blue_r  = [ 1.1*scale_scale, 0.78, 1.25 ] # Rear Hazcam
scale_v,  scale_red_v,  scale_blue_v  = [ 1.1*scale_scale, 1.12, 0.92 ] # Grayscale VCE Navcam
scale_hr, scale_red_hr, scale_blue_hr = [ 1.0*scale_scale, 0.75, 1.43 ] # Inginuity RTE
scale_hn, scale_red_hn, scale_blue_hn = [ 1.0*scale_scale, 1.1 , 0.92 ] # Inginuity Navcam
