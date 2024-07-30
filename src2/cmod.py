


# def cahvor_pds2ms( cahvor ):

def cahvor_separate( cahvor ):

    c = cahvor[ 0: 3]
    a = cahvor[ 3: 6]
    h = cahvor[ 6: 9]
    v = cahvor[ 9:12]
    o = cahvor[12:15]
    r = cahvor[15:18]

    return c,a,h,v,o,r

def cahvor_combine( c,a,h,v,o,r ):

    cahvor = np.stack([c,a,h,v,o,r])

    return cahvor

def cahvor_from_pds( cmod ):

    c = cmod['MODEL_COMPONENT_1']
    a = cmod['MODEL_COMPONENT_2']
    h = cmod['MODEL_COMPONENT_3']
    v = cmod['MODEL_COMPONENT_4']
    o = cmod['MODEL_COMPONENT_5']
    r = cmod['MODEL_COMPONENT_6']

    cahvor = cahvor_combine( c,a,h,v,o,r )

    return cahvor


# def cahvor_to_pds( cahvor ):

def cahvor_int( cahvor ):

        '''
        This function decomposes a cahvor 16-vector into its intrinsic parameters

        Input:  cahvor model in f-frame 
        Output: rotation matrix and offset vector from the camera-frame to f-frame
    '''

    return f, b1, b2, cx, cy


def cahvor_dist( cahvor ):

    return f, b1, b2, cx, cy



def cahvor_ext( cahvor_f ):

    '''
        This function decomposes a cahvor 16-vector into its extrinsic parameters in the form of a matrix and vector. These are the transformation from 

        Input:  cahvor model in f-frame 
        Output: rotation matrix and offset vector from the camera-frame to f-frame
    '''
    return R_fc, t_fc