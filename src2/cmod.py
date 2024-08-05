import numpy as np
from scipy.spatial.transform import Rotation as R
import requests
from enum import Enum

url = "https://mars.nasa.gov/mmgis-maps/M20/Layers/json/M20_waypoints.json"
m20_waypoints_data = requests.get(url).json()
m20_sitedrive2tr = {(feat["properties"]["site"], feat["properties"]["drive"]): np.array(
    [feat["properties"]["northing"], feat["properties"]["easting"], -feat["properties"]["elev_geoid"]]) for feat in m20_waypoints_data["features"]}
site2drive = {k[0]: [] for k in m20_sitedrive2tr.keys()}
for k in m20_sitedrive2tr.keys():
    site2drive[k[0]].append(k[1])


def q_wxyz2xyzw(q_wxyz):
    """
    Converts quaternion with real part at idx 0 to real part at idx 3
    """
    return np.array([q_wxyz[1], q_wxyz[2], q_wxyz[3], q_wxyz[0]])


def cahvor_from_pds(cmod):
    """
    Converts PDS label to 18 element cahvor vector
    """
    c = cmod['MODEL_COMPONENT_1']
    a = cmod['MODEL_COMPONENT_2']
    h = cmod['MODEL_COMPONENT_3']
    v = cmod['MODEL_COMPONENT_4']
    o = cmod['MODEL_COMPONENT_5']
    r = cmod['MODEL_COMPONENT_6']

    cahvor = cahvor_combine(c, a, h, v, o, r)

    return cahvor


def Rt_rm_from_pds(cmod):
    """
    Parse the rotation and translation from mast frame to rover frame
    """
    t_rm = cmod['MODEL_TRANSFORM_VECTOR']
    q_rm = q_wxyz2xyzw(cmod['MODEL_TRANSFORM_QUATERNION'])
    R_rm = R.from_quat(q_rm).as_matrix()

    return R_rm, t_rm


def cahvor_H(cahvor):
    """
    Returns a H vector given an 18 element cahvor vector
    """

    a = cahvor[3: 6]
    h = cahvor[6: 9]
    v = cahvor[9:12]

    H = np.array([[h[0], v[0], a[0]],
                  [h[1], v[1], a[1]],
                  [h[2], v[2], a[2]]])

    return H


def cahvor_Hf(cahvor):
    """
    Returns a fancy H vector given an 18 element cahvor vector
    """
    c = cahvor[0: 3]
    a = cahvor[3: 6]
    h = cahvor[6: 9]
    v = cahvor[9:12]

    Hf = np.array([[h[0], v[0], a[0], c[0]],
                   [h[1], v[1], a[1], c[1]],
                   [h[2], v[2], a[2], c[2]],
                   [0,    0,    0,   1]])

    return Hf


def cahvor_separate(cahvor):
    """
    Seperates an 18 elemetn cahvor vector to c, a, h, v, o & r vectors
    """
    c = cahvor[0: 3]
    a = cahvor[3: 6]
    h = cahvor[6: 9]
    v = cahvor[9:12]
    o = cahvor[12:15]
    r = cahvor[15:18]

    return c, a, h, v, o, r


def cahvor_combine(c, a, h, v, o, r):
    """
    Combines individual c,a,h,v,o,r vectors in a single 18 element cahvor vector
    """

    cahvor = np.stack([c, a, h, v, o, r]).flatten()

    return cahvor


def cahvor_intr(cahvor):

    c, a, h, v, o, r = cahvor_separate(cahvor)

    hs = np.linalg.norm(np.cross(h, a))
    vs = np.linalg.norm(np.cross(v, a))
    hc = np.dot(h, a)
    vc = np.dot(v, a)

    hp = (h - hc*a) / hs
    vp = (v - vc*a) / vs

    theta = np.arcsin(np.clip(np.linalg.norm(
        np.cross(vp, hp)), a_min=-1, a_max=1))
    if theta > 0:
        theta *= -1

    f = vs
    b1 = - hs * np.sin(theta) - vs
    b2 = hs * np.cos(theta)
    cx = hc
    cy = vc

    intr = np.array([f, b1, b2, cx, cy])

    return intr


def cahvor_K(cahvor):
    '''
        This function returns the calibration matrix K from 18 element cahvor vector

        Input:  cahvor model vector
        Output: K_c camera matrix
    '''

    f, b1, b2, cx, cy = cahvor_intr(cahvor)
    K_c = np.array([[f+b1,  b2,  cx],
                    [0,   f,  cy],
                    [0,   0,   1]])

    return K_c


def cahvor_xyz(cahvor_a):

    xyz_ac = cahvor_a[0: 3]

    return xyz_ac


def find_opk_from_R(R_mat):
    """
    Converts rotation matrix to omega, phi and kappa angles
    """
    # finds omega, phi, kappa from rotation matrix R_cam2site
    R_ = R.from_matrix(R_mat)
    R_enu2ned = R.from_matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    angles = (R_enu2ned * R_).as_euler('XYZ', degrees=1)
    opk = [angles[0], angles[1], angles[2]+90]

    return np.array(opk)


def cahvor_opk(R_ac):
    return find_opk_from_R(R_ac)


def cahvor_dist(cahvor):
    '''
        This function decomposes a cahvor 16-vector into its distortion parameters, in 5-vector (k1,k2,k3,p1,p2)

        Input:  cahvor model vector
        Output: distortion vector
    '''

    c, a, h, v, o, r = cahvor_separate(cahvor)

    k1 = r[1]
    k2 = r[2]
    k3 = 0
    p1 = 0
    p2 = 0

    dist = np.array([k1, k2, k3, p1, p2])

    return dist


def cahvor_Rt(cahvor_a):
    '''
        This function decomposes a cahvor 16-vector into its extrinsic parameters in the form of a matrix and vector. These are the transformation from

        Input:  cahvor model in a-frame
        Output: rotation matrix and offset vector from the camera-frame to a-frame, R_ac and t_ac
    '''
    K_c = cahvor_K(cahvor_a)
    H_c = K_c.T
    H_f = cahvor_H(cahvor_a)
    R_ac = H_f @ np.linalg.inv(H_c)
    t_ac = cahvor_a[0: 3]

    return R_ac, t_ac


def cahvor_transform(cahvor_a, R_ba, t_ba):
    """
    Applies the transformation R_ba and t_ba on cahvor vector in a
    """
    c_a, a_a, h_a, v_a, o_a, r_a = cahvor_separate(cahvor_a)
    c_b = R_ba @ c_a + t_ba
    a_b = R_ba @ a_a
    h_b = R_ba @ h_a
    v_b = R_ba @ v_a
    o_b = R_ba @ o_a
    r_b = r_a
    cahvor_b = cahvor_combine(c_b, a_b, h_b, v_b, o_b, r_b)

    return cahvor_b


def get_t_s3r(site, drive, label):
    """
      Returns the translation t_s3r corresponding to the site and drive count
    """
    # find the nearest drive count before and after

    nav_t = label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR']
    tr_curr = m20_sitedrive2tr[(site, 0)]
    tr_s3 = m20_sitedrive2tr[(3, 0)]

    t_sr = nav_t + tr_curr - tr_s3

    return t_sr


def get_t_sr(site, drive):
    """
      Returns the translation t_sr corresponding to the site and drive count
    """
    t_sr = m20_sitedrive2tr[(site, drive)] - m20_sitedrive2tr[(site, 0)]

    return t_sr


def Rt_enu(R_a, t_a):
    '''
    This takes north-east-down frame and returns a east-north-up frame
    '''

    R_e = np.array([[0, 1, 0],
                    [1, 0, 0],
                    [0, 0, -1]])

    R_ae = R_a @ R_e
    t_ae = t_a @ R_e

    return R_ae, t_ae


class Frame(Enum):
    ROVER = 'rover'
    MAST = 'mast'
    ROVER_P = 'rover prime'
    CAMERA = 'camera'
    NAV = 'navigation'
    SITE = 'site'


def create_cmod(label, frame: Frame = Frame.ROVER, cmod_version=1,
                site=3, drive=0):

    # return in a-frame, where "a"

    # model version
    if cmod_version == 1:
        # just calculate the cmod from values in label
        cahvor_r = cahvor_from_pds(label['GEOMETRIC_CAMERA_MODEL'])

    if cmod_version == 2:
        raise NotImplementedError("Create CMOD for CMOD V2 is not implemented")

    R_rc, t_rc = cahvor_Rt(cahvor_r)
    R_cr = R_rc.T
    t_cr = (-1 * R_rc.T @ np.expand_dims(t_rc, axis=1)).flatten()
    cahvor_c = cahvor_transform(cahvor_r, R_cr, t_cr)
    # find transform to requested frame
    match(frame):
        case Frame.CAMERA:
            R_, t = np.eye(3), np.zeros(3)

        case Frame.MAST:
            # R_mc, t_mc = find_Rt_mc(  )
            R_rm, t_rm = Rt_rm_from_pds(label)

            T_rm = np.vstack([np.hstack([R_rm, t_rm]), np.array([0, 0, 0, 1])])
            T_rc = np.vstack([np.hstack([R_rc, t_rc]), np.array([0, 0, 0, 1])])
            T_mc = np.linalg.inv(T_rm) @ T_rc
            R_, t = T_mc[:3, :3], T_mc[:3, 3]

        case Frame.ROVER_P:
            if cmod_version == 1:
                R_, t = R_rc, t_rc

            else:
                raise NotImplementedError(
                    "Create CMOD for CMOD V2 is not implemented")
                # R_rpc, t_rpc = find_Rt_rpc( )
                # cahvor_a

        case Frame.ROVER:
            # R_rc, t_rc = find_Rt_rc( )
            R_, t = R_rc, t_rc

        case Frame.NAV:
            t_nr = label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR']
            q_nr = q_wxyz2xyzw(
                label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'])
            R_nr = R.from_quat(q_nr).as_matrix()

            R_, t = R_nr @ R_rc, R_nr @ t_rc + t_nr

        case Frame.SITE:
            # note: here we assume R_s3r = R_nr
            q_sr = q_wxyz2xyzw(
                label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'])
            R_sr = R.from_quat(q_sr).as_matrix()
            t_sr = get_t_sr(drive, site)

            R_ = R_sr @ R_rc
            t = R_sr @  t_rc + t_sr

        case Frame.ROVER_E:
            # R_rc,  t_rc  = find_Rt_rc( )
            R_, t = Rt_enu(R_rc, t_rc)
            # cahvor_a

        case Frame.NAV_E:
            t_nr = label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR']
            q_nr = q_wxyz2xyzw(
                label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'])
            R_nr = R.from_quat(q_nr).as_matrix()

            R_, t = Rt_enu(R_nr @ R_rc, R_nr @ t_rc + t_nr)

        case Frame.SITE_E:
            q_sr = q_wxyz2xyzw(
                label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'])
            R_sr = R.from_quat(q_sr).as_matrix()
            t_sr = get_t_sr(drive, site)

            R_, t = Rt_enu(R_sr @ R_rc, R_sr @  t_rc + t_sr)

        case _:
            print("Invalid site")

    return cahvor_transform(cahvor_c, R_, t)


def create_output(label, frame: Frame = Frame.ROVER, cmod_version=1):

    # model version
    if cmod_version == 1:
        # just calculate the cmod from values in label
        cahvor_r = cahvor_from_pds(label['GEOMETRIC_CAMERA_MODEL'])

    if cmod_version == 2:
        raise NotImplementedError("Create CMOD for CMOD V2 is not implemented")

    R_rc, t_rc = cahvor_Rt(cahvor_r)
    R_cr = R_rc.T
    t_cr = (-1 * R_rc.T @ np.expand_dims(t_rc, axis=1)).flatten()
    cahvor_c = cahvor_transform(cahvor_r, R_cr, t_cr)
    # print(f"cahvor_c:{cahvor_c}")
    # find transform to requested frame
    match(frame):
        case Frame.CAMERA:
            R_, t = np.eye(3), np.zeros(3)

        case Frame.MAST:
            # R_mc, t_mc = find_Rt_mc(  )
            R_rm, t_rm = Rt_rm_from_pds(label)

            T_rm = np.vstack([np.hstack([R_rm, t_rm]), np.array([0, 0, 0, 1])])
            T_rc = np.vstack([np.hstack([R_rc, t_rc]), np.array([0, 0, 0, 1])])
            T_mc = np.linalg.inv(T_rm) @ T_rc
            R_, t = T_mc[:3, :3], T_mc[:3, 3]

        case Frame.ROVER_P:
            if cmod_version == 1:
                R_, t = R_rc, t_rc

            else:
                raise NotImplementedError(
                    "Create CMOD for CMOD V2 is not implemented")
                # R_rpc, t_rpc = find_Rt_rpc( )
                # cahvor_a

        case Frame.ROVER:
            # R_rc, t_rc = find_Rt_rc( )
            R_, t = R_rc, t_rc

        case Frame.NAV:
            t_nr = label['ROVER_COORDINATE_SYSTEM']['ORIGIN_OFFSET_VECTOR']
            q_nr = q_wxyz2xyzw(
                label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'])
            R_nr = R.from_quat(q_nr).as_matrix()

            R_, t = R_nr @ R_rc, R_nr @ t_rc + t_nr

        case Frame.SITE:
            # note: here we assume R_s3r = R_nr
            q_sr = q_wxyz2xyzw(
                label['ROVER_COORDINATE_SYSTEM']['ORIGIN_ROTATION_QUATERNION'])
            R_sr = R.from_quat(q_sr).as_matrix()
            site = label['ROVER_COORDINATE_SYSTEM']['COORDINATE_SYSTEM_INDEX'][0]
            drive = label['ROVER_COORDINATE_SYSTEM']['COORDINATE_SYSTEM_INDEX'][1]
            t_sr = get_t_s3r(site, drive, label)

            R_ = R_sr @ R_rc
            t = R_sr @  t_rc + t_sr

        case _:
            print("Invalid site")

    R_enu, t_enu = Rt_enu(R_, t)
    cahvor_a = cahvor_transform(cahvor_c, R_enu, t_enu)
    xyz_ae = cahvor_xyz(cahvor_a)
    opk_ae = cahvor_opk(R_enu)
    intr_ae = cahvor_intr(cahvor_a)
    dist_ae = cahvor_dist(cahvor_a)

    return xyz_ae, opk_ae, intr_ae, dist_ae
