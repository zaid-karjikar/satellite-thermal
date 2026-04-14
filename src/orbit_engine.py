from math import cos, sin, sqrt, pi, radians
from src.constants import EARTH_RADIUS_KM, EARTH_MU_KM3_S2

def _sun_unit_vector(raan_rad, declination_rad):
    """
    Convert spherical sky coordinates to a Cartesian unit vector.
    """
    cos_d = cos(declination_rad)
    return (cos_d * cos(raan_rad),
            cos_d * sin(raan_rad),
            sin(declination_rad))

def _position_eci(radius_m, 
                  inclination_rad, 
                  raan_rad, 
                  true_anomaly_0_rad, 
                  mean_motion_rad_s, 
                  t_s):
    """
    Returns the ECI position vector (x, y, z) in meters for a circular orbit at time t_s.
    """
    true_anomaly_rad = true_anomaly_0_rad + mean_motion_rad_s * t_s

    # step 1 - position in orbital plane
    x_p = radius_m * cos(true_anomaly_rad)
    y_p = radius_m * sin(true_anomaly_rad)
    z_p = 0.0

    # step 2 - tilt by inclination
    cos_i, sin_i = cos(inclination_rad), sin(inclination_rad)
    x_n = x_p
    y_n = y_p * cos_i - z_p * sin_i
    z_n = y_p * sin_i + z_p * cos_i

    # step 3 - rotate by RAAN
    cos_o, sin_o = cos(raan_rad), sin(raan_rad)
    x = x_n * cos_o - y_n * sin_o
    y = x_n * sin_o + y_n * cos_o
    z = z_n

    return (x, y, z)

def _in_cylindrical_eclipse(position_m, sun_unit_vector, earth_radius_m):
    """
    Returns True if the satellite position is inside Earth's cylindrical shadow.
    """
    px, py, pz = position_m
    sx, sy, sz = sun_unit_vector

    projection_onto_sun_axis = px*sx + py*sy + pz*sz
    if projection_onto_sun_axis >= 0.0:
        return False          # satellite is on the sunlit side of Earth

    position_magnitude_sq    = px*px + py*py + pz*pz
    perpendicular_distance_sq = position_magnitude_sq - projection_onto_sun_axis * projection_onto_sun_axis
    earth_radius_sq = earth_radius_m * earth_radius_m

    return perpendicular_distance_sq < earth_radius_sq