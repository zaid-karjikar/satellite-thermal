from math import cos, sin, sqrt, pi, radians
from src.constants import EARTH_MU_M3_S2, EARTH_RADIUS_M

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

def _in_cylindrical_eclipse(position_m, sun_unit_vector):
    """
    Returns True if the satellite position is inside Earth's cylindrical shadow.
    """
    px, py, pz = position_m
    sx, sy, sz = sun_unit_vector

    projection_onto_sun_axis = px*sx + py*sy + pz*sz
    if projection_onto_sun_axis >= 0.0:
        return False          # satellite is not behind Earth (sunlit or on the terminator plane)

    position_magnitude_sq    = px*px + py*py + pz*pz
    perpendicular_distance_sq = position_magnitude_sq - projection_onto_sun_axis * projection_onto_sun_axis
    earth_radius_sq = EARTH_RADIUS_M * EARTH_RADIUS_M

    return perpendicular_distance_sq < earth_radius_sq

class Orbit:
    """
    Circular orbit in ECI with a fixed sun direction.

    Accepts altitude and angles in kilometres and degrees.
    All internal state is in SI units (metres, seconds, radians).
    """
    def __init__(self, 
                 altitude_km, 
                 inclination_deg, 
                 raan_deg=0.0, 
                 true_anomaly_deg=0.0, 
                 sun_ra_deg=0.0, 
                 sun_dec_deg=0.0):
        if altitude_km <= 0:
            raise ValueError(f"altitude_km must be positive, got {altitude_km}")
        
        self.orbital_radius_m   = EARTH_RADIUS_M + altitude_km * 1000
        self.inclination_rad    = radians(inclination_deg)
        self.raan_rad           = radians(raan_deg)
        self.true_anomaly_0_rad = radians(true_anomaly_deg)

        self.orbital_period_s  = 2 * pi * sqrt(self.orbital_radius_m**3 / EARTH_MU_M3_S2)
        self.mean_motion_rad_s = 2 * pi / self.orbital_period_s

        self.sun_vector = _sun_unit_vector(radians(sun_ra_deg), radians(sun_dec_deg))