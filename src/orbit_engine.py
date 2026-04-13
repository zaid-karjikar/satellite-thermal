from math import cos, sin

def _sun_unit_vector(raan_rad, declination_rad):
    """
    Convert spherical sky coordinates to a Cartesian unit vector.
    """
    cos_d = cos(declination_rad)
    return (cos_d * cos(raan_rad),
            cos_d * sin(raan_rad),
            sin(declination_rad))
