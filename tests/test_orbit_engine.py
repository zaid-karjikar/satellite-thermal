import pytest
from math import sqrt, isclose, pi
from src import orbit_engine

# -- helpers ---------------------------------------------------------------
def magnitude(v):
    return sqrt(sum(c ** 2 for c in v))

def allclose(v, expected, tol=1e-9):
    return all(isclose(a, b, abs_tol=tol) for a, b in zip(v, expected))

# -- tests sun unit vector -------------------------------------------------
def test_vernal_equinox():
    """Raan=0, declination=0 points along +x (toward vernal equinox)."""
    assert allclose(orbit_engine._sun_unit_vector(0, 0), (1.0, 0.0, 0.0))

def test_north_pole():
    """Declination=+π/2 points along +z regardless of raan."""
    assert allclose(orbit_engine._sun_unit_vector(0, pi /2), (0.0, 0.0, 1.0))

def test_south_pole():
    """Declination=-π/2 points along -z."""
    assert allclose(orbit_engine._sun_unit_vector(0, -pi / 2), (0.0, 0.0, -1.0))

def test_90_deg_raan():
    """Raan=π/2, declination=0 points along +y."""
    assert allclose(orbit_engine._sun_unit_vector(pi / 2, 0), (0.0, 1.0, 0.0))

def test_180_deg_raan():
    """Raan=π, declination=0 points along -x."""
    assert allclose(orbit_engine._sun_unit_vector(pi, 0), (-1.0, 0.0, 0.0))

def test_unit_magnitude_equator():
    """Magnitude is 1 for a point on the equator."""
    assert isclose(magnitude(orbit_engine._sun_unit_vector(1.234, 0)), 1.0)

def test_unit_magnitude_arbitrary():
    """Magnitude is 1 for an arbitrary raan/declination."""
    assert isclose(magnitude(orbit_engine._sun_unit_vector(1.1, 0.5)), 1.0)

def test_unit_magnitude_negative_declination():
    """Magnitude is 1 for negative declination."""
    assert isclose(magnitude(orbit_engine._sun_unit_vector(3.5, -0.7)), 1.0)

def test_symmetry_declination():
    """Negating declination negates only the z component."""
    v_pos = orbit_engine._sun_unit_vector(1.0, 0.4)
    v_neg = orbit_engine._sun_unit_vector(1.0, -0.4)
    assert isclose(v_pos[0], v_neg[0])
    assert isclose(v_pos[1], v_neg[1])
    assert isclose(v_pos[2], -v_neg[2])

def test_full_raan_rotation():
    """Sweeping raan by 2π returns to the starting point."""
    v1 = orbit_engine._sun_unit_vector(0.6, 0.3)
    v2 = orbit_engine._sun_unit_vector(0.6 + 2 * pi, 0.3)
    assert allclose(v1, v2)
