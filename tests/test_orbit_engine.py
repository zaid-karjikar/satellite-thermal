import pytest
from math import sqrt, isclose, pi, radians, sin
from src import orbit_engine
from src.constants import EARTH_RADIUS_M

# -- helpers ---------------------------------------------------------------------------------------

def magnitude(v):
    return sqrt(sum(c ** 2 for c in v))

def allclose(v, expected, abs_tol=1e-6, rel_tol=1e-9):
    return all(isclose(a, b, abs_tol=abs_tol, rel_tol=rel_tol) for a, b in zip(v, expected))

R = 6_771_000.0                       # ISS-altitude circular orbit radius, meters
N = 2 * pi / 5560.0                   # approximate mean motion, rad/s

# -- tests _sun_unit_vector() ----------------------------------------------------------------------

class TestSunUnitVector:

    # cardinal directions --------------------------------------------------------------------------

    def test_vernal_equinox(self):
        """Raan=0, declination=0 points along +x (toward vernal equinox)."""
        assert allclose(orbit_engine._sun_unit_vector(0, 0), (1.0, 0.0, 0.0))

    def test_north_pole(self):
        """Declination=+π/2 points along +z regardless of raan."""
        assert allclose(orbit_engine._sun_unit_vector(0, pi /2), (0.0, 0.0, 1.0))

    def test_south_pole(self):
        """Declination=-π/2 points along -z."""
        assert allclose(orbit_engine._sun_unit_vector(0, -pi / 2), (0.0, 0.0, -1.0))

    def test_90_deg_raan(self):
        """Raan=π/2, declination=0 points along +y."""
        assert allclose(orbit_engine._sun_unit_vector(pi / 2, 0), (0.0, 1.0, 0.0))

    def test_180_deg_raan(self):
        """Raan=π, declination=0 points along -x."""
        assert allclose(orbit_engine._sun_unit_vector(pi, 0), (-1.0, 0.0, 0.0))

    # magnitude ------------------------------------------------------------------------------------
    
    def test_unit_magnitude_equator(self):
        """Magnitude is 1 for a point on the equator."""
        assert isclose(magnitude(orbit_engine._sun_unit_vector(1.234, 0)), 1.0)

    def test_unit_magnitude_arbitrary(self):
        """Magnitude is 1 for an arbitrary raan/declination."""
        assert isclose(magnitude(orbit_engine._sun_unit_vector(1.1, 0.5)), 1.0)

    def test_unit_magnitude_negative_declination(self):
        """Magnitude is 1 for negative declination."""
        assert isclose(magnitude(orbit_engine._sun_unit_vector(3.5, -0.7)), 1.0)

    # symmetry and periodicity ---------------------------------------------------------------------

    def test_symmetry_declination(self):
        """Negating declination negates only the z component."""
        v_pos = orbit_engine._sun_unit_vector(1.0, 0.4)
        v_neg = orbit_engine._sun_unit_vector(1.0, -0.4)
        assert isclose(v_pos[0], v_neg[0])
        assert isclose(v_pos[1], v_neg[1])
        assert isclose(v_pos[2], -v_neg[2])

    def test_full_raan_rotation(self):
        """Sweeping raan by 2π returns to the starting point."""
        v1 = orbit_engine._sun_unit_vector(0.6, 0.3)
        v2 = orbit_engine._sun_unit_vector(0.6 + 2 * pi, 0.3)
        assert allclose(v1, v2)

# -- tests _position_eci()--------------------------------------------------------------------------

class TestPositionECI:

    # radius conservation --------------------------------------------------------------------------

    def test_radius_constant_at_t0(self):
        """Radius equals input radius at t=0."""
        x, y, z = orbit_engine._position_eci(R, 0, 0, 0, N, 0)
        assert isclose(magnitude((x, y, z)), R)

    def test_radius_conserved_over_time(self):
        """Radius stays constant across multiple timesteps (rotation preserves length.)"""
        for t in [0, 500, 1000, 2000, 5000]:
            x, y, z = orbit_engine._position_eci(R, 0, 0, 0, N, t)
            assert isclose(magnitude((x, y, z)), R, rel_tol=1e-9)

    # zero inclination, zero RAAN ------------------------------------------------------------------
    
    def test_equatorial_t0_points_along_x(self):
        """Zero inclination, zero RAAN, true anomaly = 0 -> position along +x."""
        x, y, z = orbit_engine._position_eci(R, 0, 0, 0, N, 0)
        assert allclose((x, y, z), (R, 0.0, 0.0))

    def test_equitorial_true_anomaly_90_points_along_y(self):
        """Zero inclination, zero RAAN, true anomaly = π/2 -> position along +y."""
        x, y, z = orbit_engine._position_eci(R, 0, 0, pi / 2, N, 0)
        assert allclose((x, y, z), (0.0, R, 0.0))

    def test_equatorial_stays_in_xy_plane(self):
        """Zero inclination orbit never leaves the equitorial plane (z=0)."""
        for t in [0, 500, 1000, 2000]:
            x, y, z = orbit_engine._position_eci(R, 0, 0, 0, N, t)
            assert isclose(z, 0.0, abs_tol=1e-6)

    # inclination ----------------------------------------------------------------------------------

    def test_polar_orbit_reaches_north_pole(self):
        """inclination = 90°, true anomaly = π/2 -> satellite directly over north pole (+z)."""
        x, y, z = orbit_engine._position_eci(R, pi / 2, 0, pi / 2, N, 0)
        assert allclose((x, y, z), (0.0, 0.0, R))

    def test_polar_orbit_reaches_south_pole(self):
        """inclination = 90°, true anomaly = 3π/2 → satellite directly over south pole (-z)."""
        x, y, z = orbit_engine._position_eci(R, pi / 2, 0, 3 * pi / 2, N, 0)
        assert allclose((x, y, z), (0.0, 0.0, -R), rel_tol=1e-9)

    def test_inclination_max_z(self):
        """Maximum z reached equals radius * sin(inclination)."""
        inc = radians(51.6)
        x, y, z = orbit_engine._position_eci(R, inc, 0, pi / 2, N, 0)
        assert isclose(z, R * sin(inc), rel_tol=1e-9)

    # RAAN -----------------------------------------------------------------------------------------
    
    def test_raan_90_rotates_ascending_node(self):
        """RAAN = π/2 rotates the orbit 90° around z - ascending node moves to +y axis."""
        x, y, z = orbit_engine._position_eci(R, 0, pi / 2, 0, N, 0)
        assert allclose((x, y, z), (0.0, R, 0.0))

    def test_raan_does_not_affect_z(self):
        """RAAN rotation is around z, so z component is unchanged."""
        z_ref = orbit_engine._position_eci(R, radians(51.6), 0,       pi / 2, N, 0)[2]
        z_rot = orbit_engine._position_eci(R, radians(51.6), pi / 3,  pi / 2, N, 0)[2]
        assert isclose(z_ref, z_rot, rel_tol=1e-9)

# -- tests _position_eci()--------------------------------------------------------------------------

class TestInCylindricalEclipse:

    # sun-side early exit --------------------------------------------------------------------------

    def test_directly_toward_sun_is_sunlit(self):
        """Satellite directly between Earth and Sun is not in eclipse."""
        assert not orbit_engine._in_cylindrical_eclipse((8e6, 0, 0), (1, 0, 0))

    def test_on_sun_axis_positive_is_sunlit(self):
        """Any positive projection onto sun axis returns False immediately."""
        assert not orbit_engine._in_cylindrical_eclipse((4e6, 4e6, 0), (1, 0, 0))

    def test_at_origin_projection_zero_is_sunlit(self):
        """Zero projection (perpendicular to sun) is treated as sunlit."""
        assert not orbit_engine._in_cylindrical_eclipse((0, 8e6, 0), (1, 0, 0))

    # on shadow axis -------------------------------------------------------------------------------

    def test_directly_behind_earth_is_eclipse(self):
        """Satellite directly behind Earth on shadow axis is in eclipse."""
        assert orbit_engine._in_cylindrical_eclipse((-8e6, 0, 0), (1, 0, 0))

    def test_far_behind_earth_is_eclipse(self):
        """Satellite far behind Earth but still on axis is in eclipse."""
        assert orbit_engine._in_cylindrical_eclipse((-40e6, 0, 0), (1, 0, 0))

    # cylinder boundary ----------------------------------------------------------------------------

    def test_beside_earth_outside_cylinder_is_sunlit(self):
        """Satellite beside Earth with perpendicular distance > R is sunlit."""
        assert not orbit_engine._in_cylindrical_eclipse((0, 8e6, 0), (1, 0, 0))

    def test_just_outside_cylinder_is_sunlit(self):
        """Satellite just outside the shadow cylinder is sunlit."""
        assert not orbit_engine._in_cylindrical_eclipse((-8e6, EARTH_RADIUS_M * 1.01, 0), (1, 0, 0))

    def test_just_inside_cylinder_is_eclipse(self):
        """Satellite just inside the shadow cylinder is in eclipse."""
        assert orbit_engine._in_cylindrical_eclipse((-8e6, EARTH_RADIUS_M * 0.99, 0), (1, 0, 0))

    # sun direction independence -------------------------------------------------------------------

    def test_eclipse_with_sun_along_y_axis(self):
        """Eclipse check works regardless of sun direction — sun along +y."""
        assert orbit_engine._in_cylindrical_eclipse((0, -8e6, 0), (0, 1, 0))

    def test_eclipse_with_sun_along_z_axis(self):
        """Eclipse check works regardless of sun direction — sun along +z."""
        assert orbit_engine._in_cylindrical_eclipse((0, 0, -8e6), (0, 0, 1))

    def test_sunlit_with_sun_along_y_axis(self):
        """Sunlit check works regardless of sun direction — sun along +y."""
        assert not orbit_engine._in_cylindrical_eclipse((0, 8e6, 0), (0, 1, 0))