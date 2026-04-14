# constants.py

# Convention: the simulation engine works in SI units (meters, seconds, radians).
# Constants are provided in both km (for readability and matching published values)
# and m (for use by the engine).
# Always import the _M variant in src/orbit_engine.py

# Earth
EARTH_RADIUS_KM: float = 6378.137                  # WGS-84 equatorial radius
EARTH_RADIUS_M: float = EARTH_RADIUS_KM * 1000

EARTH_MU_KM3_S2: float = 398600.4418               # Standard gravitational parameter
EARTH_MU_M3_S2: float = EARTH_MU_KM3_S2 * 1e9
