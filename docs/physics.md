# Physics
## 1. Overview
This document describes the equations, constants, and assumptions behind `sat-thermal`.
It is aimed at thermal engineers who want to understand what the library is actually computing before they trust its outputs.
Every assumption made by the model is stated here, with a note on how violating it affects the result.

## 2. Constants
| Symbol | Value | Units | Source |
|---|---|---|---|
| $R_\oplus$ | 6378.137 | km | WGS-84 equatorial radius |
| $\mu_\oplus$ | 398600.4418 | km³/s² | Standard gravitational parameter |

All values are defined in `src/constants.py`.

## 3. Orbit
### 3.1 Sun direction

The Sun is treated as infinitely far away, so its direction is the same everywhere in Earth's vicinity and is represented as a unit vector $\hat{s}$ in ECI.

#### Frame definition

The ECI frame is right-handed:

| Axis | Points toward |
|------|--------------|
| $x$ | Vernal equinox $\gamma$ |
| $y$ | 90° east on the celestial equator |
| $z$ | Celestial north pole |

This frame is shared by all vectors in this document. It is defined here because the Sun direction is the first vector that requires it.

#### Conversion from right ascension and declination

The Sun's position on the celestial sphere is specified by right ascension $\alpha \in [0, 2\pi)$ and declination $\delta \in [-\pi/2, +\pi/2]$, both in radians, and converted to a Cartesian unit vector via:

$$\hat{s} = (\cos\delta\cos\alpha,\ \cos\delta\sin\alpha,\ \sin\delta)$$

The derivation is a two-step projection.

**Step 1 — project onto the equatorial plane.**
A unit vector tilted $\delta$ above the equator casts a shadow of length $\cos\delta$ onto the equatorial plane. That shadow points along angle $\alpha$, so it decomposes as:

$$x = \cos\delta\cos\alpha, \qquad y = \cos\delta\sin\alpha$$

**Step 2 — vertical component.**
The remaining component along $z$ is the sine of the elevation angle:

$$z = \sin\delta$$

The result is a unit vector by construction for any valid $(\alpha, \delta)$:

$$|\hat{s}|^2 = \cos^2\delta\,(\cos^2\alpha + \sin^2\alpha) + \sin^2\delta = \cos^2\delta + \sin^2\delta = 1$$

No normalisation step is required.

#### Fixed-sun assumption and validity

$\hat{s}$ is computed once at initialisation and held constant for the simulation duration.
The Sun moves approximately $360° / 365.25 \approx 0.986°$ per day relative to the stars, so a 1-day run accumulates less than 1° of directional error and a 10-day run less than 10°.
For runs longer than a few days, $\hat{s}$ should be updated at each timestep using an ephemeris.

The default `raan=0, declination=0` places the Sun on the $+x$ axis (at the vernal equinox), corresponding approximately to the March equinox. This is a convenient reference point for tests but is only physically accurate for simulations starting around that date.