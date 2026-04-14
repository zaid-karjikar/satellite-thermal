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

### 3.2 Satellite position in ECI

The satellite is assumed to follow a **circular orbit** of constant radius $r$.
Position is computed from five orbital elements plus time.

#### Orbital elements

| Symbol | Parameter | Unit |
|--------|-----------|------|
| $r$ | Orbital radius (Earth centre to satellite) | m |
| $i$ | Inclination — tilt of orbital plane from equator | rad |
| $\Omega$ | RAAN — where the orbit crosses the equator | rad |
| $\nu_0$ | True anomaly at $t = 0$ | rad |
| $n$ | Mean motion — angular speed around orbit | rad s⁻¹ |

#### Derivation

The transformation from orbital elements to ECI is a three-step rotation.

**Step 1 — position in the orbital plane.**
For a circular orbit the satellite sweeps angle $\nu$ at constant rate $n$:

$$\nu(t) = \nu_0 + n \cdot t$$

In its own plane the orbit is flat 2D circular motion:

$$\vec{r}_p = (r\cos\nu,\ r\sin\nu,\ 0)$$

**Step 2 — tilt by inclination.**
The orbital plane is tilted relative to the equator by inclination $i$.
This is a rotation around the $x$-axis:

$$R_x(i) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos i & -\sin i \\ 0 & \sin i & \cos i \end{pmatrix}$$

Since $z_p = 0$ this simplifies to:

$$\vec{r}_n = (r\cos\nu,\ r\sin\nu\cos i,\ r\sin\nu\sin i)$$

For an orbit with $i = 51.6°$ (ISS), the satellite reaches latitudes of $\pm 51.6°$.

**Step 3 — rotate by RAAN.**
RAAN $\Omega$ rotates the orbital plane around the $z$-axis, setting where the orbit
crosses the equator relative to the vernal equinox.
This is a rotation around the $z$-axis:

$$R_z(\Omega) = \begin{pmatrix} \cos\Omega & -\sin\Omega & 0 \\ \sin\Omega & \cos\Omega & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

$z$ is unchanged because the rotation axis is $z$.

**Combined transformation.**
The full rotation from orbital plane to ECI is:

$$\vec{r}_{ECI} = R_z(\Omega) \cdot R_x(i) \cdot \vec{r}_p$$

Expanding:

$$\vec{r}_{ECI} = \begin{pmatrix}
r\cos\nu\cos\Omega - r\sin\nu\cos i\sin\Omega \\
r\cos\nu\sin\Omega + r\sin\nu\cos i\cos\Omega \\
r\sin\nu\sin i
\end{pmatrix}$$

#### Magnitude check

Rotations preserve vector length, so the orbital radius must remain constant for all $t$:

$$|\vec{r}_{ECI}|^2 = x^2 + y^2 + z^2 = r^2(\cos^2\nu + \sin^2\nu) = r^2$$

$$|\vec{r}_{ECI}| = r \quad \text{for all } t$$

#### Assumptions and limitations

| Assumption | Effect if violated |
|------------|-------------------|
| Circular orbit ($e = 0$) | Elliptical orbits require true anomaly computed from Kepler's equation, not $\nu_0 + nt$ |
| Two-body gravity | Ignores $J_2$ oblateness, drag, and third-body perturbations; orbital plane drifts over days to weeks |
| Inertial frame | ECI is not truly inertial over long runs due to precession; negligible for thermal simulations |