# Satellite Thermal

Open-source thermal analysis for small satellites.

`sat-thermal` is a Python library for modelling the thermal behaviour of CubeSats and other small spacecraft in low Earth Orbit.
You define a satellite, list the components inside it, run a simulation, and get a temperature history for the spacecraft along with any component limits it violated.

## What it models
`sat-thermal` treats the spacecraft as a **single isothermal node**. The whole satellite has one temperature at any instant, and that temperature evolves under the balance of heat coming in and heat going out.
This is the standard first-order approach used in early CubeSat thermal trade studies.

**Heat in**
- Direct solar flux (Sun treated as infinitely far away)
- Earth albedo (sunlight reflected off Earth's day side)
- Earth infrared (longwave emission from Earth's surface)
- Internal power dissipation summed across all components

**Heat out**
- Radiative loss to deep space via Stefan-Boltzmann law

**Orbit**
- Circular orbit in ECI
- Cylindrical (hard-shadow) eclipse model
- Altitude-dependent view factor to Earth for albedo and IR

**Integration**
- Forward Euler, fixed time step

## What it does *not* model
The current version does not handle:

- **Internal temperature gradients:** The whole spacecraft is one node. If you need to know the temperature difference between two components, that requires a multi-node model with conduction links, which is on the roadmap but not built.
- **Elliptical orbits:** Circular orbits only. No Kepler equations solver yet.
- **Penumbra during eclipse transitions:** Hard shadow only.
- **Attitude dynamics:** The model assumes the full external area is available for absorption when in sunlight. Directional pointing and cosine losses on individual faces are not modelled.
- **Duty-cycle power dissipation:** Internal heat loads are constant throughout the simulation. Mode switching is on the roadmap.
- **J2, drag, orbit decay, or any orbit perturbation**
- **Adaptive time stepping:** Fixed dt, user's choice.

## License
MIT.

## Citing
If you use 'sat-thermal' in academic work, please cite the repository. A `CITATION.cff` will be added before the first tagged release.
