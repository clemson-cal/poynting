# Poynting

Performant Python functions for calculating electromagnetic radiation from moving charges.

## About

This project is authored by Jonathan Zrake as part of a graduate-level electromagnetic theory course at Clemson University. It can be used as a pedagogical tool to explore radiation from relativistic moving charges, and might have research applications in astrophysics and plasma physics.

## Overview

Poynting provides tools for computing the electromagnetic field and energy flow (Poynting vector) from charged particles following prescribed trajectories. It combines C++ for performance-critical calculations with a convenient Python interface.

## Features

- Calculation of Liénard-Wiechert potentials for moving charges
- Computation of electromagnetic fields via the Faraday tensor
- Evaluation of Poynting vectors (energy flux) at arbitrary points
- Time-averaging of radiation patterns
- High-order finite differencing for improved accuracy
- Relativistic calculations using proper four-vector formalism

## Example Usage

```python
from poynting import ParticleTrajectory, poynting_vector, poynting_vector_time_avg

# Create a charged particle oscillating along the z-axis
# with angular frequency omega and amplitude ell (the speed
# of light is normalized to c=1 throughout the project)
particle = ParticleTrajectory(omega=0.9, ell=1.0)

# Calculate the instantaneous Poynting vector at a specific point and time
S = poynting_vector(particle, t=0.5, x=10.0, y=0.0, z=0.0)
print(f"Energy flux: {S.t}, direction: ({S.x}, {S.y}, {S.z})")

# Calculate the time-averaged energy flux in the radial direction
avg_flux = time_averaged_energy_flux(particle, x=10.0, y=0.0, z=5.0)
print(f"Time-averaged radial energy flux: {avg_flux}")
```

## Installation

```bash
pip install poynting
```

## Technical Details

Poynting implements:
- Four-vector formalism for relativistic calculations
- Retarded time calculations for proper causality
- High-order finite difference methods (2nd, 4th, and 5th order available)
- Numerical root-finding for retarded-time equations
- Time integration for period averaging

## License

MIT License
