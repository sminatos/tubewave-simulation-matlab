# tubewave-simulation-matlab
Matlab scripts for fast and accurate forward modeling tube waves in vertical seismic profiling using a propagator matrix method.

## Overview
This repository collects scripts that calculate hydrophone VSP (vertical seismic profiling) response. The following situations are particularly considered:

- An open borehole
- A borehole embedded in a horizontally layered (1-D) poroelastic medium
- A borehole with a step-like radius change
- A plane P wave propagates vertically downward

In this case, the VSP response contains elastic waves and Stoneley waves (_tube waves_). Tube wave is a guided wave propagating along the borehole. The scripts are based on approximations to complete Biot dynamic poroelasticity equations.

The approach first solves 1-D potential fields for a plane P wave using a propagator matrix method. Then it solves potential fields for the borehole fluid using another propagator matrix method. The final waveforms include both elastic waves and tube waves.


- [Modeling_VSP_tubewave_Elastic_Layer.m](/Modeling_VSP_tubewave_Elastic_Layer.m): A water-filled borehole located at two elastic half-spaces. Tube wave is generated at the elastic impedance boundary due to an incident P wave. Analytical solutions are also provided.

- [Modeling_VSP_tubewave_Porous_Layer.m](/Modeling_VSP_tubewave_Porous_Layer.m): A water-filled borehole is located at a three-layer medium. In the medium, a poroelastic layer is sandwiched between two elastic layers. Deformation of the poroelastic layer produces the fluid flow. Then it generates a tube wave in the borehole fluid. Analytical solutions are also provided.

- [Modeling_VSP_tubewave_Caliper.m](/Modeling_VSP_tubewave_Caliper.m): A water-filled borehole is located at a homogeneous elastic half-space. Tube wave is generated at the radius change due to an interaction with a P wave. Analytical solutions are also provided.
