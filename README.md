# tubewave-simulation-matlab
Matlab scripts for fast and accurate forward modeling tube waves in vertical seismic profiling using a propagator matrix method ([Minato et al., 2021](#references)).

## Overview
This repository collects scripts that calculate hydrophone VSP (vertical seismic profiling) response. The following situations are particularly considered:

- An open borehole
- A borehole embedded in a horizontally layered (1-D) poroelastic medium
- A borehole with a step-like radius change
- A plane P wave propagates vertically downward

In this case, the VSP response contains elastic waves and Stoneley waves (_tube waves_). Tube wave is a guided wave propagating along the borehole. The theory is based on low-frequency approximations to Biot dynamic poroelasticity. The approach first solves 1-D potential fields for a plane P wave using a propagator matrix method. Then it solves potential fields for the borehole fluid using another propagator matrix method. The final output is *a complete* response including both elastic waves and tube waves. The outputs have been verified with a finite-difference modeling of Biot poroelasticity. Please see [Minato et al.(2021)](#references) for more details.


- [Modeling_VSP_tubewave_Elastic_Layer.m](/Modeling_VSP_tubewave_Elastic_Layer.m): A water-filled borehole located at two elastic half-spaces. Tube wave is generated at the elastic impedance boundary due to an incident P wave. Analytical solutions are also provided.

- [Modeling_VSP_tubewave_Porous_Layer.m](/Modeling_VSP_tubewave_Porous_Layer.m): A water-filled borehole is located at a three-layer medium. In the medium, a poroelastic layer is sandwiched between two elastic layers. Deformation of the poroelastic layer produces the fluid flow. Then it generates a tube wave in the borehole fluid. Analytical solutions are also provided.

- [Modeling_VSP_tubewave_Caliper.m](/Modeling_VSP_tubewave_Caliper.m): A water-filled borehole is located at a homogeneous elastic half-space. Tube wave is generated at the radius change due to an interaction with a P wave. Analytical solutions are also provided.

## References
- Minato et al. (2021), arXiv:2112.03410 [physics.geo-ph] (available at http://arxiv.org/abs/2112.03410).
