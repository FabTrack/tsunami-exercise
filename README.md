# Tsunami Wave Simulation

A numerical simulation of the 2D wave equation, modeling the propagation of a tsunami as it approaches the coast. Implemented in C.

This project was developed in 2020 as a learning exercise in numerical methods and computational physics, and is preserved here as a reference.

## Features 

- Simulates tsunami wave propagation using finite difference methods
- Models variable bathymetry to approximate coastal shallowing
- Visual output (gif) for viewing wave evolution over time

## Build

Requires a C compiler (e.g., `gcc`).

```bash
gcc tsunami.c -o tsunami_sim -lm
