# Lid-driven cavity problem

## Description
Consider a fluid inside a square cavity of dimension L bounded by four walls. All the walls except the top wall are fixed in space and time. The top wall is given a constant velocity 𝑢𝑇in the u direction, and we 
want to see the development of the u and v velocity profiles of the fluid in the whole domain, as it reaches steady state.
![lid_driven_cavity](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/382fa46b-ac14-42aa-8618-fbe46c894d83)

## Problem statement
The entire problem statement can be found [here](https://github.com/nilot-pal/Lid-driven-cavity/blob/main/Problem_statement.pdf).

## Governing equations
The incompressible NS equations in primitive variables are given by the equation:
![image](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/6ca38f37-3876-4d27-a556-e275df0f5f29)

If we apply **2nd order central difference scheme** for space discretization and **2nd order Adams Bashforth** for time discretization, the above equation can be represented as:
![image](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/433b7301-b720-4bbf-bf70-8c844fd90776)

## Solution method
The fractional step, or time-splitting, method solves the unsteady Navier-Stokes equations in a segregated manner. At each time step, an incomplete form of momentum equations is integrated to obtain an approximate velocity field, which is, in general, not divergence-free, then the velocity field is projected into the divergence-free field without changing vorticity. This projection step is achieved by solving the Poisson equation for pressure.

## Results
Benchmarking is done against [Ghia and Ghia](https://github.com/nilot-pal/Lid-driven-cavity/blob/main/ghia1982.pdf). The project technical report can be found [here](https://github.com/nilot-pal/Lid-driven-cavity/blob/main/Technical_report.pdf).
## 1. u velocity along vertical line through geometric centre of cavity**
![image](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/f7b1d83b-5d50-4fd2-b818-5e9318253dd8)
## 1b. v velocity along horizontal line through geometric centre of cavity**
![image](https://github.com/nilot-pal/Lid-driven-cavity/assets/72824334/867552d2-42df-46fe-9c01-5fb507c7a974)




