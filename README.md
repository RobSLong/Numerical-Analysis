# Numerical-Analysis

This repository contains projects focused around Numerical Analysis and each directory is intended to be a self-contained project with adequate enough documentation to be used as educational content. My academic background and research experience is in fluid dynamics (specifically convection - buoyancy driven flows relevant to all aspects of natural and industrial fluid flows) and most of the problems here are relevant to principles in fluid dynamics.

Each directory herein contains a self-contained project with the naming convention \<**application-method**\> e.g. DiffusionEq-FiniteDifferences is a project which describes how to solve the diffusion equation using the finite difference numerical scheme. The \<readme\> in each directory should provide enough information to understand the problem and provide a foundation for similar projects to be undertaken!

## Projects Highlights
 ### Advection
This project outlines how to discretize and numerically solve the linear advection equation in 1D using both explicit and implicit finite difference formulations. We present an analysis of how the error scales with the computing time and suggest a method for **optimising** the numerical scheme to most efficiently reach a specified level of accuracy. |
 <p align="center">
<img src = "https://github.com/RobSLong/FD_Advection_comparison/blob/main/advection.gif" width="350" />
</p>
       
### Diffusion 
This project outlines how to discretize and numerically solve the diffusion equation in both 1D and 2D with different boundary conditions (Dirichlet and Neumann). We investigate how dimensionality affects the diffusive process and present different methods of visualising the numerical solutions.

### Thermal convection
This project uses Dedalus (a spectral solver written in Python) to formulate and numerically solve the Navier Stokes equations to simulate 2D convection. The emphasis is on how to verify and validate numerical models and we do this through comparison with published solutions, comparing with linear theory and finally internal consistency checks.

