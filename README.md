# SS-Div-Conforming-B-splines-Incompressible-NS
Code references for the numerical examples in the paper: 

> [*Skeleton-stabilized divergence-conforming B-spline discretizations for incompressible flow problems of high Reynolds number, Guoxiang Grayson Tong, David Kamensky and John A. Evans, Computers & Fluids, 2022*](https://www-sciencedirect-com.proxy.library.nd.edu/science/article/pii/S0045793022002602)

Usage of the code requires the [tIGAr](https://github.com/david-kamensky/tIGAr) package, constructed via the [FEniCS](https://fenicsproject.org/) background. For more detailed information, please contact: guto1826@colorado.edu

### MFS.py
Code for the steady Navier-Stokes solver, using method of manufactured solutions to test convergence

### Lid2d.py
Code for the benchmark 2d Lid-driven cavity problem

### TG3D.py
Code for the benchmark 3D Taylor-Green vortex problem
