# Fortran_Numerical_Modelling

This is the homework of ETH course.
T: temperature
x,y: direction


## day4 ##

2D diffusion with initialization of spike and random. 

![day4_1](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

## day5 ##

2D advection-diffusion with a fixed flow field and zero flux at the boundary x and 0 and 1 at the boundary y. 

![day5_1](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=-v_x\frac{\partial%20T}{\partial%20x}%20-%20v_y\frac{\partial%20T}{\partial%20y}%20+%20\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

![day5_2](https://latex.codecogs.com/svg.latex?(v_x,%20v_y)%20=%20(\frac{\partial%20\psi}{\partial%20y},%20-\frac{\partial%20\psi}{\partial%20x})) 

![day5_3](https://latex.codecogs.com/svg.latex?\psi=Bsin(\frac{\pi%20x}{x_{max}})sin(\frac{\pi%20y}{y_{max}})) 

loop:
ψ -> v -> T

## day6 ##

Iterative Poisson solver.

![day6_1](https://latex.codecogs.com/svg.latex?\frac{\partial^2%20u}{\partial%20x^2}%20+%20\frac{\partial^2%20u}{\partial%20y^2}=f)

## day7 ##

Use Poisson solver for two times for 2D convection-diffusion along time.
With initialization of spike or random, then calculate Poisson-eqs every time step. 
Boundary condition same as in day5.

![day7_1](https://latex.codecogs.com/svg.latex?-\frac{\partial%20p}{\partial%20x}+\nabla^2v_x%20=%200)

![day7_2](https://latex.codecogs.com/svg.latex?-\frac{\partial%20p}{\partial%20y}+\nabla^2v_y%20=%20-RaT)

combine and cancel p related term:

![day7_3](https://latex.codecogs.com/svg.latex?\nabla^2(\frac{\partial%20v_x}{\partial%20y}%20-%20\frac{\partial%20v_y}{\partial%20x})=%20Ra\frac{\partial%20T}{\partial%20x})

Use streamfunction：

![day7_4](https://latex.codecogs.com/svg.latex?\nabla^2(\frac{\partial^2%20\psi}{\partial%20y^2}%20+%20\frac{\partial^2%20\psi}{\partial%20x^2})=%20Ra\frac{\partial%20T}{\partial%20x})

Then: 

![day7_5](https://latex.codecogs.com/svg.latex?\nabla^2\psi=%20w)

![day7_6](https://latex.codecogs.com/svg.latex?\nabla^2%20w%20=%20Ra\frac{\partial%20T}{\partial%20x})

![day7_7](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=-v_x\frac{\partial%20T}{\partial%20x}%20-%20v_y\frac{\partial%20T}{\partial%20y}%20+%20\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

loop:
ω(T) -> ψ -> v -> T

