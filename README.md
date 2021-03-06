# Fortran_Numerical_Modelling

This is the homework of ETH course.
T: temperature
x,y: direction


## day1-3 ##

Basic Fortran knowledge.

## day4 ##

2D diffusion with initialization of spike and random. 

![day4_1](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

<img src="day4/ex2-4/delta_result.png" width="320" height="240"> <img src="day4/ex2-4/random_result.png" width="320" height="240">


## day5 ##

2D advection-diffusion with a fixed flow field and zero flux at the boundary x and 0 and 1 at the boundary y. 

![day5_1](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=-v_x\frac{\partial%20T}{\partial%20x}%20-%20v_y\frac{\partial%20T}{\partial%20y}%20+%20\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

Streamfunction：

![day5_2](https://latex.codecogs.com/svg.latex?(v_x,%20v_y)%20=%20(\frac{\partial%20\psi}{\partial%20y},%20-\frac{\partial%20\psi}{\partial%20x})) 

![day5_3](https://latex.codecogs.com/svg.latex?\psi=Bsin(\frac{\pi%20x}{x_{max}})sin(\frac{\pi%20y}{y_{max}})) 

loop t:\
T(t) -> ψ -> v -> T(t+1)

<img src="day5/01_image.png" width="320" height="240"> <img src="day5/1_image.png" width="320" height="240"> 

<img src="day5/10_image.png" width="320" height="240"> <img src="day5/100_image.png" width="320" height="240">

## day6 ##

Iterative Poisson solver.

![day6_1](https://latex.codecogs.com/svg.latex?\frac{\partial^2%20u}{\partial%20x^2}%20+%20\frac{\partial^2%20u}{\partial%20y^2}=f)

<img src="day6/V_257_image.png" width="320" height="240">

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

![day7_6](https://latex.codecogs.com/svg.latex?\nabla^2%20w%20=%20Ra\frac{\partial%20T}{\partial%20x})

![day7_5](https://latex.codecogs.com/svg.latex?\nabla^2\psi=%20w)

![day7_7](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=-v_x\frac{\partial%20T}{\partial%20x}%20-%20v_y\frac{\partial%20T}{\partial%20y}%20+%20\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

loop t:\
T(t) -> ω(T(t)) -> ψ -> v -> T(t+1)

Ra = 0.1

ω, ψ:\
<img src="day7/Pr-None_time-0.1_Ra-100000.0/w_result.png" width="320" height="240"> <img src="day7/Pr-None_time-0.1_Ra-100000.0/S_result.png" width="320" height="240"> 

T:\
<img src="day7/Pr-None_time-0.1_Ra-100000.0/T_result.png" width="320" height="240">



## day8 ##

2D convection-diffusion along time with vorticity w depending on Prandtl number Pr. \
This is for no so high viscous fluid. When Pr = inf., it returns Poisson equation as in day7.

![day8_1](https://latex.codecogs.com/svg.latex?\frac{1}{Pr}(\frac{\partial%20w}{\partial%20t}+v_x\frac{\partial%20w}{\partial%20x}+v_y\frac{\partial%20w}{\partial%20y})=\nabla^2%20w%20-%20Ra\frac{\partial%20T}{\partial%20x})

![day8_2](https://latex.codecogs.com/svg.latex?\frac{\partial%20T}{\partial%20t}=-v_x\frac{\partial%20T}{\partial%20x}%20-%20v_y\frac{\partial%20T}{\partial%20y}%20+%20\kappa(\frac{\partial^2%20T}{\partial%20x^2}%20+%20\frac{\partial^2%20T}{\partial%20y^2})) 

loop t:\
T(t) -> ω(T(t),t) -> ψ -> v -> T(t+1), w(t+1)

Ra = 0.1
Pr = 0.01

ω, ψ:\
<img src="day8/Pr-0.01_time-0.1_Ra-100000.0/w_result.png" width="320" height="240"> <img src="day8/Pr-0.01_time-0.1_Ra-100000.0/S_result.png" width="320" height="240"> 

T:\
<img src="day8/Pr-0.01_time-0.1_Ra-100000.0/T_result.png" width="320" height="240">

Ra = 0.1
Pr = 1

ω, ψ:\
<img src="day8/Pr-1_time-0.1_Ra-100000.0/w_result.png" width="320" height="240"> <img src="day8/Pr-1_time-0.1_Ra-100000.0/S_result.png" width="320" height="240"> 

T:\
<img src="day8/Pr-1_time-0.1_Ra-100000.0/T_result.png" width="320" height="240">

## day9 ##

Implicit, semi-implicit diffusion.
