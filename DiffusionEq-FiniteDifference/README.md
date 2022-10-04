# Solving the diffusion equation using Finite Differences

## Summary



## Numerically solving the diffusion equation
Diffusion describes motion from high concenrtation to low concentration; think of how the amount of heat you feel from a radiator decreases with increasing distance. The time-dependent diffusion equation for a scalar field, $T$, is given by
$$\dfrac{\partial T}{\partial t}= \nu \dfrac{\partial ^2 T}{\partial x^2}.$$
Here $alpha$ is the diffusivity with units $m^2/s$.

For a given initial condition, $u_0=u(x,0)$, this equation has an analytical solution,
$$u(x,t)=u_0(x-ut).$$
This analytical solution allows us to compute the error of different implementations and directly compare finite difference schemes against one another.

## Solving the Diffusion Equation in 1D
The diffusion equation can be discretised on to a grid of points in space and time (indexed by $j$ and $n$, respectively). Here, we implement the forward in time and centred in space (FTCS) scheme which approximates temporal and spatial derivatives, respectively, as
$$\dfrac{\partial T}{\partial t} = \dfrac{T^{n+1}_j-T^n_j}{\Delta t}, \dfrac{\partial^2 T}{\partial x^2} = \dfrac{T^n\_{j+1}-2T^n_j+T^n\_{j-1}}{2 \Delta x},$$

The discretised equation then reads
$$ \dfrac{T^{n+1}_j-T^n_j}{\Delta t}= \nu \dfrac{T^n_{j+1}-2T^n_j+T^n\_{j-1}}{2 \Delta x},$$

$$T^{n+1}_j = \dfrac{\nu \Delta t}{\Delta x^2} \left(T^n\_{j+1}-2T^n_j+T^n\_{j-1} \right).$$
where $c=v\Delta t/ \Delta x$, known as the Courant number. The Courant number indicates how much the information travels, if greater than one, it means that the information propagates through more than one grid cell at each time step, making the solution inaccurate and potentially leading to nonphysical results or divergence of the solution.
 
To close the system, two boundary conditions are needed. We choose to apply the same discretisation scheme to the first node
$$u^{n+1}_0=u^n_0-\dfrac{c}{2}\left(u^n_1-u^n\_{N-2} \right).$$
The periodic boundary conditions dictate that the two boundary grid points have the same value,
$$u^{n+1}\_{N-1}=u^{n+1}_0.$$

In practice, the discretised evolution equation and boundary conditions are solved over all space $(j)$ and then looped over the total number of timesteps.

Other explicit schemes exist and can be derived from different approximations of the gradient terms in an equation. Using the forward in time and forward in space approximations,
$$\dfrac{\partial u}{\partial t} = \dfrac{u^{n+1}_j-u^n_j}{\Delta t}, \dfrac{\partial u}{\partial x} = \dfrac{u^n\_{j+1}-u^n\_{j}}{\Delta x},$$
respectively.
The FTFS scheme 
$$u^{n+1}_j = u^n_j - c\left(u^n\_{j+1}-u^n\_{j}\right).$$
This looks similar to the scheme discussed before but notice that the spatial derivative is now only approximated from the grid cell to the left. 
Alternative schemes can be derived from different approximations, I suggest the reader to use Tannehill as a starting point.

## Implicit Scheme
Unlike explicit schemes which use only values at the current and previous timesteps to predict the next-step, implicit schemes also use values at the time-step $(n+1)$ to calculate values at time-step, $n+1$. When correctly formulated, implicit schemes can be stable for arbitrarily large time-step, although they will not necessarily give high accuracy in time and they can be computationally expensive per time-step. 

The backward in time, centred in space scheme is:
$$u^{n+1}_j=u^n_j-\dfrac{c}{2}\left(u^{n+1}\_{j+1}-u^{n+1}\_{j-1}\right).$$
In order to evaluate all the values of $u_j^{n+1}$ we solve the set of simultaneous equations given in matrix form. This matrix is dependent on the boundary conditions, our case of periodic boundary conditions gives $u_0=u\_{N}$:
$$u^{n+1}_0=u^n_0-\dfrac{c}{2}\left(u^{n+1}\_{1}-u^{n+1}\_{N-1}\right),$$
$$u^{n+1}\_{N-1}=u^n\_{N-1}-\dfrac{c}{2}\left(u^{0}\_{1}-u^{n+1}\_{N-2}\right).$$
The matrix equation,

<p align="center">
  <img src="linearsystem.png" />
</p>

is solved using Gaussian elimination giving $\boldsymbol{u}^{n+1}=M^{-1}\boldsymbol{u}^n$. Other computational tools to solves sets of linear simultaneous equations include the Jacobi method, Gauss-Seidel method or the conjugate gradient.

## A comparison
