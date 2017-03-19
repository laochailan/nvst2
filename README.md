![example solution](https://github.com/laochailan/nvst2/blob/master/screenshots/colormode2.png?raw=true)
# nvst2

This is a simulation of the 2D [incompressible Navier-Stokes equation][1]. There is
some graphical output and an FFMPEG backend to record the results in a video.

[1]: https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations#Incompressible_flow

## How to build

```
$ mkdir build; cd build
$ cmake ..
$ make
```

run with

```
$ src/nvst2
```

## Method

In case you wonder about how I did this:


For an incompressible fluid of constant density, the pressure is like a constraint force enforcing

```div v = 0```

In the code I handled the pressure contribution by projecting the solution on
the subspace of divergence free functions each time step (by subtracting a
gradient field with the right divergence).
This requires solving a Poisson equation.

The convection term is handled using an explicit differencing scheme. Even if
you add diffusion terms this does not improve as they themselves aren’t stable
in an explicit scheme.

The key is to use an implicit scheme for the diffusion term. This requires
solving another poisson-like equation, but is unconditionally stable. For high
enough diffusion constants, the overall scheme is then stable.
The implicit scheme requires the solution of a Poisson-like equation

Additionally, I put a cap on the velocity so it cannot exceed Δx/Δt anywhere
(where Δx, and Δt are the space and time discretization).

The linear systems of equations are solved using the Fast Sine Transform.

The boundary conditions are some cheaty periodic boundary conditions. I didn’t do that very cleanly, but it works kind of.

It’s probably not very accurate, but accurate enough to produce something that looks like a fluid so I’m happy with it.

![Animated output](https://github.com/laochailan/nvst2/blob/master/screenshots/example.gif?raw=true)
