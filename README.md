# maxwell_fdfd

### Two point sources
![](https://github.com/kjabon/maxwell_fdfd/blob/main/demo6.gif)

### Gaussian source
![](https://github.com/kjabon/maxwell_fdfd/blob/main/demo8.gif)

### Summary
A general matrix equation solver for 2D Yee grids, supporting perfectly matched layers (PMLs)
In principle this could be applied to any differential matrix equation, e.g. the heat equation.
Here we apply it to Maxwell's equations in the frequency domain (hence, finite difference frequency domain, or FDFD).

You can add arbitrary materials to the grid, e.g. a Si waveguide or lens (with the correct permittivity $\eps_r$ and permeability $\mu_r$ at their respective grid points).
You can then solve for the optical eigenmodes of the system (e.g., a waveguide). E.g., see ``/demos/scalarModes2D.m``.

You can add, also, arbitrary light sources to the grid, solve the corresponding matrix equation, and step the phase to see the light propagating through your system. E.g., see ``/demos/unidirectionalGaussianSourceDemo.m``
An astute observer may consider solving for the modes of a waveguide and then using these as light sources for a larger system...! 
As this is simply a general solver, the limit is your time and creativity.

### Installation
Just install Matlab (after getting a license...) and you're good to clone! Get started in the ``/demos`` folder.
