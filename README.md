# PBD-2d-sim
Position-Based Dynamics 2d simulation

A small personal project I did to play with the Position-Based Dynamics (PBD) method for simulating particles with friction (e.g. sand).

![Quick demo](https://media.giphy.com/media/UkyPhN6lSwJjNIPhfV/giphy.gif)


It relies on [Magnum](https://magnum.graphics/), a cross-platform  graphics library, and [Intel's TBB](https://github.com/oneapi-src/oneTBB) for multi-threading. Some of the code is derived from the APIC 2d Fluid simulation from the [Magnum example gallery](https://magnum.graphics/showcase/). 

An excellent introduction to PBD is available in this paper: http://mmacklin.com/EG2015PBD.pdf.

While the hybrid Material Point Method (MPM) is gaining popularity for simulating particulate effects such as snow or sand in the Computer Graphics community, PBD remains an amazingly fast simulation tool that offers a simple trade off between precision and speed by controlling the number of constraint iterations.

I have used extensively [Houdini's implementation of (X)PBD](https://www.sidefx.com/docs/houdini/grains/about.html) on the movie How to Train your Dragon: the Hidden World. See [some sand examples here](https://vimeo.com/156511737#t=35s).

