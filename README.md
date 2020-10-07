# PBD-2d-sim
Position-Based Dynamics 2d simulation

This is an implementation of the Position-Based Dynamics (PBD) method for simulating particles with friction (e.g. sand). 
It relies on [Magnum](https://magnum.graphics/), a cross-platform  graphics library and some of the code is derived from the APIC 2d Fluid simulation from the [Magnum example gallery](https://magnum.graphics/showcase/). This implementation employs [Intel's TBB](https://github.com/oneapi-src/oneTBB) for multi-threading.

An excellent introduction to PBD is available in this paper: http://mmacklin.com/EG2015PBD.pdf. The more recent extension known as XPBD is described here: http://mmacklin.com/xpbd.pdf

While the hybrid Material Point Method (MPM) is gaining popularity for simulating particulate effects such as snow or sand in the Computer Graphics community, the PBD method remains an amazingly fast simulation tool that offers a simple trade off between precision and speed by controlling the number of constraint iterations.

I have used extensively Houdini's implementation of (X)PBD, a.k.a the [grain solver](https://www.sidefx.com/docs/houdini/grains/about.html) on the movie How to Train your Dragon: the Hidden World. See [some sand examples in my reel](https://vimeo.com/156511737#t=35s).

