# PBD-2d-sim
Position-Based Dynamics 2d simulation

A small personal project in C++ I did to play with the Position-Based Dynamics (PBD) method for simulating particles with friction (e.g. sand).

![Quick demo](https://media.giphy.com/media/UkyPhN6lSwJjNIPhfV/giphy.gif)![Another quick demo](https://media.giphy.com/media/WnRcLGYGB2TKOzxWz2/giphy.gif)

It relies on [Magnum](https://magnum.graphics/), a cross-platform  graphics library, and [Intel's TBB](https://github.com/oneapi-src/oneTBB) for multi-threading. Except for the physics part, a good chunk of the code is derived from the APIC Fluid simulation from the [Magnum example gallery](https://magnum.graphics/showcase/). 

An excellent introduction to PBD is available in this paper: http://mmacklin.com/EG2015PBD.pdf.

While the hybrid Material Point Method (MPM) is gaining popularity for simulating particulate effects such as snow or sand in the Computer Graphics community, PBD remains an amazingly fast simulation tool that offers a simple trade off between precision and speed by controlling the number of constraint iterations.

I have used extensively [Houdini's implementation of PBD](https://www.sidefx.com/docs/houdini/grains/about.html) on the movie How to Train your Dragon: the Hidden World. See [some sand examples here](https://vimeo.com/156511737#t=35s). Here is another example from the movie Abominable where I simulated [an avalanche](https://vimeo.com/156511737#t=94s) using 120 million simulated particles.

Having fun with the washing machine:

![Washer](https://media.giphy.com/media/N9l1VG8Yl08TzS8tuu/giphy.gif)
