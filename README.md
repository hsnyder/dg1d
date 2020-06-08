# dg1d
Arbitrary-order solver for the 1D advection equation using the discontinuous Galerkin method.
Time discretisation is by the 4th order Runge-Kutta method
Written as self-study. There are a number of inefficiencies. Not recommended for "real" use.

Run make to build it, execute ./dg1d, then run ./vis.py to create an .mp4 file with the animated plot.
Animation requires that ffmpeg is installed (along with python3, numpy, matplotlib).


