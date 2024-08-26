# Sympy Projection Matrices

This book is the more extensive documentation of the sympy projection matrices
library. The purpose of the library is threefold: speed up deriving projection
matrices; explore using sympy as a tool for assisting in developing algorithms
in more mathematical engineering domains; develop the documentation and exposition
of a library in interactive book form using Jupyter. 

The library exposes the perspective projection, perspective fov
projection, and orthographic projection matrices using a canonically chosen coordinate
system. In this case, the canonically chose coordinate system is OpenGL's clip space.
The library documentation derives the generalized projection matrices using sympy, and 
then it shows how to apply it to different platform APIs. The docs contain cases for Vulkan,
OpenGL, DirectX, and Metal. The basic workflow for other situations is more or less 
identical.
