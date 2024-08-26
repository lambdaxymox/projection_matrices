# Sympy Projection Matrices
**projection_matrices** is a sympy module that allows one to work with perspective projection
and orthographic projection matrices symbolically using the **sympy** computer algebra system.

## Getting Started
To build the project, run 

```bash
poetry build
```

in the root of the source tree. To run the test suite, run

```bash
poetry run pytest
```

from the text shell. The test suite uses the `pytest` framework. To use 
**projection_matrices** in your project, simply import the module

```python
import projection_matrices as pm
```

and you're done.

## Building The Documentation

The documentation found in the `docs` folder is a jupyter book. It covers the 
derivation of the matrices provided by the library. Enter

```bash
poetry run jupyter-book build docs/
```

to build the documentation for the project. The documentation is a static site 
that can be read in a web browser.
