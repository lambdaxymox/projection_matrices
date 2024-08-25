import sympy

from dataclasses import dataclass


@dataclass(frozen=True)
class FrustumBounds:
    """
    A data class describing the shape of the viewing frustum for a projection.
    """
    left: sympy.Symbol
    right: sympy.Symbol
    bottom: sympy.Symbol
    top: sympy.Symbol
    near: sympy.Symbol
    far: sympy.Symbol


@dataclass(frozen=True)
class NDCBounds:
    """
    a data class describing the bounds of the canonical view volume.
    """
    horizontal_min: sympy.Symbol
    horizontal_max: sympy.Symbol
    vertical_min: sympy.Symbol
    vertical_max: sympy.Symbol
    depth_min: sympy.Symbol
    depth_max: sympy.Symbol

    def __str__(self):
        h_min = self.horizontal_min
        h_max = self.horizontal_max
        v_min = self.vertical_min
        v_max = self.vertical_max
        d_min = self.depth_min
        d_max = self.depth_max

        return f'[{h_min}, {h_max}] x [{v_min}, {v_max}] x [{d_min} {d_max}]'


@dataclass(frozen=True)
class FrustumFovBounds:
    """
    A data class describing the shape of the viewing frustum for a perspective projection.
    """
    aspect_ratio: sympy.Symbol
    vfov: sympy.Symbol
    near: sympy.Symbol
    far: sympy.Symbol


def perspective(frustum_bounds: FrustumBounds, ndc_bounds: NDCBounds) -> sympy.Matrix:
    """
    Generate an instance of a perspective projection in the canonical orthonormal frames.

    The **canonical view space** is a vector space with a left-handed orthonormal frame
    defined as follows.

    * The **origin** of the coordinate system is `[0, 0, 0]^T`.
    * The **positive x-axis** is the horizontal direction points to the right.
    * The **positive y-axis** is the vertical direction and points up.
    * The **positive z-axis** is the depth direction and points into the viewing volume.

    The **canonical clip space** is a vector space with a left-handed orthonormal frame
    defined as follows.

    * The **origin** of the coordinate system is `[0, 0, 0]^T`.
    * The **positive x-axis** is the horizontal direction points to the right.
    * The **positive y-axis** is the vertical direction and points up.
    * The **positive z-axis** is the depth direction and points into the viewing volume.

    The **canonical view volume** is a vector space with a left-handed orthonormal frame
    identical to the clip space frame with bounds specified by `ndc_bounds`.

    Parameter Specification:
    The fundamental parametrization of the perspective projection transformation
    is the specification based on defining the placement of the frustum bounds.
    We represent the frustum bounds by defining the placements with respect to the
    **view space** orthonormal frame vectors. More precisely, the fundamental
    parametrization is given by the parameters `left`, `right`, `bottom`, `top`,
    `near`, and `far` such that

    * left   > 0
    * right  > 0
    * bottom > 0
    * top    > 0
    * far    > near > 0

    where the parameters define the placement of the planes. The plane placement
    definitions follow.

    * `left` defines the location of the **left plane** by its distance along
      the **negative x-axis** from the origin of the coordinate frame.
      The **left plane** is a plane parallel to the **yz-plane**.
    * `right` defines the location of the **right plane** by its distance along
      the **positive x-axis** from the origin of the coordinate frame.
      The **right plane** is a plane parallel to the **yz-plane**.
    * `bottom` defines the location of the **bottom plane** by its distance along
      the **negative y-axis** from the origin of the coordinate frame.
      The **bottom plane** is a plane parallel to the **zx-plane**.
    * `top` defines the location of the **top plane** by its distance along
      the **positive y-axis** from the origin of the coordinate frame.
      The **top plane** is a plane parallel to the **zx-plane**.
    * `near` defines the location of the **near plane** by its distance along
      the **negative z-axis** from the origin of the coordinate frame.
      The **near plane** is a plane parallel to the **xy-plane**.
    * `far` defines the location of the **far plane** by its distance along
      the **negative z-axis** from the origin of the coordinate frame.
      The **far plane** is a plane parallel to the **xy-plane**.

    Parameters:
    - frustum_bounds: The bounds of the frustum defined in terms of relative displacements
    along the coordinate axes.
    - ndc_bounds: The bounds of the viewing volume in normalized device coordinates.

    Returns:
    - A 4x4 perspective projection matrix.
    """
    l = frustum_bounds.left
    r = frustum_bounds.right
    b = frustum_bounds.bottom
    t = frustum_bounds.top
    n = frustum_bounds.near
    f = frustum_bounds.far

    h_min = ndc_bounds.horizontal_min
    h_max = ndc_bounds.horizontal_max
    v_min = ndc_bounds.vertical_min
    v_max = ndc_bounds.vertical_max
    d_min = ndc_bounds.depth_min
    d_max = ndc_bounds.depth_max

    c0r0 = ((h_max - h_min) * n) / (r - (-l))
    c1r0 = 0
    c2r0 = (h_min * r - h_max * (-l)) / (r - (-l))
    c3r0 = 0

    c0r1 = 0
    c1r1 = ((v_max - v_min) * n) / (t - (-b))
    c2r1 = (v_min * t - v_max * (-b)) / (t - (-b))
    c3r1 = 0

    c0r2 = 0
    c1r2 = 0
    c2r2 = (d_max * f - d_min * n) / (f - n)
    c3r2 = -((d_max - d_min) * f * n) / (f - n)

    c0r3 = 0
    c1r3 = 0
    c2r3 = 1
    c3r3 = 0

    matrix = sympy.Matrix([
        [c0r0, c1r0, c2r0, c3r0],
        [c0r1, c1r1, c2r1, c3r1],
        [c0r2, c1r2, c2r2, c3r2],
        [c0r3, c1r3, c2r3, c3r3]
    ])

    return matrix


def orthographic(frustum_bounds: FrustumBounds, ndc_bounds: NDCBounds) -> sympy.Matrix:
    """
    Generate an instance of a orthographic projection in the canonical orthonormal frames.

    The **canonical view space** is a vector space with a left-handed orthonormal frame
    defined as follows.

    * The **origin** of the coordinate system is `[0, 0, 0]^T`.
    * The **positive x-axis** is the horizontal direction points to the right.
    * The **positive y-axis** is the vertical direction and points up.
    * The **positive z-axis** is the depth direction and points into the viewing volume.

    The **canonical clip space** is a vector space with a left-handed orthonormal frame
    defined as follows.

    * The **origin** of the coordinate system is `[0, 0, 0]^T`.
    * The **positive x-axis** is the horizontal direction points to the right.
    * The **positive y-axis** is the vertical direction and points up.
    * The **positive z-axis** is the depth direction and points into the viewing volume.

    The **canonical view volume** is a vector space with a left-handed orthonormal frame
    identical to the clip space frame with bounds specified by `ndc_bounds`.

    Parameter Specification:
    The fundamental parametrization of the orthographic projection transformation
    is the specification based on defining the placement of the frustum bounds.
    We represent the frustum bounds by defining the placements with respect to the
    **view space** orthonormal frame vectors. More precisely, the fundamental
    parametrization is given by the parameters `left`, `right`, `bottom`, `top`,
    `near`, and `far` such that

    * left   > 0
    * right  > 0
    * bottom > 0
    * top    > 0
    * far    > near > 0

    where the parameters define the placement of the planes. The plane placement
    definitions follow.

    * `left` defines the location of the **left plane** by its distance along
      the **negative x-axis** from the origin of the coordinate frame.
      The **left plane** is a plane parallel to the **yz-plane**.
    * `right` defines the location of the **right plane** by its distance along
      the **positive x-axis** from the origin of the coordinate frame.
      The **right plane** is a plane parallel to the **yz-plane**.
    * `bottom` defines the location of the **bottom plane** by its distance along
      the **negative y-axis** from the origin of the coordinate frame.
      The **bottom plane** is a plane parallel to the **zx-plane**.
    * `top` defines the location of the **top plane** by its distance along
      the **positive y-axis** from the origin of the coordinate frame.
      The **top plane** is a plane parallel to the **zx-plane**.
    * `near` defines the location of the **near plane** by its distance along
      the **negative z-axis** from the origin of the coordinate frame.
      The **near plane** is a plane parallel to the **xy-plane**.
    * `far` defines the location of the **far plane** by its distance along
      the **negative z-axis** from the origin of the coordinate frame.
      The **far plane** is a plane parallel to the **xy-plane**.

    Parameters:
    - frustum_bounds: The bounds of the frustum defined in terms of relative displacements
    along the coordinate axes.
    - ndc_bounds: The bounds of the viewing volume in normalized device coordinates.

    Returns:
    - An 4x4 orthographic projection matrix.
    """
    l = frustum_bounds.left
    r = frustum_bounds.right
    b = frustum_bounds.bottom
    t = frustum_bounds.top
    n = frustum_bounds.near
    f = frustum_bounds.far

    h_min = ndc_bounds.horizontal_min
    h_max = ndc_bounds.horizontal_max
    v_min = ndc_bounds.vertical_min
    v_max = ndc_bounds.vertical_max
    d_min = ndc_bounds.depth_min
    d_max = ndc_bounds.depth_max

    c0r0 = (h_max - h_min) / (r - (-l))
    c1r0 = 0
    c2r0 = 0
    c3r0 = (h_min * r - h_max * (-l)) / (r - (-l))

    c0r1 = 0
    c1r1 = (v_max - v_min) / (t - (-b))
    c2r1 = 0
    c3r1 = (v_min * t - v_max * (-b)) / (t - (-b))

    c0r2 = 0
    c1r2 = 0
    c2r2 = (d_max - d_min) / (f - n)
    c3r2 = (d_min * f - d_max * n) / (f - n)

    c0r3 = 0
    c1r3 = 0
    c2r3 = 0
    c3r3 = 1

    matrix = sympy.Matrix([
        [c0r0, c1r0, c2r0, c3r0],
        [c0r1, c1r1, c2r1, c3r1],
        [c0r2, c1r2, c2r2, c3r2],
        [c0r3, c1r3, c2r3, c3r3]
    ])

    return matrix


def perspective_fov(frustum_fov_bounds: FrustumFovBounds, ndc_bounds: NDCBounds) -> sympy.Matrix:
    """
    Generate an instance of a perspective projection in the canonical orthonormal frames.

    The **canonical view space** is a vector space with a left-handed orthonormal frame
    defined as follows.

    * The **origin** of the coordinate system is `[0, 0, 0]^T`.
    * The **positive x-axis** is the horizontal direction points to the right.
    * The **positive y-axis** is the vertical direction and points up.
    * The **positive z-axis** is the depth direction and points into the viewing volume.

    The **canonical clip space** is a vector space with a left-handed orthonormal frame
    defined as follows.

    * The **origin** of the coordinate system is `[0, 0, 0]^T`.
    * The **positive x-axis** is the horizontal direction points to the right.
    * The **positive y-axis** is the vertical direction and points up.
    * The **positive z-axis** is the depth direction and points into the viewing volume.

    The **canonical view volume** is a vector space with a left-handed orthonormal frame
    identical to the clip space frame with bounds specified by `ndc_bounds`.

    Parameter Specification:
    The fundamental parametrization of the perspective projection transformation
    is the specification based on defining the placement of the frustum bounds.
    We represent the frustum bounds by defining the placements with respect to the
    **view space** orthonormal frame vectors. More precisely, the fundamental
    parametrization is given by the parameters `aspect_ratio`, `vfov`, `near`, and
    `far` such that

    * aspect_ratio > 0
    * vfov         > 0
    * vfov         < pi
    * far          > near > 0

    where the parameters define the placement of the planes. The plane placement
    definitions follow.

    * `left` defines the location of the **left plane** by its distance along
      the **negative x-axis** from the origin of the coordinate frame.
      The **left plane** is a plane parallel to the **yz-plane**.
    * `right` defines the location of the **right plane** by its distance along
      the **positive x-axis** from the origin of the coordinate frame.
      The **right plane** is a plane parallel to the **yz-plane**.
    * `bottom` defines the location of the **bottom plane** by its distance along
      the **negative y-axis** from the origin of the coordinate frame.
      The **bottom plane** is a plane parallel to the **zx-plane**.
    * `top` defines the location of the **top plane** by its distance along
      the **positive y-axis** from the origin of the coordinate frame.
      The **top plane** is a plane parallel to the **zx-plane**.
    * `near` defines the location of the **near plane** by its distance along
      the **negative z-axis** from the origin of the coordinate frame.
      The **near plane** is a plane parallel to the **xy-plane**.
    * `far` defines the location of the **far plane** by its distance along
      the **negative z-axis** from the origin of the coordinate frame.
      The **far plane** is a plane parallel to the **xy-plane**.

    In the case of the field of view specification, we have

    * left   == aspect_ratio * near * tan(vfov / 2)
    * right  == aspect_ratio * near * tan(vfov / 2)
    * bottom == near * tan(vfov / 2)
    * top    == near * tan(vfov / 2)

    Parameters:
    - frustum_fov_bounds: The bounds of the frustum defined in terms of the vertical field
      of view and aspect ratio.
    - ndc_bounds: The bounds of the viewing volume in normalized device coordinates.

    Returns:
    - A 4x4 perspective projection matrix.
    """
    aspect_ratio = frustum_fov_bounds.aspect_ratio
    theta_vfov = frustum_fov_bounds.vfov
    n = frustum_fov_bounds.near
    f = frustum_fov_bounds.far

    h_min = ndc_bounds.horizontal_min
    h_max = ndc_bounds.horizontal_max
    v_min = ndc_bounds.vertical_min
    v_max = ndc_bounds.vertical_max
    d_min = ndc_bounds.depth_min
    d_max = ndc_bounds.depth_max

    c0r0 = ((h_max - h_min) / (sympy.Rational(2, 1) * aspect_ratio)) * (sympy.S.One / sympy.tan(theta_vfov / 2))
    c1r0 = 0
    c2r0 = sympy.Rational(1, 2) * (h_max + h_min)
    c3r0 = 0

    c0r1 = 0
    c1r1 = (sympy.Rational(1, 2) * (v_max - v_min)) * (sympy.S.One / sympy.tan(sympy.Rational(1, 2) * theta_vfov))
    c2r1 = sympy.Rational(1, 2) * (v_max + v_min)
    c3r1 = 0

    c0r2 = 0
    c1r2 = 0
    c2r2 = (d_max * f - d_min * n) / (f - n)
    c3r2 = -(d_max - d_min) * ((f * n) / (f - n))

    c0r3 = 0
    c1r3 = 0
    c2r3 = 1
    c3r3 = 0

    matrix = sympy.Matrix([
        [c0r0, c1r0, c2r0, c3r0],
        [c0r1, c1r1, c2r1, c3r1],
        [c0r2, c1r2, c2r2, c3r2],
        [c0r3, c1r3, c2r3, c3r3]
    ])

    return matrix
