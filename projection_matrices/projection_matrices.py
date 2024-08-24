import sympy

from dataclasses import dataclass


@dataclass(frozen=True)
class FrustumBounds:
    left: sympy.Symbol
    right: sympy.Symbol
    bottom: sympy.Symbol
    top: sympy.Symbol
    near: sympy.Symbol
    far: sympy.Symbol


@dataclass(frozen=True)
class NDCBounds:
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
    aspect_ratio: sympy.Symbol
    vfov: sympy.Symbol
    near: sympy.Symbol
    far: sympy.Symbol
    

def perspective(frustum_bounds: FrustumBounds, ndc_bounds: NDCBounds) -> sympy.Matrix:
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

