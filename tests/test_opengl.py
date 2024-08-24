import projection_matrices as pm
import sympy
import pytest


def change_of_orientation_lh_to_rh() -> sympy.Matrix:
    return sympy.Matrix([
        [1, 0,  0, 0],
        [0, 1,  0, 0],
        [0, 0, -1, 0],
        [0, 0,  0, 1]
    ])

def change_of_orientation_rh_to_lh() -> sympy.Matrix:
    return sympy.Matrix([
        [1, 0,  0, 0],
        [0, 1,  0, 0],
        [0, 0, -1, 0],
        [0, 0,  0, 1]
    ])


class TestOpenGLLeftHanded:
    def test_perspective_projection_asymmetric(self):
        l, r, b, t, n, f = sympy.symbols('l r b t n f')
        frustum_bounds = pm.FrustumBounds(l, r, b, t, n, f)
        ndc_bounds = pm.NDCBounds(-1, 1, -1, 1, -1, 1)
        expected = sympy.Matrix([
            [ (2 * n) / (r - (-l)), 0,                    -(r + (-l)) / (r - (-l)),  0                     ],
            [ 0,                    (2 * n) / (t - (-b)), -(t + (-b)) / (t - (-b)),  0                     ],
            [ 0,                    0,                     (f + n) / (f - n),       -(2 * f * n) / (f - n) ],
            [ 0,                    0,                     1,                        0                     ]
        ])
        x_lh_lh = sympy.Matrix.eye(4)
        m_coord = sympy.Matrix.eye(4)
        m_coord_inv = sympy.Matrix.eye(4)
        m_canonical_lh_lh = pm.perspective(frustum_bounds, ndc_bounds)
        result = (x_lh_lh * m_coord_inv) * m_canonical_lh_lh * (m_coord * x_lh_lh)
        
        assert result.equals(expected)
    

    def test_perspective_fov_projection_symmetric(self):
        aspect, theta_vfov, n, f = sympy.symbols('aspect theta_vfov n f')
        frustum_bounds = pm.FrustumFovBounds(aspect, theta_vfov, n, f)
        ndc_bounds = pm.NDCBounds(-1, 1, -1, 1, -1, 1)

        c0r0 = 1 / (aspect * sympy.tan(sympy.Rational(1, 2) * theta_vfov))
        c1r1 = 1 / sympy.tan(sympy.Rational(1, 2) * theta_vfov)
        c2r2 = (f + n) / (f - n)
        c3r2 = -(2 * f * n) / (f - n)

        expected = sympy.Matrix([
            [ c0r0, 0,     0,    0    ],
            [ 0,    c1r1,  0,    0    ],
            [ 0,    0,     c2r2, c3r2 ],
            [ 0,    0,     1,    0    ]
        ])

        x_lh_lh = sympy.Matrix.eye(4)
        m_coord = sympy.Matrix.eye(4)
        m_coord_inv = sympy.Matrix.eye(4)
        m_canonical_lh_lh = pm.perspective_fov(frustum_bounds, ndc_bounds)
        result = (x_lh_lh * m_coord_inv) * m_canonical_lh_lh * (m_coord * x_lh_lh)
        
        assert result.equals(expected)
    
    def test_orthographic_projection(self):
        l, r, b, t, n, f = sympy.symbols('l r b t n f')
        frustum_bounds = pm.FrustumBounds(l, r, b, t, n, f)
        ndc_bounds = pm.NDCBounds(-1, 1, -1, 1, -1, 1)
        expected = sympy.Matrix([
            [ 2 / (r - (-l)), 0,              0,           -(r + (-l)) / (r - (-l)) ],
            [ 0,              2 / (t - (-b)), 0,           -(t + (-b)) / (t - (-b)) ],
            [ 0,              0,              2 / (f - n), -(f + n) / (f - n)       ],
            [ 0,              0,              0,            1                       ]
        ])
        x_lh_lh = sympy.Matrix.eye(4)
        m_coord = sympy.Matrix.eye(4)
        m_coord_inv = sympy.Matrix.eye(4)
        m_canonical_lh_lh = pm.orthographic(frustum_bounds, ndc_bounds)
        result = (x_lh_lh * m_coord_inv) * m_canonical_lh_lh * (m_coord * x_lh_lh)
        
        assert result.equals(expected)


class TestOpenGLRightHanded:
    def test_perspective_projection_asymmetric(self):
        l, r, b, t, n, f = sympy.symbols('l r b t n f')
        frustum_bounds = pm.FrustumBounds(l, r, b, t, n, f)
        ndc_bounds = pm.NDCBounds(-1, 1, -1, 1, -1, 1)
        expected = sympy.Matrix([
            [ (2 * n) / (r - (-l)), 0,                     (r + (-l)) / (r - (-l)),  0                     ],
            [ 0,                    (2 * n) / (t - (-b)),  (t + (-b)) / (t - (-b)),  0                     ],
            [ 0,                    0,                    -(f + n) / (f - n),       -(2 * f * n) / (f - n) ],
            [ 0,                    0,                    -1,                        0                     ]
        ])
        x_lh_lh = sympy.Matrix.eye(4)
        x_rh_lh = change_of_orientation_rh_to_lh()
        m_coord = sympy.Matrix.eye(4)
        m_coord_inv = sympy.Matrix.eye(4)
        m_canonical_lh_lh = pm.perspective(frustum_bounds, ndc_bounds)
        result = (x_lh_lh * m_coord_inv) * m_canonical_lh_lh * (m_coord * x_rh_lh)

        assert result.equals(expected)
    
    def test_perspective_fov_projection_symmetric(self):
        aspect, theta_vfov, n, f = sympy.symbols('aspect theta_vfov n f')
        frustum_bounds = pm.FrustumFovBounds(aspect, theta_vfov, n, f)
        ndc_bounds = pm.NDCBounds(-1, 1, -1, 1, -1, 1)

        c0r0 = 1 / (aspect * sympy.tan(sympy.Rational(1, 2) * theta_vfov))
        c1r1 = 1 / sympy.tan(sympy.Rational(1, 2) * theta_vfov)
        c2r2 = -(f + n) / (f - n)
        c3r2 = -(2 * f * n) / (f - n)

        expected = sympy.Matrix([
            [ c0r0, 0,      0,    0    ],
            [ 0,    c1r1,   0,    0    ],
            [ 0,    0,      c2r2, c3r2 ],
            [ 0,    0,     -1,    0    ]
        ])
        
        x_lh_lh = sympy.Matrix.eye(4)
        x_rh_lh = change_of_orientation_rh_to_lh()
        m_coord = sympy.Matrix.eye(4)
        m_coord_inv = sympy.Matrix.eye(4)
        m_canonical_lh_lh = pm.perspective_fov(frustum_bounds, ndc_bounds)
        result = (x_lh_lh * m_coord_inv) * m_canonical_lh_lh * (m_coord * x_rh_lh)
        
        assert result.equals(expected)
    
    def test_orthographic_projection(self):
        l, r, b, t, n, f = sympy.symbols('l r b t n f')
        frustum_bounds = pm.FrustumBounds(l, r, b, t, n, f)
        ndc_bounds = pm.NDCBounds(-1, 1, -1, 1, -1, 1)
        expected = sympy.Matrix([
            [ 2 / (r - (-l)), 0,               0,           -(r + (-l)) / (r - (-l)) ],
            [ 0,              2 / (t - (-b)),  0,           -(t + (-b)) / (t - (-b)) ],
            [ 0,              0,              -2 / (f - n), -(f + n) / (f - n)       ],
            [ 0,              0,               0,            1                       ]
        ])
        x_lh_lh = sympy.Matrix.eye(4)
        x_rh_lh = change_of_orientation_rh_to_lh()
        m_coord = sympy.Matrix.eye(4)
        m_coord_inv = sympy.Matrix.eye(4)
        m_canonical_lh_lh = pm.orthographic(frustum_bounds, ndc_bounds)
        result = (x_lh_lh * m_coord_inv) * m_canonical_lh_lh * (m_coord * x_rh_lh)
        
        assert result.equals(expected)