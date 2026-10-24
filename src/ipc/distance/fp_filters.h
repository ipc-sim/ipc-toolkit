#pragma once
/* Automatically generated code, do not edit the functions! */

/*
The first filter evaluates the sign of this expression:
dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p0)) =
= dot(p1-p0, p3-p0) * dot(p1-p0, p2-p0) - dot(p1-p0, p1-p0) * dot(p2-p0, p3-p0)

Explanation:
let t0,t1,t2,p be 3D points.
e0 := t1 - t0;
e1 := t2 - t1;
e2 := t0 - t2;
r0 := p - t0;
r1 := p - t1;
r2 := p - t2;
n := cross(e0, -e1) = cross(e1, -e2) = cross(e2, -e0) up to rescaling

for any i=0,1,2, call j=i+1 mod 3 and k=i+2 mod 3.

dot(n, cross(ri, ei))
= dot(cross(ej, -ek), cross(ri, ei))
= dot(cross(tj-ti, tk-ti), cross(p-ti, tj-ti))
= dot(tj-ti, p-ti) * dot(tk-ti, tj-ti) - dot(tj-ti, tj-ti) * dot(tk-ti, p-ti)
so we evaluate the predicate with
p0 = ti
p1 = tj
p2 = tk
p3 = p
*/

constexpr int FPG_UNCERTAIN_VALUE = 0;

inline int cross_dot_cross_1_3d_filter(
    const double* p0, const double* p1, const double* p2, const double* p3)
{
    double d1_0;
    d1_0 = (p1[0] - p0[0]);
    double d1_1;
    d1_1 = (p1[1] - p0[1]);
    double d1_2;
    d1_2 = (p1[2] - p0[2]);
    double d2_0;
    d2_0 = (p2[0] - p0[0]);
    double d2_1;
    d2_1 = (p2[1] - p0[1]);
    double d2_2;
    d2_2 = (p2[2] - p0[2]);
    double d3_0;
    d3_0 = (p3[0] - p0[0]);
    double d3_1;
    d3_1 = (p3[1] - p0[1]);
    double d3_2;
    d3_2 = (p3[2] - p0[2]);
    double m11;
    m11 = (((d1_0 * d1_0) + (d1_1 * d1_1)) + (d1_2 * d1_2));
    double m12;
    m12 = (((d1_0 * d2_0) + (d1_1 * d2_1)) + (d1_2 * d2_2));
    double m13;
    m13 = (((d1_0 * d3_0) + (d1_1 * d3_1)) + (d1_2 * d3_2));
    double m23;
    m23 = (((d2_0 * d3_0) + (d2_1 * d3_1)) + (d2_2 * d3_2));
    double r;
    r = ((m12 * m13) - (m11 * m23));
    int int_tmp_result;
    double eps;
    double max1 = fabs(d1_0);
    if ((max1 < fabs(d1_1))) {
        max1 = fabs(d1_1);
    }
    if ((max1 < fabs(d1_2))) {
        max1 = fabs(d1_2);
    }
    double max2 = fabs(d2_0);
    if ((max2 < fabs(d2_1))) {
        max2 = fabs(d2_1);
    }
    if ((max2 < fabs(d2_2))) {
        max2 = fabs(d2_2);
    }
    double max3 = fabs(d3_0);
    if ((max3 < fabs(d3_1))) {
        max3 = fabs(d3_1);
    }
    if ((max3 < fabs(d3_2))) {
        max3 = fabs(d3_2);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if ((max2 < lower_bound_1)) {
        lower_bound_1 = max2;
    } else {
        if ((max2 > upper_bound_1)) {
            upper_bound_1 = max2;
        }
    }
    if ((max3 < lower_bound_1)) {
        lower_bound_1 = max3;
    } else {
        if ((max3 > upper_bound_1)) {
            upper_bound_1 = max3;
        }
    }
    if ((lower_bound_1 < 3.14773426688569445494e-74)) {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if ((upper_bound_1 > 7.23700557733225900010e+75)) {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.26648152760393650857e-14 * (((max1 * max2) * max1) * max3));
        if ((r > eps)) {
            int_tmp_result = 1;
        } else {
            if ((r < -eps)) {
                int_tmp_result = -1;
            } else {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

/*
The second filter evaluates the sign of this expression:
dot(cross(p1-p0, p2-p0), cross(p3-p0, p1-p2)) =
= dot(p1-p0, p3-p0) * dot(p1-p2, p2-p0) - dot(p1-p0, p1-p2) * dot(p2-p0, p3-p0)
Explanation:
We separately determine if the closest point on A to B is at endpoint
0, endpoint 1, or the interior. Given one of the endpoints, we can
assess whether it is a local minimum of distance by checking if the other
endpoint is behind the plane defined by the closest direction from the first
endpoint to segment B, that is, the gradient of the distance function.
This predicate can be tested with dot_3. If neither endpoint of A is the
local minimum, then the minimum must be in the interior of A.
To determine the gradient of the distance function at a point, it is enough to
determine the point-edge distance type for that point, for which we have
another function.
If the distance type is point-point, the gradient points in the direction
connecting the two points. Otherwise, the gradient points in the direction
connecting the endpoint to the interior of B; up to changes in magnitude,
this is g=cross(cross(b1-ai, b0-ai), b1-b0). We must compute dot(g, aj-ai),
which is dot(cross(cross(b1-ai, b0-ai), b1-b0), aj-ai) =
= dot(cross(b1-ai, b0-ai), cross(b1-b0, aj-ai)) =
= dot(cross(b1-ai, b0-ai), cross(aj-ai, b0-b1)) =
= dot(b1-ai, aj-ai) * dot(b0-ai, b0-b1) - dot(b1-ai, b0-b1) * dot(b0-ai, aj-ai)
*/
inline int cross_dot_cross_2_3d_filter(
    const double* p0, const double* p1, const double* p2, const double* p3)
{
    double d1_0;
    d1_0 = (p1[0] - p0[0]);
    double d1_1;
    d1_1 = (p1[1] - p0[1]);
    double d1_2;
    d1_2 = (p1[2] - p0[2]);
    double d2_0;
    d2_0 = (p2[0] - p0[0]);
    double d2_1;
    d2_1 = (p2[1] - p0[1]);
    double d2_2;
    d2_2 = (p2[2] - p0[2]);
    double d3_0;
    d3_0 = (p3[0] - p0[0]);
    double d3_1;
    d3_1 = (p3[1] - p0[1]);
    double d3_2;
    d3_2 = (p3[2] - p0[2]);
    double dX_0;
    dX_0 = (p1[0] - p2[0]);
    double dX_1;
    dX_1 = (p1[1] - p2[1]);
    double dX_2;
    dX_2 = (p1[2] - p2[2]);
    double mX1;
    mX1 = (((dX_0 * d1_0) + (dX_1 * d1_1)) + (dX_2 * d1_2));
    double mX2;
    mX2 = (((dX_0 * d2_0) + (dX_1 * d2_1)) + (dX_2 * d2_2));
    double m13;
    m13 = (((d1_0 * d3_0) + (d1_1 * d3_1)) + (d1_2 * d3_2));
    double m23;
    m23 = (((d2_0 * d3_0) + (d2_1 * d3_1)) + (d2_2 * d3_2));
    double r;
    r = ((mX2 * m13) - (mX1 * m23));
    int int_tmp_result;
    double eps;
    double max1 = fabs(d1_0);
    if ((max1 < fabs(d1_1))) {
        max1 = fabs(d1_1);
    }
    if ((max1 < fabs(d1_2))) {
        max1 = fabs(d1_2);
    }
    double max2 = fabs(d2_0);
    if ((max2 < fabs(d2_1))) {
        max2 = fabs(d2_1);
    }
    if ((max2 < fabs(d2_2))) {
        max2 = fabs(d2_2);
    }
    double max3 = fabs(d3_0);
    if ((max3 < fabs(d3_1))) {
        max3 = fabs(d3_1);
    }
    if ((max3 < fabs(d3_2))) {
        max3 = fabs(d3_2);
    }
    double max4 = fabs(dX_0);
    if ((max4 < fabs(dX_1))) {
        max4 = fabs(dX_1);
    }
    if ((max4 < fabs(dX_2))) {
        max4 = fabs(dX_2);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if ((max2 < lower_bound_1)) {
        lower_bound_1 = max2;
    } else {
        if ((max2 > upper_bound_1)) {
            upper_bound_1 = max2;
        }
    }
    if ((max3 < lower_bound_1)) {
        lower_bound_1 = max3;
    } else {
        if ((max3 > upper_bound_1)) {
            upper_bound_1 = max3;
        }
    }
    if ((max4 < lower_bound_1)) {
        lower_bound_1 = max4;
    } else {
        if ((max4 > upper_bound_1)) {
            upper_bound_1 = max4;
        }
    }
    if ((lower_bound_1 < 3.14773426688569445494e-74)) {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if ((upper_bound_1 > 7.23700557733225900010e+75)) {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (2.26648152760393650857e-14 * (((max4 * max2) * max1) * max3));
        if ((r > eps)) {
            int_tmp_result = 1;
        } else {
            if ((r < -eps)) {
                int_tmp_result = -1;
            } else {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

inline int dot3_2d_filter(const double* p0, const double* p1, const double* p2)
{
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double Delta;
    Delta = ((a11 * a21) + (a12 * a22));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if ((max1 < fabs(a12))) {
        max1 = fabs(a12);
    }
    double max2 = fabs(a21);
    if ((max2 < fabs(a22))) {
        max2 = fabs(a22);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if ((max2 < lower_bound_1)) {
        lower_bound_1 = max2;
    } else {
        if ((max2 > upper_bound_1)) {
            upper_bound_1 = max2;
        }
    }
    if ((lower_bound_1 < 5.00368081960964802120e-147)) {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if ((upper_bound_1 > 1.67597599124282389316e+153)) {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (8.88720573725927779595e-16 * (max1 * max2));
        if ((Delta > eps)) {
            int_tmp_result = 1;
        } else {
            if ((Delta < -eps)) {
                int_tmp_result = -1;
            } else {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

inline int dot3_3d_filter(const double* p0, const double* p1, const double* p2)
{
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double Delta;
    Delta = (((a11 * a21) + (a12 * a22)) + (a13 * a23));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if ((max1 < fabs(a12))) {
        max1 = fabs(a12);
    }
    if ((max1 < fabs(a13))) {
        max1 = fabs(a13);
    }
    double max2 = fabs(a21);
    if ((max2 < fabs(a22))) {
        max2 = fabs(a22);
    }
    if ((max2 < fabs(a23))) {
        max2 = fabs(a23);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if ((max2 < lower_bound_1)) {
        lower_bound_1 = max2;
    } else {
        if ((max2 > upper_bound_1)) {
            upper_bound_1 = max2;
        }
    }
    if ((lower_bound_1 < 3.78232824369468580207e-147)) {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if ((upper_bound_1 > 1.67597599124282389316e+153)) {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.55534235888797938037e-15 * (max1 * max2));
        if ((Delta > eps)) {
            int_tmp_result = 1;
        } else {
            if ((Delta < -eps)) {
                int_tmp_result = -1;
            } else {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}

inline int cross_null_3d_filter(
    const double* p0, const double* p1, const double* q0, const double* q1)
{
    double v_0;
    v_0 = (p1[0] - p0[0]);
    double v_1;
    v_1 = (p1[1] - p0[1]);
    double v_2;
    v_2 = (p1[2] - p0[2]);
    double w_0;
    w_0 = (q1[0] - q0[0]);
    double w_1;
    w_1 = (q1[1] - q0[1]);
    double w_2;
    w_2 = (q1[2] - q0[2]);
    double c_i;
    c_i = ((v_1 * w_2) - (v_2 * w_1));
    double c_j;
    c_j = ((v_2 * w_0) - (v_0 * w_2));
    double c_k;
    c_k = ((v_0 * w_1) - (v_1 * w_0));
    double cross;
    cross = (((c_i * c_i) + (c_j * c_j)) + (c_k * c_k));
    int int_tmp_result;
    double eps;
    double max1 = fabs(v_0);
    if ((max1 < fabs(v_1))) {
        max1 = fabs(v_1);
    }
    if ((max1 < fabs(v_2))) {
        max1 = fabs(v_2);
    }
    double max2 = fabs(w_0);
    if ((max2 < fabs(w_1))) {
        max2 = fabs(w_1);
    }
    if ((max2 < fabs(w_2))) {
        max2 = fabs(w_2);
    }
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if ((max2 < lower_bound_1)) {
        lower_bound_1 = max2;
    } else {
        if ((max2 > upper_bound_1)) {
            upper_bound_1 = max2;
        }
    }
    if ((lower_bound_1 < 3.53675409555095779245e-74)) {
        return FPG_UNCERTAIN_VALUE;
    } else {
        if ((upper_bound_1 > 1.44740111546645180002e+76)) {
            return FPG_UNCERTAIN_VALUE;
        }
        eps = (1.42208308574965596077e-14 * (((max1 * max2) * max1) * max2));
        if ((cross > eps)) {
            int_tmp_result = 1;
        } else {
            if ((cross < -eps)) {
                int_tmp_result = -1;
            } else {
                return FPG_UNCERTAIN_VALUE;
            }
        }
    }
    return int_tmp_result;
}