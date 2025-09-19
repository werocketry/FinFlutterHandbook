"""
Fin Flutter Velocity Calculator
Based on: 
- NACA TN 4197
- https://www.apogeerockets.com/education/downloads/Newsletter291.pdf
- https://www.apogeerockets.com/education/downloads/Newsletter411.pdf
- https://github.com/jkb-git/Fin-Flutter-Velocity-Calculator
"""

import numpy as np

# ------------------------------
# Helper function
# ------------------------------
def speed_of_sound(T: float, κ: float = 1.4, R: float = 8.3144598, M: float = 0.0289644) -> float:
    """Returns speed of sound in air at temperature T (Kelvin)."""
    return np.sqrt(κ * R * T / M)


# ------------------------------
# Flutter velocity functions
# ------------------------------
def flutter_velocity_trapezoidal_fin(input_dict: dict) -> float:
    """Returns the flutter velocity of a trapezoidal fin."""
    G = input_dict['G']
    T2T = input_dict['T2T']
    c_r = input_dict['c_r']
    c_t = input_dict['c_t']
    b = input_dict['b']
    m = input_dict['m']
    t = input_dict['t']
    P = input_dict['P']
    T = input_dict['T']

    κ = 1.4

    t_to_c_r = t / c_r
    λ = c_t / c_r
    S = (c_r + c_t) * b / 2
    AR = b**2 / S
    C_x = (2 * c_t * m + c_t**2 + m * c_r + c_r * c_t + c_r**2) / (3 * (c_t + c_r))
    ε = C_x / c_r - 0.25

    denom_const = 24 * ε / np.pi * (λ + 1)/2 * (AR**3 / (t_to_c_r**3 * (AR + 2)))
    fin_const = G / denom_const

    if T2T:
        fin_const *= 2

    a = speed_of_sound(T, κ=κ)
    return a * np.sqrt(fin_const / (P * κ))


def flutter_velocity_polygon_fin(geometry: list, conditions: dict, verbose=False) -> float:
    """Returns flutter velocity of a fin defined by polygon points."""
    import triangle as tr
    verts = np.array(geometry)
    segs = [[i, (i + 1) % len(verts)] for i in range(len(verts))]
    tr_input = dict(vertices=verts, segments=segs)
    tr_output = tr.triangulate(tr_input, 'p')
    tris = tr_output['triangles'].tolist()

    # Compute areas and centroid x positions
    tri_areas = []
    tri_Cxs = []
    for tri in tris:
        pts = verts[tri]
        area = 0.5 * abs(
            pts[0,0]*(pts[1,1]-pts[2,1]) +
            pts[1,0]*(pts[2,1]-pts[0,1]) +
            pts[2,0]*(pts[0,1]-pts[1,1])
        )
        cx = np.mean(pts[:,0])
        tri_areas.append(area)
        tri_Cxs.append(cx)

    fin_area = np.sum(tri_areas)
    C_x = np.sum(np.array(tri_areas) * np.array(tri_Cxs)) / fin_area
    b = np.max(verts[:,1])
    c_r = np.max(verts[:,0]) - np.min(verts[:,0])

    # Inputs
    G = conditions['G']
    T2T = conditions['T2T']
    t = conditions['t']
    P = conditions['P']
    T = conditions['T']

    # Aerodynamic terms
    κ = 1.4
    t_to_c_r = t / c_r
    AR = b**2 / fin_area
    λ = (fin_area / b) / c_r
    ε = C_x / c_r - 0.25

    denom_const = 24 * ε / np.pi * (λ + 1)/2 * (AR**3 / (t_to_c_r**3 * (AR + 2)))
    fin_const = G / denom_const

    if T2T:
        fin_const *= 2

    a = speed_of_sound(T, κ=κ)
    Vf = a * np.sqrt(fin_const / (P * κ))

    if verbose:
        import matplotlib.pyplot as plt
        tr.compare(plt, tr_input, tr_output)
        plt.show()
        print(f"Fin area: {fin_area:.3f}")
        print(f"Aspect ratio (AR): {AR:.3f}")
        print(f"Centroid offset (ε): {ε:.3f}")
        print(f"Flutter velocity: {Vf:.2f} m/s")

    return Vf


def flutter_velocity_elliptical_fin(input_dict: dict) -> float:
    """Returns flutter velocity of an elliptical fin."""
    G = input_dict['G']
    T2T = input_dict['T2T']
    c_r = input_dict['c_r']
    b = input_dict['b']
    m = input_dict['m']
    t = input_dict['t']
    P = input_dict['P']
    T = input_dict['T']

    κ = 1.4
    S = ((c_r / 2) * np.pi) * b / 2
    bc = (S / b * 2) - c_r
    t_to_c_r = t / c_r
    λ = bc / c_r
    AR = b**2 / S
    C_x = (2 * bc * m + bc**2 + m * c_r + c_r * bc + c_r**2) / (3 * (bc + c_r))
    ε = C_x / c_r - 0.25

    denom_const = 24 * ε / np.pi * (λ + 1)/2 * (AR**3 / (t_to_c_r**3 * (AR + 2)))
    fin_const = G / denom_const

    a = speed_of_sound(T, κ=κ)
    return a * np.sqrt(fin_const / (P * κ))


# ------------------------------
# Main script
# ------------------------------
if __name__ == '__main__':
    # Trapezoidal fin example
    trapezoid_result = flutter_velocity_trapezoidal_fin({
        'G': 4136854000,
        'T2T': False,
        'c_r': 7.5,
        'c_t': 2.5,
        'b': 3,
        'm': 4.285,
        't': 0.125,
        'P': 49633,
        'T': 251.56
    }) * 3.28084  # ft/s

    # Polygon fin example
    polygon_result = flutter_velocity_polygon_fin(
        [(0, 0), (2, 3), (4, 5), (7, 7), (1, 0)],
        {
            'G': 4136854000,
            'T2T': False,
            't': 0.125,
            'P': 49633,
            'T': 251.56
        },
        verbose=True
    )

    # Elliptical fin example
    elliptical_result = flutter_velocity_elliptical_fin({
        'G': 4136854000,
        'T2T': False,
        'c_r': 7.5,
        'b': 3,
        'm': 4.285,
        't': 0.125,
        'P': 49633,
        'T': 251.56
    })

    # Print results
    print("\n--- Flutter Velocity Results ---")
    print(f"Trapezoidal fin: {trapezoid_result:.2f} ft/s")
    print(f"Polygon fin:     {polygon_result:.2f} m/s")
    print(f"Elliptical fin:  {elliptical_result:.2f} m/s")
