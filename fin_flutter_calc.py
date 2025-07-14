# based on https://github.com/jkb-git/Fin-Flutter-Velocity-Calculator, which took inspiration from:
    # NACA TN 4197
    # https://www.apogeerockets.com/education/downloads/Newsletter291.pdf
    # https://www.apogeerockets.com/education/downloads/Newsletter411.pdf


# TODO take speed of sound calc out of the function and make it a separate function?

import numpy as np

def flutter_velocity_trapezoidal_fin(input_dict : dict):
    # what kinds of assumptions are made in this analysis?
        # from POF 291:
            # materials are assumed to be isotropic
        # from POF 615:
            # air is an ideal gas
            # dry air
                # but humid air only gives a slight difference (0.5% in warm air), and makes the actual flutter velocity higher, so making this assumption yields a very slightly conservative flutter velocity
            # no tip to tip, which would increase flutter velocity by a factor of √2
            # trapezoidal fin

    """ Returns the flutter velocity of a trapezoidal fin.

    Args
    ----
    input_dict : dict
        Dictionary containing the following key-value pairs:
            G (float): shear modulus, Pa
            T2T (bool): tip to tip, True if tip to tip, False if not
            c_r (float): root chord, unit of length
            c_t (float): tip chord, unit of length
            b (float): semi-span, unit of length
            m (float): fin sweep length, unit of length
            t (float): thickness, unit of length
            P (float): pressure, Pa
            T (float): temperature, K

    Returns
    -------
    float
        flutter velocity, m/s

    Notes
    -----
    Ensure to use the same unit for all length inputs, and use Pa for the shear modulus and air pressure. All length units cancel out except for the speed of sound's length unit, m, which is what final flutter velocity is returned as (m/s) regardless of the unit of length used for the inputs.

    This function also works for triangular fins by setting c_t = 0 (as a triangle can be considered a trapezoid with a tip chord of 0).
    """

    # Inputs
    G = input_dict['G'] # shear modulus, Pa
    T2T = input_dict['T2T'] # tip to tip, True if tip to tip, False if not

    c_r = input_dict['c_r'] # root chord, unit of length
    c_t = input_dict['c_t'] # tip chord, unit of length
    b = input_dict['b'] # semi-span aka height, unit of length
    m = input_dict['m'] # fin sweep length, unit of length
    t = input_dict['t'] # thickness, unit of length

    P = input_dict['P'] # pressure, Pa
    T = input_dict['T'] # temperature, K

    # Constants
    κ = 1.4 # adiabatic index aka ratio of specific heats for air, unitless

    # Calculations
    t_to_c_r = t / c_r # thickness ratio, unitless
    λ = c_t / c_r # taper ratio, unitless
    S = (c_r + c_t) * b / 2 # planform (fin) area, unit of area (unit of length squared)
    AR = b**2 / S # aspect ratio, unitless
    C_x = (2 * c_t * m + c_t ** 2 + m * c_r + c_r * c_t + c_r ** 2)/(3 * (c_t + c_r)) # location of the centroid of the fin in the axis along the fin's chord, unit of length
    ε = C_x / c_r - 0.25 # distance of the centroid behind the quarter chord, unitless

    denom_const = 24 * ε / np.pi * (λ + 1)/2 * (AR ** 3 / (t_to_c_r ** 3 * (AR + 2))) # values of the denominator inside the radical that depend on fin geometry, unitless
    fin_const = G / denom_const # values in the radical that depend on the fin (geometry and material), Pa

    if T2T:
        fin_const *= 2 # tip to tip, which would increase flutter velocity by a factor of √2

    a = np.sqrt(κ * T * 8.3144598 / 0.0289644) # speed of sound, m/s

    return a * np.sqrt(fin_const / (P * κ)) # flutter velocity, m/s

# https://github.com/jkb-git/Fin-Flutter-Velocity-Calculator/blob/main/Calculating_Fin_Flutter_Velocity_For_Complex_Fin_Shapes.pdf

def flutter_velocity_polygon_fin(geometry : list, conditions : dict, verbose = False):
    """ Returns the flutter velocity of a fin defined by a series of points.

    Args
    ----
    geometry : list
        List of tuples containing the x and y coordinates of the fin's points. Order of points matters, as the points are connected in the order they are listed. y = 0 is the outer diameter of the rocket.
    conditions : dict
        Dictionary containing the following key-value pairs:
            G (float): shear modulus, Pa
            T2T (bool): tip to tip, True if tip to tip reinforcement, False if not
            t (float): thickness, unit of length
            P (float): pressure, Pa
            T (float): temperature, K
    verbose : bool, optional
        True to display/print a plot of the fin, the results of intermediate calculations, and the final flutter velocity. False to not display/print anything. Default is False.
    """
    import triangle as tr
    verts = np.array(geometry) # vertices of the fin
    segs = []                       # must specify segments as well as verts to create a PSLG 
    for i in range (len(verts)):    # create segments of a closed PSLG from list of vertices 
        seg = []    # create a new empty segment 
        seg.append(i)   # add the current vertex index as the start of this segment 
        if (i == (len(verts)-1)): # if we are on the last vertex, loop back to vertex zero 
            seg.append(0)  # to close the polygon 
        else:    # otherwise, make the next vertex end this segment 
            seg.append(i+1) 
        segs.append(seg)  # add the newly created segment to our list of segments 
    tr_input = dict(vertices=verts, segments=segs)  # triangle expects Python dict as input 
    tr_output = tr.triangulate(tr_input,'p')        # make triangles. 'p' = input is a PSLG 
    tris = tr_output['triangles'].tolist()          # convert output of triangle to Python list
    if verbose: # compare the input shape to the output triangles
        import matplotlib.pyplot as plt
        tr.compare(plt, tr_input, tr_output)
        plt.show()

    # For each triangle, calculate the area and the Cx
    tri_areas = []      # this list will hold the area of each triangle
    tri_Cxs = []        # this list will hold the Cx of each triangle
    fin_area = 0
    fin_Cx = 0




if __name__ == '__main__':
    if 0: # example in POF 615
        print(flutter_velocity_trapezoidal_fin({
            'G': 4136854000, # shear modulus, Pa = 600000 psi
            'T2T': False, # tip to tip, True if tip to tip, False if not
            'c_r': 7.5, # root chord, in
            'c_t': 2.5, # tip chord, in
            'b': 3, # semi-span, in
            'm': 4.285, # fin sweep length, in
            't': 0.125, # thickness, in
            'P': 49633, # pressure, Pa = 7.1986 psi
            'T': 251.56 # temperature, K = -6.86 F
        }) * 3.28084) # flutter velocity, ft/s
    if 1:
        flutter_velocity_polygon_fin([(0,0), (2,3), (4,5), (7,7), (1,0)], {
            'G': 4136854000, # shear modulus, Pa = 600000 psi
            'T2T': False, # tip to tip, True if tip to tip, False if not
            't': 0.125, # thickness, in
            'P': 49633, # pressure, Pa = 7.1986 psi
            'T': 251.56 # temperature, K = -6.86 F
        }, verbose=True)