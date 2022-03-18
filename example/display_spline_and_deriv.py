import gustav as gus
import splinepy as sp
import numpy as np
import vedo
import sys

if __name__ == "__main__":

    # Read file names
    raw_spline_list = sp.load_splines(sys.argv[1])
    raw_spline_deriv_list = sp.load_splines(sys.argv[2])
    # load. I know, there should be a better way
    spline_list = [gus.NURBS(**a.todict()) for a in raw_spline_list]
    spline_derivative_list = [gus.NURBS(**a.todict()) for a in raw_spline_deriv_list]

    # derivative and its plot origin sample
    sam_res = [10, 10]
    origin = spline_list[0].sample(sam_res)
    der = spline_derivative_list[0].sample(sam_res)

    for i in range(1, len(spline_list)):
      origin = np.vstack((origin,spline_list[i].sample(sam_res)))
      der = np.vstack((der,spline_derivative_list[i].sample(sam_res)))
    length = der.view()

    # scale w.r.t. spline cp bounds
    max_length = np.max(np.linalg.norm(der, axis=1))
    min_length = np.min(np.linalg.norm(der, axis=1))
    b_diagonal = max_length - min_length
    length *= b_diagonal / 10

    arrows = vedo.shapes.Arrows(origin, origin + length, c="r")
    # arrows = vedo.shapes.Arrows(origin, origin + length)

    # gus.show.show_vedo([spline, arrows])
    # gus.show.show_vedo([spline, arrows])
    gusSplines = [a.show(return_discrete=True)["spline"] for a in spline_list]
    gusKnots = [a.show(return_discrete=True)["knots"] for a in spline_list]
    plt = gus.show.show_vedo([*gusSplines,*gusKnots, arrows])
