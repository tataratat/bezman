import gustav as gus
import splinepy as sp
import numpy as np
import vedo
import sys

if __name__ == "__main__":
    # load. I know, there should be a better way
    s = gus.NURBS(**sp.load_splines(sys.argv[1])[0].todict())
    ds = gus.NURBS(**sp.load_splines(sys.argv[2])[0].todict())

    # derivative and its plot origin sample
    sam_res = [10, 10]
    origin = s.sample(sam_res)
    der = ds.sample(sam_res)
    length = der.view()

    # scale w.r.t. spline cp bounds
    bounds = ds.control_point_bounds
    b_diagonal = bounds[1] - bounds[0]
    length *= b_diagonal / 10

    arrows = vedo.shapes.Arrows(origin, origin + length, c="r")
    # arrows = vedo.shapes.Arrows(origin, origin + length)

    # gus.show.show_vedo([s, arrows])
    # gus.show.show_vedo([s, arrows])
    gusobjs = s.show(return_discrete=True)
    gus.show.show_vedo([gusobjs["spline"], gusobjs["knots"], arrows])
