# my file name is viewjacspline.py
import sys

import gustav as gus
import splinepy as sp

print(sys.argv)

spls = sp.load_splines(sys.argv[1])
gs = [gus.NURBS(**a.todict()) for a in spls]

print(gs)

gus.show.show_vedo(gs)
