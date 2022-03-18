from optparse import OptionParser
import gustav as gus
import splinepy as sp


# Quick parse options
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="filename",\
                  help="File of the input name",type="string", default=None)

parser.add_option("-k", "--knots",\
                  action="store_true", dest="show_knots", default=False,\
                  help="Show Knotlines on spline")
parser.add_option("-c", "--controlpoints",\
                  action="store_true", dest="show_ctps", default=False,\
                  help="Show Control Points")
parser.add_option("-m", "--controlpointmesh",\
                  action="store_true", dest="show_ctps_mesh", default=False,\
                  help="Show Control Points Mesh")
parser.add_option("-s", "--Spline",\
                  action="store_true", dest="show_spline", default=False,\
                  help="Show spline")

(options, args) = parser.parse_args()

if options.filename == None:
  print("No filename specified. , set with -i <filename>")
  quit()

spls = sp.load_splines(options.filename)
gs = [gus.NURBS(**a.todict()) for a in spls]

# Add objects
gsList = [g.show(return_discrete=True) for g in gs]
gsObjects = []
if options.show_spline :
  gsObjects += [iGsList["spline"] for iGsList in gsList]
if options.show_ctps :
  gsObjects += [iGsList["control_points"] for iGsList in gsList]
if options.show_knots:
  gsObjects += [iGsList["knots"] for iGsList in gsList]
if options.show_ctps_mesh :
  gsObjects += [iGsList["control_mesh"] for iGsList in gsList]
if len(gsObjects)==0:
  print("no option chosen - see -h for help")

gus.show.show_vedo(gsObjects)
