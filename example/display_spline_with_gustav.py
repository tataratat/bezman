from optparse import OptionParser, OptionGroup
import gustav as gus
import vedo
import splinepy as sp
import numpy as np

if __name__ == "__main__":
  # Quick parse options
  parser = OptionParser()
  parser.add_option("-i", "--inputfile", dest="filename",\
                    help="File of the input name",type="string", default=None)
  parser.add_option("-r", "--sampleresolution", dest="q_res",\
                    help="Samples per parametric dimension",type="int", default=None)
  
  group = OptionGroup(parser, "Plotting Options")
  group.add_option("-k", "--knots",\
                    action="store_true", dest="show_knots", default=False,\
                    help="Show Knotlines on spline")
  group.add_option("-c", "--controlpoints",\
                    action="store_true", dest="show_ctps", default=False,\
                    help="Show Control Points")
  group.add_option("-m", "--controlpointmesh",\
                    action="store_true", dest="show_ctps_mesh", default=False,\
                    help="Show Control Points Mesh")
  group.add_option("-s", "--Spline",\
                    action="store_true", dest="show_spline", default=False,\
                    help="Show spline")
  parser.add_option_group(group)

  group = OptionGroup(parser, "Derivative Options", "If the Derivative is calculated with "
                            "together with the spline, an additional file can be provided")
  group.add_option("-d", "--derivativefile", dest="derivfilename",\
                    help="File (name) that stores derivative",type="string", default=None)
  group.add_option("-D", "--derivativesamples", dest="derv_q_res",\
                    help="Samples for the spline derivative", type="int", default=10)
  group.add_option("-A", "--arrowlength", dest="max_length_arrows",\
                    help="Maximum Length of the derivative arrows", type="float", default=1.)
  parser.add_option_group(group)
  
  (options, args) = parser.parse_args()

  if options.filename == None:
    print("No filename specified. , set with -i <filename>")
    quit()
  
  spls = sp.load_splines(options.filename)
  gs = [gus.NURBS(**a.todict()) for a in spls]
  spls = None
  
  if options.q_res == None:
    if len(gs[0].degrees) == 3 :
      options.q_res = 10
    else :
      options.q_res = 100
  
    print("Using default resolution : " + str(options.q_res))
  
  # Add objects
  gsList = [g.show(return_discrete=True, resolutions=options.q_res) for g in gs]
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

  # Now check the derivative if required
  if options.derivfilename != None:
    print("Start Sampling and reading the derivative file")
    spls = sp.load_splines(options.derivfilename)
    gsDeriv = [gus.NURBS(**a.todict()) for a in spls]

    # Free memory
    spls = None

    # Retrieve parametric dimension
    par_dim = len(gs[0].degrees)

    # Transform points
    sam_res = [options.derv_q_res] * par_dim
    origin = gs[0].sample(sam_res)
    der = gsDeriv[0].sample(sam_res)
    for i in range(1, len(gs)):
      origin = np.vstack((origin, gs[i].sample(sam_res)))
      der = np.vstack((der, gsDeriv[i].sample(sam_res)))
    length = der.view()

    # scale w.r.t. spline cp bounds
    max_length = np.max(np.linalg.norm(der, axis=1))
    length /= max_length / options.max_length_arrows

    # Add arrows to the list
    gsObjects += [vedo.shapes.Arrows(origin, origin + length, c="r")]
    


  
  gus.show.show_vedo(gsObjects)

