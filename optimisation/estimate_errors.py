import numpy
import lsystem_temp as lsys
import calc_growth_space as cgs

def estimate_error_1(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    # Generate L-system points for params
    ls_pts = lsys.generate_lsystem_tree_points(params)
    ls_pts = numpy.asarray(ls_pts)
    ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = cgs.cgs_ls(ls_pts)
    # Calculate error
    Ebbox = numpy.sqrt((((ls_bbox[1]-ls_bbox[0])-(targ_bbox[1]-targ_bbox[0]))/(targ_bbox[1]-targ_bbox[0]))**2 + \
        (((ls_bbox[3]-ls_bbox[2])-(targ_bbox[3]-targ_bbox[2]))/(targ_bbox[3]-targ_bbox[2]))**2) #BBox error
    Eth = numpy.sqrt(((ls_th-targ_th)/targ_th)**2) # Trunk height error
    Ech = numpy.sqrt(((ls_ch-targ_ch)/targ_ch)**2) # Crown height error
    #Eel_theta = (ls_el[0]-targ_el[0])/(2*numpy.pi)
    Eel_radii = numpy.mean(((ls_el[1:2]-targ_el[1:2])/(targ_el[1:2]))**2)
    #Eeu_theta = (ls_eu[0]-targ_eu[0])/(2*numpy.pi)
    Eeu_radii = numpy.mean(((ls_eu[1:2]-targ_eu[1:2])/(targ_eu[1:2]))**2)
    E_total = numpy.sum(numpy.array(\
        [Eth,Ech,Eel_radii,Eeu_radii]))
    #print params, E_total
    return E_total

def estimate_error_2(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    # Generate L-system points for params
    ls_pts = lsys.generate_lsystem_tree_points(params)
    ls_pts = numpy.asarray(ls_pts)
    ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = cgs.cgs_ls(ls_pts)
    # Calculate error
    Ebbox = numpy.sqrt((((ls_bbox[1]-ls_bbox[0])-(targ_bbox[1]-targ_bbox[0]))/(targ_bbox[1]-targ_bbox[0]))**2 + \
        (((ls_bbox[3]-ls_bbox[2])-(targ_bbox[3]-targ_bbox[2]))/(targ_bbox[3]-targ_bbox[2]))**2) #BBox error
    Eth = numpy.sqrt(((ls_th-targ_th)/targ_th)**2) # Trunk height error
    Ech = numpy.sqrt(((ls_ch-targ_ch)/targ_ch)**2) # Crown height error
    #Eel_theta = (ls_el[0]-targ_el[0])/(2*numpy.pi)
    Eel_radii = numpy.mean(((ls_el[1:2]-targ_el[1:2])/(targ_el[1:2]))**2)
    #Eeu_theta = (ls_eu[0]-targ_eu[0])/(2*numpy.pi)
    Eeu_radii = numpy.mean(((ls_eu[1:2]-targ_eu[1:2])/(targ_eu[1:2]))**2)

    E_mu = numpy.mean((ls_mu-targ_mu)**2)

    E_total = numpy.sum(numpy.array(\
        [Eth,Ech,Eel_radii,Eeu_radii,E_mu]))
    #print params, E_total
    return E_total
