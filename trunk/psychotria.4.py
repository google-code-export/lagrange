import sys
sys.path.append("..")
from pprint import pprint
from model09 import *
import scipy
from scipy import optimize
LARGE = 10e5
PMAX = 10.0

## (1) Kauai/Niihau
## (2) Oahu
## (3) Molokai, Lanai, Maui, and Kahoolawe (together Maui Nui)
## (4) Hawaii

data = {
    "P_mariniana_Kokee2":        (1,0,0,0),
##  combine following two as P_mariniana_Oahu
##     "P_mariniana_WN2":           (0,1,0,0),
##     "P_mariniana_Manoa2":        (0,1,0,0),
    "P_mariniana_Oahu":          (0,1,0,0),
## combine these as sister to hawaiiensis
##     "P_mariniana_Kipahulu2":     (0,0,1,0),
##     "P_mariniana_Kaw2":          (0,0,1,0),
##     "P_mariniana_Lanai2":        (0,0,1,0),
    "P_mariniana_MauiNui":       (0,0,1,0),
    "P_hawaiiensis_Makaopuhi":   (0,0,0,1),
    "P_wawraeDL7428":            (1,0,0,0),
## group following two as sisters
    "P_kaduana_PuuKukuiAS":      (0,0,1,0),
    "P_mauiensis_PepeAS":        (0,0,1,0),
    "P_hawaiiensis_WaikamoiL1":  (0,0,1,0),
    "P_mauiensis_Eke":           (0,0,1,0),
    "P_fauriei2":                (0,1,0,0),
    "P_hathewayi_1":             (0,1,0,0),
    "P_kaduana_HawaiiLoa":       (0,1,0,0),
## combine these
##     "P_greenwelliae07":          (1,0,0,0),
##     "P_greenwelliae907":         (1,0,0,0),
    "P_greenwelliae_Kauai":      (1,0,0,0),
    "P_grandiflora_Kal2":        (1,0,0,0),
    "P_hobdyi_Kuia":             (1,0,0,0),
    "P_hexandra_K1":             (1,0,0,0),
    "P_hexandra_M":              (1,0,0,0),
    "P_hexandra_Oahu":           (0,1,0,0),
    }

#newickstr = "(((((P_mariniana_Kokee2:0.033539,(P_mariniana_Oahu:0.030646,(P_mariniana_MauiNui:0.02241,P_hawaiiensis_Makaopuhi:0.02241):0.008236):0.002893):0.005171,P_wawraeDL7428:0.03871):0.008237,((((P_kaduana_PuuKukuiAS:0.020803,P_mauiensis_PepeAS:0.020803):0.00001,((P_hawaiiensis_WaikamoiL1:0.010853,P_mauiensis_Eke:0.010853):0.007964,(P_fauriei2:0.013826,P_hathewayi_1:0.013826):0.004991):0.001986):0.003762,P_kaduana_HawaiiLoa:0.024565):0.003398,P_greenwelliae_Kauai:0.027963):0.018984):0.008255,(P_grandiflora_Kal2:0.027864,P_hobdyi_Kuia:0.027864):0.027338):0.003229,((P_hexandra_K1:0.026568,P_hexandra_M:0.026568):0.005204,P_hexandra_Oahu:0.031771):0.026659);"

# branches rotated to match figures in Nepokroeff paper
newickstr = "((((((((P_hawaiiensis_WaikamoiL1:0.010853,P_mauiensis_Eke:0.010853):0.007964,(P_fauriei2:0.013826,P_hathewayi_1:0.013826):0.004991):0.001986,(P_kaduana_PuuKukuiAS:0.020803,P_mauiensis_PepeAS:0.020803):1e-05):0.003762,P_kaduana_HawaiiLoa:0.024565):0.003398,P_greenwelliae_Kauai:0.027963):0.018984,((((P_mariniana_MauiNui:0.02241,P_hawaiiensis_Makaopuhi:0.02241):0.008236,P_mariniana_Oahu:0.030646):0.002893,P_mariniana_Kokee2:0.033539):0.005171,P_wawraeDL7428:0.03871):0.008237):0.008255,(P_grandiflora_Kal2:0.027864,P_hobdyi_Kuia:0.027864):0.027338):0.003229,((P_hexandra_K1:0.026568,P_hexandra_M:0.026568):0.005204,P_hexandra_Oahu:0.031771):0.026659);"

periods = [10.0]
model = RateModelGE(4, periods=periods)
model.setup_Q()
tree = Tree(newickstr, periods, root_age=5.2)
tree.set_default_model(model)
tree.set_tip_conditionals(data)

## d = 0.1
## e = 0.1
## model.setup_D(d)
## model.setup_E(e)
## model.setup_Q()
## ancdist_conditional_lh(tree.root)
## print -(math.log(sum(tree.root.dist_conditionals)))
## sys.exit()

#print tree.root.render_ascii(scaled=1, minwidth=80); sys.exit()

def f(params, model, tree):
    for p in params:
        if (p < 0) or (p > PMAX):
            return LARGE
    d, e = params
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    try:
        ancdist_conditional_lh(tree.root)
        return -(math.log(sum(tree.root.dist_conditionals)))
    except:
        #return 100
        return LARGE

for ai in model.arange:
    dist = [0]*model.nareas
    dist[ai] = 1
    dist = tuple(dist)
    for c in tree.root.children():
        c.segments[-1].startdist = dist
    v = optimize.fmin_powell(f, [0.1, 0.1], args=(model, tree))
    print dist
    print v
    print

sys.exit()
