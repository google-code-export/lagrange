import sys
import newick
from pprint import pprint
from model import *
import model_optimize
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
    "P_greenwelliae07":          (1,0,0,0),
    "P_greenwelliae907":         (1,0,0,0),
    "P_grandiflora_Kal2":        (1,0,0,0),
    "P_hobdyi_Kuia":             (1,0,0,0),
    "P_hexandra_K1":             (1,0,0,0),
    "P_hexandra_M":              (1,0,0,0),
    "P_hexandra_Oahu":           (0,1,0,0),
    }

# newick tree with all taxa from Nepokroeff paper
# (((((P_mariniana_Kokee2:0.033539,(P_mariniana_WN2:0.030646,P_mariniana_Manoa2:0.030646,(P_mariniana_Kipahulu2:0.02241,P_mariniana_Kaw2:0.02241,P_mariniana_Lanai2:0.02241,P_hawaiiensis_Makaopuhi:0.02241):0.008236):0.002893):0.005171,P_wawraeDL7428:0.03871):0.008237,(((P_kaduana_PuuKukuiAS:0.020803,((P_hawaiiensis_WaikamoiL1:0.010853,P_mauiensis_Eke:0.010853):0.007964,(P_fauriei2:0.013826,P_hathewayi_1:0.013826):0.004991):0.001986,P_mauiensis_PepeAS:0.020803):0.003762,P_kaduana_HawaiiLoa:0.024565):0.003398,(P_greenwelliae07:0.012715,P_greenwelliae907:0.012715):0.015248):0.018984):0.008255,(P_grandiflora_Kal2:0.027864,P_hobdyi_Kuia:0.027864):0.027338):0.003229,((P_hexandra_K1:0.026568,P_hexandra_M:0.026568):0.005204,P_hexandra_Oahu:0.031771):0.026659):0.021119;


newickstr = "((((((((P_hawaiiensis_WaikamoiL1:0.010853,P_mauiensis_Eke:0.010853):0.007964,(P_fauriei2:0.013826,P_hathewayi_1:0.013826):0.004991):0.001986,(P_kaduana_PuuKukuiAS:0.020803,P_mauiensis_PepeAS:0.020803):1e-05):0.003762,P_kaduana_HawaiiLoa:0.024565):0.003398,(P_greenwelliae07:0.01271500,P_greenwelliae907:0.01271500):0.01524800):0.018984,((((P_mariniana_MauiNui:0.02241,P_hawaiiensis_Makaopuhi:0.02241):0.008236,P_mariniana_Oahu:0.030646):0.002893,P_mariniana_Kokee2:0.033539):0.005171,P_wawraeDL7428:0.03871):0.008237):0.008255,(P_grandiflora_Kal2:0.027864,P_hobdyi_Kuia:0.027864):0.027338):0.003229,((P_hexandra_K1:0.026568,P_hexandra_M:0.026568):0.005204,P_hexandra_Oahu:0.031771):0.026659);"

#dists = [
#    (1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),
#    (1,1,0,0),(0,1,1,0),(0,0,1,1),
#    (1,1,1,0),(0,1,1,1),(1,1,1,1),
#    ]
dists = None

#periods = [10.0]
# geological connections
periods = [0.5,1.9,3.9,5.6]

model = RateModel(4, periods=periods, dists=dists)

# geological connections
## default all but the 1 to be zero until they form
#island formations
#5.6 auai/Niihau
model.Dmask[3,0,1] = 0.0
model.Dmask[3,0,2] = 0.0
model.Dmask[3,0,3] = 0.0
model.Dmask[3,1,0] = 0.0
model.Dmask[3,1,2] = 0.0
model.Dmask[3,1,3] = 0.0
model.Dmask[3,2,0] = 0.0
model.Dmask[3,2,1] = 0.0
model.Dmask[3,2,3] = 0.0
model.Dmask[3,3,0] = 0.0
model.Dmask[3,3,1] = 0.0
model.Dmask[3,3,2] = 0.0
#3.9 Oahu
model.Dmask[2,0,1] = 1.0
model.Dmask[2,0,2] = 0.0
model.Dmask[2,0,3] = 0.0
model.Dmask[2,1,0] = 1.0
model.Dmask[2,1,2] = 0.0
model.Dmask[2,1,3] = 0.0
model.Dmask[2,2,0] = 0.0
model.Dmask[2,2,1] = 0.0
model.Dmask[2,2,3] = 0.0
model.Dmask[2,3,0] = 0.0
model.Dmask[2,3,1] = 0.0
model.Dmask[2,3,2] = 0.0
#1.9 Molokai, Lanai, Maui, and Kahoolawe (together Maui Nui)
model.Dmask[1,0,1] = 1.0
model.Dmask[1,0,2] = 0.5
model.Dmask[1,0,3] = 0.0
model.Dmask[1,1,0] = 1.0
model.Dmask[1,1,2] = 1.0
model.Dmask[1,1,3] = 0.0
model.Dmask[1,2,0] = 0.5
model.Dmask[1,2,1] = 1.0
model.Dmask[1,2,3] = 0.0
model.Dmask[1,3,0] = 0.0
model.Dmask[1,3,1] = 0.0
model.Dmask[1,3,2] = 0.0
#0.5 Hawaii
model.Dmask[0,0,1] = 1.0
model.Dmask[0,0,2] = 0.5
model.Dmask[0,0,3] = 0.25
model.Dmask[0,1,0] = 1.0
model.Dmask[0,1,2] = 1.0
model.Dmask[0,1,3] = 0.5
model.Dmask[0,2,0] = 0.5
model.Dmask[0,2,1] = 1.0
model.Dmask[0,2,3] = 1.0
model.Dmask[0,3,0] = 0.25
model.Dmask[0,3,1] = 0.5
model.Dmask[0,3,2] = 1.0

# allow only dispersal to adjacent islands
## model.Dmask[0,0,:] = [1,1,0,0]
## model.Dmask[0,1,:] = [1,1,1,0]
## model.Dmask[0,2,:] = [0,1,1,1]
## model.Dmask[0,3,:] = [0,0,1,1]

# allow only dispersal to adjacent islands,
# old->new
#model.Dmask[0,0,:] = [1,1,0,0]
#model.Dmask[0,1,:] = [0,1,1,0]
#model.Dmask[0,2,:] = [0,0,1,1]
#model.Dmask[0,3,:] = [0,0,0,1]

#print model.Q_repr(0); sys.exit()
#print model.P_repr(0, 1.0); sys.exit()
#print "\n".join(map(str, model.iter_dist_splits((1,1,1,0)))); sys.exit()
tree = Tree(newickstr, periods, root_age=5.2)
tree.set_default_model(model)
tree.set_tip_conditionals(data)

## print newick.to_string(tree.root); sys.exit()
## print tree.root.render_ascii(scaled=1, minwidth=80); sys.exit()

## for i, dist in model.enumerate_dists():
##     if sum(dist) == 1:
##         v = model_optimize.ancsplit_optimize_de(
##             tree.root, Ancsplit(dist, (dist, dist), weight=1.0), model
##             )
##         print dist
##         print v[0], v[1]
##         print
## sys.exit()
## # result:
## (1, 0, 0, 0)
## [ 0.08549811  0.00748089] 33.8841200491
## (0, 1, 0, 0)
## [ 0.08740711  0.006893  ] 37.8725989252
## (0, 0, 1, 0)
## [ 0.11688822  0.06155934] 46.423931664
## (0, 0, 0, 1)
## [ 0.16278906  0.11076712] 55.0916741519

# ml dispersal and extinction for root dist = 1000
d, e = 0.29707871, 10
c1, c2 = tree.root.children()

#map
#x = 0.1
#for i in range(100):
#    x = x + 0.1
#    d, e = 0.29707871, x
#    model.setup_D(d)
#    model.setup_E(e)
#    model.setup_Q()
#    tree.set_default_model(model)
#    tree.set_tip_conditionals(data)
#    print str(x)+"\t"+str(math.log(tree.eval_likelihood()));
#end map
#
#optimize before
#
c1.segments[-1].startdist = (1,0,0,0)
c2.segments[-1].startdist = (1,0,0,0)
v = optimize.fmin_powell(
    model_optimize.likelihood_de, [d,e], args=(model, tree),
    full_output=True,disp=0
    )
print v[:2]
#
#end optimize
#
results = {}
for node in (c1, c2):
    results.update(model_optimize.calculate_global_for_all_nodes(
        node, model, d, e
        ))
# results
## 2, P_hawaiiensis_WaikamoiL1, P_mauiensis_Eke
## 33.8877967403 0010 0010
## 5, P_fauriei2, P_hathewayi_1
## 33.8894728129 0100 0100
## 6, 2, 5
## 33.9034886929 0010 0100
## 9, P_kaduana_PuuKukuiAS, P_mauiensis_PepeAS
## 33.8842840368 0010 0010
## 10, 6, 9
## 33.9324254898 0110 0010
## 12, 10, P_kaduana_HawaiiLoa
## 34.2987766191 0100 0100
## 35.1084022602 0110 0100
## 15, P_greenwelliae07, P_greenwelliae907
## 33.8874848097 1000 1000
## 16, 12, 15
## 34.1608360966 0100 1000
## 35.4025529457 0110 1000
## 19, P_mariniana_MauiNui, P_hawaiiensis_Makaopuhi
## 34.0364389271 0010 0001
## 21, 19, P_mariniana_Oahu
## 34.1966309612 0010 0100
## 35.5943671521 0011 0100
## 23, 21, P_mariniana_Kokee2
## 34.6414925177 0100 1000
## 34.7484798877 0110 1000
## 36.4672621948 0111 1000
## 25, 23, P_wawraeDL7428
## 34.4793322113 1000 1000
## 35.1312158518 1100 1000
## 35.9886883588 1110 1000
## 26, 16, 25
## 34.1779232143 1000 1000
## 29, P_grandiflora_Kal2, P_hobdyi_Kuia
## 33.885580068 1000 1000
## 30, 26, 29
## 33.9042573072 1000 1000
## 33, P_hexandra_K1, P_hexandra_M
## 33.890110262 1000 1000
## 35, 33, P_hexandra_Oahu
## 34.000757579 1000 0100

## for label, ds1, ds2 in (
##     ("25", '1000', '1000'),
##     #("25", '1100', '1000'),
##     #("25", '1110', '1000'),
##     ("26", '1000', '1000'),
##     ("21", '0010', '0100'),
##     #("21", '0011', '0100'),
##     ("23", '0100', '1000'),
##     #("23", '0110', '1000'),
##     #("23", '0111', '1000'),
##     ("29", '1000', '1000'), ("2", '0010', '0010'),
##     ("5", '0100', '0100'), ("6", '0010', '0100'), ("9", '0010', '0010'),
##     ("10", '0110', '0010'),
##     ("12", '0100', '0100'),
##     #("12", '0110', '0100'),
##     ("15", '1000', '1000'),
##     ("16", '0100', '1000'),
##     #("16", '0110', '1000'),
##     ("19", '0010', '0001'), ("30", '1000', '1000'),
##     ("35", '1000', '0100'), ("33", '1000', '1000'),
##     ):
##     node = tree.label2node[label]
##     c1, c2 = node.children()
##     c1.segments[-1].startdist = tuple(map(int, ds1))
##     c2.segments[-1].startdist = tuple(map(int, ds2))

for k, v in results.items():
    node = tree.label2node[k]
    if v:
        ds1, ds2 = v[0][1:]

##         if len(v) > 1:
##             ds1, ds2 = v[1][1:]
##             print k, v[1]
##         else:
##             print k, v[0]

        print k, v[0]

        c1, c2 = node.children()
        c1.segments[-1].startdist = tuple(map(int, ds1))
        c2.segments[-1].startdist = tuple(map(int, ds2))

v = optimize.fmin_powell(
    model_optimize.likelihood_de, [d,e], args=(model, tree),
    full_output=True,disp=0
    )
print v[:2]

## print model_optimize.likelihood_de([d,e], model, tree)
