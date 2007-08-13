from lagrange import *

#comments are denoted as such
"""
    multi-line
    comments are
    like this
"""
anames = ["1","2","3"]
data = {
    "a":    (1,0,0),
    "b":    (1,0,0),
    "c":    (1,0,0),
    "d":    (1,0,1),
    "e":    (0,1,0),
    "f":    (0,1,0),    
    }

newickstr = "((((a:0.25,b:0.25):0.25,c:0.5):0.25,(d:0.65,e:0.65):0.1):0.1,f:0.85);"

dists = None

#from recent to old
periods = [0.5,1.0]

model = RateModelGE(3,labels=anames, periods=periods)

#period, from, to = value (0-1)
model.Dmask[1,2,0] = 0.0
model.Dmask[1,0,2] = 0.5

tree = Tree(newickstr,periods)

tree.set_default_model(model)
tree.set_tip_conditionals(data)

d = 0.1
e = 0.1
c1, c2 = tree.root.children()
results = {}
for node in (c1, c2):
    results.update(optimize.calculate_global_for_all_nodes(
        node, model, d, e
        ))

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
        

