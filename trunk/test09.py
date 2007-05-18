#!/usr/bin/env python
import sys
from pprint import pprint
from model09 import *
import scipy
from scipy import optimize

LARGE = 10e5
PMAX = 10.0

def f(params, model, tree):
    """
    demo function for optimizing dispersal and extinction rates
    """
    d, e = params
    for p in (d, e):
        if (p < 0) or (p > PMAX):
            return LARGE
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    try:
        lh = tree.eval_likelihood()
        return -(math.log(lh))
    except:
        return LARGE

def calculate_local_for_all_nodes(node,model,tree):
    if not node.istip:
        c1, c2 = node.children()
        calculate_local_for_all_nodes(c1, model, tree)
        calculate_local_for_all_nodes(c2, model, tree)
        print c1.label+"\t"+ c2.label
        results = []
        clear_startdist(tree.root)
        for dist in model.dists:
            for split, weight in iter_dist_splits_weighted(dist):
                d1, d2 = split
                c1.segments[-1].startdist = d1
                c2.segments[-1].startdist = d2
                v = optimize.fmin_powell(
                    f, [0.1]*2, args=(model, tree),
                    full_output=True,disp=0
                    )
                params, opt = v[:2]
                #opt = scipy.log(exp(-opt) * weight) * -1
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                results.append((float(opt), tuple(params), ds1, ds2, root))
        results.sort()
        for opt, params, d1, d2, root in results:
            if opt > (results[0][0] + 2):
                break
            print opt, d2, d1, params

def calculate_global_for_all_nodes(node,model,tree,d,e):
    if not node.istip:
        c1, c2 = node.children()
        calculate_global_for_all_nodes(c1, model, tree,d,e)
        calculate_global_for_all_nodes(c2, model, tree,d,e)
        print c1.label+"\t"+ c2.label
        results = []
        clear_startdist(tree.root)
        for dist in model.dists:
            for split, weight in iter_dist_splits_weighted(dist):
                d1, d2 = split
                c1.segments[-1].startdist = d1
                c2.segments[-1].startdist = d2
                model.setup_D(d)
                model.setup_E(e)
                model.setup_Q()
                lh = tree.eval_likelihood()
                #opt = scipy.log(exp(-opt) * weight) * -1
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                results.append((-(math.log(lh)), ds1, ds2, root))
        results.sort()
        for opt, d1, d2, root in results:
            if opt > (results[0][0] + 2):
                break
            print opt, d2, d1

def clear_startdist(node):
    if not node.istip:
        c1, c2 = node.children()
        clear_startdist(c1)
        clear_startdist(c2)
        c1.segments[-1].startdist = []
        c2.segments[-1].startdist = []

if __name__ == "__main__":

    #periods = [0.1,1.5]
    periods = [2.0,4.0]

    model = RateModel(3, periods=periods)
    model.Dmask[0,1,0] = 0.5
    model.Dmask[0,0,1] = 0.5
    #model.setup_D(0.01)
    #model.setup_E(0.06)
    d = 0.1
    e = 1.1
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    newickstr = "((((s1:1.5, s2:1.5)X:0.5, s3:2.0)Y:1.0, s4:3.0)Z:1.0, (s5:2.5, s6:2.5)M:1.5)N;"
##                             -----------------+ s1  
##                        ----X+                     
##             ----------Y+    -----------------+ s2  
##             :          :                          
##  ----------Z+          ----------------------+ s3 
##  :          :                                     
## N+          ---------------------------------+ s4 
##  :                                                
##  :               ----------------------------+ s5 
##  ---------------M+                                
##                  ----------------------------+ s6 


    tree = Tree(newickstr, periods, root_age=4.0)
    #print newick.to_string(tree.root); sys.exit()
    #print tree.root.render_ascii(scaled=1); sys.exit()

    data = {
        "s1": tuple(map(int, "010")),
        "s2": tuple(map(int, "100")),
        "s3": tuple(map(int, "010")),
        "s4": tuple(map(int, "001")),
        "s5": tuple(map(int, "001")),
        "s6": tuple(map(int, "001")),
        }

##     data = {
##         "s1": tuple(map(int, "010")),
##         "s2": tuple(map(int, "011")),
##         "s3": tuple(map(int, "001")),
##         "s4": tuple(map(int, "011")),
##         "s5": tuple(map(int, "001")),
##         "s6": tuple(map(int, "001")),
##         }
    
    tree.set_default_model(model)
    tree.set_tip_conditionals(data)

    # optimize dispersal and extinction
    v = optimize.fmin_powell(f, [0.1, 0.1], args=(model, tree))
    print v
    lh = tree.eval_likelihood()
    print -(math.log(lh))
    tree.print_dist_conds()
    #calculate_local_for_all_nodes(tree.root, model, tree)
    calculate_global_for_all_nodes(tree.root, model, tree,v[0],v[1])
    sys.exit()

    # iterate over ancestral ranges at node X, optimizing dispersal
    # and extinction rates, and recording most likely states at root
    results = []
    for dist in model.dists:
        for split, weight in iter_dist_splits_weighted(dist):
            d1, d2 = split
            clear_startdist(tree.root)
            tree.label2node["Y"].segments[-1].startdist = d1
            tree.label2node["s4"].segments[-1].startdist = d2
            #v = optimize.fmin_powell(
            #    f, [0.1]*2, args=(model, tree),
            #    full_output=True
            #    )
            #params, opt = v[:2]
            #opt = scipy.log(exp(-opt) * weight) * -1
            lh = tree.eval_likelihood()
            ds1 = model.diststrings[model.dist2i[d1]]
            ds2 = model.diststrings[model.dist2i[d2]]
            root = zip(scipy.log(tree.root.dist_conditionals),
                       model.diststrings)
            root.sort(); root.reverse()
            for i, r in enumerate(root):
                dc, ds = r
                if dc < (root[0][0] - 2):
                    break
            root = root[:i]
            results.append((-math.log(lh), ds1, ds2, root))
            print str(-math.log(lh)) +"\t" +ds1+"\t"+ str(ds2)
            #results.append((float(opt), tuple(params), ds1, ds2, root))
            
    #results.sort()
    #for opt, d1, d2, root in results:
    #    if opt > (results[0][0] + 2):
    #        break
    #    print opt, d2, d1
        #print "\t", "\n\t".join(map(str, root))



