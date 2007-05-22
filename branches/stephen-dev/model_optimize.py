#!/usr/bin/env python
import sys
from pprint import pprint
from marginal_range_reconstructor import *
import scipy
from scipy import optimize

LARGE = 10e5
PMAX = 10.0

def likelihood_de(params, model, rec):
    """
    demo function for optimizing dispersal and extinction rates, for
    use with scipy.optimize

    * returns negative log-likelihood of the tree, model, and data for
      dispersal and extinction rates given in params

    * assumes that the tree has only one model assigned to all branch
      segments
    """
    d, e = params
    for p in (d, e):
        if (p < 0) or (p > PMAX):
            return LARGE
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    try:
        lh = rec.eval_likelihood()
        return -(math.log(lh))
    except:
        return LARGE

def ancsplit_likelihood_de(node, ancsplit, model, d, e):
    """
    calculate likelihood of ancsplit at node with dispersal and
    extinction rates d, e

    """
    c1, c2 = node.children()
    c1.segments[-1].startdist = ancsplit.descdists[0]
    c2.segments[-1].startdist = ancsplit.descdists[1]
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    lh = node.tree.eval_likelihood()
    c1.segments[-1].startdist = None
    c2.segments[-1].startdist = None
    return lh

def ancdist_likelihood_de(node, dist, model, d, e):
    """
    calculate likelihood of dist at node with dispersal and
    extinction rates d, e (weighted average of split likelihoods)
    """
    c1, c2 = node.children()
    model.setup_D(d)
    model.setup_E(e)
    model.setup_Q()
    v = []
    for ancsplit in model.iter_ancsplits(dist):
        c1.segments[-1].startdist = ancsplit.descdists[0]
        c2.segments[-1].startdist = ancsplit.descdists[1]
        lh = node.tree.eval_likelihood()
        v.append(lh * ancsplit.weight)
    return sum(v)

def ancsplit_optimize_de(node, ancsplit, model):
    """
    optimize likelihood of ancsplit at node with dispersal and
    extinction rates

    """
    c1, c2 = node.children()
    c1.segments[-1].startdist = ancsplit.descdists[0]
    c2.segments[-1].startdist = ancsplit.descdists[1]
    v = optimize.fmin_powell(
        likelihood_de, [0.1]*2, args=(model, node.tree),
        full_output=True,disp=0
        )
    params, negloglikelihood = v[:2]
    return params, negloglikelihood


def calculate_local_for_all_nodes(node, model):
    """
    optimize dispersal and extinction rates on tree with only one model

    * for each internal node, calculates likelihood and optimal
      dispersal and extinction rates for all split scenarios
    """
    tree = node.tree
    if not node.istip:
        c1, c2 = node.children()
        calculate_local_for_all_nodes(c1, model, tree)
        calculate_local_for_all_nodes(c2, model, tree)
        print ",".join([node.label, c1.label, c2.label])
        results = []
        clear_startdist(tree.root)
        for dist in model.dists:
            for split, weight in iter_dist_splits_weighted(dist):
                d1, d2 = split
                c1.segments[-1].startdist = d1
                c2.segments[-1].startdist = d2
                v = optimize.fmin_powell(
                    likelihood_de, [0.1]*2, args=(model, tree),
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
            print opt, d1, d2, params

def calculate_global_for_all_nodes(node, model, d, e, skip=[], _rv={}):
    if (not node.istip) and (not node.label in skip):
        c1, c2 = node.children()
        calculate_global_for_all_nodes(c1, model, d, e)
        calculate_global_for_all_nodes(c2, model, d, e)
        print ", ".join([node.label, c1.label, c2.label])
        results = []
        #node.tree.clear_startdist()
        for disti, dist in model.enumerate_dists():
            S = set()
            for ancsplit in model.iter_ancsplits(dist):
                d1, d2 = ancsplit.descdists
                # calcuate likelihood of this split at this node
                lh = ancsplit_likelihood_de(node, ancsplit, model, d, e)
                
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(node.tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                try:
                    results.append((-(math.log(lh)), ds1, ds2, root))
                except:
                    pass

        results.sort()
        v = []
        for opt, ds1, ds2, root in results:
            if opt > (results[0][0] + 2):
                break
            v.append((opt, ds1, ds2))
            print "  -lnL %g, %s, %s" % (opt, ds1, ds2)
        _rv[node.label] = v

    return _rv

def comb(ds1, ds2):
    ret = ""
    for i in range(len(ds1)):
        if ds1[i] == '1' or ds2[i] == '1':
            ret += '1'
        else:
            ret += '0'
    return ret

def calculate_global_for_all_nodes(node, model, d, e, skip=[], _rv={}):
    if (not node.istip) and (not node.label in skip):
        c1, c2 = node.children()
        calculate_global_for_all_nodes(c1, model, d, e)
        calculate_global_for_all_nodes(c2, model, d, e)
        print ", ".join([node.label, c1.label, c2.label])
        results = []
        #node.tree.clear_startdist()
        for disti, dist in model.enumerate_dists():
            S = set()
            for ancsplit in model.iter_ancsplits(dist):
                d1, d2 = ancsplit.descdists
                # calcuate likelihood of this split at this node
                lh = ancsplit_likelihood_de(node, ancsplit, model, d, e)
                
                ds1 = model.diststrings[model.dist2i[d1]]
                ds2 = model.diststrings[model.dist2i[d2]]
                root = zip(scipy.log(node.tree.root.dist_conditionals),
                           model.diststrings)
                root.sort(); root.reverse()
                for i, r in enumerate(root):
                    dc, ds = r
                    if dc < (root[0][0] - 2):
                        break
                root = root[:i]
                try:
                    results.append((-(math.log(lh)), ds1, ds2, root))
                except:
                    pass

        results.sort()
        v = []
        for opt, ds1, ds2, root in results:
            if opt > (results[0][0] + 2):
                break
            v.append((opt, ds1, ds2))
            print "  -lnL %g, %s, %s" % (opt, ds1, ds2)
        _rv[node.label] = v

    return _rv
