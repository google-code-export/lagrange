"""
should combine biogeotree and model_utils
delete once checked
"""

import sys, math, random, sets
import scipy
import n_node, tree_reader_n
import rate_matrix as rates
from rate_model import *
import nchoosem
from branch_segment import *
from ancsplit import *
import tree_printer_n

rand = random.Random()

class MarginalRangeReconstructor:
    def __init__(self, ratemodel, newickstr, periods=None, root_age=None):
        self.ratemodel = ratemodel
        self.root = tree_reader_n.parse(newickstr)
        n_node.polarize(self.root)
        self.periods = periods

        # initialize nodes (label interiors, etc)
        # and collect leaves and postorder sequence
        self.postorder_nodes = []
        self.leaves = []
        for i, node in enumerate(self.root.descendants(n_node.POSTORDER)):
            node.tree = self
            node.number = i
            node.segments = []
            if (not node.istip) and (not node.label):
                node.label = str(node.number)
            node.age = None
            if node.istip:
                node.age = 0.0
                self.leaves.append(node)
            self.postorder_nodes.append(node)

        self.root_age = root_age
        if root_age:
            self.calibrate(root_age)
            self.root_age = root_age

        self.label2node = dict([(n.label, n) for n in self.postorder_nodes ])

        for node in self.postorder_nodes:
            if node.parent and (node.parent.age is None):
                node.parent.age = node.age + node.length

        # initialize branch segments
        for node in self.postorder_nodes:
            if node.parent:
                periods = self.periods
                anc = node.parent.age
                des = node.age
                t = des
                for i, p in enumerate(periods):
                    s = sum(periods[:i+1])
                    if t < s:
                        duration = min((s - t, anc - t))
                        if duration > 0:
                            seg = BranchSegment(duration, i)
                            node.segments.append(seg)
                        t += p
                    if t > anc:
                        break
                #print node.label, anc, des, node.length, [ s.duration for s in node.segments ]
        self.ancdist_conditional_lh = self.ancdist_conditional_lh1
        
    def iter_dist_splits(self, dist):
        if sum(dist) == 1:
            yield (dist, dist)
        else:
            for i in scipy.nonzero(dist)[0]:
                x = scipy.zeros((len(dist),), dtype="i")
                x[i] = 1
                x = tuple(x)
                yield (x, dist)
                yield (dist, x)
                y = tuple(scipy.array(scipy.logical_xor(dist, x), dtype="i"))
                yield (x, y)
                if sum(y) > 1:
                    yield (y, x)
    
    ## d = (1,1,0,0)
    ## for x in iter_dist_splits(d):
    ##     print x
    ## sys.exit()
    
    def dist_splits(dist):
        return sets.Set([ s for s in iter_dist_splits(dist) ])
    
    def iter_ancsplits(self, dist):
        splits = [ s for s in iter_dist_splits(dist) ]
        nsplits = len(splits)
        weight = (1.0/nsplits)
        for split in splits:
            as = Ancsplit(dist, split, weight=weight)
            yield as
    
    ## d = (1,1,0,0)
    ## for x in iter_ancsplits(d):
    ##     print x
    ## sys.exit()
    
    def iter_dist_splits_weighted(self,dist):
        s = sum(dist)
        if s == 1:
            yield (dist, dist), 1.0
        else:
            wt = 1.0/(s*4)
            for sp in self.iter_dist_splits(dist):
                yield sp, wt
    
    ## d = (1,1,0)
    ## for x in iter_dist_splits_weighted(d):
    ##     print x
    ## print
    ## sys.exit()
    
    def test_conditionals(distconds, seglens, model):
        """
        small test function to make sure evaluating likelihoods along
        segmented branches works properly
        """
        distrange = model.distrange
        for p, seglen in enumerate(seglens):
            print "distconds", distconds
            P = model.P(p, seglen)
            v = scipy.zeros((model.ndists,))
            for i in distrange:
                # P[i] is the vector of probabilities of going from dist i
                # to all other dists
                v[i] = sum(distconds * P[i])
            distconds = v
        return distconds
    
    def calculate_conditionals(self,node):
        """
        calculate the conditional likelihoods of ranges at the start of a
        branch given the likelihoods of states at its end
        """
    
        # conditional likelihoods of dists at the end of the branch: for
        # tips, this is a zeroed vector with 1.0 at the observed range
        distconds = node.segments[0].dist_conditionals
        #print node.label, node.segments, node.length
        #for model, time, startdist in node.segments:
        for seg in node.segments:
            seg.dist_conditionals = distconds
            model = seg.model
            P = model.get_P(seg.period, seg.duration)
            v = scipy.zeros((model.ndists,))
            if seg.startdist:
                distrange = (model.dist2i[seg.startdist],)
            else:
                distrange = model.distrange
    
            for i in distrange:
                # P[i] is the vector of probabilities of going from dist i
                # to all other dists
                v[i] = sum(distconds * P[i])
            distconds = v
        return distconds
    
    def test_conditionals(model, period, duration, distconds):
        P = model.P(period, duration)
        v = scipy.zeros((model.ndists,))
        distrange = model.distrange
        for i in distrange:
            # P[i] is the vector of probabilities of going from dist i
            # to all other dists
            v[i] = sum(distconds * P[i])
        return v
    
    def ancdist_conditional_lh1(self,node):
        """
        recursive calculation of fractional likelihoods for dists at
        internal nodes in a tree
        """
        if not node.istip:
            c1, c2 = node.children()
            if node.parent:
                model = node.segments[0].model
            else:
                model = node.model
            dist2i = model.dist2i
    
            self.ancdist_conditional_lh1(c1)
            self.ancdist_conditional_lh1(c2)
    
            v1 = self.calculate_conditionals(c1)
            v2 = self.calculate_conditionals(c2)
    
            distconds = scipy.zeros((model.ndists,))
            
            ancsplits = []
            for distidx, dist in model.enumerate_dists():
                lh = 0.0
    
                split_db = []
                for split, wt in self.iter_dist_splits_weighted(dist):
                    if split in split_db:
                        continue
                    else:
                        split_db.append(split)
                    d1, d2 = split
    
                    lh_part = (v1[dist2i[d1]] * v2[dist2i[d2]])
    
                    lh += (lh_part * wt)
    
                    ancsplits.append(Ancsplit(dist, split,
                                              weight=wt, likelihood=lh_part))
                distconds[distidx] = lh
            node.ancsplits = ancsplits
    
        else:
            distconds = node.segments[0].dist_conditionals
        
        if node.parent:
            node.segments[0].dist_conditionals = distconds
        else:
            node.dist_conditionals = distconds
    
    def ancdist_conditional_lh2(self, node):
        """
        recursive calculation of fractional likelihoods for dists at
        internal nodes in a tree
        """
        if not node.istip:
            c1, c2 = node.children()
            if node.parent:
                model = node.segments[0].model
            else:
                model = node.model
            dist2i = model.dist2i
    
            self.ancdist_conditional_lh2(c1)
            self.ancdist_conditional_lh2(c2)
    
            v1 = self.calculate_conditionals(c1)
            v2 = self.calculate_conditionals(c2)
    
            distconds = scipy.zeros((model.ndists,))
            ancsplits = []
            for distidx, dist in model.enumerate_dists():
                lh = 0.0
                for ancsplit in model.iter_ancsplits(dist):
                    d1, d2 = ancsplit.descdists
                    lh_part = (v1[dist2i[d1]] * v2[dist2i[d2]])
                    lh += (lh_part * ancsplit.weight)
                    ancsplit.likelihood = lh_part
                    ancsplits.append(ancsplit)
                distconds[distidx] = lh
            node.ancsplits = ancsplits
    
        else:
            distconds = node.segments[0].dist_conditionals
        
        if node.parent:
            node.segments[0].dist_conditionals = distconds
        else:
            node.dist_conditionals = distconds
    
    
    def nondiag_indices(m):
        """
        m is a square array - iterate over the (i,j) pairs indexing the
        non-diagonal cells
        """
        ind = 0
        N = len(m)
        R = range(N)
        for i in R:
            for j in R:
                if i != j:
                    yield (i,j)
                ind += 1
    
    def set_default_model(self, model):
        "assigns the same model to all segments of all branches"
        for node in self.postorder_nodes:
            for seg in node.segments:
                seg.model = model
            if not node.parent:
                node.model = model

    def set_tip_conditionals(self, data):
        for leaf in self.leaves:
            segments = leaf.segments
            model = segments[0].model
            cond = scipy.zeros((model.ndists,))
            dist = data[leaf.label]
            cond[model.dist2i[dist]] = 1.0
            segments[0].dist_conditionals = cond

    def calibrate(self, depth):
        "scale an ultrametric tree to given root-to-tip depth"
        len2tip = 0.0
        node = self.leaves[0]
        while 1:
            len2tip += (node.length or 0.0)
            node = node.parent
            if not node: break

        scale = depth/len2tip

        for node in self.postorder_nodes:
            if node.parent:
                node.length *= scale
            else:
                node.length = 0.0

    def eval_likelihood(self):
        """
        evaluate fractional likelihoods at root node
        """
        self.ancdist_conditional_lh(self.root)
        return scipy.sum(self.root.dist_conditionals)

    def print_dist_conds(self):
        "report the fractional likelihoods"
        print "Likelihoods at root"
        model = self.root.model
        for i, s in enumerate(model.diststrings):
            x = self.root.dist_conditionals[i]
            if x:
                try:
                    print s, math.log(x)
                except:
                    print s, "Undefined"
            else:
                print s, "NaN"

    def clear_startdist(self, node=None):
        "recursively remove startdists from all branches"
        if node is None:
            node = self.root
        if not node.istip:
            c1, c2 = node.children()
            self.clear_startdist(c1)
            self.clear_startdist(c2)
            c1.segments[-1].startdist = []
            c2.segments[-1].startdist = []
