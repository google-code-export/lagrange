#!/usr/bin/env python
"""
formerly model09.py, changed to model.py since it's now under version
control

conventions:

dist = shorthand for distribution = binary vector indicating
       presence/absence in a vector of areas

"""
import sys, math, random, sets
import scipy
import phylo, newick
import rates_scipy as rates
import nchoosem

rand = random.Random()

class RateModel:
    """
    Model of dispersal and local extinction between discrete
    geographic areas.
    """
    def __init__(self,
                 nareas,         # number of areas
                 anames=None,    # ordered list of labels for areas
                 periods=[1.0],  # list of durations, ordered from
                                 # present to past, corresponding to
                                 # "strata" of dispersal and
                                 # extinction parameters; sum(periods)
                                 # must exceed the root age of any
                                 # tree that uses the model
                 dists=None
                 ):
        """
        initialize variables
        """
        self.nareas = nareas
        self.arange = range(nareas)
        if not anames:
            self.anames = [ "a%d" % i for i in self.arange ]
        else:
            assert len(anames) == nareas, \
                   "Mismatch between number of areas and number of area names"
        self.periods = periods
        self.nperiods = len(periods)
        self.prange = range(self.nperiods)
        # set up self.dists, the list of dists corresponding to the
        # number of areas; does not include the empty dist (vector of
        # all zeros) (cf. RateModelGE)
        self.setup_dists(dists)
        # make self.dist look like
        # ['100', '010', '001', '110', '101', '011', '111']
        # for three areas
        self.diststrings = [ "".join(map(str, d)) for d in self.dists ]
        self.ndists = len(self.dists)
        # make self.dist a dictionary 
        # looking like 
        # {(1, 1, 0): 3, (0, 1, 1): 5, (1, 0, 0): 0, (0, 0, 1): 2, (1, 0, 1): 4, (0, 1, 0): 1, (1, 1, 1): 6}
        # for three areas
        self.dist2i = dict([
            (dist, i) for i, dist in enumerate(self.dists)
            ])
        #self.dist_priors = [1.0/float(self.ndists)] * self.ndists
        self.distrange = range(self.ndists)
        # instantiate a matrix for each period
        # of nareas x nareas
        self.Dmask = scipy.ones((self.nperiods, nareas, nareas))
        default_d = 0.01
        default_e = 0.01
        self.setup_D(default_d)
        self.setup_E(default_e)
        self.setup_Q()

    def setup_dists(self, dists=None):
        # results in somthing like 
        # [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        # for three areas
        if dists is None:
            self.dists = [ dist for dist in \
                           nchoosem.iterate_all_bv(self.nareas) ]
        else:
            self.dists = dists[:]

    def setup_D(self, d):
        nareas = self.nareas
        self.D = scipy.ones((self.nperiods, nareas, nareas)) * d
        #for period, duration in enumerate(self.periods):
        for period in self.prange:
            for i in self.arange:
                self.D[period,i,i] = 0.0
        self.D *= self.Dmask

    def set_D_cell(self, period, fromdist, todist, d):
        i = self.dist2i[fromdist]
        j = self.dist2i[todist]
        self.D[period, i,j] = d

    def setup_E(self, e):
        self.E = scipy.ones((self.nperiods, self.nareas)) * e

    def setup_Q(self):
        self.Q = scipy.zeros((self.nperiods, self.ndists, self.ndists))
        #D = self.D * self.Dmask
        for p in range(self.nperiods):
            for i, d1 in enumerate(self.dists):
                s1 = sum(d1)
                if s1 > 0:
                    for j, d2 in enumerate(self.dists):
                        s2 = sum(d2)
                        xor = scipy.logical_xor(d1, d2)
                        # only consider transitions between dists that are
                        # 1 step (dispersal or extinction) apart
                        if sum(xor) == 1:
                            dest = scipy.nonzero(xor)[0]
                            #prior = self.dist_priors[i]
                            if s1 < s2: # dispersal
                                rate = 0.0
                                for src in scipy.nonzero(d1)[0]:
                                    rate += self.D[p,src,dest] * \
                                            self.Dmask[p,src,dest]
                                # for each area in d1, add rate of
                                # dispersal to dest
                            else: # extinction
                                rate = self.E[p,dest]
                            #self.Q[i][j] = (prior * rate)
                            self.Q[p,i,j] = rate
            self.set_Qdiag(p)

    def set_Qdiag(self, period):
        for i in self.distrange:
            self.Q[period,i,i] = (sum(self.Q[period,i,:]) - \
                                  self.Q[period,i,i]) * -1.0

    def P(self, period, t):
        """
        return P, the matrix of dist-to-dist transition probabilities,
        from the model's rate matrix (Q) over a time duration (t)
        """
        p = rates.Q2P(self.Q[period], t)
        # filter out impossible dists
        for i, dist in self.enumerate_dists():
            for area, present in enumerate(dist):
                if present and (sum(self.Dmask[period,area,:])+\
                                sum(self.Dmask[period,:,area]) == 0):
                    p[i,:] *= 0.0
                    break
        return p

    def Q_repr(self, period):
        lines = []
        widths = [ max([ max((len("%g" % x), self.nareas)) \
                         for x in self.Q[period,:,col] ]) \
                   for col in self.distrange ]
        lines.append("  ".join(
            [ ("".join([ str(a) for a in d ])).rjust(widths[col]) \
              for col, d in enumerate(self.dists) ]))
        lines.append("".join(["-"] * len(lines[0])))
                     
        for row in self.distrange:
            lines.append(
                "  ".join([ ("%g" % x).rjust(widths[col]) \
                            for col, x in enumerate(self.Q[period,row,:]) ])
                )

        w = max(len(r"From\To"), self.nareas)
        for i, line in enumerate(lines):
            if i == 0:
                lines[i] = "%s %s" % (r"From\To".rjust(w), line)
            elif i == 1:
                lines[i] = "%s%s" % ("".join(["-"]*(w+1)), line)
            else:
                d = self.dists[i-2]
                lines[i] = "%s|%s" % ("".join([ str(a) for a in d ]).ljust(w),
                                      line)

        return "\n".join(lines)

    def P_repr(self, period, t):
        P = self.P(period, t)#rates.Q2P(self.Q[period], t)
        lines = []
        widths = [ max([ max((len("%g" % x), self.nareas)) \
                         for x in P[:,col] ]) \
                   for col in self.distrange ]
        lines.append("  ".join(
            [ ("".join([ str(a) for a in d ])).rjust(widths[col]) \
              for col, d in enumerate(self.dists) ]))
        lines.append("".join(["-"] * len(lines[0])))
                     
        for row in self.distrange:
            lines.append(
                "  ".join([ ("%g" % x).rjust(widths[col]) \
                            for col, x in enumerate(P[row,:]) ])
                )

        w = max(len(r"From\To"), self.nareas)
        for i, line in enumerate(lines):
            if i == 0:
                lines[i] = "%s %s" % (r"From\To".rjust(w), line)
            elif i == 1:
                lines[i] = "%s%s" % ("".join(["-"]*(w+1)), line)
            else:
                d = self.dists[i-2]
                lines[i] = "%s|%s" % ("".join([ str(a) for a in d ]).ljust(w),
                                      line)

        return "\n".join(lines)

    def enumerate_dists(self):
        "enumerate non-empty dists"
        return [ (i, d) for i, d in zip(self.distrange, self.dists) \
                 if sum(d) ]

    def iter_dist_splits(self, dist):
        assert dist in self.dists
        if sum(dist) == 1:
            yield (dist, dist)
        else:
            for i in scipy.nonzero(dist)[0]:
                x = scipy.zeros((len(dist),), dtype="i")
                x[i] = 1
                x = tuple(x)
                if x in self.dists:
                    yield (x, dist)
                    yield (dist, x)
                    y = tuple(scipy.array(scipy.logical_xor(dist, x),
                                          dtype="i"))
                    if y in self.dists:
                        yield (x, y)
                        if sum(y) > 1:
                            yield (y, x)
                    
    def iter_ancsplits(self, dist):
        splits = [ s for s in self.iter_dist_splits(dist) ]
        nsplits = len(splits)
        weight = (1.0/nsplits)
        for split in splits:
            as = Ancsplit(dist, split, weight=weight)
            yield as

class RateModelGE(RateModel):
    """
    Subclass of RateModel that allows global extinction
    """
        
    def setup_dists(self, dists=None):
        if dists is None:
            self.dists = [ dist for dist in \
                           nchoosem.iterate_all_bv2(self.nareas) ]
        else:
            self.dists = dists[:]
            emptydist = tuple([0]*self.nareas)
            if emptydist not in dists:
                self.dists.insert(0, emptydist)

## print RateModelGE(3).enumerate_dists()
## sys.exit()

class Tree:
    def __init__(self, newickstr, periods=None, root_age=None):
        self.root = newick.parse(newickstr)
        phylo.polarize(self.root)
        self.periods = periods

        # initialize nodes (label interiors, etc)
        # and collect leaves and postorder sequence
        self.postorder_nodes = []
        self.leaves = []
        for i, node in enumerate(self.root.descendants(phylo.POSTORDER)):
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
        ancdist_conditional_lh(self.root)
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

class BranchSegment:
    def __init__(self, duration, period, model=None, startdist=None):
        self.duration = duration
        self.period = period
        self.model = model
        self.startdist = startdist


class Ancsplit:
    """
    convenience class for encapsulating an ancestor range splitting
    into descendant ranges
    """
    def __init__(self, ancdist, descdists, weight=None, likelihood=None):
        self.ancdist = ancdist
        self.descdists = descdists
        self.weight = weight
        self.likelihood = likelihood

    def __repr__(self):
        d1, d2 = [ "".join(map(str, x)) for x in self.descdists ]
        lh = self.likelihood
        if lh: lh = "%g" % lh
        w = self.weight
        if w: w = "%g" % w
        return "Ancsplit(%s, %s, w=%s, lh=%s)" % (d1, d2, w, lh)

#
# utility functions
#

def iter_dist_splits(dist):
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

def iter_ancsplits(dist):
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

def iter_dist_splits_weighted(dist):
    s = sum(dist)
    if s == 1:
        yield (dist, dist), 1.0
    else:
        wt = 1.0/(s*4)
        for sp in iter_dist_splits(dist):
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

def conditionals(node):
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
        P = model.P(seg.period, seg.duration)
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

def ancdist_conditional_lh1(node):
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

        ancdist_conditional_lh1(c1)
        ancdist_conditional_lh1(c2)

        v1 = conditionals(c1)
        v2 = conditionals(c2)

        distconds = scipy.zeros((model.ndists,))
        
        ancsplits = []
        for distidx, dist in model.enumerate_dists():
            lh = 0.0

            split_db = []
            for split, wt in iter_dist_splits_weighted(dist):
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

def ancdist_conditional_lh2(node):
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

        ancdist_conditional_lh2(c1)
        ancdist_conditional_lh2(c2)

        v1 = conditionals(c1)
        v2 = conditionals(c2)

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

ancdist_conditional_lh = ancdist_conditional_lh2

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
## for x in nondiag_indices(scipy.zeros([6,6])):
##     print x
## sys.exit()

## m = RateModelGE(4)
## s = set()
## for i, d in m.enumerate_dists():
##     for as in iter_ancsplits(d):
##         if as.descdists not in s:
##             s.add(as.descdists)
##         else:
##             print "!", as.descdists
