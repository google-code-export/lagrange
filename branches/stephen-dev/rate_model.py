 #!/usr/bin/env python
"""
CHECKED AND GOOD - SAS

conventions:

dist = shorthand for distribution = binary vector indicating
       presence/absence in a vector of areas

"""
#import psyco
#psyco.full()
import sys, math, random, sets
import scipy
import n_node, tree_reader_n
import rate_matrix as rates
from ancsplit import *
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

    def get_P(self, period, t):
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
