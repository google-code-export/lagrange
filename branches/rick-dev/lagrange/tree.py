import phylo, newick, ratemodel
import scipy

class Tree:
    def __init__(self, newickstr, periods=None, root_age=None):
        """
        newickstr: newick tree with branch lengths
        periods: list of durations ordered from present to past
        root_age: age of root node for branch length scaling
        """
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
                #node.age = 0.0
                self.leaves.append(node)
            self.postorder_nodes.append(node)

        self.root_age = root_age
        if root_age:
            self.calibrate(root_age)

        self.label2node = dict([(n.label, n) for n in self.postorder_nodes ])

        self.assign_node_ages()
##         for node in self.postorder_nodes:
##             if node.parent and (node.parent.age is None):
##                 node.parent.age = node.age + node.length

        # initialize branch segments
        for node in self.postorder_nodes:
            if node.parent:
                periods = self.periods
                anc = node.parent.age
                des = node.age
                assert anc > des
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

    def length2root(self, node):
        n = node
        v = 0.0
        while n.parent:
            assert n.length is not None
            v += n.length
            n = n.parent
        return v

    def assign_node_ages(self):
        maxlen = max([ self.length2root(lf) for lf in self.leaves ])
        for node in self.postorder_nodes:
            if node.parent:
                node.age = maxlen - self.length2root(node)
            else:
                node.age = maxlen

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
        self.root_age = depth

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

def ancdist_conditional_lh(node):
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

        ancdist_conditional_lh(c1)
        ancdist_conditional_lh(c2)

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
