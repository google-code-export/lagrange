import sys, scipy, random
#sys.path.append("/home/rree/src/som")
RND = random.Random()
expovariate = RND.expovariate
take = scipy.take; array = scipy.array
PREORDER = 0; POSTORDER = 1
BRANCHLENGTH = 0; INTERNODES = 1
PSETS = ((0,3,5), (1,4,6), (0,1,2,3,4))

def newickstr(node, length_fmt=":%s"):
    if not node.istip:
        node_str = "(%s)%s" % \
                   (",".join([ newickstr(child, length_fmt) \
                               for child in node.children ]),
                    node.label or ""
                    )
    else:
        node_str = "%s" % node.label

    if node.length is not None:
        length_str = length_fmt % node.length
    else:
        length_str = ""

    s = "%s%s" % (node_str, length_str)
    return s


class Node:
    def __init__(self):
        self.state = None
        self.i = None
        self.pi = [] # parameter index list
        self.pv = []
        self.psum = None
        self.isroot = False
        self.istip = False
        self.label = None
        self.length = None
        self.parent = None
        self.children = []
        self.nchildren = 0

    def add_child(self, child):
        assert child not in self.children
        self.children.append(child)
        child.parent = self
        self.nchildren += 1

    def remove_child(self, child):
        assert child in self.children
        self.children.remove(child)
        child.parent = None
        self.nchildren -= 1

    def leaves(self):
        return [ n for n in self.iternodes() if n.istip ]

    def iternodes(self, order=PREORDER, v=None):
        """
        returns a list of nodes descendant from self - including self
        """
        if order == PREORDER:
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order == POSTORDER:
            yield self

    def descendants(self, order=PREORDER, v=None):
        """
        returns a list of nodes descendant from self - not including self!
        """
        if v is None:
            v = []
        assert order in (PREORDER, POSTORDER)
        for child in self.children:
            if order == PREORDER:
                v.append(child)
            else:
                v.insert(0, child)
            if child.children:
                child.descendants(order, v)
        return v

    def speciate(self, params, event, nodenum):
        c1 = Node()
        c1.i = nodenum
        c1.label = "N%s" % nodenum; nodenum += 1
        self.add_child(c1)
        c2 = Node()
        c2.i = nodenum
        c2.label = "N%s" % nodenum; nodenum += 1
        self.add_child(c2)

        if event != 2: #self.state in (0, 1):
            c1.state = self.state; c2.state = self.state
        else:
            lam_a = params[0]; lam_b = params[1]; lam_ab = params[2]
            psum = lam_a + lam_b + lam_ab
            h = [lam_a/psum, lam_b/psum, lam_ab/psum]
            e = choose_weighted(h)
            if e == 0:
                c1.state = 0; c2.state = 2
            elif e == 1:
                c1.state = 1; c2.state = 2
            else:
                c1.state = 0; c2.state = 1

        for c in (c1, c2):
            c.istip = True
            c.length = 0.0
            c.pi = PSETS[c.state]
            c.pv = take(params, c.pi)
            c.psum = sum(c.pv)

        self.istip = False

        return nodenum

    def rootpath(self):
        n = self
        while 1:
            yield n
            if n.parent:
                n = n.parent
            else:
                break
            
    def subtree_mapping(self, labels, clean=False):
        """
        find the set of nodes in 'labels', and create a new tree
        representing the subtree connecting them.  nodes are assumed to be
        non-nested.

        return value is a mapping of old nodes to new nodes and vice versa.
        """
        d = {}
        oldtips = [ x for x in self.leaves() if x.label in labels ]
        for tip in oldtips:
            path = list(tip.rootpath())
            for node in path:
                if node not in d:
                    newnode = Node()
                    newnode.istip = node.istip
                    newnode.length = node.length
                    newnode.label = node.label
                    d[node] = newnode
                    d[newnode] = node
                else:
                    newnode = d[node]

                for child in node.children:
                    if child in d:
                        newchild = d[child]
                        if newchild not in newnode.children:
                            newnode.add_child(newchild)
        d["oldroot"] = self
        d["newroot"] = d[self]
        if clean:
            n = d["newroot"]
            while 1:
                if n.nchildren == 1:
                    oldnode = d[n]
                    del d[oldnode]; del d[n]
                    child = n.children[0]
                    child.parent = None
                    child.isroot = True
                    d["newroot"] = child
                    d["oldroot"] = d[child]
                    n = child
                else:
                    break
                    
            for newnode in d["newroot"].leaves():
                path = list(newnode.rootpath())
                for n in path:
                    if len(n.children) == 1:
                        child = n.children[0]
                        if n.length:
                            child.length += n.length
                        n.remove_child(child)
                        if n.parent:
                            parent = n.parent
                            parent.remove_child(n)
                            parent.add_child(child)
                        del d[d[n]]; del d[n]
            
        return d

def choose_weighted(h):
    """
    h is a list of probabilities where sum(h) == 1.0; return a random
    index into h weighted by its probability
    """
    #assert sum(h) - 1.0 < 0.00001
    rv = RND.random()
    j = 0.0
    for i, c in enumerate(h):
        j += c
        if rv < j: return i
    raise "error in choose_weighted: rv = %s, j = %s" % (rv, j)

def step(nodes, params):
    "choose a node"
    nodev = array([ n.psum for n in nodes ])
    nodev_sum = sum(nodev)
    dt = expovariate(nodev_sum)
    i = choose_weighted(nodev/nodev_sum)
    return (dt, i)

def rndevent(node):
    i = choose_weighted(node.pv/node.psum)
    return node.pi[i]

def sim(lam_a, lam_b, lam_ab, mu_a, mu_b, q_ab, q_ba, ntips):
    pnames = "lam_a lam_b lam_ab mu_a mu_b q_ab q_ba".split()
    params = [lam_a, lam_b, lam_ab, mu_a, mu_b, q_ab, q_ba]
    nodes = []
    tries = 0
    while 1:
        if not nodes:
            root = Node()
            root.label = "root"
            root.state = 0
            nodenum = root.speciate(params, 0, 1)
            nodes = list(root.children)
            t = 0
            tries += 1
            if tries > 100: return None
        else:
            dt, i = step(nodes, params)
            t += dt
            for node in nodes:
                node.length += dt
            node = nodes[i]; event = rndevent(node)
            if event in (0,1,2): # speciation
                #print "speciation"
                if len(nodes) >= ntips:
                    break
                nodenum = node.speciate(params, event, nodenum)
                nodes.remove(node); nodes.extend(node.children)
            elif event == 3: # extinction
                #print "extinction a"
                if node.state == 2:
                    node.state = 1
                    node.pi = PSETS[1]
                    node.pv = take(params, node.pi)
                    node.psum = sum(node.pv)
                else:
                    node.label = "E" + node.label[1:]
                    nodes.remove(node) # global extinction
            elif event == 4: # extinction
                #print "extinction b"
                if node.state == 2:
                    node.state = 0
                    node.pi = PSETS[1]
                    node.pv = take(params, node.pi)
                    node.psum = sum(node.pv)
                else:
                    node.label = "E" + node.label[1:]
                    nodes.remove(node) # global extinction
            else: # dispersal
                #print "dispersal"
                node.state = 2
                node.pi = PSETS[2]
                node.pv = take(params, node.pi)
                node.psum = sum(node.pv)

    return root
    

if __name__ == "__main__":
    from scipy import optimize
    from model import GSE2_Parameters, evaluate_node_downpass
    sys.path.append("..")
    import newick, ascii

    while 1:
        simtree = sim(0.1, 0.1, 0.3, 0.003, 0.003, 0.015, 0.015, 100)
        data = {}
        for lf in simtree.leaves():
            if lf.label.startswith("N"):
                x = [0,0,0]; x[lf.state] = 1
                data[lf.label] = x

        labels = [ lf.label for lf in simtree.leaves() \
                   if lf.label.startswith("N") ]
        stmap = simtree.subtree_mapping(labels, clean=True)
        extant = stmap["newroot"]
        if simtree == stmap[extant]:
            break
    
    for node in extant.descendants():
        if len(node.children) == 1:
            print node.label, stmap[node]

    tipstates = [0,0,0]
    for leaf in [ lf for lf in simtree.leaves() if lf.label.startswith("N") ]:
        tipstates[leaf.state] += 1
    print "tipstates", tipstates

    nstr = newickstr(extant)+";"
    tree = newick.parse(nstr)
    #print ascii.render(tree, data=data, scaled=1)
    
    def f(params):
        lam_w, lam_b, mu, q = params
        for x in params:
            if (x <= 0) or x >= 10:
                return 10e10
##         if (lam_w + lam_b) > 1:
##             return 10e10
        p = GSE2_Parameters((lam_w, lam_w, lam_b, mu, mu, q, q))
        v = evaluate_node_downpass(tree, data, p)
        rv = -sum(scipy.log(v))
        #rv = -(scipy.log(v[0]))
        if (rv < 0) or rv == scipy.nan:
            return 10e10
        return rv

    p = GSE2_Parameters((0.1, 0.1, 0.3, 0.003, 0.003, 0.02, 0.02))
    v = evaluate_node_downpass(tree, data, p)
    print v
    p = GSE2_Parameters((1.3, 1.3, 2.2, 2.3, 1.3, 1.2, 1.2))
    v = evaluate_node_downpass(tree, data, p)
    print v

##     print optimize.fmin_powell(f, [0.1, 0.3, 0.03, 0.01],
##                                full_output=True, disp=1)
