import sets, time
from gluon.storage import Storage
PREORDER = 0; POSTORDER = 1
BRANCHLENGTH = 0; INTERNODES = 1

class Node:
    def __init__(self):
        self.data = Storage()
        self.isroot = False
        self.istip = False
        self.label = None
        self.common_prefix = ""
        self.length = None
        self.parent = None
        self.children = []
        self.nchildren = 0
        self.hbranch_style = {}
        self.vbranch_style = {}
        self.label_style = {}
        self.ref_style = {}
        self.length_style = {}
        self.support_style = {}
        self.id = None
        self.tree_id = None

    def order_subtrees_by_size(self, n2s=None, recurse=False, reverse=False):
        if n2s is None:
            n2s = node2size(self)
        if not self.istip:
            v = [ (n2s[c], c.label, c) for c in self.children ]
            v.sort()
            if reverse:
                v.reverse()
            self.children = [ x[-1] for x in v ]
            if recurse:
                for c in self.children:
                    c.order_subtrees_by_size(n2s, recurse=True, reverse=reverse)

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

##     def leaves(self, v=None):
##         if v is None:
##             v = []
##         if not self.children:
##             v.append(self)
##         else:
##             for child in self.children:
##                 child.leaves(v)
##         return v

    def leaves(self, order=PREORDER, collapsed={}):
        return [ n for n in self.iternodes(order=order, collapsed=collapsed) \
                 if (n.istip or (n.id in collapsed)) ]

    def iternodes(self, order=PREORDER, v=None, collapsed={}):
        """
        returns a list of nodes descendant from self - including self
        """
        if collapsed is None:
            collapsed = {}
        if order == PREORDER:
            yield self
        if self.id not in collapsed:
            for child in self.children:
                for d in child.iternodes(order=order, v=v, collapsed=collapsed):
                    yield d
        if order == POSTORDER:
            yield self

    def descendants(self, order=PREORDER, v=None, collapsed={}):
        """
        returns a list of nodes descendant from self - not including self!
        """
        if v is None:
            v = []
        assert order in (PREORDER, POSTORDER)
        if self.id not in collapsed:
            for child in self.children:
                if order == PREORDER:
                    v.append(child)
                else:
                    v.insert(0, child)
                if child.children:
                    child.descendants(order=order, v=v, collapsed=collapsed)
        return v

    def find_descendant(self, label):
        if label == self.label:
            return self
        else:
            for child in self.children:
                n = child.find_descendant(label)
                if n:
                    return n
        return None

    def prune(self):
        p = self.parent
        if p:
            p.remove_child(self)
        return p

    def graft(self, node):
        parent = self.parent
        parent.remove_child(self)
        n = Node()
        n.add_child(self)
        n.add_child(node)
        parent.add_child(n)

    def leaf_distances(self, store=None, measure=BRANCHLENGTH,
                       collapsed={}):
        """
        for each internal node, calculate the distance to each leaf,
        measured in branch length or internodes
        """
        if store is None:
            store = {}
        leaf2len = {}
        if self.children and self.id not in collapsed:
            for child in self.children:
                if measure == BRANCHLENGTH:
                    assert child.length is not None
                    dist = child.length
                elif measure == INTERNODES:
                    dist = 1
                else:
                    raise "InvalidMeasure"
                child.leaf_distances(store=store, measure=measure,
                                     collapsed=collapsed)
                if child.istip or child.id in collapsed:
                    leaf2len[child.id] = dist
                else:
                    for k, v in store[child].items():
                        leaf2len[k] = v + dist
        else:
            #leaf2len[self] = {self.id: 0}
            leaf2len[self.id] = 0
        store[self] = leaf2len
        return store

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
                    
            for tip in oldtips:
                newnode = d[tip]
                while 1:
                    newnode = newnode.parent
                    oldnode = d[newnode]
                    if newnode.nchildren == 1:
                        child = newnode.children[0]
                        if newnode.length:
                            child.length += newnode.length
                        newnode.remove_child(child)
                        if newnode.parent:
                            parent = newnode.parent
                            parent.remove_child(newnode)
                            parent.add_child(child)
                        del d[oldnode]; del d[newnode]
                    if not newnode.parent:
                        break
            
        return d

def node2size(node, d=None):
    "recursively map node.id to number of descendant tips"
    if d is None:
        d = {}
    size = int(node.istip)
    if not node.istip:
        for child in node.children:
            node2size(child, d)
            size += d[child.id]
    d[node.id] = size
    return d

def node2visible(node, collapsed, d=None):
    "recursively map node.id to number of visible tips"
    if d is None:
        d = {}
    x = node.istip or (node.id in collapsed)
    size = int(x)
    if not x:
        for child in node.children:
            node2visible(child, collapsed, d)
            size += d[child.id]
    d[node.id] = size
    return d

def common_prefix(labels):
    words = [ (x or "").split("_") for x in labels ]
    common = []
    while words and (min([ len(w) for w in words ]) > 0):
        s = set([ x.pop(0) for x in words ])
        if len(s) == 1:
            common.append(list(s)[0])
        else:
            break
    return "_".join(common)

def autolabel(node, store=None):
    v = []
    for child in node.children:
        if child.istip:
            v.append(child.label)
        else:
            autolabel(child, store)
            v.append(child.common_prefix)
    node.common_prefix = common_prefix(v)
    if store:
        store[node.id] = node.common_prefix
    for child in node.children:
        if child.common_prefix == node.common_prefix:
            child.common_prefix = ""
            if store and child.id in store:
                del store[child.id]

def auto_collapse_info(node, collapsed, visible=True):
    if visible and node.id in collapsed:
        visible = False
    nnodes = 1 # total number of nodes, including node
    # number of visible leaves
    nvisible = int((visible and node.istip) or node.id in collapsed)
    ntips = int(node.istip)
    ntips_visible = int(node.istip and visible)
    node.has_labeled_descendant = False

    node.data.depth = 1

    if node.common_prefix and (not node.label):
        node.label = node.common_prefix
    for child in node.children:
        auto_collapse_info(child, collapsed, visible)
        nnodes += child.nnodes
        nvisible += child.nvisible
        ntips += child.ntips
        ntips_visible += child.ntips_visible
        if (child.label and (not child.istip)) \
           or (child.has_labeled_descendant):
            node.has_labeled_descendant = True
        if child.data.depth >= node.data.depth:
            node.data.depth = child.data.depth+1
    node.nnodes = nnodes
    node.nvisible = nvisible
    node.ntips = ntips
    node.ntips_visible = ntips_visible

def auto_collapse(root, collapsed, keep_visible, max_tips_visible=125):
    ntries = 0
    while True:
        if ntries > 10:
            return
        ntries += 1
        auto_collapse_info(root, collapsed)
        nvisible = root.nvisible
        #print "nvisible", nvisible
        if nvisible <= max_tips_visible:
            return
        
        v = []
        for node in root.iternodes():
            if node.label and (not node.istip) and node.parent and \
                   (node.id not in keep_visible):
                w = node.nvisible/float(node.data.depth)
                if node.has_labeled_descendant:
                    w *= 0.25
                v.append((w, node))
        v.sort(); v.reverse()
        for w, node in v:
            if node.id not in keep_visible and node.nvisible < (nvisible-1):
                collapsed[node.id] = 1
                nvisible -= node.nvisible
            if nvisible <= max_tips_visible:
                break

def reroot(oldroot, newroot):
    oldroot.isroot = False
    newroot.isroot = True
    v = []
    n = newroot
    while 1:
        v.append(n)
        if not n.parent: break
        n = n.parent
    #print [ x.label for x in v ]
    v.reverse()
    for i, cp in enumerate(v[:-1]):
        node = v[i+1]
        # node is current node; cp is current parent
        #print node.label, cp.label
        cp.remove_child(node)
        node.add_child(cp)
        cp.length = node.length
    return newroot


if __name__ == "__main__":
    import newick, ascii, os
    from numpy import array
    #tree = newick.parse("(a,(b,(c,(d,e))));")
    f = os.path.expanduser("~/Projects/pedic-sympatry/matrices/")
    tree = eval(file(f+"garli-ml.tree").read())
    treespp = tree["species"]
    root = newick.parse(tree["newick"])
    spp = ['alaschanica', 'cheilanthifolia', 'dichotoma', 'kansuensis',
           'oederi', 'plicata', 'przewalskii', 'remotiloba',
           'rhinanthoides', 'roylei', 'rupicola', 'scolopax']
    print root.subtree_mapping(spp, clean=1)
