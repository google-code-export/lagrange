# Mavric -- a module for manipulating and visualizing phylogenies

# Copyright (C) 2000 Rick Ree
# Email : rree@post.harvard.edu
# 	   
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version.
#   
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

import ascii

PREORDER = 0
POSTORDER = 1

class Fnode:
    """

    Each internal 'node' in a tree is represented as a linked list of
    Fnodes (the prefix F is in acknowledgement of Joe Felsenstein,
    author of the PHYLIP package, from which I got the idea.)  The
    Fnodes form a circular linked list, i.e., each Fnode points to the
    next Fnode via its 'next' attribute.  Each Fnode is connected to a
    branch (descendant or parent node) via its 'back' attribute.
    E.g., consider a node X that has parent node P and descendant
    nodes A and B.  Because X connects 3 nodes, it will consist of 3
    Fnodes: X1, X2, and X3.  Node X is constructed as follows: X1.next
    == X2, X2.next == X3, and X3.next == X1.  X1.back == (some Fnode
    of P), X2.back == (some Fnode of A), and X3.back == (some Fnode of
    B).  In this case, because P is the ancestor to X, X1 is the
    'root' Fnode of node X.

    Note that all Fnodes in a node can (and do) point to a common
    'data' member, which in this case is a dictionary.

    Leaf nodes only consist of a single Fnode, which does not have the
    'next' attribute; all it has is the 'back' attribute, which points
    back to the ancestral Fnode.  E.g., if node A is a leaf (with
    single Fnode A1), then X2.back == A1; A1.back == X2; and, thus,
    X2.back.back == X2.  The same goes for the root node of a tree.

    The advantage of this kind of tree structure is that rerooting
    nodes (and trees) becomes trivial: in the example above, you could
    call node B the root, and P and A become descendants of X (and X3
    becomes the 'root' Fnode of X).  Also, it becomes possible to
    traverse the entire tree non-recursively, just by following
    'fnode.next.back', or just 'fnode.back' in the case of leaf nodes.

    """
    
    render_ascii = ascii.tree2ascii

    def __init__(self, next=None, back=None, data=None,
                 isroot=0, istip=0, label=None, length=None):
        self.next = next
        self.back = back
        if not (data == None): self.data = data
        else: self.data = {}
        self.isroot = isroot
        self.istip = istip
        self.length = length
        self.label = label
        self.excluded_dists = set()

    def __setitem__(self, item, value):
        self.data[item] = value

    def __delitem__(self, item):
        del self.data[item]

    def __getitem__(self, item):
        try:
            return self.data[item]
        except KeyError:
            raise KeyError, 'item: %s not in %s' % \
                  (item, self.data.keys())

    def get(self, item):
        return self.data.get(item)

    def unlink(self):
        """Remove references to next, back, and data, to let refcounts
        go to zero"""
        del self.next; del self.back; del self.data
        #self.next = None; self.back = None; del self.data

    def fnodes(self):
        """Returns a list of the Fnodes linked to self, including self"""
        nodes = []
        n = self
        while 1:
            nodes.append(n)
            if n.next == self: break
            else: n = n.next
        return nodes

    def iterchildren(self):
        """Returns the immediate descendants of this node"""
        if not self.istip:
            n = self.next
            while 1:
                yield n.back
                if n.next == self: break
                else: n = n.next
                
    def children(self):
        """Returns the immediate descendants of this node"""
        return list(self.iterchildren())
##         if self.istip: return None
##         children = []
##         n = self.next
##         while 1:
##             children.append(n.back)
##             if n.next == self: break
##             else: n = n.next
##         return children

    def iterdescendants(self, order=PREORDER):
        """Returns a list of all descendants of this node."""
        if order == PREORDER:
            yield self
        for child in self.children():
            for d in child.descendants(order):
                yield d
        if order == POSTORDER:
            yield self

    def descendants(self, order=PREORDER):
        """Returns a list of all descendants of this node."""
        return list(self.iterdescendants(order))
##         d = []
##         if not self.istip:
##             d.append(self)
##             for child in self.children():
##                 if child.istip: d.append(child)
##                 else: d.extend(child.descendants())
##         return d
                
    def iterleaves(self):
        """Returns a list of leaf nodes that are descendant from this
        node."""
        for child in self.children():
            if child.istip:
                yield child
            else:
                for leaf in child.leaves():
                    yield leaf

    def leaves(self):
        """Returns a list of leaf nodes that are descendant from this
        node."""
        return list(self.iterleaves())
##         lvs = []
##         if not self.istip:
##             for child in self.children():
##                 if child.istip: lvs.append(child)
##                 else: lvs.extend(child.leaves())
##         return lvs

    def insert_fnode(self, node):
        """Add node to self's linked list of fnodes"""
        assert self.next
        old_next = self.next
        self.next = node
        node.next = old_next
        node.data = self.data

    def add_child(self, node):
        """Adds node to self's children"""
        assert self.next

        if self.next.back:
            self.insert_fnode(Fnode(data=self.data))

        self.next.back = node
        node.back = self.next

    def prune(self):
        """Remove self from its linked list of fnodes"""
        last = self; next = self.next
        while last.next != self:
            last = last.next

        last.next = next
        self.next = None
        del self.data; self.data = {}

        return last, next

    def bisect(self):
        """Bisect the branch between self and self.back with a new
        internal node"""
        assert self.back
        back = self.back
        n = InternalNode()

        # join the new internal node to self and self.back
        n.back = back; back.back = n
        self.back = n.next; n.next.back = self

        length = self.length
        if length:
            half = length*0.5
            self.length = half; n.length = half

        return n

    def make_polytomy(self):
        if self.istip: return
        for c in self.children():
            d = c.descendants()
            for n in d:
                if not n.istip:
                    n.collapse()

    def make_pectinate(self):
        """
        Order descendant branches according to their size, so largest
        subtrees are first, etc.
        """
        children = self.children()
        first = children[0]
        def sort_func(n1, n2):
            if n1.istip and n2.istip:
                lab1 = n1.label; lab2 = n2.label
                if lab1 < lab2: return -1
                if lab1 == lab2: return 0
                if lab1 > lab2: return 1
            else:
                n1lvs = len(n1.leaves()); n2lvs = len(n2.leaves())
                if n1lvs > n2lvs:  return -1
                if n1lvs == n2lvs:
                    lab1 = [x.label for x in n1.leaves()]; lab1.sort()
                    lab2 = [x.label for x in n2.leaves()]; lab2.sort()
                    if lab1 < lab2: return -1
                    if lab1 == lab2: return 0
                    if lab1 > lab2: return 1
                if n1lvs < n2lvs:  return 1
        for child in children:
            if not child.istip:
                child.make_pectinate()
        children.sort(sort_func)
        nptr = self.next
        for child in children:
            nptr.back = child
            child.back = nptr
            nptr = nptr.next

    def rotate(self):
        """shifts the position of node's children by one"""
        if self.istip: return
        children = self.children()
        nptr = children[-1]
        del children[-1]
        children.insert(0,nptr)
        nptr = self.next
        for child in children:
            nptr.back = child
            child.back = nptr
            nptr = nptr.next

    def swivel_180_degrees(self):
        """reverses the order of the node's children"""
        if self.istip: return
        children = self.children()
        children.reverse()
        nptr = self.next
        for child in children:
            nptr.back = child
            child.back = nptr
            nptr = nptr.next

    def swivel_180_recursive(self):
        if self.istip: return
        self.swivel_180_degrees()
        for child in self.children():
            child.swivel_180_recursive()

    def collapse(self):
        """collapse node"""
        assert (not self.istip)
        back = self.back
        children = self.children()
        back.back = children.pop(0)
        back.back.back = back
        for child in children:
            back.add_child(child)
        for node in self.fnodes():
            node.unlink()

    def has_descendant(self, node):
        """check of node is descendant of self"""
        if self.istip: return 0
        flag = 0
        for child in self.children():
            if child.data == node.data:
                return 1
            else:
                if child.has_descendant(node):
                    flag = 1
        return flag

    def print_leaves(self):
        for leaf in self.leaves():
            print leaf.label,
        print

    def collapse_descendants(self, threshold):
        """collapses descendant branches if length < threshold"""
        if not self.istip:
            for child in self.children():
                child.collapse_descendants(threshold)
            if (not self.istip) and (self.length < threshold):
                self.collapse()

def InternalNode(isroot=0):
    """return an internal node"""
    node = Fnode(isroot=isroot)
    node.next = Fnode(data=node.data)
    node.next.next = node
    return node

def connect(node1, node2):
    node1.back = node2
    node2.back = node1
    if node1.length is not None and node2.length is None:
        node2.length = node1.length
    if node2.length is not None and node1.length is None:
        node1.length = node2.length

def polarize(node, parent=None):
    node.parent = parent
    if not node.istip:
        for child in node.children():
            polarize(child, node)

def depolarize(node):
    for n in node.descendants():
        n.parent = None

def nodes_to_tips(node, vect=None, n=0):
    """return a list of how many internodes are between node and its
    leaves """

    if vect == None: vect = []
    
    if node.istip:
        vect.append(n)
    else:
        for child in node.children():
            nodes_to_tips(child, vect, n+1)
    return vect

def length_to_tips(node, vect=None, length=0.0):
    """return a list of total lengths between node and its leaves"""
    if vect == None:
        vect = []
        node_length = 0.0
    else:
        node_length = node.length or 1.0

    if node.istip:
        vect.append(length+node_length)
    else:
        for child in node.children():
            length_to_tips(child, vect, length+node_length)
    return vect

def reroot(oldroot, newroot):
    ofn = oldroot.fnodes()
    rooted = 1
    if oldroot.back != None:
        rooted = 0

    #if len(ofn) == 3:
    if rooted:
        oldroot.prune()
        if len(ofn) == 3:
            connect(ofn[1].back, ofn[2].back)
        for n in ofn:
            n.unlink()
    #else:
    #    ofn[-1].next = ofn[1]

    if newroot.istip:
        newroot = newroot.back
    else:
        nfn = newroot.fnodes()
        if rooted:
            nr = Fnode(isroot=1)
            nfn[-1].insert_fnode(nr)
            newroot = nr
        else:
            newroot.isroot = 1

    oldroot.isroot = 0
    #oldroot.unlink()

    print newroot.length

    return newroot

def split(node, numchildren=2):
    """split a node into multiple nodes"""
    if node.istip:
        node.istip = 0
        node.next = Fnode(data=node.data)
        node.next.next = node
    for c in range(numchildren):
        child = Fnode(istip=1)
        node.add_child(child)

def balanced_tree(depth):
    """tree will have 2**depth taxa!!"""
    root = InternalNode(isroot=1)

    nodes = [root]
    for i in range(depth):
        for n in nodes:
            split(n, 2)
        nodes = root.leaves()
    
    return root
            
def pectinate_tree(numtaxa):
    assert numtaxa >= 2
    leaves = []
    for i in range(numtaxa):
        leaves.append(Fnode(istip=1))
    node = InternalNode()
    node.add_child(leaves.pop())
    node.add_child(leaves.pop())
    while leaves:
        leaf = leaves.pop()
        intnode = InternalNode()
        intnode.add_child(leaf)
        intnode.add_child(node)
        node = intnode
    node.isroot = 1
    return node

def __determine_ancstatus(n, descendant_list):
    """private recursive function called by most_recent_common_ancestor"""
    if n.istip:
        n['labels'] = [n.label,]
        n['is_anc'] = 0
    else:
        labels = []
        children = n.children()
        for c in children:
            mrca = __determine_ancstatus(c, descendant_list)
            if mrca != None:
                return mrca
            labels.extend(c['labels'])
        
        n['labels'] = labels
        #n['is_anc'] = (label1 in labels and label2 in labels)
        n['is_anc'] = 1
        for d in descendant_list:
            if not d in labels:
                n['is_anc'] = 0;
                break

        if n['is_anc'] and not \
               filter(lambda x: x['is_anc'], children):
            #print n['is_anc'], map(lambda x: x['is_anc'], children)
            return n
    return None

def clades(n, label2index=None, _d=None):
    """
    return a mapping of taxon label tuples to nodes in a tree
    n: root node of tree
    label2index: optional mapping of labels to index (of an array)
    _d: dummy variable, do not pass in
    """
    return_d_flag = (_d is None)
    if return_d_flag: _d = {}
    labels = []
    if n.istip:
        if label2index:
            labels.append(label2index[n.label])
        else:
            labels.append(n.label)
    else:
        children = n.children()
        for c in children:
            labels.extend(clades(c, label2index, _d))
        labels.sort()
        _d[tuple(labels)] = n
    if return_d_flag: return _d
    else: return labels

def most_recent_common_ancestor(node, descendant_list):
    """find most recent common ancestor of two taxa
       parameters:
         node: root of tree to search
         descendant_list: list of descendants"""
    def clean(n):
        try:
            del n['labels']
            del n['is_anc']
        except KeyError:
            pass

    mrca = __determine_ancstatus(node, descendant_list)
    map(clean, [node,]+node.descendants())
    return mrca
    
def mrca(node, descendant_list):
    return most_recent_common_ancestor(node, descendant_list)

def random_tree(labels):
    """
    Given a list of labels, create a list of leaf nodes, and then one
    by one pop them off, randomly grafting them on to the growing tree.

    Return the root node.
    """
    assert len(labels) > 2
    import RandomArray; RandomArray.seed()
    leaves = []
    for label in labels:
        leaves.append(Fnode(istip=1, label=label))

    leaf_indices = list(RandomArray.permutation(len(leaves)))

    joined = [leaves[leaf_indices.pop()]]
    remaining = leaf_indices
    while remaining:
        i = RandomArray.randint(0, len(joined)-1)
        c1 = joined[i]
        if c1.back:
            n = c1.bisect()
        else:
            n = InternalNode()
            n.add_child(c1)
        c = leaves[remaining.pop()]
        n.add_child(c)
        joined.append(c)
        joined.append(n)

    for node in joined:
        if not node.back:
            node.isroot = 1
            return node


if __name__ == "__main__":
    import newick
    tree = newick.parse("(a,(b,(c,(d,e))));")
    for i, n in enumerate(tree.descendants(order=PREORDER)):
        print i, n.label or n
    #print tree.draw_ascii
