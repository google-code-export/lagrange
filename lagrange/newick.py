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

import string, sys
from shlex import shlex
import phylo
from types import StringType
from cStringIO import StringIO

class Tokenizer(shlex):
    """Provides tokens for parsing Newick-format trees"""
    def __init__(self, infile):
        shlex.__init__(self, infile)
        self.commenters = ''
        self.wordchars = self.wordchars+'-.'
        self.quotes = "'"

    def parse_comment(self):
        while 1:
            token = self.get_token()
            if token == '':
                sys.stdout.write('EOF encountered mid-comment!\n')
                break
            elif token == ']':
                break
            elif token == '[':
                self.parse_comment()
            else:
                pass

def parse(input, ttable=None):
    """
    Parse a Newick-formatted tree description
    input is any file-like object that can be coerced into shlex,
    or a string (converted to StringIO)
    """
    if type(input) is StringType:
        input = StringIO(input)
    
    start_pos = input.tell()
    tokens = Tokenizer(input)

    node = None; root = None
    lp=0; rp=0; rooted=1

    prev_tok = None

    while 1:
        token = tokens.get_token()
        #print token,
        if token == ';' or token == '':
            assert lp == rp, \
                   'unbalanced parentheses in tree description'
            break

        # internal node
        elif token == '(':
            lp = lp+1
            newnode = phylo.InternalNode()
            if node:
                if node.istip:
                    if not node.back:
                        node.back = newnode
                        newnode.back = node
                    else:
                        node.back.add_child(newnode)
                else:
                    node.add_child(newnode)
            node = newnode

        elif token == ')':
            rp = rp+1
            node = traverse(node)
            
        elif token == ',':
            if lp == rp:
                if node.back == None:
                    pass
                else:
                    rooted = 0
                    if node.next:
                        node.insert_fnode(phylo.Fnode())
                    node.next.isroot = 1
                    root = node.next
            else:
                node = traverse(node)
            
        # branch length
        elif token == ':':
            token = tokens.get_token()

            if not (token == ''):
                try:
                    brlen = float(token)
                except ValueError:
                    raise 'NewickError', \
                          "invalid literal for branch length, '%s'" % token
            else:
                raise 'NewickError', \
                      'unexpected end-of-file (expecting branch length)'

            if node.istip: node.length = brlen
            else: node.next.length = brlen

        # comment
        elif token == '[':
            tokens.parse_comment()

        # leaf node or label
        else:
            if prev_tok != ')':
                if ttable is not None:
                    token = ttable[token]
                newnode = phylo.Fnode(label=token, istip=1)
                if node:
                    if node.istip:
                        if not node.back:
                            intnode = phylo.InternalNode()
                            intnode.back = node
                            node.back = intnode
                        node.back.add_child(newnode)
                        newnode = node.back
                    else: node.add_child(newnode)
                node = newnode
            else:
                node.next.label = token
                #print "label %s for %s" % (token, node.next)

        prev_tok = token
        #print node

    input.seek(start_pos)

    if rooted:
        root = phylo.Fnode(isroot=1)
        root.label = node.next.label; node.next.label = None
        root.length = node.next.length; node.next.length = None
        node.insert_fnode(root)
        

    for fn in root.fnodes()[1:]:
        if not fn.back:
            fn.prune()

    return root

def traverse(node):
    if node.istip: return node.back
    else: return node.next.back
        
def to_string(node, lengths = 1, length_fmt=":%s"):
    if not node.istip:
        node_str = "(%s)%s" % \
                   (",".join([to_string(child, lengths, length_fmt) \
                              for child in node.children()]),
                    node.label or "")
    else:
        node_str = "%s" % node.label

    length_str = ""
    if lengths:
        if node.length is not None:
            length_str = length_fmt % node.length
        else:
            length_str = ""

    s = "%s%s" % (node_str, length_str)
    return s
        
             

## def to_string(node, lengths = 1):
##     nstr = ''
##     if node.isroot and node.back:
##         nstr = to_string(node.back)+','

##     if not node.istip:
##         nstr = nstr+'('
##         children = node.children()
##         for child in children:
##             nstr = nstr+to_string(child, lengths)
##             if child != children[-1]: nstr = nstr+','
##         nstr = nstr+')'
##     else:
##         label = node.label
##         nstr = nstr+label
##     if (node.length is not None) and lengths:
##         nstr = nstr+':'+'%f' % node.length
##     return nstr

def parse_from_file(filename):
    if filename == '-':
        file = sys.stdin
    else:
        file = open(filename, 'r')
    content = string.strip(file.read())
    treedescs = string.split(content, ';')
    tree = parse(treedescs[0])
    file.close()
    return tree

if __name__ == "__main__":
    import ascii
    s = "(a,(b,c)int)root;"
    n = parse(s)
    print
    print ascii.render(n)
    print to_string(n)
    print n.next.back.label
