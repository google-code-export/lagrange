
import string, sys
from shlex import shlex
import n_node
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

    curnode = None; root = None
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
            newnode = n_node.Node(istip = 0)
            if curnode:
                curnode.add_child(newnode)
            else:
                newnode.isroot = 1
                root = newnode
            curnode = newnode

        elif token == ')':
            rp = rp+1
            curnode = curnode.parent           
            
        elif token == ',':
            curnode = curnode.parent
            
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
            curnode.length = brlen
        # comment
        elif token == '[':
            tokens.parse_comment()
        # leaf node or label
        else:
            if prev_tok != ')':
                if ttable is not None:
                    token = ttable[token]
                newnode = n_node.Node(label=token, istip=1)
                if curnode:
                    curnode.add_child(newnode)
                curnode = newnode
            else:
                curnode.label = token

        prev_tok = token
        #print node

    input.seek(start_pos)
    """
    if rooted:
        root = Node(isroot=1)
        root.label = node.label #node.next.label; node.next.label = None
        root.length = node.length # node.next.length; node.next.length = None
        #node.insert_fnode(root)
        
    for fn in root.fnodes()[1:]:
        if not fn.parent:
            fn.prune()
    """
    return root
        
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
    import tree_printer_n
    s = "(a,(b,c)int)root;"
    n = parse(s)
    print tree_printer_n.render(n)
    print to_string(n)
    #print n.next.back.label
