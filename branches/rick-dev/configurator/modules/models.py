from gluon.sqlhtml import *
from gluon.html import URL
from gluon.sql import *
from gluon.storage import Storage
import cPickle, pprint, treeview, string, newick
import lagrange, scipy

nbsp = XML("&nbsp;")

class DECModel:
    def __init__(self, session):
        self.name = "Untitled%s" % len(session.models or [])
        self.datamatrix = DataMatrix(self)
        self.adjacencymatrix = AdjacencyMatrix(self)
        self.dm = DispersalMatrixMP(self)
        self.treelist = TreeList(self)
        self.rangelist = RangeList(self)
        self.max_range_size = 0
        self.base_rates = "__estimate__"

    def render_base_rates(self, r, s):
        cb = URL(r=r,f="base_rates_radio_clicked")
        target = "base_rates"
        est = INPUT(
            _type="radio", _name="base_rates", _value="estimate",
            _onclick="this.form.submit();",
            value="estimate" if self.base_rates=="__estimate__" else ""
            )
        fix = INPUT(
            _type="radio", _name="base_rates", _value="fix",
            _onclick="$('#base_rates_inputs').show('normal');",
            value="fix" if self.base_rates!="__estimate__" else ""
            )
        disp = INPUT(
            _name="dispersal", _size=4,
            _value=self.base_rates.get("dispersal") \
            if self.base_rates!="__estimate__" else ""
            )
        ext = INPUT(
            _name="extinction", _size=4,
            _value=self.base_rates.get("extinction") \
            if self.base_rates!="__estimate__" else ""
            )
        
        inputs = TABLE(
            TR(TD("dispersal:"), TD(disp),
               TD(nbsp, nbsp, "extinction:"), TD(ext),
               TD(INPUT(_type="submit", _value="Submit"))),
            _id="base_rates_inputs",
            _style="border-spacing:0.25em;margin-left:1em;"
            )
        if self.base_rates == "__estimate__":
            inputs.update(
                _style="border-spacing:0.25em;margin-left:1em;display:none;"
                )
        return FORM(
            DIV(
            DIV("Baseline rates of dispersal and local extinction",
                _class="widgetheader", _style="padding:0.5em;"),
            DIV(
            DIV(est, " estimate", _style="margin-bottom:0.5em;"),
            DIV(fix, " fix", inputs),
            _class="widgetcontent"
            ),
            _class="widgetwrapper"
            ),
            _action="base_rates_form"
            )

    def rangeconstraints_toolbar(self, request, session):
        session.rci = int(request.vars.rci or session.rci or 0)
        v = []
        for i, x in enumerate(
            ("Adjacency matrix",
             "Maximum range size",
             "Ranges allowed in analysis")
            ):
            u = URL(r=request,f="index",vars=dict(rci=i))
            a = A(x, _href=u)
            v.append(TD(a))
        v[session.rci].update(_class="selected")
        return TABLE(TR(*v), _id="rangeconstraints_toolbar",
                     _class="toolbar")

    def rangeconstraints_render(self, request, session):
        rci = int(request.vars.rci or session.rci or 0)
        if rci == 0:
            return self.adjacencymatrix.render(request, session)
        elif rci == 1:
            return self.render_max_range_size(request, session)
        elif rci == 2:
            return self.rangelist.render(request, session)

    def errors(self):
        v = [ x for x in [
            self.datamatrix.errors(),
            self.dm.errors(),
            self.treelist.errors(),
            self.rangelist.errors()
            ] if x ]
        return v

    def create_script(self):
        s = self.serialize()
        return """\
#!/usr/bin/env python
import os
import lagrange
data = \"\"\"\\
### begin data
%s
### end data
\"\"\"

i = 0
while 1:
    if not i:
        outfname = "%s.results.txt"
    else:
        outfname = "%s.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = file(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)
""" % (s, self.name, self.name)

    def calc_ranges(self):
        g = self.adjacencymatrix.graph()
        ranges = [ tuple(x) for x in self.datamatrix.data.values() ]
        nareas = len(self.datamatrix.labels)
        maxareas = self.max_range_size
        dists = set([ x for x in
                      lagrange.nchoosem.dists_by_maxsize_idx(nareas, maxareas)\
                      if x and g.are_connected(x) ] + ranges)
        dists.add(tuple())
        self.rangelist.ranges = dists
        self.rangelist.excluded = set()

    def serialize(self):
        d = {
            "model_name": self.name,
            "area_labels": self.datamatrix.labels,
            "taxon_range_data": self.datamatrix.serialize_data(),
            "taxa": self.datamatrix.taxa,
            "area_adjacency": self.adjacencymatrix.data,
            "area_dispersal": self.dm.data,
            "dispersal_durations": self.dm.ages2durations(),
            "max_range_size": self.max_range_size,
            "dm_symmetric_entry": self.dm.symmetric,
            "newick_trees": [ t.serialize() for t in self.treelist.trees ],
            "ranges": self.rangelist.ranges_sorted(),
            "excluded_ranges": self.rangelist.excluded_sorted(),
            "base_rates": self.base_rates,
            "lagrange_version": lagrange.VERSION
            }
        return pprint.pformat(d)

    def restore(self, d):
        self.name = d["model_name"]
        labels = d["area_labels"]
        self.datamatrix.labels = labels
        self.datamatrix.data = d["taxon_range_data"]
        self.datamatrix.taxa = d["taxa"]
        self.adjacencymatrix.data = d["area_adjacency"]
        self.adjacencymatrix.labels = labels
        self.dm.data = d["area_dispersal"]
        self.dm.labels = labels
        self.dm.durations2ages(d["dispersal_durations"])
        self.max_range_size = d["max_range_size"]
        self.dm.symmetric = d["dm_symmetric_entry"]
        self.treelist.trees = [ treeview.restore(t, self)
                                for t in d["newick_trees"] ]
        self.rangelist.ranges = set(d["ranges"])
        self.rangelist.excluded = set(d["excluded_ranges"])
        self.base_rates = d.get("base_rates") or "__estimate__"

    def render_name(self, request, session):
        cb = URL(r=request, f="update_modelname")
        ajax = "textinput_enter('%s', 'model:name:input', 'model:name');" % cb
        inp = INPUT(_value=self.name, _id="model:name:input",
                    _onchange=ajax)
        return inp

    def render_max_range_size(self, request, session):
        nareas = len(self.datamatrix.labels)
        values = self.datamatrix.data.values()
        if values:
            maxobs = max([ len(x) for x in values ])
        else:
            return ""
        v = self.max_range_size or max(maxobs, 2)
        self.max_range_size = v
        cb = URL(r=request, f="update_max_range_size")
        ajax = "selection('%s', 'model:max_range_size:input', "\
               "'model:max_range_size');" % cb
        opts = [ OPTION(x, _value=x) for x in \
                 range(2, nareas+1) ]
        return TABLE(
            TR(TD("Maximum number of areas in ancestral ranges:"),
               TD(SELECT(_name="max_range_size_select", value=v,
                         _id="model:max_range_size:input",
                         _onchange=ajax, *opts))),
            _style="border-spacing:0.25em;"
            )

class RangeList:
    def __init__(self, model):
        self.model = model
        self.ranges = set()
        self.excluded = set()
        self.visible = True

    def graph(self):
        m = [ x for x in self.ranges_sorted() if x ]
        edges = []
        for i, v in enumerate(m[:-1]):
            for i2, k in enumerate(m[i+1:]):
                j = i2+i+1
                intersect = set(v).intersection(set(k))
                if intersect in (set(v),set(k)):
                    edges.append((i,j))
                    edges.append((j,i))
                    print i, j, self.rng2str(v), self.rng2str(k)
        g = lagrange.graph.AreaGraph(edges=edges)
        ## g.write_svg("/tmp/tmp.svg", layout="circle")
        return g

    def errors(self):
        nareas = len(self.model.datamatrix.labels)
        if not self.ranges:
            return "no ranges have been specifed"
        v = []
        for r in self.ranges:
            for x in r:
                if x >= nareas:
                    v.append(r)
        if v:
            return "ranges in model mismatched to species data"

    def ranges_sorted(self):
        v = list(self.ranges)
        v.sort()
        return v

    def excluded_sorted(self):
        v = list(self.excluded)
        v.sort()
        return v

    def empty(self):
        self.ranges = set()
        self.excluded = set()

    def rng2str(self, v):
        labelsep = ""
        labels = self.model.datamatrix.labels
        for x in labels:
            if len(x) > 1:
                labelsep = "+"
                break
        return labelsep.join([ str(labels[i]) for i in v ])

    def render_ranges(self, request, session):
        ranges = self.ranges_sorted()
        if ranges:
            target = "rangelist_cell"
            opts = [ OPTION(self.rng2str(x), _value=str(x))
                     for i, x in enumerate(ranges)
                     if len(x) > 1 ]
            ninc = len(opts)
            inc_id = "model:rangelist:included"
            inc = SELECT(_name="included_ranges", _multiple=ON,
                         _id=inc_id, *opts)
            excluded = self.excluded_sorted()
            opts = [ OPTION(self.rng2str(x), _value=str(x))
                     for i, x in enumerate(excluded)
                     if len(x) > 1 ]
            nexc = len(opts)
            exc_id = "model:rangelist:excluded"
            exc = SELECT(_name="excluded_ranges", _multiple=ON,
                         _id=exc_id, *opts)
            sub = INPUT(_type="submit", _value="<-->")
            st = "padding-bottom:0.2em;"
            return FORM(
                INPUT(_type="hidden", _name="f", _value="index"),
                TABLE(
                TR(TH("Included (%s)" % ninc, _style=st),
                   TH(),
                   TH("Excluded (%s)" % nexc, _style=st)),
                TR(TD(inc, _style="text-align:center;"),
                   TD(sub, _style="vertical-align:middle;"),
                   TD(exc, _style="text-align:center;")),
                _style="margin-left:auto;margin-right:auto;"
                ),
                _action=URL(r=request, f="rangelist_exchange")
                )
        else:
            return "No ranges calculated"

    def render(self, request, session):
        if self.visible:
            return self.render_ranges(request, session)
        else:
            return ""
    

class AdjacencyMatrix:
    def __init__(self, model):
        self.model = model
        self.labels = []
        self.data = None
        self.view = "html"
        self.visible = True

    def check_contiguity(self, ranges):
        g = self.graph()
        return [ r for r in ranges if not g.are_connected(r) ]
            
    def graph(self):
        m = self.data
        rows, cols = scipy.array(m).shape
        assert rows == cols, "not equal: rows = %s, cols = %s" % (rows, cols)
        edges = []
        for i, v in enumerate(m):
            for j in range(i+1, cols):
                if v[j]:
                    edges.append((i,j))
        return lagrange.graph.AreaGraph(edges=edges)

    def toggle(self, request, session, i, j):
        current = self.data[i][j]
        new = int(not current)
        print "toggle before:", [ self.data[x][j] \
                                  for x in range(len(self.labels)) ]
        self.data[i][j] = new
        self.data[j][i] = new
        print "toggle after:", [ self.data[x][j] \
                                 for x in range(len(self.labels)) ]
        return new

    def default(self):
        self.labels = self.model.datamatrix.labels
        self.nlabels = len(self.labels)
        self.data = [ [1]*(self.nlabels) for x in self.labels ]
        ## for i in range(self.nlabels):
        ##     for j in range(0, i):
        ##         self.data[i][j] = 0

    def render_cell(self, request, session, i, j):
        checkbox_id = "model:adjacencymatrix:checkbox:%s:%s" % (i,j)
        target = "adjacencymatrix_cell"
        cb = URL(r=request, f="adjacencymatrix_click")
        onclick="adjacencymatrix_checkbox_clicked('%s', '%s:%s', '%s')" % \
                 (cb, i, j, target)
        v = bool(self.data[i][j])
        t = "%s--%s" % (self.labels[i],self.labels[j])
        inp = INPUT(_type="checkbox", value=v, _id=checkbox_id,
                    _onclick=onclick, _title=t)
        return inp

    def html(self, request, session):
        if (not self.labels) and self.model.datamatrix.labels:
            self.labels = self.model.datamatrix.labels
            self.default()
        labels = self.labels
        nareas = len(labels)
        R = range(nareas)
        rows = [TR([TH()]+[ TH(x) for x in labels ])]
        for i, l1 in enumerate(labels):
            r = [TH(l1, _id="model:adjacencymatrix:area:%s" % i)]
            for j in R:
                cell_id = "model:adjacencymatrix:cell:%s:%s" % (i,j)
                if i < j:
                    inp = self.render_cell(request, session, i, j)
                elif i == j:
                    inp = XML("&#8212;")
                else:
                    inp = ""
                r.append(TD(inp, _id=cell_id))
            st = "background-color:white; color:black"
            if i % 2:
                st = "background-color:lightgray; color:black"
            rows.append(TR(_style=st, *r))

        msgs = []
        if self.model.datamatrix.data:
            ranges = [ x for x in self.model.datamatrix.data.values()
                       if len(x) > 1 ]
            if ranges:
                try:
                    v = self.check_contiguity(ranges)
                    if v:
                        msg = " Discontiguous range(s) in data: %s" % \
                              ", ".join([ self.model.rangelist.rng2str(x)
                                          for x in v ])
                        alert = SPAN(
                            "Note:", _style="font-weight:bold;color:red;"
                            )
                        msgs.append(DIV(alert, msg))
                except InternalError:
                    alert = SPAN(
                        "Note:", _style="font-weight:bold;color:red;"
                        )
                    msgs.append(DIV(alert, "Observed area not in model ranges?"))
                    
        return DIV(TABLE(_class="adjacencymatrix", *rows), DIV(*msgs))

    def textarea(self, request, session):
        nareas = len(self.labels)
        lines = ["\t".join([""]+self.labels)]
        for i, lab in enumerate(self.labels):
            v = [ str(x) for x in self.data[i] ]
            lines.append("\t".join([lab]+v))
        inp = TEXTAREA(value="\n".join(lines), _rows=nareas,
                       _cols=max([ len(x)+(x.count("\t")*6) for x in lines ]),
                       _readonly=ON)
        return inp

    def render(self, request, session):
        if self.visible:
            if self.labels:
                if self.view == "html":
                    return self.html(request, session)
                else:
                    return self.textarea(request, session)
            else:
                #return self.upload_form(request, session)
                return ""
        else:
            return ""

    def upload_form(self, request, session):
        f = form_factory(
            SQLField("adjacencymatrix", "upload", requires=IS_NOT_EMPTY(),
                     label="Upload area adjacency matrix"),
            _action="upload_adjacencymatrix"
            )
        return f

class DispersalMatrixMP:
    "multi-period version"
    def __init__(self, model):
        self.model = model
        self.labels = []
        self.data = []
        self.ages = []
        self.view = "html"
        self.visible = True
        self.symmetric = True
        self.name = ""

    def errors(self):
        labels = self.model.datamatrix.labels
        nareas = len(labels)
        if not self.ages:
            return "no dispersal constraints have been specified"
        if self.labels != labels:
            return "labels mismatched with data matrix"
        for m in self.data:
            if len(m) != nareas:
                return "dimensions of matrix do not match number of areas"
            for r in m:
                if len(r) != nareas:
                    return "dimensions of matrix do not match number of areas"
        a = self.ages[-1]
        if self.model.treelist.trees:
            t = self.model.treelist.current_tree()
            age = t.root_age or \
                  max([ sum([ n.length for n in lf.rootpath() if n.parent ])
                        for lf in t.root.leaves() ])
            if a < age:
                return "the time span of dispersal constraints does not bracket the age of the root node"

    def durations2ages(self, durations):
        ages = [durations[0]]
        for v in durations[1:]:
            ages.append(v+ages[-1])
        self.ages = ages

    def ages2durations(self):
        n = 0.0
        v = []
        for i, age in enumerate(self.ages):
            v.append(age - n)
            n = age
        return v

    def set(self, request, session, i, j, v, sym=True):
        self.data[i][j] = v
        if sym:
            self.data[j][i] = v

    def default(self):
        self.labels = self.model.datamatrix.labels
        self.nlabels = len(self.labels)
        self.data = [
            [ [1.0]*(self.nlabels) for x in self.labels ]
            ]
        if self.ages:
            self.ages = self.ages[:1]
        else:
            self.ages = [10,]

    def append_clone(self, clone=-1, age=None):
        if age and age <= self.ages[clone]:
            session.flash = "Age cannot be younger than previous age"
            return
        m = [ x[:] for x in self.data[clone] ]
        self.data.append(m)
        self.ages.append(age or self.ages[clone]+10)

    def delete(self, i):
        del self.data[i]
        del self.ages[i]

    def render_symmetric(self, request, session):
        checkbox_id = "model:dm:symmetric"
        cb = URL(r=request, f="dm_symmetric_click")
        ajax="ajaxdirect('%s', '', '%s')" % (cb, checkbox_id)
        v = self.symmetric
        inp = INPUT(_type="checkbox", value=v, _id=checkbox_id,
                    _onchange=ajax)
        return inp

    def render_cell(self, request, session, p, i, j):
        inp_id = "model:dm:input:%s:%s:%s" % (p,i,j)
        target = "dm_cell"
        cb = URL(r=request, f="dm_edit_cell")
        ajax="dm_edit('%s', '%s', '%s', %s, %s, %s)" % \
              (cb, inp_id, target, p, i, j)
        v = self.data[p][i][j]
        labels = self.model.datamatrix.labels
        t = "Period %s: %s to %s" % (p, labels[i], labels[j])
        style = ""
        if not v:
            style = "background-color:lightgray;"
        inp = INPUT(_type="text", value=v, _id=inp_id,
                    _onchange=ajax, _size=4, _title=t, _style=style)
        return inp

    def add_form(self, request, session):
        ## url = URL(r=request, f="dm_add")
        ## target = "dm_cell"
        ## ajax = "ajaxdirect('%s','p=%s','%s');" % (url, p, target)
        f = FORM("Add dispersal matrix period with age: ",
                 BR(),
                 INPUT(_name="age", _size=3, _value=self.ages[-1]+10),
                 INPUT(_type="hidden", _name="f", _value="index"),
                 INPUT(_type="submit", _value="Submit"),
                 _action=URL(r=request, f="dm_add", vars=dict(f="index")))
        return f

    def html(self, request, session):
        if (not self.labels) and self.model.datamatrix.labels:
            self.labels = self.model.datamatrix.labels
            self.default()
        labels = self.labels
        nareas = len(labels)
        R = range(nareas)
        np = len(self.data)
        if np == 0:
            if self.model.datamatrix.data:
                u = URL(r=request,f='dm_default',vars=dict(f="index"))
                return A("[Create default matrix]", _href=u)
            else:
                return "No areas have been specified (suggestion: upload species range data)"
        rows = []
        for p in range(np):
            if p == 0: t0 = 0
            else: t0 = self.ages[p-1]
            duration = self.ages[p] - t0 

            cb = URL(r=request,f="dm_delete")
            query = "p=%s" % p
            target = "dm_cell"
            ajax="ajaxdirect('%s','%s','%s');" % (cb, query, target)
            d = A("[Delete]", _onclick=ajax, _style="cursor:pointer;",
                  _title="Delete time period %s" % p)

            rows.append(
                TR(TD("Time period %s: from %s to " % (p, t0),
                      self.render_age(request, session, p),
                      XML("&nbsp;&nbsp;"), d,
                      _colspan=nareas+1,
                      _style="padding:0.5em;text-align:left;"),
                   _style="background-color:lightgray;color:black;")
                )
            rows.append(TR([TH()]+[ TH(x) for x in labels ]))
            for i, l1 in enumerate(labels):
                r = [TH(l1, _id="model:dm:area:%s" % i)]
                for j in R:
                    cell_id = "model:dm:cell:%s:%s:%s"%(p,i,j)
                    if i == j:
                        inp = "--"
                    else:
                        inp = self.render_cell(request, session, p, i, j)
                    r.append(TD(inp, _id=cell_id))
                st = "background-color:white;color:black;"
                if not (i%2):
                    st = "background-color:lightgray;color:black;"
                rows.append(TR(_style=st, *r))
            rows.append(TR(TD(HR(), _colspan=nareas+1,
                              _style="padding:0.25em;")))
            if p == np-1:
                rows.append(
                    TR(TD(self.add_form(request, session),
                          _colspan=nareas+1,
                          _style="text-align:left;"))
                    )
            ## else:
            ##     rows.append(
            ##         TR(TD(self.render_age(request, session, p),
            ##               _colspan=nareas+1))
            ##         )
                
        return TABLE(_class="adjacencymatrix", *rows)

    def render_age(self, request, session, p):
        cb = URL(r=request, f="dm_update_age")
        source = "dm_age_input_%s" % p
        target = "dm_age_%s" % p
        ajax = "update_textinput_by_index('%s','%s','%s',%s)" %\
               (cb, source, target, p)
        inp = INPUT(_id=source, _value=self.ages[p], _size=3,
                    _onchange=ajax,
                _title="The lower bound of this time period")
        #return SPAN("Age: ", inp, _id=target)
        return SPAN(inp, _id=target)

    def render(self, request, session):
        return self.html(request, session)


class TreeList:
    def __init__(self, model):
        self.model = model
        self.trees = []
        self.view = "html"
        self.visible = True
        self.vi = 0

    def toolbar(self, request, session):
        session.tli = int(request.vars.tli or session.tli or 0)
        v = []
        for i, x in enumerate(
            ("Tree",
             "Node constraints",
             "Fossils")
            ):
            u = URL(r=request,f="index",vars=dict(tli=i))
            a = A(x, _href=u)
            v.append(TD(a))
        v[session.tli].update(_class="selected")
        return TABLE(TR(*v), _id="treelist_toolbar",
                     _class="toolbar")

    def current_tree(self):
        return self.trees[self.vi]

    def errors(self):
        try:
            s = self.trees[0].newick
        except:
            return "no tree specified"
        t = self.trees[0]
        return t.errors()

    def append(self, s, name=None, age=None):
        s = s.strip()
        t = treeview.Tree(self.model)
        t.newick = s
        t.root_age = age
        nt = len(self.trees)
        if not name:
            while 1:
                name = "Tree%s" % nt
                if name in [ t.name for t in self.trees ]:
                    nt += 1
                else:
                    break
        t.name = name
        self.trees.append(t)

    def remove(self, i):
        try:
            self.trees.remove(i)
        except:
            pass

    def create_datamatrix_form(self, request, session, i=0):
        url = URL(r=request,f="datamatrix_default_from_tree")
        nareas = INPUT(_size=2,_name="nareas",_value="2")
        form = FORM(
            INPUT(_type="hidden",_name="t",_value=i),
            INPUT(_type="hidden",_name="f",_value="index"),
            "Create new matrix of range data with ", nareas, " areas: ",
            INPUT(_type="submit", _value="Submit"),
            _action=url
            )
        return form
        
    def html(self, request, session):
        v = []
        for i, t in enumerate(self.trees):
            delete = A("[Delete tree]",
                       _href=URL(r=request,f="delete_tree?t=%s"%i),
                       _title="Delete this tree")
            d = TR(
                TD("Root age: %g" % (t.root_age or t.calibrate())),
                TD(delete, _style="text-align:right;")
                #TD(delete)
                )
            v.append(d)
        if not self.model.datamatrix.data:
            v.append(
                TR(TD(self.create_datamatrix_form(request, session)))
                )
        return DIV(
            self.tree_instructions(request, session),
            TABLE(_id="model:treelist", _style="width:100%;", *v),
            BR(),
            DIV(self.trees[self.vi].render(request, session),
                _id="treecss")
            )

    def tree_instructions(self, request, session):
        s = (
            B("Instructions:"),
            """
            select the nodes at which you wish to estimate ancestral
            ranges.  The default is to estimate them for all nodes.
            """,
            BR(_style="margin-bottom:0.5em;"),
            """
            Note: the tree below may not display nicely in Internet
            Explorer (try Firefox).
            """
            )
        return DIV(_class="instructions", *s)

    def render_nodeconstraints(self, request, session):
        return DIV("Range constraints for individual nodes are possible in Lagrange, but configuring them here is not yet implemented.",
                   _class="instructions")

    def render_fossils(self, request, session):
        return DIV("Branch-specific fossil constraints on ancestral ranges are possible in Lagrange, but configuring them here is not yet implemented.",
                   _class="instructions")

    def render(self, request, session):
        tli = session.tli or request.tli or 0
        if self.visible:
            if self.trees:
                if tli == 0:
                    return self.html(request, session)
                elif tli == 1:
                    return self.render_nodeconstraints(request, session)
                elif tli == 2:
                    return self.render_fossils(request, session)
                else:
                    return "if you see this, please file a bug report"
            else:
                return self.upload_form(request, session)
        else:
            return ""

    def upload_form(self, request, session):
        #name = SPAN("Name (optional): ", INPUT(_name="name", _size=16))
        age = SPAN("Age of the root node (if branch lengths are not in units of time): ",
                   INPUT(_name="age", _size=8, _id="newick_age"))
        ta = TEXTAREA(_type="text", _rows=6, _cols=70, _name="newick_string",
                      _id="newick_textarea")
        textarea = DIV("Paste the tree here:", BR(), ta)
        submit = INPUT(_type="submit", _value="Add newick tree")
        action = URL(r=request, f="upload_newick")
        f = FORM(#name,
                 #BR(_style="margin-bottom:0.5em;"),
                 textarea,
                 BR(_style="margin-bottom:0em;"),
                 age,
                 BR(_style="margin-bottom:0.5em;"),
                 submit, _action=action,
                 _id="treelist_upload")
        return f

class DataMatrix:
    def __init__(self, model):
        self.model = model
        self.labels = []
        self.taxa = []
        self.data = {}
        self.view = "html"
        self.visible = True

    def verify_rangelist(self):
        rangelist = [ x for x in self.model.rangelist.ranges_sorted() if x ]
        v = list(set([ d for d in self.data.values() if d ]))
        return [ self.model.rangelist.rng2str(rng)
                 for rng in v if rng not in rangelist ]

    def verify_connectivity(self):
        unconnected = []
        g = self.model.rangelist.graph()
        rangelist = [ x for x in self.model.rangelist.ranges_sorted() if x ]
        v = list(set([ d for d in self.data.values() if d ]))
        v.sort()
        for i, rng1 in enumerate(v[:-1]):
            assert rng1 in rangelist
            rng1_i = rangelist.index(rng1)
            for rng2 in v[i+1:]:
                assert rng2 in rangelist
                rng2_i = rangelist.index(rng2)
                if not g.edge_connectivity(rng1_i, rng2_i):
                    unconnected.append((rng1_i, rng2_i))
        return [ "%s-%s" % (self.model.rangelist.rng2str(rangelist[i]),
                            self.model.rangelist.rng2str(rangelist[j]))
                 for i, j in unconnected ]

    def errors(self):
        if not self.labels:
            return "no area labels in species range matrix"
        if not self.data:
            return "no range data for species"
        for x in self.labels:
            if not x:
                return "empty area label in species range matrix"
        for k, v in self.data.items():
            if not v:
                return "no range data for taxon %s" % k
            for x in v:
                try:
                    assert self.labels[x]
                except:
                    return "problem with range for taxon %s: value is %s" % (k, x)
        not_in_rangelist = self.verify_rangelist()
        if not_in_rangelist:
            return "observed species ranges in data matrix have been excluded: %s" % ", ".join(not_in_rangelist)

        unconnected = self.verify_connectivity()
        if unconnected:
            s = ", ".join(unconnected)
            return "the following species range pairs are disconnected: %s" % s

    def serialize_data(self):
        return dict([ (k, tuple(v)) for k, v in self.data.items() ])

    def default_from_tree(self, i, nareas):
        t = newick.parse(self.model.treelist.trees[i].newick)
        self.taxa = [ n.label for n in t.leaves() ]
        self.labels = string.uppercase[:nareas]
        self.nlabels = len(self.labels)
        self.data = dict(
            [ (x, [False]*nareas) for x in self.taxa ]
            )
        self.data = dict([ (x, []) for x in self.taxa ])
        self.model.adjacencymatrix.default()
        self.model.dm.default()

    def toggle(self, request, session, i, j):
        taxa = self.taxa
        taxon = taxa[i]
        v = list(self.data[taxon])
        if j in v:
            v.remove(j)
            flag = False
        else:
            v.append(j); v.sort()
            flag = True
        self.data[taxon] = tuple(v)
        return flag

    def render_cell(self, request, session, i, j):
        checkbox_id = "model:datamatrix:checkbox:%s:%s" % (i,j)
        taxa = self.taxa
        taxon = taxa[i]
        dist = self.data[taxon]
        v = (j in dist)
        cb = URL(r=request, f="datamatrix_click")
        onclick="cell_checkbox_clicked('%s', '%s:%s', '%s')" % \
                 (cb, i, j, checkbox_id)
        t = "%s: %s" % (taxon, self.labels[j])
        inp = INPUT(_type="checkbox", value=v, _id=checkbox_id,
                    _onclick=onclick, _title=t)
        return inp

    def render_taxon(self, request, session, i, edit=True):
        taxa = self.taxa
        w = max([ len(x) for x in taxa ])
        t = taxa[i]
        if not edit:
            return t
        else:
            cb = URL(r=request, f="update_taxon")
            target = "model:datamatrix:taxon:%s" % i
            source = "model:datamatrix:taxon:input:%s" % i
            ajax = "update_textinput_by_index('%s', '%s', '%s', %s);"\
                   % (cb, source, target, i)
            inp = INPUT(_value=t, _id=source, _onchange=ajax, _size=w)
            return inp

    def render_area_label(self, request, session, i, edit=True):
        labels = self.labels
        w = max([ len(x) for x in labels ])
        lab = labels[i]
        if not edit:
            return lab
        else:
            cb = URL(r=request, f="update_area_label")
            target = "model:datamatrix:arealabel:%s" % i
            source = "model:datamatrix:arealabel:input:%s" % i
            ajax = "update_textinput_by_index('%s', '%s', '%s', %s);"\
                   % (cb, source, target, i)
            inp = INPUT(_value=lab, _id=source, _onchange=ajax, _size=w)
            return inp

    def render(self, request, session):
        if self.visible:
            if self.labels:
                if self.view == "html":
                    return self.html(request, session)
                else:
                    return self.textarea(request, session)
            else:
                return self.upload_form(request, session)
        else:
            return ""

    def textarea(self, request, session):
        nareas = len(self.labels)
        taxa = self.taxa
        lines = ["\t".join([""]+self.labels)]
        for t in taxa:
            dist = self.data[t]
            v = ["0"]*nareas
            for x in dist:
                v[x] = "1"
            lines.append("\t".join([t]+v))
        #width = max([ len(x) for x in lines ])
        inp = TEXTAREA(value="\n".join(lines), _rows=30,_cols=80,
                       _readonly=ON)
        return inp
        
    def html(self, request, session):
        cb = URL(r=request, f="datamatrix_click")
        labels = self.labels; data = self.data
        nareas = len(labels)
        taxa = self.taxa
        ntax = len(taxa)
        labelrange = range(nareas)
        #rows = [TR([TH()]+[ TH(x) for x in labels ])]
        d = URL(r=request,f="delete_datamatrix?f=index")
        a = A("[Delete matrix]", _href=d, _title="Delete species range data")
        rows = [TR([TH(a)]+
                   [ TH(self.render_area_label(request, session, i))
                     for i in labelrange ])]
        for i, taxon in enumerate(taxa):
            #r = [TD(taxon, _id="model:datamatrix:taxon:%s" % i)]
            r = [TD(self.render_taxon(request, session, i),
                    _id="model:datamatrix:taxon:%s" % i)]
            for j in range(nareas):
                cell_id = "model:datamatrix:cell:%s:%s" % (i,j)
                inp = self.render_cell(request, session, i, j)
                r.append(TD(inp, _id=cell_id, _style="text-align:center;"))
            st = "background-color:white; color:black"
            if i % 2:
                st = "background-color:lightgray; color:black"
            rows.append(TR(_style=st, *r))
        return TABLE(*rows)

    def upload_form(self, request, session):
        ## f = form_factory(
        ##     SQLField("datamatrix", "upload", requires=IS_NOT_EMPTY(),
        ##              label="Upload data matrix of taxon ranges"),
        ##     _action="upload_datamatrix"
        ##     )
        up = INPUT(_type="file", _name="datamatrix")
        sub = INPUT(_type="submit", _value="Submit")
        t = TABLE(
            TR(
            TD("Upload your text file of species range data:"), TD(up)
            ),
            TR(
            TD(), TD(sub)
            ),
            _style="border-spacing:0.25em;"
            )
        f = FORM(t, _action="upload_datamatrix")
        return f
