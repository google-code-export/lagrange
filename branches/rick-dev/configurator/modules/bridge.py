from gluon.sqlhtml import *
from gluon.sql import *
from gluon.storage import Storage
import cPickle, pprint, treeview, string, newick
import lagrange, scipy

class DECModel:
    def __init__(self, session):
        self.name = "Untitled%s" % len(session.models or [])
        self.datamatrix = DataMatrix(self)
        self.adjacencymatrix = AdjacencyMatrix(self)
        #self.dispersalmatrix = DispersalMatrix(self)
        self.dm = DispersalMatrixMP(self)
        self.treelist = TreeList(self)
        self.rangelist = RangeList(self)
        self.max_range_size = 0

    def render(self, request, session):
        ti = session.ti
        if ti == 0:
            return self.render_instructions(request, session)
        elif ti == 1:
            #data matrix
            pass
        elif ti == 2:
            # tree
            pass
        elif ti == 3:
            # dispersal model
            pass
        
        return ""

    def render_instructions(self, request, session):
        e = self.datamatrix.errors()

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
%s
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
model, tree, data = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
lagrange.output.ancsplits(outfile, tree, model, d, e, tee=True)
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
            "newick_trees": [ dict(t) for t in self.treelist.trees ],
            "ranges": self.rangelist.ranges_sorted(),
            "excluded_ranges": self.rangelist.excluded_sorted(),
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
        self.treelist.trees = [ Storage(t) for t in d["newick_trees"] ]
        self.rangelist.ranges = set(d["ranges"])
        self.rangelist.excluded = set(d["excluded_ranges"])

    def render_name(self, request, session):
        cb = URL(r=request, f="update_modelname")
        ajax = "textinput_enter('%s', 'model:name:input', 'model:name');" % cb
        ## inp = INPUT(_value=self.name, _id="model:name:input",
        ##             _onchange=ajax)
        ## inp = FORM(INPUT(_value=self.name, _name="s",
        ##                  _onchange="$('#modelnamesubmitbutton').show();"),
        ##            INPUT(_type="submit", _value="Update name",
        ##                  _id="modelnamesubmitbutton",
        ##                  _style="display:none"),
        ##            _action="update_modelname_form")
        inp = FORM(INPUT(_value=self.name, _name="s"),
                   INPUT(_type="hidden", _name="f", _value=request.function),
                   _action="update_modelname_form",
                   _title="Edit name and hit Enter")
        return inp

    ## def render_maxareas(self, request, session):
    ##     cb = URL(r=request, f="update_maxareas")
    ##     ajax = "textinput_enter('%s', 'model:maxareas:input', "\
    ##            "'model:maxareas');" % cb
    ##     inp = INPUT(_value=self.maxareas, _id="model:maxareas:input",
    ##                 _onchange=ajax)
    ##     return inp

    def render_max_range_size(self, request, session):
        nareas = len(self.datamatrix.labels)
        maxobs = max([ len(x) for x in self.datamatrix.data.values() ])
        v = self.max_range_size or max(maxobs, 2)
        self.max_range_size = v
        cb = URL(r=request, f="update_max_range_size")
        ajax = "selection('%s', 'model:max_range_size:input', "\
               "'model:max_range_size');" % cb
        opts = [ OPTION(x, _value=x) for x in \
                 range(2, nareas+1) ]
        return SELECT(
            _name="max_range_size_select", value=v,
            _id="model:max_range_size:input", _onchange=ajax, *opts
            )

class RangeList:
    def __init__(self, model):
        self.model = model
        self.ranges = set()
        self.excluded = set()
        self.visible = True

    def errors(self):
        nareas = len(self.model.datamatrix.labels)
        if not self.ranges:
            return "ranges missing"
        v = []
        for r in self.ranges:
            for x in r:
                if x >= nareas:
                    v.append(r)
        if v:
            return "range values mismatched to data matrix"

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
            inc_id = "model:rangelist:included"
            inc = SELECT(_name="included_ranges", _multiple=ON,
                         _id=inc_id, *opts)
            excluded = self.excluded_sorted()
            opts = [ OPTION(self.rng2str(x), _value=str(x))
                     for i, x in enumerate(excluded)
                     if len(x) > 1 ]
            exc_id = "model:rangelist:excluded"
            exc = SELECT(_name="excluded_ranges", _multiple=ON,
                         _id=exc_id, *opts)
            sub = INPUT(_type="submit", _value="<-->")
            return FORM(
                INPUT(_type="hidden", _name="f", _value=request.function),
                TABLE(TR(TH("Included"), TH(), TH("Excluded")),
                      TR(TD(inc),
                         TD(sub, _style="vertical-align:middle;"),
                         TD(exc))),
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
        v = self.data[i][j]
        self.data[i][j] = int(not v)
        self.data[j][i] = int(not v)
        return bool(v)

    def default(self):
        if not self.labels:
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
        inp = INPUT(_type="checkbox", value=v, _id=checkbox_id,
                    _onclick=onclick)
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
            r = [TD(l1, _id="model:adjacencymatrix:area:%s" % i)]
            for j in R:
                cell_id = "model:adjacencymatrix:cell:%s:%s" % (i,j)
                if i < j:
                    inp = self.render_cell(request, session, i, j)
                elif i == j:
                    inp = "--"
                else:
                    inp = ""
                r.append(TD(inp, _id=cell_id))
            st = "background-color:white; color:black"
            if i % 2:
                st = "background-color:lightgray; color:black"
            rows.append(TR(_style=st, *r))

        if self.model.datamatrix.data:
            ranges = [ x for x in self.model.datamatrix.data.values()
                       if len(x) > 1 ]
            if ranges:
                v = self.check_contiguity(ranges)
                if v:
                    msg = " Discontiguous ranges in data matrix: %s" % \
                          ", ".join([ self.model.rangelist.rng2str(x)
                                      for x in v ])
                    alert = SPAN(
                        "Note:", _style="font-weight:bold;color:red;"
                        )
                    rows.append(
                        TR(TD(alert, msg, _colspan=nareas+1))
                        )
        return TABLE(*rows)

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

class DispersalMatrix:
    def __init__(self, model):
        self.model = model
        self.labels = []
        self.data = None
        self.view = "html"
        self.visible = True
        self.symmetric = True
        self.age = 0
        self.name = ""

    def set(self, request, session, i, j, v, sym=True):
        self.data[i][j] = v
        if sym:
            self.data[j][i] = v

    def default(self):
        if not self.labels:
            self.labels = self.model.datamatrix.labels
        self.nlabels = len(self.labels)
        self.data = [ [1.0]*(self.nlabels) for x in self.labels ]

    def render_symmetric(self, request, session):
        checkbox_id = "model:dispersalmatrix:symmetric"
        cb = URL(r=request, f="dispersalmatrix_symmetric_click")
        ajax="ajaxdirect('%s', '', '%s')" % (cb, checkbox_id)
        v = self.symmetric
        inp = INPUT(_type="checkbox", value=v, _id=checkbox_id,
                    _onchange=ajax)
        return inp

    def render_cell(self, request, session, i, j):
        inp_id = "model:dispersalmatrix:input:%s:%s" % (i,j)
        target = "dispersalmatrix_cell"
        cb = URL(r=request, f="dispersalmatrix_edit_cell")
        ajax="matrix_edit('%s', '%s', '%s', %s, %s)" % \
              (cb, inp_id, target, i, j)
        v = self.data[i][j]
        inp = INPUT(_type="text", value=v, _id=inp_id,
                    _onchange=ajax, _size=4)
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
            r = [TD(l1, _id="model:dispersalmatrix:area:%s" % i)]
            for j in R:
                cell_id = "model:dispersalmatrix:cell:%s:%s" % (i,j)
                if i == j:
                    inp = "--"
                else:
                    inp = self.render_cell(request, session, i, j)
                r.append(TD(inp, _id=cell_id))
            st = "background-color:white; color:black"
            if i % 2:
                st = "background-color:lightgray; color:black"
            rows.append(TR(_style=st, *r))
        rows.append(TR(TD(self.render_age(request, session),
                          _colspan=nareas+1)))
        return TABLE(*rows)

    def render_age(self, request, session):
        inp = INPUT(_name="dispersalmatrix_age", _value=self.age, _size=3)
        return SPAN("Age: ", inp)

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
        if self.visible and self.data:
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
            SQLField("dispersalmatrix", "upload", requires=IS_NOT_EMPTY(),
                     label="Upload area dispersal matrix"),
            _action="upload_dispersalmatrix"
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
        if self.labels != labels:
            return "labels mismatched with data matrix"
        for m in self.data:
            if len(m) != nareas:
                return "dimensions of matrix do not match number of areas"
            for r in m:
                if len(r) != nareas:
                    return "dimensions of matrix do not match number of areas"

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
        if not self.labels:
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
        inp = INPUT(_type="text", value=v, _id=inp_id,
                    _onchange=ajax, _size=4)
        return inp

    def add_form(self, request, session):
        ## url = URL(r=request, f="dm_add")
        ## target = "dm_cell"
        ## ajax = "ajaxdirect('%s','p=%s','%s');" % (url, p, target)
        f = FORM("Add dispersal matrix period with age: ",
                 INPUT(_name="age", _size=3, _value=self.ages[-1]+10),
                 INPUT(_type="hidden", _name="f",
                       _value=request.function),                 
                 _action=URL(r=request, f="dm_add"))
        return f

    def html(self, request, session):
        if (not self.labels) and self.model.datamatrix.labels:
            self.labels = self.model.datamatrix.labels
            self.default()
        labels = self.labels
        nareas = len(labels)
        R = range(nareas)
        rows = []
        np = len(self.data)
        for p in range(np):
            if p == 0: t0 = 0
            else: t0 = self.ages[p-1]
            duration = self.ages[p] - t0 
            rows.append(
                TR(TD("Period %s: duration %s" % (p, duration),
                      _colspan=nareas+1),
                   _style="background-color:lightgray;color:black;")
                )
            rows.append(TR([TH()]+[ TH(x) for x in labels ]))
            for i, l1 in enumerate(labels):
                r = [TD(l1, _id="model:dm:area:%s" % i)]
                for j in R:
                    cell_id = "model:dm:cell:%s:%s:%s"%(p,i,j)
                    if i == j:
                        inp = "--"
                    else:
                        inp = self.render_cell(request, session, p, i, j)
                    r.append(TD(inp, _id=cell_id))
                st = "background-color:white; color:black"
                if i % 2:
                    st = "background-color:lightgray; color:black"
                rows.append(TR(_style=st, *r))
            if p == np-1:
                rows.append(
                    TR(TD(self.render_age(request, session, p),
                          self.add_form(request, session),
                          _colspan=nareas+1))
                    )
            else:
                rows.append(
                    TR(TD(self.render_age(request, session, p),
                          _colspan=nareas+1))
                    )
                
        return TABLE(*rows)

    def render_age(self, request, session, p):
        cb = URL(r=request, f="dm_update_age")
        source = "dm_age_input_%s" % p
        target = "dm_age_%s" % p
        ajax = "update_textinput_by_index('%s','%s','%s',%s)" %\
               (cb, source, target, p)
        inp = INPUT(_id=source, _value=self.ages[p], _size=3,
                    _onchange=ajax)
        cb = URL(r=request,f="dm_delete")
        query = "p=%s" % p
        target = "dm_cell"
        ajax="ajaxdirect('%s','%s','%s');" % (cb, query, target)
        d = A("[Delete]", _onclick=ajax, _style="cursor:pointer;")
        return SPAN("Age: ", inp, XML("&nbsp;"), d, _id=target)

    def render(self, request, session):
        if self.visible and self.data and self.labels:
            return self.html(request, session)
        return ""


class TreeList:
    def __init__(self, model):
        self.model = model
        self.trees = []
        self.view = "html"
        self.visible = True

    def errors(self):
        try:
            s = self.trees[0].newick
        except:
            return "no tree specified"
        n = newick.parse(s)
        leaves = [ lf.label for lf in n.leaves() ]
        for x in self.model.datamatrix.data.keys():
            if x not in leaves:
                return "%s not found among tree tips" % x

    def append(self, s, name=None, age=None):
        s = s.strip()
        t = Storage()
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

    def html(self, request, session):
        v = []
        for i, t in enumerate(self.trees):
            ## name = SPAN("Name: ", INPUT(_value=t.name))
            ## age = SPAN(
            ##     "Root age: ", INPUT(_name="age", _size=8, _value=t.age),
            ##     _style="padding-left:1em;"
            ##     )
            ## inp = TEXTAREA(value=t.newick, _rows=6, _cols=70,
            ##                _readonly=ON)
            ## v.append(DIV(name, age, BR(), inp))
            cb = URL(r=request,c="treeview",f="index?t=%s"%i)
            onclick = "popup_winsize('%s',600,800);" % cb
            name = A(t.name, _onclick=onclick, _style="cursor:pointer;")
            cb = URL(r=request,f=request.function)
            delete = A("[Delete]",
                       _href=URL(r=request,f="delete_tree?t=%s"%i))
            ## opts = [ OPTION(x, _value=x) for x in range(2, 9) ]
            ## nareas = SELECT(_name="nareas",*opts)
            url = URL(r=request,f="datamatrix_default_from_tree")
            nareas = INPUT(_size=2,_name="nareas",_value="2")
            form = FORM(
                INPUT(_type="hidden",_name="t",_value=i),
                INPUT(_type="hidden",_name="f",_value=request.function),
                "[Create empty data matrix with ", nareas, " areas]",
                _action=url
                )
            st = "padding-left:0.5em;"
            st2 = "vertical-align:bottom;"
            d = TR(TD("%s: "%(i+1), _style=st+st2),
                   TD(name, _style=st2),
                   TD(delete, _style=st+st2),
                   TD(form, _style=st+st2))
            v.append(d)
        return TABLE(_id="model:treelist", *v)

    def render(self, request, session):
        if self.visible:
            if self.trees:
                return self.html(request, session)
            else:
                return self.upload_form(request, session)
        else:
            return ""

    def upload_form(self, request, session):
        name = SPAN("Name (optional): ", INPUT(_name="name", _size=16))
        age = SPAN("Root age (if branch lengths not in units of time): ",
                   INPUT(_name="age", _size=8))
        textarea = DIV(
            "Paste newick tree below (note: should be ultrametric):", BR(),
            TEXTAREA(_type="text", _rows=6, _cols=70, _name="newick_string")
            )
        f = INPUT(_type="hidden", _name="f",
                  _value=request.function)
        submit = INPUT(_type="submit", _value="Add newick tree")
        action = URL(r=request, f="upload_newick")
        f = FORM(name, BR(), textarea, age, BR(), submit, _action=action)
        return f

class DataMatrix:
    def __init__(self, model):
        self.model = model
        self.labels = []
        self.taxa = []
        self.data = {}
        self.view = "html"
        self.visible = True

    def errors(self):
        for x in self.labels:
            if not x:
                return "empty area label"
        for k, v in self.data.items():
            if not v:
                return "no data for taxon %s" % k
            for x in v:
                try:
                    assert self.labels[x]
                except:
                    return "problem with range: value %s" % x

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
        self.data[taxon] = v
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
        inp = INPUT(_type="checkbox", value=v, _id=checkbox_id,
                    _onclick=onclick)
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
        rows = [TR([TH()]+
                   [ TH(self.render_area_label(request, session, i))
                     for i in labelrange ])]
        for i, taxon in enumerate(taxa):
            #r = [TD(taxon, _id="model:datamatrix:taxon:%s" % i)]
            r = [TD(self.render_taxon(request, session, i),
                    _id="model:datamatrix:taxon:%s" % i)]
            for j in range(nareas):
                cell_id = "model:datamatrix:cell:%s:%s" % (i,j)
                inp = self.render_cell(request, session, i, j)
                r.append(TD(inp, _id=cell_id))
            st = "background-color:white; color:black"
            if i % 2:
                st = "background-color:lightgray; color:black"
            rows.append(TR(_style=st, *r))
        return TABLE(*rows)

    def upload_form(self, request, session):
        f = form_factory(
            SQLField("datamatrix", "upload", requires=IS_NOT_EMPTY(),
                     label="Upload data matrix of taxon ranges"),
            _action="upload_datamatrix"
            )
        return f
