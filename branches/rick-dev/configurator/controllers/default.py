import string, cPickle, sys, pprint, string
#lagrange = local_import("lagrange", reload=True)
import lagrange
import models
import treeview
import newick
## models = local_import("models", reload=True)
## treeview = local_import("treeview", reload=True)
## newick = local_import("newick", reload=True)
from gluon.storage import Storage
from gluon.sqlhtml import form_factory
from cStringIO import StringIO
from types import StringTypes

def toolbar():
    session.ti = int(request.vars.ti or session.ti or 0)
    m = current_model()
    v = []
    for i, x in enumerate(
        (("Species ranges","Upload/configure species range data"),
         ("Phylogenetic tree", "Upload the phylogenetic tree and specify nodes at which you want to reconstruct ancestral ranges"),
         ("Range constraints", "Specify constraints on ancestral ranges, based on spatial layout of areas and other considerations"),
         ("Dispersal constraints", "Scale the base dispersal rate between areas by fixed values, to reflect prior beliefs about the probability of lineage movement"),
         ("Rate parameters", "Configure how rates of dispersal and local extinction are estimated/fixed"),
         ("Save/Download", "Create an executable Python script that performs the analysis"))
        ):
        x, instr = x
        u = URL(r=request,f="index",vars=dict(ti=i))
        a = A(x, _href=u, _title=instr)
        v.append(TD(a))
    v[session.ti].update(_class="selected")
    return TABLE(TR(*v), _id="toolbar")

def rangelist_exchange():
    data = current_model().datamatrix.data
    m = current_model().rangelist
    inc = m.ranges_sorted()
    exc = m.excluded_sorted()
    ## inc2exc = [ inc[int(i)] for i in request.vars.included_ranges or [] ]
    ## exc2inc = [ exc[int(i)] for i in request.vars.excluded_ranges or [] ]
    var_inc = request.vars.included_ranges or []
    var_exc = request.vars.excluded_ranges or []
    if type(var_inc) in StringTypes and "import" not in var_inc:
        inc2exc = [eval(var_inc, {}, {})]
    else:
        inc2exc = [ eval(x) for x in var_inc ]
    if type(var_exc) in StringTypes and "import" not in var_exc:
        exc2inc = [eval(var_exc, {}, {})]
    else:
        exc2inc = [ eval(x, {}, {}) for x in var_exc if "import" not in x ]

    for x in inc2exc:
        print x, data
        if data and (x in data.values()):
            session.flash = "Cannot exclude ranges occupied by species"
            continue
        try:
            m.ranges.remove(x); m.excluded.add(x)
        except:
            pass
    for x in exc2inc:
        try:
            m.excluded.remove(x); m.ranges.add(x)
        except:
            pass
    redirect(URL(r=request, f=request.vars.f or "index"))

def delete_current_model():
    m = current_model()
    session.models.remove(m)
    u = URL(r=request, f=request.vars.f)
    redirect(u)

def delete_treelist():
    m = current_model()
    m.treelist.trees = []
    redirect(URL(r=request, f=request.vars.f or "index"))

def delete_tree():
    m = current_model()
    try:
        i = int(request.vars.t)
        del m.treelist.trees[i]
    except:
        session.flash = "Could not delete tree '%s'" % request.vars.i
    redirect(URL(r=request, f=request.vars.f or "index"))

def delete_dm_all():
    m = current_model()
    m.dm = models.DispersalMatrixMP(m)
    redirect(URL(r=request, f=request.vars.f or "index"))

def dispersalmatrix_symmetric_click():
    m = current_model().dm
    m.symmetric = bool(not m.symmetric)
    return m.render_symmetric(request, session)

def tree_view_stretch():
    t = current_model().treelist.current_tree()
    t.wscale += 0.2
    return t.render(request, session)

def tree_view_squeeze():
    t = current_model().treelist.current_tree()
    t.wscale -= 0.2
    return t.render(request, session)

def tree_view_select_all():
    t = current_model().treelist.current_tree()
    t.included = "__all__"
    return t.render(request, session)

def tree_view_select_none():
    t = current_model().treelist.current_tree()
    t.included = []
    return t.render(request, session)

def nodecheck_clicked():
    m = current_model()
    t = m.treelist.current_tree()
    for k, v in request.vars.items():
        if k.startswith("nodecheck"):
            label = k.split("_")[-1]
            v = not t.is_included(label)
            t.set_included(label, v)
    return m.treelist.render(request, session)

def tree_scaled_checkbox_clicked():
    t = current_model().treelist.current_tree()
    t.scaled_view = bool(not t.scaled_view)
    return t.render(request, session)

def dm_update_age():
    age = float(request.vars.s)
    p = int(request.vars.i)
    m = current_model().dm
    np = len(m.ages)
    if p == 0:
        upper = 0
    else:
        upper = m.ages[p-1]
    if p == np-1:
        lower = None
    else:
        lower = m.ages[p+1]
    if (age > upper) and (lower and (age < lower)) \
           or (lower is None):
        m.ages[p] = age
    return m.render_age(request, session, p)

def dm_delete():
    m = current_model().dm
    p = int(request.vars.p or 0)
    m.delete(p)
    return m.render(request, session)

def dm_add():
    m = current_model().dm
    age = float(request.vars.age or m.ages[-1]+10)
    m.append_clone(age=age)
    redirect(URL(r=request, f=request.vars.f or "index"))

def dm_symmetric_click():
    m = current_model().dm
    m.symmetric = bool(not m.symmetric)
    return m.render_symmetric(request, session)

def dm_edit_cell():
    try:
        v = float(request.vars.s)
    except:
        v = 1
    try:
        p = int(request.vars.p or 0)
        i = int(request.vars.i)
        j = int(request.vars.j)
        m = current_model().dm
        m.data[p][i][j] = v
        if m.symmetric:
            m.data[p][j][i] = v
        return m.render(request, session)
    except:
        redirect("index")

def upload_newick():
    #n = "_".join(request.vars.name.strip().split())
    n = None
    age = request.vars.age.strip()
    if age:
        try:
            age = float(age)
        except:
            session.flash = "Age must be a number"
    s = request.vars.newick_string.strip()
    if s:
        r = None
        try:
            r = newick.parse(s)
        except:
            session.flash = "Problem with newick string"
        if r:
            m = current_model()
            m.treelist.append(s, name=n, age=age)
            t = m.treelist.trees[-1]
            t.parse()
            if not t.scaled:
                for n in t.root.iternodes():
                    if n.parent and (n.length is None):
                        n.length = 1
                t.calibrate()
                t.scaled = True
                t.newick = newick.tostring(t.root)+";"
                session.flash = "Branches without lengths assigned a value of 1"
            if not t.is_ultrametric():
                session.flash = "Tree is not ultrametric"
    redirect(URL(r=request, f=request.vars.f or "index"))

def upload_pickle():
    f = request.vars.pickle
    m = cPickle.loads(f.file.read())
    session.models.append(m)
    redirect(URL(r=request, f=request.vars.f or "index"))

def upload_serial():
    #d = dict(eval(request.vars.d.file.read()))
    m = models.DECModel(session)
    m.restore(d)
    session.models.append(m)
    redirect(URL(r=request, f=request.vars.f or "index"))

def upload_analysis():
    v = []
    flag = False
    try:
        for line in [ x for x in request.vars.d.file.readlines()
                      if "import" not in x ]:
            if line.startswith("### begin data"):
                flag = True
                continue
            if line.startswith("### end data"):
                break
            if flag:
                v.append(line)
        s = "".join(v)
        d = dict(eval(s, {}, {}))
    except:
        session.flash = "Cannot parse data"
        redirect(URL(r=request, f=request.vars.f or "index"))
    version = d["lagrange_version"]
    if version < lagrange.VERSION:
        session.flash = "Scripts created in earler versions of this web site may not upload correctly"
    m = models.DECModel(session)
    m.restore(d)
    session.models.append(m)
    redirect(URL(r=request, f=request.vars.f or "index"))

## def download():
##     m = current_model()
##     response.headers["Content-Type"] = "text/plain"
##     s = "attachment; filename=%s.pickle" % m.name
##     response.headers["Content-Disposition"] = s
##     return cPickle.dumps(m)

def download():
    m = current_model()
    response.headers["Content-Type"] = "text/plain"
    s = "attachment; filename=%s.config.txt" % m.name
    response.headers["Content-Disposition"] = s
    return m.serialize()

def create_script():
    m = current_model()
    ## errors = m.errors()
    response.headers["Content-Type"] = "text/plain"
    s = "attachment; filename=%s.lagrange.py" % m.name
    response.headers["Content-Disposition"] = s
    return m.create_script()
    ## if not errors:
    ##     response.headers["Content-Type"] = "text/plain"
    ##     s = "attachment; filename=%s.lagrange.py" % m.name
    ##     response.headers["Content-Disposition"] = s
    ##     return m.create_script()
    ## else:
    ##     session.flash = "errors in configuration"
    ##     session.model_errors = errors
    ##     redirect(URL(r=request,f="index"))

def errors():
    return dict(errors = session.model_errors,
                linkback = A("Return", _href=URL(r=request,f="index")))

def export():
    m = current_model()
    response.headers["Content-Type"] = "text/plain"
    s = "attachment; filename=%s.params.py" % m.name
    response.headers["Content-Disposition"] = s
    return m.serialize()

def adjacencymatrix_default():
    m = current_model()
    m.adjacencymatrix.labels = m.datamatrix.labels
    m.adjacencymatrix.default()
    u = URL(r=request, f=request.vars.f)
    redirect(u)
    #return m.adjacencymatrix.render(request, session)
    
## def dispersalmatrix_default():
##     m = current_model()
##     m.dispersalmatrix.labels = m.datamatrix.labels
##     m.dispersalmatrix.default()
##     u = URL(r=request, f=request.vars.f)
##     redirect(u)

def dm_default():
    m = current_model()
    m.dm.labels = m.datamatrix.labels
    m.dm.default()
    redirect(URL(r=request, f=request.vars.f or "index"))

def datamatrix_default_from_tree():
    try:
        i = int(request.vars.t)
    except:
        session.flash = "no tree '%s'" % i
        redirect(URL(r=request,f=request.vars.f or "index"))
    try:
        nareas = int(request.vars.nareas)
        if nareas > 8:
            session.flash = "that is a lot of areas! (more than 8)"
    except:
        session.flash = "areas must be an integer"
        redirect(URL(r=request,f=request.vars.f or "index"))
        
    m = current_model()
    m.datamatrix = models.DataMatrix(m)
    t = newick.parse(m.treelist.trees[i].newick)
    m.datamatrix.default_from_tree(i, nareas)
    u = URL(r=request, f=request.vars.f or "index")
    session.ti = 0
    redirect(u)

def datamatrix_toggle_view():
    m = current_model()
    d = m.datamatrix
    if d.view == "html":
        d.view = "textarea"
    else:
        d.view = "html"
    return d.render(request, session)

def datamatrix_toggle_visible():
    m = current_model()
    d = m.datamatrix
    d.visible = not d.visible
    return d.render(request, session)

def treelist_toggle_visible():
    m = current_model().treelist
    m.visible = not m.visible
    return m.render(request, session)

def adjacencymatrix_toggle_visible():
    m = current_model()
    d = m.adjacencymatrix
    d.visible = not d.visible
    return d.render(request, session)

def adjacencymatrix_toggle_view():
    m = current_model()
    d = m.adjacencymatrix
    if d.view == "html":
        d.view = "textarea"
    else:
        d.view = "html"
    return d.render(request, session)

## def dispersalmatrix_toggle_visible():
##     m = current_model()
##     d = m.dispersalmatrix
##     d.visible = not d.visible
##     return d.render(request, session)

def dm_toggle_visible():
    m = current_model()
    m.dm.visible = not m.dm.visible
    return m.dm.render(request, session)

def rangelist_toggle_visible():
    m = current_model()
    d = m.rangelist
    d.visible = not d.visible
    return d.render(request, session)

## def dispersalmatrix_toggle_view():
##     m = current_model()
##     d = m.dispersalmatrix
##     if d.view == "html":
##         d.view = "textarea"
##     else:
##         d.view = "html"
##     return d.render(request, session)

def update_taxon():
    m = current_model().datamatrix
    s = "_".join(request.vars.s.strip().split())
    i = int(request.vars.i)
    assert s
    taxa = m.taxa
    t = taxa[i]
    taxa[i] = s
    for k, v in m.data.items():
        if k == t:
            m.data[s] = v
            del m.data[k]
            break
    return m.render_taxon(request, session, i)

def update_area_label():
    m = current_model().datamatrix
    s = "".join(request.vars.s.strip().split())
    i = int(request.vars.i)
    assert s
    labels = m.labels
    lab = labels[i]
    labels[i] = s
    return m.render_area_label(request, session, i)

def update_modelname():
    m = current_model()
    m.name = "_".join(request.vars.s.split())
    return m.render_name(request, session)

def update_modelname_form():
    m = current_model()
    m.name = "_".join(request.vars.s.split())
    redirect(URL(r=request, f=request.vars.f or "index"))

def update_max_range_size():
    m = current_model()
    m.max_range_size = int(request.vars.s)
    return m.render_max_range_size(request, session)

def delete_datamatrix():
    m = current_model()
    m.datamatrix = models.DataMatrix(m)
    m.rangelist.empty()
    u = URL(r=request, f=request.vars.f)
    redirect(u)

def delete_adjacencymatrix():
    m = current_model()
    m.adjacencymatrix = models.AdjacencyMatrix(m)
    u = URL(r=request, f=request.vars.f)
    redirect(u)

def datamatrix_click():
    i, j = map(int, request.vars.cell.split(":"))
    mod = current_model()
    dm = mod.datamatrix
    v = dm.toggle(request, session, i,j)
    return dm.render_cell(request, session, i, j)

def adjacencymatrix_click():
    i, j = map(int, request.vars.cell.split(":"))
    m = current_model().adjacencymatrix
    v = m.toggle(request, session, i,j)
    print "adj mat click:", request.vars.cell, i, j, v
    return m.render(request, session)

def edit():
    if request.args:
        n = int(request.args[0])
        m = session.models[n]
        if isinstance(m.labels, str):
            m.labels = list(m.labels)
        e = [SQLField("nareas", "integer", default=m.nareas,
                      label="Number of areas",
                      requires=IS_INT_IN_RANGE(2,10))]
        for i, x in enumerate(m.labels):
            f = SQLField("area_%s" % i, default=x, label=i,
                         requires=IS_NOT_EMPTY())
            e.append(f)
        form = form_factory(*e)
        if form.accepts(request.vars, session):
            for i in range(m.nareas):
                v = form.vars["area_%s" % i]
                #print i, v
                m.labels[i] = v
                e[i+1].default = v
            form = form_factory(*e)
    else:
        form=None

    return dict(form=form)

def pickle_upload_form():
    f = form_factory(
        SQLField("pickle", "upload", requires=IS_NOT_EMPTY(),
                 label="Upload model from pickle file"),
        _action="upload_pickle"
        )
    return f

def serial_upload_form():
    f = FORM(INPUT(_name="d", _type="file"), BR(),
             INPUT(_type="submit", _value="Submit"),
             _action="upload_serial",
             _title="Upload a previously saved configuration from {name}.config.txt file"
             )
    return f

def analysis_upload_form():
    f = FORM(INPUT(_name="d", _type="file"), BR(),
             INPUT(_type="submit", _value="Submit"),
             _action="upload_analysis",
             _title="Upload a previously saved analysis file (python script)"
             )
    return f

def model_selector():
    if not session.models:
        session.models = []
    session.mi = min(session.mi or 0, len(session.models))
    opts = []
    for i, m in enumerate(session.models):
        opts.append(OPTION(m.name, _value=i))
    if opts:
        f = FORM(
            SELECT(_name="mi", value=session.mi,
                   _onchange="this.form.submit()",
                   _id="model:selector:input",
                   *opts),
            _action=URL(r=request,f="select_mi")
            )
        ## return TABLE(
        ##     TR(TD("Current analysis:", _style="vertical-align:middle;"),
        ##        TD(f, _style="vertical-align:middle;")),
        ##     _style="padding:0px;"
        ##     )
        return f
    else:
        return ""

def select_mi():
    session.mi = max(0, int(request.vars.mi))
    session.mi = min(session.mi, len(session.models))
    redirect("index")

def delete_all_models():
    uid = response.session_id+"-model"
    cache.ram.clear(uid)
    session.ti = 0
    session.rci = 0
    session.models = []
    redirect("index")

def index():
    session.ti = int(request.vars.ti or session.ti or 0)
    ti = session.ti
    m = current_model()
    #session.models = []
    return dict(
        model = m,
        model_selector = model_selector(),
        analysis_upload_form = analysis_upload_form(),
        toolbar=toolbar()
        )

def create_range_list_ajax():
    m = current_model()
    m.calc_ranges()
    return m.rangelist.render(request, session)

def create_range_list():
    session.rci = 2
    m = current_model()
    m.calc_ranges()
    redirect(URL(r=request, f=request.vars.f or "index"))

def create_new_model():
    session.ti = int(request.vars.ti or session.ti or 0)

    ## if session.models is None:
    ##     session.models = []
    ## session.models.append(models.DECModel(session))

    uid = response.session_id+"-model"
    m = cache.ram(uid, lambda: models.DECModel(session), time_expire=30000)
    ## session.models = [m]
    session.mi = 0
    
    ## session.mi = len(session.models)-1

    u = URL(r=request, f=request.vars.f)
    redirect(u)

def current_model():
    uid = response.session_id+"-model"
    m = cache.ram(uid, lambda: models.DECModel(session), time_expire=30000)
    
    ## if session.models:
    ##     session.mi = max(0, session.mi)
    ##     session.mi = min(session.mi, len(session.models)-1)
    ##     m = session.models[session.mi]
    ##     session.model = m
    ## else:
    ##     m = None
    ##     session.mi = 0
    return m

def upload_datamatrix():
    f = request.vars.datamatrix
    m = current_model()
    infile = None
    fname = None
    if hasattr(f, "file"):
        infile = f.file
        fname = f.filename.split(".")[0]
    elif request.vars.datamatrix=="psychotria":
        infile = open("applications/%s/static/psychotria.matrix.txt"
                      % request.application)
        fname = "psychotria_demo"
    else:
        session.flash = "No file specified"
        pass
    if infile:
        try:
            labels, taxa, data = lagrange.input.parse_matrix3(infile)
        except:
            session.flash = "Unable to parse data"
            redirect(URL(r=request, f="index"))
        m.datamatrix.labels = labels
        m.datamatrix.taxa = taxa 
        m.datamatrix.data = data
        m.adjacencymatrix.default()
        nareas = len(m.datamatrix.labels)
        maxobs = max([ len(x) for x in m.datamatrix.data.values() ])
        v = m.max_range_size or max(maxobs, 2)
        m.max_range_size = v
        ## if nareas*v > 100:
        ##     session.flash = "Too many areas!"
        ##     redirect(URL(r=request, f="index"))
        m.calc_ranges()
        m.dm.default()
        if fname:
            m.name = fname

        if not m.treelist.trees:
            session.flash = "Species ranges uploaded. Next, upload a tree"
        
    redirect(URL(r=request, f="index"))

def base_rates_radio_clicked():
    print request.vars
    return current_model().render_base_rates(request, session)

def base_rates_form():
    m = current_model()
    br = request.vars.base_rates
    if br == "estimate":
        m.base_rates = "__estimate__"
    elif br == "fix":
        try:
            d = float(request.vars.dispersal)
            e = float(request.vars.extinction)
            m.base_rates = dict(dispersal=d, extinction=e)
        except:
            session.flash = "Dispersal and extinction must be numbers!"
    else:
        pass
    redirect(URL(r=request, f="index"))


def instructions():
    return dict()
