import applications.lagrange.modules.layout as layout
import applications.lagrange.modules.phylo as phylo
import applications.lagrange.modules.newick as newick
import applications.lagrange.modules.bridge as bridge

def current_model():
    if session.models:
        session.mi = max(0, session.mi)
        session.mi = min(session.mi, len(session.models)-1)
        m = session.models[session.mi]
        session.model = m
    else:
        m = None
        session.mi = 0
    return m

def index():
    i = request.vars.i or 0
    #s = "((((((((((T._windsori_oli117:0.10632847,T._windsori_oli116:0.10632847):0.04722967,(T._montislewisi_oli86:0.11755736,T._montislewisi_oli81:0.11755736):0.03600078):0.07131713,(T._sp._nov._oli122:0.02829776,T._sp._nov._oli119:0.02829776):0.19657751):0.01112654,T._carbinensis_oli76:0.23600181):0.0298459,(((T._nashi_OLI83:0.02125559,T._nashi_OLI29:0.02125559):0.02378355,T._nashi_oli23:0.04503914):0.10589362,(T._nashi_oli67:0.1104632,(T._liber_Paris:0.03137278,T._liber_oli88:0.03137278):0.07909041):0.04046956):0.11491495):0.015917,T._sp._nov._oli85:0.28176471):0.00699509,((((T._terrareginae_oli79:0.12986896,(T._canaliculatus_oli82:0.05718215,T._canaliculatus_oli10:0.05718215):0.07268681):0.02825186,(T._canaliculatus_oli92:0.08427881,(T._canaliculatus_OLI37:0.01226744,T._canaliculatus_OLI31:0.01226744):0.07201137):0.07384202):0.03424563,T._kuranda_oli55:0.19236646):0.06071168,T._miseriae_oli77:0.25307813):0.03568166):0.00806363,T._mcilwraithi_oli75:0.29682342):0.01063206,(T._grandis_oli126:0.0,T._grandis_oli125:0.0):0.30745548):0.03417052,(((T._phalacrus_oli90:0.09751735,T._phalacrus_oli44:0.09751735):0.16894837,(T._millamilla_oli48:0.0,T._millamilla_oli39:0.0):0.26646572):0.0285159,((T._erici_oli91:0.11421808,T._erici_oli89:0.11421808):0.06572013,(T._erici_oli57:0.1279707,(T._erici_OLI35:0.03200512,T._erici_43:0.03200512):0.09596559):0.05196751):0.11504341):0.04664438):0.01293956;"
    m = current_model()
    s = m.treelist.trees[i].newick
    root = newick.parse(s)
    return dict(treediv=nodes2css(root))

def newick2css():
    i = request.vars.i or 0
    m = current_model()
    s = m.treelist.trees[i].newick
    root = newick.parse(s)
    return nodes2css(root)

def nodes2css(root, interactive=False, wscale=1.0):
    width, height = style_nodes(root, wscale=wscale)
    divs = []
    for node in root.iternodes(phylo.PREORDER):
        if interactive is True:
            onclick = "ajax3(['%s','%s'], ['t=%s&n=%s','t=%s&n=%s'],"\
                      " ['treecss', 'clipboard'], 0);" % \
                      (URL(r=request,f="render"), URL(r=request,f="clipboard"),
                       node.tree_id, node.id, node.tree_id, node.id)
        else:
            onclick = ""
        divs.append(DIV(
            _style="position:absolute; "\
            "width:%(width)sem; height:%(height)sex; "\
            "background-color:%(background-color)s; "\
            "top:%(top)sex; left:%(left)sem" % node.hbranch_style,
            _onclick=onclick, _id="hbranch%s" % node.id,
            _title="node_id = %s" % node.id
            ))
        divs.append(DIV(
            _style="position:absolute; height:%(height)sex; "\
            "width:%(width)sem; background-color:%(background-color)s; "\
            "top:%(top)sex; left:%(left)sem" % node.vbranch_style,
            _onclick=onclick, _id="vbranch%s" % node.id,
            _title="node_id = %s" % node.id
            ))

        if len(node.children) == 1 and node.parent:
            style = node.hbranch_style.copy()
            style["left"] = style["left"]+style["width"]-0.5
            style["width"] = 0.75
            style["top"] -= style["height"]*0.5
            style["height"] *= 1.75
            #style["background-color"] = "blue"
            d = DIV(
                _style="position:absolute; "\
                "width:%(width)sem; height:%(height)sex; "\
                "background-color:%(background-color)s; "\
                "top:%(top)sex; left:%(left)sem" % style,
                _onclick=onclick, _id="singlechild%s" % node.id
                )
            divs.append(d)

        if node.istip:
            if session.snode_clickaction == "edit_label" and \
               node.id == session.selected_snode_id:
                divs.append(
                    DIV(FORM(INPUT(_type="hidden", _name="edit_node_id",
                                   _value=str(node.id)),
                             INPUT(_type="text", _name="edit_node_label",
                                   _size=8, value=node.label or "")),
                        _style="position:absolute; top:%(top)sex; "\
                        "left:%(left)sem" % node.label_style,
                        _id="label%s" % node.id)
                    )
            else:
                style=""
                if node.id in {}:#session.snode_collapsed:
                    style="background-color:yellow"
                if node.istip:
                    span = SPAN(node.label or "[collapsed]", _style=style,
                                _onclick=onclick)
                else:
                    span = SPAN(node.label or "[collapsed]", _style=style,
                                _onclick=onclick,
                                _title="%s leaves" % node.ntips)
                divs.append(
                    DIV(span, _style="position:absolute; top:%(top)sex; "\
                        "left:%(left)sem" % node.label_style,
                        _id="label%s" % node.id)
                    )
        else:
            if session.snode_clickaction == "edit_label" and \
               node.id == session.selected_snode_id:
                style = node.label_style.copy()
                if not node.label:
                    style["left"] = style["left"] - 6
                divs.append(
                    DIV(FORM(INPUT(_type="hidden", _name="edit_node_id",
                                   _value=str(node.id)),
                             INPUT(_type="text", _name="edit_node_label",
                                   _size=8, value=node.label or "")),
                        _style="position:absolute; width:%(width)sem; "\
                        "text-align:%(text-align)s; top:%(top)sex; "\
                        "left:%(left)sem" % style,
                        _id="label%s" % node.id)
                    )
            else:
                divs.append(
                    DIV(SPAN(node.label or "", _style="background-color:yellow",
                             _onclick=onclick),
                        _style="position:absolute; width:%(width)sem; "\
                        "text-align:%(text-align)s; top:%(top)sex; "\
                        "left:%(left)sem" % node.label_style,
                        _id="label%s" % node.id)
                    )
##             divs.append(
##                 DIV(SPAN(node.label or "", _style="background-color:yellow",
##                          _onclick="ajax2('%s','n=%s','label%s')" % \
##                          (URL(r=request,f="editlabel"), node.id, node.id)),
##                     _style="position:absolute; width:%(width)sem; "\
##                     "text-align:%(text-align)s; top:%(top)sex; "\
##                     "left:%(left)sem" % node.label_style,
##                     _id="label%s" % node.id)
##                 )
    d = DIV(_style="width:%sem; height:%sex;"\
            % (width, height),
            *divs)
    return d


def style_nodes(root, collapsed={}, selected_node_id=None, wscale=1.0,
                scaled=True):
    bgcolor = "black"
    selcolor = "red"
    leaves = root.leaves(collapsed=collapsed)
    l2d = root.leaf_distances(measure=phylo.INTERNODES,
                              collapsed=collapsed)[root]
    height = len(leaves)*3
    unit = 3.0
    width = max([ l2d[lf.id] for lf in leaves ]) * unit * wscale
    width = min(width, 65)
    rpad = max([ len(lf.label or "") for lf in leaves ]) * 0.7
    lpad = max(1, len(root.label or []) * 0.7)
    width += rpad+2 + lpad
    branchwidth = 0.75
    n2c = layout.calc_node_positions(
        root, width, height,
        lpad=lpad+1, tpad=1, bpad=2.5, rpad=rpad+1,
        collapsed=collapsed,
        scaled=scaled
        )
    n2c[root].px = 1
    for node in root.iternodes(collapsed=collapsed):
        coords = n2c[node]
        w = coords.x-coords.px
        node.hbranch_style["width"] = w
        node.hbranch_style["height"] = branchwidth
        node.hbranch_style["top"] = coords.y+0.5
        node.hbranch_style["left"] = coords.px
        node.hbranch_style["background-color"] = bgcolor
        if coords.py is None:
            coords.py = coords.y
        if coords.py < coords.y:
            h = coords.y-coords.py
            y = coords.py
        else:
            h = coords.py-coords.y
            y = coords.y
        node.vbranch_style["width"] = 0.5 * branchwidth
        node.vbranch_style["height"] = h
        node.vbranch_style["top"] = y+0.5
        node.vbranch_style["left"] = coords.px
        node.vbranch_style["background-color"] = bgcolor
        if node.istip or node.id in collapsed:
            node.label_style["top"] = coords.y-0.5
            node.label_style["left"] = coords.x+0.25
            node.label_style["width"] = len(node.label or "")
            node.label_style["text-align"] = "left"
        else:
            node.label_style["text-align"] = "right"
            node.label_style["width"] = len(node.label or "")
            node.label_style["top"] = coords.y-0.5
            node.label_style["left"] = coords.x-len(node.label or "")

        node.ref_style["top"] = coords.y-0.5
        node.ref_style["left"] = coords.x+0.25
        
    if selected_node_id:
        for n in root.iternodes(collapsed=collapsed):
            if n.id == selected_node_id:
                for m in n.iternodes():
                    m.vbranch_style["background-color"] = selcolor
                    m.hbranch_style["background-color"] = selcolor
                break
                
    return width, height
