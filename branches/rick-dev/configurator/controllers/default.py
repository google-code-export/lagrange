import random
from gluon.storage import Storage
import applications.lagrange.modules.lagrange as lagrange
import applications.lagrange.modules.bridge as bridge

reload(bridge)

def decmodel():
    name = request.vars.name
    if not session.counter:
        session.counter = 0
    ## if not session.models:
    ##     session.models = []
    ##     session.name2model = {}
    ## if name:
    ##     print name, session.name2model
    ##     assert name in session.name2model
    ##     model = session.name2model[name]
    ## else:

    if not session.model:
        model = Storage()
        model.name = "Untitled" + str(session.counter)
        ## session.models.append(model)
        ## session.name2model[model.name] = model
        session.counter += 1
        #session.model = model
    else:
        model = session.model

    f = FORM(TABLE(
        TR("Name:", INPUT(_type="text", _name="name",
                          _value=model.name or "",
                          requires=IS_NOT_EMPTY())),
        TR("Number of areas:", INPUT(_type="text", _name="nareas",
                                     _value=model.nareas,
                                     requires=IS_INT_IN_RANGE(2,10))),
        TR("",INPUT(_type="submit",_value="SUBMIT"))
        ))

    if f.accepts(request.vars, session):
        print "here!"

    return dict(f=f)

def testassign():
    session.foo = "bar"
    session.x = random.randint(0,100)
    session.model = bridge.X()
    redirect(URL(r=request,f="index"))

def index():
    #return dict()
    redirect(URL(r=request,c="configurator",f="index"))

def counter():
    if not session.counter: session.counter=0
    if not session.xd: session.xd = {"model": lagrange.DECModel(4)}
    session.counter+=1
    return dict(counter=session.counter, xd=id(session.xd["model"]))
