from gluon.sqlhtml import *
from gluon.sql import *
from gluon.storage import Storage

def d2s(d):
    return "".join([ "%s:%s;" % (k, v) for k, v in d.items() ])

class ViewPane:
    def __init__(self, obj):
        self.color = "lightgray"
        self.outer_style = {
            "-moz-border-radius": "1em",
            "border-style": "solid",
            "border-color": self.color,
            "display": "inline-block",
            "padding": "0em",
            }
        self.header_style = {
            "-moz-border-radius": "0.75em 0.75em 0 0",
            "background-color": "lightgray",
            "border-style": "solid",
            "border-color": self.color,
            }
        self.obj = obj

    def render_header(self):
        return DIV(_style=d2s(self.header_style))
            
    def render(self):
        return DIV(_style=d2s(self.outer_style))
