class Ancsplit:
    """
    convenience class for encapsulating an ancestor range splitting
    into descendant ranges
    """
    def __init__(self, ancdist, descdists, weight=None, likelihood=None):
        self.ancdist = ancdist
        self.descdists = descdists
        self.weight = weight
        self.likelihood = likelihood

    def __repr__(self):
        d1, d2 = [ "".join(map(str, x)) for x in self.descdists ]
        lh = self.likelihood
        if lh: lh = "%g" % lh
        w = self.weight
        if w: w = "%g" % w
        return "Ancsplit(%s, %s, w=%s, lh=%s)" % (d1, d2, w, lh)
