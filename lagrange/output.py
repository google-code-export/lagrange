import scipy

def summarize_splits(splits, weighted=True):
    rows = []
    if weighted:
        v = [ (x.likelihood * x.weight, x) for x in splits ]
    else:
        v = [ (x.likelihood, x) for x in splits ]
    ptot = sum([ x[0] for x in v ])
    v.sort(); v.reverse()
    opt = scipy.log(v[0][0])

    rows.append(["split", "lnL", "Rel.Prob"])
    sumprob = 0.0
    for L, split in v:
        lnL = scipy.log(L)
        relprob = (L/ptot)
        #if sumprob < 0.95:
        if (opt - lnL) < 2:
            rows.append([str(split), "%.4g" % lnL, "%.4g" % relprob])
        sumprob += relprob
    widths = []
    for i in range(3):
        w = max([ len(x[i]) for x in rows ])
        for x in rows:
            x[i] = x[i].ljust(w)
    return [ "  ".join(x) for x in rows ]
