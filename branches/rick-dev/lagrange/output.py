import scipy

def summarize_splits(splits, weighted=True):
    lines = []
    if weighted:
        v = [ (x.likelihood * x.weight, x) for x in splits ]
    else:
        v = [ (x.likelihood, x) for x in splits ]
    ptot = sum([ x[0] for x in v ])
    v.sort(); v.reverse()
    opt = scipy.log(v[0][0])

    lines.append("".join(["split".ljust(8), "lnL".ljust(8), "Rel.Prob"]))
    for L, split in v:
        lnL = scipy.log(L)
        if (opt - lnL) < 2:
            lines.append("".join(
                [str(split).ljust(8), ("%.4g" % lnL).ljust(8),
                 "%.4g" % (L/ptot)]
                ))
    return lines
