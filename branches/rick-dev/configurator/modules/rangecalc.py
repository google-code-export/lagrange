import lagrange
scipy = lagrange.scipy

def test(model):
    labels = model.datamatrix.labels
    data = model.datamatrix.data
    nareas = len(labels)
    edges = []
    dm = scipy.ones((nareas, nareas))
    adj = model.adjacencymatrix
    for i in range(nareas):
        for j in range(i+1, nareas+1):
            edges.append((i,j))
    ag = lagrange.graph.AreaGraph(edges=edges)
    ranges = data.values()
