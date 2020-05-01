def is_colourable(G, n):
    """Return True if the graph G is n-colourable, False otherwise."""
    # easy naming of the variables
    lst = ['x%s' %v for v in G.vertices()]
    R = PolynomialRing(QQ, lst)
    # ah, now we want easy access to the variables as elements of R
    # (rather than just strings)
    x = dict()
    for v in G.vertices():
        x[v] = R('x%s' %v)
    # now we're ready for the vertex colours
    vlst = [x[v]^n - 1 for v in G.vertices()]
    # and the constraints from the edges
    # slightly hacky
    elst = [(x[v]^n - x[w]^n).quo_rem(x[v] - x[w])[0] for (v, w, _) in G.edges()]
    I = R.ideal(vlst + elst)
    # just get Sage to decide whether I contains 1
    # (uses Groebner bases behind the scenes)
    return (1 not in I)
