# one variable for each vertex
R.<x0, x1, x2, x3> = PolynomialRing(QQ)

# equations for colouring the vertices
vlst = [x0^3 - 1, x1^3 - 1, x2^3 - 1, x3^3 - 1]
# constraints given by the edges
elst = [x0^2 + x0*x1 + x1^2,
        x0^2 + x0*x2 + x2^2,
        x0^2 + x0*x3 + x3^2,
        x1^2 + x1*x3 + x3^2,
        x2^2 + x2*x3 + x3^2]

I = R.ideal(vlst + elst)

G = I.groebner_basis()

# output:
# [x3^3 - 1,
#  x2^2 + x2*x3 + x3^2,
#  x0 + x2 + x3,
#  x1 - x2]
# pretty cool!
# in particular we see that 1 does not appear to be in G (and hence in I)
# there is a shortcut for checking the latter:
1 in I

# output: False
# so 1 is not in I, hence the system is solvable, hence the graph is
# 3-colourable
# in fact, the Groebner basis gives us some hints toward finding a colouring,
# e.g. that vertices 1 and 2 must have the same colour
