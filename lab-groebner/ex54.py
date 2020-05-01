# surely the complete graph on 4 vertices is not 3-colourable
# let's check this
R.<x0, x1, x2, x3> = PolynomialRing(QQ)

vlst = [x0^3 - 1, x1^3 - 1, x2^3 - 1, x3^3 - 1]
elst = [x0^2 + x0*x1 + x1^2,
        x0^2 + x0*x2 + x2^2,
        x0^2 + x0*x3 + x3^2,
        x1^2 + x1*x2 + x2^2,
        x1^2 + x1*x3 + x3^2,
        x2^2 + x2*x3 + x3^2]

I = R.ideal(vlst + elst)

G = I.groebner_basis()

# output: [1]
# aha!
# just to double-check:
1 in I

# output: True

# Sage also knows some graph theory:
K4 = graphs.CompleteGraph(4)
K4.chromatic_number()

# output: 4
# which means that the smallest number of colours necessary for colouring
# the vertices of this is 4
# (in particular, not 3-colourable)
