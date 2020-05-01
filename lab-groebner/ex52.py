# same overall approach as in 5.1
R.<x, y> = QQ[]

I = R.ideal([6*x^2*y - x + 4*y^3 - 1, 2*x*y + y^3])

# get a Groebner basis
G = I.groebner_basis()

G

# really scary output:
# [x^3 - 4/3*x^2 + 1/12*x*y - 7/6*y^2 - 7/3*x + 1/12*y,
#  x^2*y - 4/3*x*y - 1/6*x - 1/6,
#  x*y^2 + 2*x^2 + y^2 + 2*x,
#  y^3 + 2*x*y]
#
# can't really see what to conclude from this
# but remember that an ideal has many Groebner bases, depending
# (in particular) on the monomial order we choose
# looking at PolynomialRing? suggests that the default order is
# degrevlex (Sage speak for grevlex), and that one can choose other
# orders with the 'order' argument
R.<x, y> = PolynomialRing(QQ, order='lex')

I = R.ideal([6*x^2*y - x + 4*y^3 - 1, 2*x*y + y^3])

G = I.groebner_basis()

G

# output:
# [x - 3/2*y^5 - 4*y^3 + 1,
#  y^6 + 8/3*y^4 + 1/3*y^3 - 2/3*y]
# aha, now the second polynomial is in y only

fy = G[1].univariate_polynomial()

fy.roots()

# output: [(0, 1)]
# 0 is the only rational root
# more of them, symbolically?
fy.roots(SR)

# output: [(0, 1)]
# bummer.
# but there must be more real roots (because of the Intermediate Value Theorem)
fy.roots(RR)

# output: [(0.000000000000000, 1), (0.571237275081126, 1)]
# so these are two real roots
# we can ask for the complex roots to convince ourselves we haven't
# missed anything
fy.roots(CC)

# output:
# [(0.000000000000000, 1),
#  (0.571237275081126, 1),
#  (-0.390690658252738 - 0.524320823326721*I, 1),
#  (-0.390690658252738 + 0.524320823326721*I, 1),
#  (0.105072020712175 - 1.64881462666250*I, 1),
#  (0.105072020712175 + 1.64881462666250*I, 1)]
