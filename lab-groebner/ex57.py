# define alpha
K.<alpha> = NumberField(x^5 - x - 2)

K

# output:
# Number Field in alpha with defining polynomial x^5 - x - 2

alpha^5 - alpha - 2

# output: 0

beta = (1 - alpha - 2*alpha^3)/alpha

# let's get Theorem 5.6 working
R.<x, y> = QQ[]

p = x^5 - x - 2

f = 1 - x - 2*x^3

g = x

I = R.ideal(p, g*y - f)

I

I.groebner_basis()

# output:
# [x*y^2 + 61/30*x*y + 2/15*y^2 - 73/30*x - 26/15*y - 313/30,
#  y^3 + 719/120*x*y + 52/15*y^2 + 2863/120*x + 11/60*y - 497/120,
#  x^2 - 1/60*x*y - 1/15*y^2 - 17/60*x + 11/30*y + 43/60]

# oh no, what happened?  we were supposed to have a polynomial of y only!
# ah yes, lexicographic order
# here we go again
R.<x, y> = PolynomialRing(QQ, order='lex')
p = x^5 - x - 2
f = 1 - x - 2*x^3
g = x
I = R.ideal(p, g*y - f)

G = I.groebner_basis()

G

# output:
# [x - 1438/45887*y^4 - 2183/45887*y^3 + 10599/45887*y^2 - 8465/45887*y - 101499/45887,
#  y^5 + 11/2*y^4 + 4*y^3 - 5*y^2 + 95*y + 259]
# there it is
mp = G[1].univariate_polynomial()

mp(beta)

# output: 0

mp.is_irreducible()

# output: True

# by the way:
beta.minpoly()

# output:
# x^5 + 11/2*x^4 + 4*x^3 - 5*x^2 + 95*x + 259
