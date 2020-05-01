# declare the polynomial ring
R.<x, y> = QQ[]

# and define the ideal
I = R.ideal([x^2 + y^2 - 1, x^2 - y])

# get a Groebner basis
G = I.groebner_basis()

G

# output: [x^2 - y, y^2 + y - 1]
# note that the second polynomial is in the single variable y

f = G[1]

# but Sage still thinks it's in two variables
parent(f)

# so we force Sage to think of it in one variable
fy = f.univariate_polynomial()

parent(fy)

# now we can ask for its roots
fy.roots()

# output: []
# no roots? yes, no roots in QQ
# let's ask for roots in RR
fy.roots(RR)

# output: [(-1.61803398874989, 1), (0.618033988749895, 1)]
# you can now take these values of y and plug them into
# the first polynomial of the Groebner basis to get the corresponding
# values of x

# or you can be more demanding and ask Sage for the symbolic (exact) roots
fy.roots(SR)

# output: [(-1/2*sqrt(5) - 1/2, 1), (1/2*sqrt(5) - 1/2, 1)]

