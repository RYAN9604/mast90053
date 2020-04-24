# declare that x is a polynomial variable
R.<x> = QQ[]

f1 = 3*x^4 - x^3 + x^2 - x - 2

f2 = 3*x^4 + 2*x^3 + 6*x + 4

# long division of f1 by f2
(q1, r1) = f1.quo_rem(f2)

r1

# output: -3*x^3 + x^2 - 7*x - 6

(q2, r2) = f2.quo_rem(r1)

r2

# output: -6*x^2 - 7*x - 2

(q3, r3) = r1.quo_rem(r2)

r3

# output: -45/4*x - 15/2

(q4, r4) = r1.quo_rem(r3)

r4

# output: 0

# so the previous remainder r4 is the gcd
# this is only defined up to multiplication
# by a nonzero scalar

gcd(f1, f2)

# output: x + 2/3

# getting the coefficients u1 and u2 is known as the
# extended gcd
xgcd(f1, f2)

# output: (x + 2/3, 2/45*x^2 - 1/45*x - 7/45, -2/45*x^2 + 1/15*x + 4/45)

f3 = x^4 - x^3 + 2*x - 2

gcd([f1, f2, f3])

# output: 1
