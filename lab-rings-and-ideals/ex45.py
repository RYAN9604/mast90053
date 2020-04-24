I = ZZ.ideal([52500, 10725])

I

# output: Principal ideal (75) of Integer Ring

I.is_principal()

# output: True

R.<x> = QQ[]

I = R.ideal([3*x^4 - x^3 + x^2 - x - 2, 3*x^4 + 2*x^3 + 6*x + 4])

I

# output: Principal ideal (x + 2/3) of Univariate Polynomial Ring in x over Rational Field

I.is_principal()

# output: True

K.<alpha> = QuadraticField(-5)

R = K.ring_of_integers()

alpha^2

# output: -5

I = R.ideal([2, 1 + alpha])

I

# output: Fractional ideal (2, alpha + 1)

I.is_principal()

# output: False

R.<x, y> = QQ[]

I = R.ideal([x, y])

I

# output: Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field

I.is_principal()

# output: sadly, NotImplementedError
