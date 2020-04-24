n1 = 52500

n2 = 10725

# long division of n1 by n2
# q1 is the quotient, r1 the remainder
(q1, r1) = n1.quo_rem(n2)

r1

# output: 9600

# continue with n2 and r1
(q2, r2) = n2.quo_rem(r1)

r2

# output: 1125

(q3, r3) = r1.quo_rem(r2)

r3

# output: 600

(q4, r4) = r2.quo_rem(r3)

r4

# output: 525

(q5, r5) = r3.quo_rem(r4)

r5

# output: 75

(q6, r6) = r4.quo_rem(r5)

r6

# output: 0

# so the previous remainder, r5 = 75, is the gcd

# this is the Euclidean algorithm for computing the gcd
# of course it's implemented in Sage so we can just do
gcd(n1, n2)

# output: 75

# getting the coefficients u1 and u2 is known as the
# extended gcd; they are obtained by rolling back the
# steps of the Euclidean algorithm.
# Or you can just ask Sage:
xgcd(n1, n2)

# output: (75, 19, -93)

75 == 19*n1 - 93*n2

# output: True
