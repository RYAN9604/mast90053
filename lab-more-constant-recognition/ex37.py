def intrel(xlst, A):
    n = len(xlst)
    mlst = identity_matrix(n).rows()
    ylst = [(A*x).floor() for x in xlst]
    mlst.append(ylst)
    M = matrix(ZZ, mlst)
    MR = M.transpose().LLL().transpose()
    return MR.column(0)[:-1]

eta(1+I)

# output 0.7420487758365647 + 0.1988313702299107*I

beta = eta(sqrt(-5)/2)^2/(2*eta(2*sqrt(-5))^2)

beta

# output 2.8900536382639648 - 3.130036836413252e-16*I

# it looks like a complex number, but the imaginary part is really small
# (actually just numerical noise)
# so we work with the real part
alpha = beta.real_part()

alpha

# output 2.8900536382639648

xlst = [alpha^j for j in range(4, -1, -1)]

intrel(xlst, 10^16)

# output (1, -2, -2, -2, 1)

# guess: alpha is a root of the polynomial x^4 - 2x^3 - 2x^2 - 2x + 1
# this was proved by Abel long before computers
