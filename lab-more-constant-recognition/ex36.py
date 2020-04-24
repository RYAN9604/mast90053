def intrel(xlst, A):
    n = len(xlst)
    mlst = identity_matrix(n).rows()
    ylst = [(A*x).floor() for x in xlst]
    mlst.append(ylst)
    M = matrix(ZZ, mlst)
    MR = M.transpose().LLL().transpose()
    return MR.column(0)[:-1]

alpha = RR(1+sqrt(1+sqrt(1+sqrt(2))))

xlst = [alpha^j for j in range(8, -1, -1)]

intrel(xlst, 10^3)

# output (0, -1, 2, 1, 0, 3, 1, 2, 1)

intrel(xlst, 10^4)

# output (-1, 1, 3, 1, 4, 3, 0, 0, 1)

intrel(xlst, 10^8)

# output (-1, 3, 0, -1, -3, -7, 11, -6, -2)

intrel(xlst, 10^10)

# output (4, -12, 8, -4, -15, 5, -15, -2, 0)

intrel(xlst, 10^25)

# output (-36, 58, 35, 115, 91, 8, -73, 77, 18)

# let's increase the precision a lot
R = RealField(5000)

alpha = R(1+sqrt(1+sqrt(1+sqrt(2))))

xlst = [alpha^j for j in range(8, -1, -1)]

intrel(xlst, 10^50)

# output (-1, 8, -24, 32, -14, -8, 8, 0, 1)


# in fact, Sage can solve this problem exactly (without approximations)
# we define our number symbolically
beta = 1+sqrt(1+sqrt(1+sqrt(2)))

beta

# output sqrt(sqrt(sqrt(2) + 1) + 1) + 1

beta.minpoly()

# output x^8 - 8*x^7 + 24*x^6 - 32*x^5 + 14*x^4 + 8*x^3 - 8*x^2 - 1

