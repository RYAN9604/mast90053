def intrel(xlst, A):
    n = len(xlst)
    mlst = identity_matrix(n).rows()
    ylst = [(A*x).floor() for x in xlst]
    mlst.append(ylst)
    M = matrix(ZZ, mlst)
    MR = M.transpose().LLL().transpose()
    return MR.column(0)[:-1]

R = RealField(10000)

alpha = R(pi)

alpha

# output: lots of decimals

# let's try with degree 10
xlst = [alpha^j for j in range(10, -1, -1)]

# try with increasing multiplier, seem to just be getting random
# numbers; indication that pi is not algebraic of degree 10
# same phenomenon with other degrees
intrel(xlst, 10^10)

# output (-1, 3, -1, 0, 16, -3, -8, -3, 10, 5, 3)

intrel(xlst, 10^20)

# output (23, -43, -80, -30, -9, -42, -35, 63, 46, 5, 69)

intrel(xlst, 10^40)

# output (754, -1895, 1016, -7849, 1388, -4681, 811, -1425, -1121, 6564, -2753)

intrel(xlst, 10^100)

# output (-231778308, 539954899, 98037843, 1219530518, 786398747, 814478433, 457822548, -1623804602, -837409217, 1407541980, 755233122)

