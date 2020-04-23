alpha = RR(sqrt(2))

x1 = alpha^2
x2 = alpha
x3 = 1

def intrel(xlst, A):
    n = len(xlst)
    mlst = identity_matrix(n).rows()
    ylst = [(A*x).floor() for x in xlst]
    mlst.append(ylst)
    M = matrix(ZZ, mlst)
    MR = M.transpose().LLL().transpose()
    return MR.column(0)[:-1]

intrel([x1, x2, x3], 10^15)

# outputs (-1, 0, 2)
