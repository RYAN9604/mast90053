def is_hypergeometric_term(f, v):
    ratio = f.subs({v: v+1})/f.subs({v: v})
    ratio = ratio.simplify_full()
    num, denom = ratio.numerator_denominator()
    return num.is_polynomial(v) and denom.is_polynomial(v)


def hypergeometric_term_ratio(f, v):
    if not is_hypergeometric_term(f, v):
        raise ValueError("%s is not a hypergeometric term in %s" % (f, v))
    ratio = f.subs({v: v+1})/f.subs({v: v})
    ratio = ratio.simplify_full()
    num, denom = ratio.numerator_denominator()
    # R = PolynomialRing(SR, str(v))
    # return R(num), R(denom)
    return num, denom


def hypergeometric_term_recurrence(f, n, k, I=1, J=1):
    """
    Return the coefficients of a k-free recurrence relation for the
    hypergeometric summand f(n, k).

    Implements Celine Fasenmyer's algorithm.

    If a relation does not exist with the given box parameters I and J,
    raises a RuntimeError.
    """
    num   = [[0] * (J + 1) for _ in range(I + 1)]
    denom = [[1] * (J + 1) for _ in range(I + 1)]
    for i in range(I + 1):
        for j in range(J + 1):
            g = f.subs({n: n + j, k: k + i}) / f.subs({n: n, k: k})
            g = g.simplify_full()
            num[i][j], denom[i][j] = g.numerator_denominator()
    LCM = lcm(flatten(denom))
    poly = [[(num[i][j] / denom[i][j] * LCM).simplify_full() for j in range(J + 1)] for i in range(I + 1)]
    coeff = [[SR.var("a%s_%s" % (i, j)) for j in range(J + 1)] for i in range(I + 1)]
    lhs = sum(flatten([[coeff[i][j] * poly[i][j] for j in range(J + 1)] for i in range(I + 1)]))
    try:
        v = solve_for_coefficients_homog(lhs, k, flatten(coeff))
    except IndexError:
        raise RuntimeError("bounds I = %s and J = %s are not large enough" % (I, J))
    return [SR(a).factor() for a in v]


def dot_product(a, b):
    if len(a) != len(b):
        raise ValueError("%s and %s must have the same length" %(a, b))
    return sum([a[j] * b[j] for j in range(len(a))])


def solve_for_coefficients_homog(f, k, coeff):
    """
    Given a symbolic expression f, a distinguished variable k, and a
    set of variables coeff such that f is polynomial in k and
    (homogeneous) linear in coeff, return a nonzero solution for coeff
    in the equation f = 0.

    If the only solution is zero this raises an IndexError.
    """
    rels = [r for [r, _] in f.coefficients(k)]
    rows = [[r.coefficient(c).simplify_full() for c in coeff] for r in rels]
    lst = list(f.variables())
    lst.remove(k)
    for c in coeff:
        lst.remove(c)
    S = PolynomialRing(QQ, lst)
    m = matrix(S.fraction_field(), rows)
    K = m.right_kernel()
    return K.basis()[0]


def _test_all():
    _test_hypergeometric_term_recurrence()


def _test_hypergeometric_term_recurrence():
    k, n = SR.var("k, n")
    # first test
    f = k * binomial(n, k)
    res = hypergeometric_term_recurrence(f, n, k, I=1, J=1)
    correct = [1, 0, 1, -n / (n+1)]
    if res != correct:
        raise RuntimeError("hypergeometric_term_recurrence(%s, %s, %s, %s, %s) = %s but should be %s" % (f, n, k, 1, 1, res, correct))
    # second test, part 1
    x = SR.var("x")
    f = binomial(n, k) * binomial(-n-1, k) * ((1-x) / 2)^k
    passed = False
    I = 1
    J = 1
    try:
        res = hypergeometric_term_recurrence(f, n, k, I, J)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("hypergeometric_term_recurrence(%s, %s, %s, %s, %s) = %s but should have raised error" % (f, n, k, I, J, res))
    # second test, part 2
    f = binomial(n, k) * binomial(-n-1, k) * ((1-x) / 2)^k
    I = 1
    J = 2
    res = hypergeometric_term_recurrence(f, n, k, I, J)
    correct = [0, 1, 0, -(n+1) / ((x-1) * (2*n+3)), 1 / (x-1), -(n+2) / ((x-1)*(2*n+3))]
    if res != correct:
        raise RuntimeError("hypergeometric_term_recurrence(%s, %s, %s, %s, %s) = %s but should be %s" % (f, n, k, I, J, res, correct))
    # third test
    # commented out because it takes a long time: 6.5 minutes on my laptop
    # x, y = var("x, y")
    # f = binomial(n, k) * factorial(k) * 3^k / factorial(3*k) * x^(n-k) * y^k
    # I = 1
    # J = 3
    # res = hypergeometric_term_recurrence(f, n, k, I, J)
    # correct = [0, 0, 1, 0, 9*(n + 2)*(n + 1)*x^3/y, -9*(3*n + 5)*(n + 2)*x^2/y, (27*n^2 + 117*n + 128)*x/y, -(3*n + 8)*(3*n + 7)/y]
    # if res != correct:
    #     raise RuntimeError("hypergeometric_term_recurrence(%s, %s, %s, %s, %s) = %s but should be %s" % (f, n, k, I, J, res, correct))
