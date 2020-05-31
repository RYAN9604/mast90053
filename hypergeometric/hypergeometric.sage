"""
Collection of tools for hypergeometric summation in Sage.

Includes:
* fasenmyer_kfree(): Celine Fasenmyer's algorithm for computing
  recurrence relations
* gosper_sum(): Gosper's algorithm for hypergeometric indefinite
  summation
* wz_certificate(): Wilf-Zeilberger's algorithm for definite
  summation
* zeilberger(): Zeilberger's algorithm for computing recurrence
  relations
"""


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


def fasenmyer_kfree(f, n, k, I=1, J=1):
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


def _test_fasenmyer_kfree():
    k, n = SR.var("k, n")
    # first test
    f = k * binomial(n, k)
    res = fasenmyer_kfree(f, n, k, I=1, J=1)
    correct = [1, 0, 1, -n / (n+1)]
    if res != correct:
        raise RuntimeError("fasenmyer_kfree(%s, %s, %s, %s, %s) = %s but should be %s" % (f, n, k, 1, 1, res, correct))
    # second test, part 1
    x = SR.var("x")
    f = binomial(n, k) * binomial(-n-1, k) * ((1-x) / 2)^k
    passed = False
    I = 1
    J = 1
    try:
        res = fasenmyer_kfree(f, n, k, I, J)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("fasenmyer_kfree(%s, %s, %s, %s, %s) = %s but should have raised error" % (f, n, k, I, J, res))
    # second test, part 2
    f = binomial(n, k) * binomial(-n-1, k) * ((1-x) / 2)^k
    I = 1
    J = 2
    res = fasenmyer_kfree(f, n, k, I, J)
    correct = [0, 1, 0, -(n+1) / ((x-1) * (2*n+3)), 1 / (x-1), -(n+2) / ((x-1)*(2*n+3))]
    if res != correct:
        raise RuntimeError("fasenmyer_kfree(%s, %s, %s, %s, %s) = %s but should be %s" % (f, n, k, I, J, res, correct))
    # third test
    # commented out because it takes a long time: 6.5 minutes on my laptop
    # x, y = var("x, y")
    # f = binomial(n, k) * factorial(k) * 3^k / factorial(3*k) * x^(n-k) * y^k
    # I = 1
    # J = 3
    # res = fasenmyer_kfree(f, n, k, I, J)
    # correct = [0, 0, 1, 0, 9*(n + 2)*(n + 1)*x^3/y, -9*(3*n + 5)*(n + 2)*x^2/y, (27*n^2 + 117*n + 128)*x/y, -(3*n + 8)*(3*n + 7)/y]
    # if res != correct:
    #     raise RuntimeError("fasenmyer_kfree(%s, %s, %s, %s, %s) = %s but should be %s" % (f, n, k, I, J, res, correct))


def dot_product(a, b):
    if len(a) != len(b):
        raise ValueError("%s and %s must have the same length" %(a, b))
    return sum([a[j] * b[j] for j in range(len(a))])


def _test_dot_product():
    # first test
    a = [1, 2, 3]
    b = [3, 2, 0]
    res = dot_product(a, b)
    correct = 7
    if res != correct:
        raise RuntimeError("dot_product(%s, %s) = %s but should be %s" % (a, b, res, correct))
    # second test
    a = [1, 2, 3]
    b = [3, 2]
    passed = False
    try:
        res = dot_product(a, b)
    except ValueError:
        passed = True
    if not passed:
        raise RuntimeError("dot_product(%s, %s) = %s but should have raised error" % (a, b, res))
    # third test: symbolic case
    x, y, z = SR.var("x, y, z")
    a = [1, 2, 3]
    b = [x, y, z]
    res = dot_product(a, b)
    correct = x + 2*y + 3*z
    if res != correct:
        raise RuntimeError("dot_product(%s, %s) = %s but should be %s" % (a, b, res, correct))


def solve_for_coefficients_homog(f, k, coeff):
    """
    Given a symbolic expression f, a distinguished variable k, and a
    set of variables coeff such that f is polynomial in k and
    (homogeneous) linear in coeff, return a nonzero solution for coeff
    in the equation f = 0.

    If the only solution is zero this raises an IndexError.

    Used by fasenmyer_kfree().
    """
    rels = [r for [r, _] in f.coefficients(k)]
    rows = [[r.coefficient(c).simplify_full() for c in coeff] for r in rels]
    m = matrix(SR, rows)
    K = m.right_kernel()
    return K.basis()[0]


def irreducible_dispersion(s, t, k):
    """
    Used by gosper_certificate().
    """
    m = s.degree(k)
    n = t.degree(k)
    if m != n:
        return []
    if n == 0:
        return []
    a = s.coefficient(k, n)
    b = s.coefficient(k, n-1)
    c = t.coefficient(k, n)
    d = t.coefficient(k, n-1)
    j = (b*c - a*d) / (a*c*n)
    if j < 0 or not j in ZZ:
        return []
    if c * s - a * t.subs({k: k+j}) == 0:
        return [j]
    return []


def _test_irreducible_dispersion():
    k = SR.var("k")
    # first test
    s = k
    t = k - 97
    res = irreducible_dispersion(s, t, k)
    correct = [97]
    if res != correct:
        raise RuntimeError("irreducible_dispersion(%s, %s, %s) = %s but should be %s" % (s, t, k, res, correct))
    # second test
    s = k^2 + 5
    t = k^2 + 2
    res = irreducible_dispersion(s, t, k)
    correct = []
    if res != correct:
        raise RuntimeError("irreducible_dispersion(%s, %s, %s) = %s but should be %s" % (s, t, k, res, correct))


def dispersion_set(q, r, k):
    """
    Used by gosper_certificate().
    """
    J = []
    for s, _ in q.factor_list():
        for t, _ in r.factor_list():
            J = J + irreducible_dispersion(s, t, k)
    return J


def _test_dispersion_set():
    k = SR.var("k")
    # first test
    q = (k+5)^2 * (2*k+7)
    r = k^2 * (2*k+1)
    res = dispersion_set(q, r, k)
    correct = [3, 5]
    if res != correct:
        raise RuntimeError("dispersion_set(%s, %s, %s) = %s but should be %s" % (q, r, k, res, correct))
    # second test
    q = (k+2) * (k+3)
    r = k^2+1
    res = dispersion_set(q, r, k)
    correct = []
    if res != correct:
        raise RuntimeError("dispersion_set(%s, %s, %s) = %s but should be %s" % (q, r, k, res, correct))


def find_polys(n, d, k):
    """
    Used by gosper_certificate().
    """
    R = parent(n)
    p = R(1)
    q = n.subs({k: k-1}).simplify_full()
    r = d.subs({k: k-1}).simplify_full()
    J = dispersion_set(q, r, k)
    while len(J) > 0:
        j = J[0]
        g = gcd([q, r.subs({k: k+j})])
        p = p * prod([g.subs({k: k-i}) for i in range(j)])
        q = (q / g).simplify_full()
        r = (r / g.subs({k: k-j})).simplify_full()
        J = dispersion_set(q, r, k)
    return (p, q, r)


def _test_find_polys():
    k = SR.var("k")
    # first test
    n = (k+1)^3
    d = k^3
    res = find_polys(n, d, k)
    correct = (k^3, 1, 1)
    if res != correct:
        raise RuntimeError("find_polys(%s, %s, %s) = %s but should be %s" % (n, d, k, res, correct))
    # second test
    n = k
    d = k+1
    res = find_polys(n, d, k)
    correct = (1, k-1, k)
    if res != correct:
        raise RuntimeError("find_polys(%s, %s, %s) = %s but should be %s" % (n, d, k, res, correct))
    # third test
    x = SR.var("x")
    n = x-k
    d = k+1
    res = find_polys(n, d, k)
    correct = (1, x-k+1, k)
    if res != correct:
        raise RuntimeError("find_polys(%s, %s, %s) = %s but should be %s" % (n, d, k, res, correct))


def degree_bound(p, q, r, k):
    """
    Used by gosper_certificate().
    """
    sigma = (q.subs({k: k+1}) + r).simplify_full()
    delta = (q.subs({k: k+1}) - r).simplify_full()
    if sigma == 0:
        s = -Infinity
    else:
        s = sigma.degree(k)
    if delta == 0:
        d = -Infinity
    else:
        d = delta.degree(k)
    if s <= d:
        return p.degree(k) - d
    a = sigma.coefficient(k, s)
    b = delta.coefficient(k, s-1)
    n = -2 * b / a
    if n < 0 or not n in ZZ:
        return p.degree(k) - s + 1
    return max([n, p.degree(k) - s + 1])


def _test_degree_bound():
    k = SR.var("k")
    # first test
    p = k^3
    q = SR(1)
    r = SR(1)
    res = degree_bound(p, q, r, k)
    correct = 4
    if res != correct:
        raise RuntimeError("degree_bound(%s, %s, %s, %s) = %s but should be %s" % (p, q, r, k, res, correct))
    # second test
    x = SR.var("x")
    p = SR(1)
    q = x-k+1
    r = k
    res = degree_bound(p, q, r, k)
    correct = -1
    if res != correct:
        raise RuntimeError("degree_bound(%s, %s, %s, %s) = %s but should be %s" % (p, q, r, k, res, correct))
    # third test
    x = SR.var("x")
    p = SR(1)
    q = k-x-1
    r = k
    res = degree_bound(p, q, r, k)
    correct = 0
    if res != correct:
        raise RuntimeError("degree_bound(%s, %s, %s, %s) = %s but should be %s" % (p, q, r, k, res, correct))


def gosper_sum(a, k):
    R = gosper_certificate(a, k)
    return R * a


def _test_gosper_sum():
    # first test
    a = 1/k - 1/(k+1)
    res = gosper_sum(a, k)
    res = res.simplify_full().factor()
    correct = -1/k
    if res != correct:
        raise RuntimeError("gosper_sum(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # second test
    a = factorial(k)
    passed = False
    try:
        res = gosper_sum(a, k)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("gosper_sum(%s, %s) = %s but should have raised error" % (a, k, res))


def gosper_certificate(a, k):
    """
    Main ingredient for gosper_sum().
    """
    n, d = (a.subs({k: k+1}) / a).simplify_factorial().numerator_denominator()
    p, q, r = find_polys(n, d, k)
    N = degree_bound(p, q, r, k)
    if N < 0:
        raise RuntimeError("%s is not Gosper-summable in %s" % (a, k))
    # declare the variables c[0], ..., c[N]
    coeff = [SR.var("c_%s" % j) for j in range(N + 1)]
    f = sum([coeff[j] * k^j for j in range(N+1)])
    eq = p - q.subs({k: k+1}) * f + r * f.subs({k: k-1})
    try:
        sol = solve_for_coefficients(eq, k, coeff)
    except ValueError:
        raise RuntimeError("%s is not Gosper-summable in %s" % (a, k))
    f_sol = sum([sol[j] * k^j for j in range(N+1)])
    R = r / p * f_sol.subs({k: k-1})
    return R


def _test_gosper_certificate():
    k = SR.var("k")
    # first test
    a = k^3
    res = gosper_certificate(a, k)
    res = res.simplify_full().factor()
    correct = (k-1)^2 / (4 * k)
    if res != correct:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # second test
    a = 1/k
    passed = False
    try:
        res = gosper_certificate(a, k)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should have raised error" % (a, k, res))
    # third test
    a = 1 / (k * (k+6))
    res = gosper_certificate(a, k)
    res = res.simplify_full().factor()
    correct = -(3*k^4+30*k^3+95*k^2+100*k+24) * (2*k+5) * (k+6) / (6 * (k+5) * (k+4) * (k+3) * (k+2) * (k+1))
    if res != correct:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # fourth test
    n = var("n")
    a = binomial(n, k)
    passed = False
    try:
        res = gosper_certificate(a, k)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should have raised error" % (a, k, res))
    # fifth test
    n = var("n")
    a = (-1)^k * binomial(n, k)
    res = gosper_certificate(a, k)
    correct = -k/n
    if res != correct:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # sixth test
    a = k * factorial(k)
    res = gosper_certificate(a, k)
    correct = 1/k
    if res != correct:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # seventh test
    a = (4*k+1) * factorial(k) / factorial(2*k+1)
    res = gosper_certificate(a, k)
    correct = -2 * (2*k+1) / (4*k+1)
    if res != correct:
        raise RuntimeError("gosper_certificate(%s, %s) = %s but should be %s" % (a, k, res, correct))


def gosper_verify(a, R, k):
    """
    Verify that R is indeed a rational certificate for Gosper's algorithm
    on the term a in the variable k.
    """
    check = ((R + 1) * a - a.subs({k: k+1}) * R.subs({k: k+1})).simplify_full()
    if check == 0:
        return True
    return False


def _test_gosper_verify():
    k = SR.var("k")
    # first test
    a = k
    R = (k - 1) / 2
    res = gosper_verify(a, R, k)
    correct = True
    if res != correct:
        raise RuntimeError("gosper_verify(%s, %s, %s) = %s but should be %s" % (a, R, k, res, correct))
    # second test
    a = k^2
    R = k + 1
    res = gosper_verify(a, R, k)
    correct = False
    if res != correct:
        raise RuntimeError("gosper_verify(%s, %s, %s) = %s but should be %s" % (a, R, k, res, correct))
    # third test
    a = 1/k - 1/(k + 1)
    R = -k - 1
    res = gosper_verify(a, R, k)
    correct = True
    if res != correct:
        raise RuntimeError("gosper_verify(%s, %s, %s) = %s but should be %s" % (a, R, k, res, correct))
    # fourth test
    n = var("n")
    a = (-1)^k * binomial(n, k)
    R = -k / n
    res = gosper_verify(a, R, k)
    correct = True
    if res != correct:
        raise RuntimeError("gosper_verify(%s, %s, %s) = %s but should be %s" % (a, R, k, res, correct))
    # fifth test
    a = k * factorial(k)
    R = 1 / k
    res = gosper_verify(a, R, k)
    correct = True
    if res != correct:
        raise RuntimeError("gosper_verify(%s, %s, %s) = %s but should be %s" % (a, R, k, res, correct))


def wz_certificate(a, n, k):
    d = a.subs({n: n+1, k: k}) - a
    Rg = gosper_certificate(d, k)
    R = Rg * (a.subs({n: n+1, k: k}) / a - 1)
    return R


def _test_wz_certificate():
    n, k = SR.var("n, k")
    # first test
    a = binomial(n, k) / 2^n
    res = wz_certificate(a, n, k).simplify_full().factor()
    correct = k / (2*(k - n - 1))
    if res != correct:
        raise RuntimeError("wz_certificate(%s, %s, %s) = %s but should be %s" % (a, n, k, res, correct))
    # second test
    a = binomial(n, k)
    passed = False
    try:
        res = wz_certificate(a, n, k)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("wz_certificate(%s, %s, %s) = %s but should have raised error" % (a, n, k, res))
    # third test
    # sadly, fails at the moment
    # x, y = SR.var("x, y")
    # a = binomial(n, k) * x^k * y^(n-k) / (x+y)^n
    # res = wz_certificate(a, n, k)
    # correct = k / (2*(k - n - 1))
    # if res != correct:
    #     raise RuntimeError("wz_certificate(%s, %s, %s) = %s but should be %s" % (a, n, k, res, correct))


def wz_verify(a, R, n, k):
    """
    Verify that R is a Wilf-Zeilberger certificate for the term a in variables n and k.
    """
    check = a * (R - 1) + a.subs({n: n+1}) - R.subs({k: k+1}) * a.subs({k: k+1})
    if check.simplify_full() == 0:
        return True
    return False


def _test_wz_verify():
    n, k = SR.var("n, k")
    # first test
    a = binomial(n, k) / 2^n
    R = k / (k - n - 1)
    res = wz_verify(a, R, n, k)
    correct = False
    if res != correct:
        raise RuntimeError("wz_verify(%s, %s, %s, %s) = %s but should be %s" % (a, R, n, k, res, correct))
    # second test
    a = binomial(n, k) / 2^n
    R = k / (2*(k - n - 1))
    res = wz_verify(a, R, n, k)
    correct = True
    if res != correct:
        raise RuntimeError("wz_verify(%s, %s, %s, %s) = %s but should be %s" % (a, R, n, k, res, correct))


def zeilberger(F, n, k, J=1):
    sigma = [SR.var("sigma_%s" %j) for j in range(1, J+1)]
    a = F + sum([sigma[j-1] * F.subs({n: n+j, k: k}) for j in range(1, J+1)])
    num, denom = (a.subs({k: k+1}) / a).simplify_factorial().numerator_denominator()
    p, q, r = find_polys(num, denom, k)
    N = degree_bound(p, q, r, k)
    if N < 0:
        raise RuntimeError("%s is not Gosper-summable in %s" % (a, k))
    # declare the variables c[0], ..., c[N]
    coeff = [SR.var("c_%s" % j) for j in range(N + 1)]
    f = sum([coeff[j] * k^j for j in range(N+1)])
    eq = p - q.subs({k: k+1}) * f + r * f.subs({k: k-1})
    eq = eq.simplify_full()
    try:
        sol = solve_for_coefficients(eq, k, coeff+sigma)
    except ValueError:
        raise RuntimeError("%s is not Gosper-summable in %s" % (a, k))
    return [SR(1)] + list(sol[N+1:])


def _test_zeilberger():
    n, k = SR.var("n, k")
    # first test
    f = binomial(n, k)
    J = 1
    res = zeilberger(f, n, k, J)
    correct = [1, -1/2]
    if res != correct:
        raise RuntimeError("zeilberger(%s, %s, %s, %s) = %s but should be %s" % (f, n, k, J, res, correct))
    # second test
    f = binomial(n, k)^2
    J = 1
    res = zeilberger(f, n, k, J)
    correct = [1, -(n+1)/(2*(2*n+1))]
    if res != correct:
        raise RuntimeError("zeilberger(%s, %s, %s, %s) = %s but should be %s" % (f, n, k, J, res, correct))


def solve_for_coefficients(f, k, coeff):
    rels = [r for [r, _] in f.coefficients(k)]
    rows = [[r.coefficient(c).simplify_full() for c in coeff] for r in rels]
    lst = list(f.variables())
    lst.remove(k)
    for c in coeff:
        if c in lst:
            lst.remove(c)
    if len(lst) == 0:
        # annoying special case
        S = QQ
    else:
        S = PolynomialRing(QQ, lst).fraction_field()
    m = matrix(S, rows)
    consts = matrix(S, [-(rels[j] - dot_product(rows[j], coeff)).simplify_full() for j in range(len(rels))]).transpose()
    sol = m.solve_right(consts)
    return vector(sol)


def _test_solve_for_coefficients():
    k = SR.var("k")
    # first test
    c = SR.var("c")
    coeff = [c]
    f = SR(1)
    passed = False
    try:
        res = solve_for_coefficients(f, k, coeff)
    except ValueError:
        passed = True
    if not passed:
        raise RuntimeError("solve_for_coefficients(%s, %s, %s) = %s but should have raised error" % (f, k, coeff, res))


def _test_all():
    _test_fasenmyer_kfree()
    _test_irreducible_dispersion()
    _test_dispersion_set()
    _test_find_polys()
    _test_degree_bound()
    _test_dot_product()
    _test_solve_for_coefficients()
    _test_gosper_sum()
    _test_gosper_certificate()
    _test_gosper_verify()
    #_test_wz_certificate()
    _test_wz_verify()
    _test_zeilberger()

