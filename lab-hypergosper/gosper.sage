def irreducible_dispersion(s, t, k):
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


def dispersion_set(q, r, k):
    J = []
    for s, _ in q.factor_list():
        for t, _ in r.factor_list():
            J = J + irreducible_dispersion(s, t, k)
    return J


def find_polys(n, d, k):
    R = parent(n)
    p = R(1)
    q = n.subs({k: k-1})
    r = d.subs({k: k-1})
    J = dispersion_set(q, r, k)
    while len(J) > 0:
        j = J[0]
        g = gcd([q, r.subs({k: k+j})])
        p = p * prod([g.subs({k: k-i}) for i in range(j)])
        q = q / g
        r = r / g.subs({k: k-j})
        J = dispersion_set(q, r, k)
    return (p, q, r)


def degree_bound(p, q, r, k):
    sigma = q.subs({k: k+1}) + r
    delta = q.subs({k: k+1}) - r
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


def gosper_sum(a, k):
    n, d = (a.subs({k: k+1}) / a).numerator_denominator()
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
    A = r / p * f_sol.subs({k: k-1}) * a
    return A


def dot_product(a, b):
    if len(a) != len(b):
        raise ValueError("%s and %s must have the same length" %(a, b))
    return sum([a[j] * b[j] for j in range(len(a))])


def solve_for_coefficients(f, k, coeff):
    rels = [r for [r, _] in f.coefficients(k)]
    rows = [[r.coefficient(c).simplify_full() for c in coeff] for r in rels]
    consts = [-(rels[j] - dot_product(rows[j], coeff)).simplify_full() for j in range(len(rels))]
    m = matrix(SR, rows)
    return m.solve_right(vector(consts))


def _test_all():
    _test_irreducible_dispersion()
    _test_dispersion_set()
    _test_find_polys()
    _test_degree_bound()
    _test_dot_product()
    _test_solve_for_coefficients()
    _test_gosper_sum()


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


def _test_solve_for_coefficients():
    pass


def _test_gosper_sum():
    k = SR.var("k")
    # first test
    a = k^3
    res = gosper_sum(a, k)
    res = res.simplify_full().factor()
    correct = k^2 * (k-1)^2 / 4
    if res != correct:
        raise RuntimeError("gosper_sum(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # second test
    a = 1/k
    passed = False
    try:
        res = gosper_sum(a, k)
    except RuntimeError:
        passed = True
    if not passed:
        raise RuntimeError("gosper_sum(%s, %s) = %s but should have raised error" % (a, k, res))
    # third test
    a = 1/k - 1/(k+1)
    res = gosper_sum(a, k)
    res = res.simplify_full().factor()
    correct = -1/k
    if res != correct:
        raise RuntimeError("gosper_sum(%s, %s) = %s but should be %s" % (a, k, res, correct))
    # fourth test
    a = 1 / (k * (k+6))
    res = gosper_sum(a, k)
    res = res.simplify_full().factor()
    correct = -(3*k^4+30*k^3+95*k^2+100*k+24) * (2*k+5) / (6 * (k+5) * (k+4) * (k+3) * (k+2) * (k+1) * k)
    if res != correct:
        raise RuntimeError("gosper_sum(%s, %s) = %s but should be %s" % (a, k, res, correct))
