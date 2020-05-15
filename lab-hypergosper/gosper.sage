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
