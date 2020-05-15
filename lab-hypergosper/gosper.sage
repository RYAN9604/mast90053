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
    return set(J)


def find_polys(n, d, k):
    R = parent(n)
    p = R(1)
    q = n.subs({k: k-1})
    r = d.subs({k: k-1})
    J = dispersion_set(q, r, k)
    while len(J) > 0:
        j = list(J)[0]
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



