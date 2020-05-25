def gosper_verify(a, R, k):
    check = ((R + 1) * a - a.subs({k: k+1}) * R.subs({k: k+1})).simplify_full()
    if check == 0:
        return True
    return False


def wz_verify(a, R, n, k):
    check = a * (R - 1) + a.subs({n: n+1}) - R.subs({k: k+1}) * a.subs({k: k+1})
    if check.simplify_full() == 0:
        return True
    return False


def wz_certificate(a, n, k):
    d = a.subs({n: n+1, k: k}) - a
    Rg = gosper_certificate(d, k)
    R = Rg * (a.subs({n: n+1, k: k}) / a - 1)
    return R


def _test_all():
    _test_gosper_verify()
    _test_wz_verify()
    _test_wz_certificate()


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
    x, y = SR.var("x, y")
    a = binomial(n, k) * x^k * y^(n-k) / (x+y)^n
    res = wz_certificate(a, n, k)
    correct = k / (2*(k - n - 1))
    if res != correct:
        raise RuntimeError("wz_certificate(%s, %s, %s) = %s but should be %s" % (a, n, k, res, correct))
