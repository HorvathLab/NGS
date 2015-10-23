from scipy.stats.distributions import hypergeom, binom
from math import log


def memoize(f):
    cache = {}

    def decorated_function(*args, **kw):
        key = (args, tuple(sorted(kw.iteritems())))
        if key not in cache:
            cache[key] = f(*args, **kw)
        return cache[key]
    return decorated_function


@memoize
def binom_test_low(n, N, p):
    return binom.cdf(n, N, p)


@memoize
def binom_test_high(n, N, p):
    if n == 0:
        return 1.0
    return min(max(0.0, binom.sf(n - 1, N, p)), 1.0)


@memoize
def binom_test(n, N, p, direction=None):
    if direction in ('high', 1):
        return binom_test_high(n, N, p)
    elif direction in ('low', -1):
        return binom_test_low(n, N, p)
    return min(binom_test_low(n, N, p), binom_test_high(n, N, p))


@memoize
def pvscore(x):
    if x < 1e-10:
        x = 1e-10
    return max(0.0, -10.0 * log(x, 10.0))


@memoize
def lod(x, N, n, M, pseudocount=0.5, base=2.0):
    n00 = x
    n01 = (N - x)
    n10 = (n - x)
    n11 = (M - n00 - n01 - n10)
    if min(n00, n01, n10, n11) < 0 or \
            (pseudocount == 0.0 and min(n00, n01, n10, n11) <= 0):
        return None
    n00 += pseudocount
    n01 += pseudocount
    n10 += pseudocount
    n11 += pseudocount
    return (log(n00, base) - log(n01, base)) - (log(n10, base) - log(n11, base))


@memoize
def fisher_exact_high(x, N, n, M):
    return min(1, max(0, sum(hypergeom.pmf(x1, M, n, N) for x1 in range(x, N + 1))))
    # return min(1,max(0,hypergeom.sf(x-1,M,n,n)))


@memoize
def fisher_exact_low(x, N, n, M):
    return min(1, max(0, sum(hypergeom.pmf(x1, M, n, N) for x1 in range(0, x + 1))))
    # return min(1,max(0,hypergeom.cdf(x,M,n,N)))


@memoize
def fisher_exact(x, N, n, M, direction=None):
    if direction in ('high', 1):
        return fisher_exact_high(x, N, n, M)
    elif direction in ('low', -1):
        return fisher_exact_low(x, N, n, M)
    return min(fisher_exact_low(x, N, n, M), fisher_exact_high(x, N, n, M))


def bonferroni(pvs):
    n = len(pvs)
    return map(lambda pv: min(n * pv, 1), pvs)


def fdr(pvs):
    n = len(pvs)
    ind = sorted(range(n), key=lambda i: pvs[i])
    fdr = [-1.0] * n
    for i in range(n):
        fdr[ind[i]] = min(pvs[ind[i]] * n / (i + 1), 1)
    if n > 2:
        for i in range(n - 2, -1, -1):
            fdr[ind[i]] = min(fdr[ind[i]], fdr[ind[i + 1]])
    return fdr

if __name__ == '__main__':
    import sys
    tests = """
	9 1 0 10
        1 9 6 4
        1 9 0 10
        7 1 0 5
        2000 1500 6711 5323
    """
    for l in tests.splitlines():
        sl = map(int, l.split())
        if len(sl) != 4:
            continue
        x = sl[0]
        N = sl[0] + sl[1]
        n = sl[0] + sl[2]
        M = sum(sl)
        v = (lod(x, N, n, M), fisher_exact_low(x, N, n, M),
             fisher_exact_high(x, N, n, M), fisher_exact(x, N, n, M))
        print v[0], v[1], v[2], v[3]
        if min(v[1:3]) == 0:
            print "Problem with hypergeometic precision..."
