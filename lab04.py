from numpy.polynomial.polynomial import Polynomial as poly
# import math
import numpy as np
from lab01t2 import func
from typing import Union

def PolyInverseModOverQ(f: poly, r):  # f*g=1modx^r
    if r < 1:
        raise ValueError
    l = int(np.ceil(np.log2(r)))
    if f.coef[0] != 0:
        g=poly([1 / f.coef[0]])
    else:
        return None
    x = 2
    if g==None:
        return None
    for number in range(l):
        g = (2 * g - f * g ** 2).truncate(x)
        x <<= 1
    return g


def PolyInverseModOverZn(f: poly, r: int, n: int) :
    if r < 1 or n < 1:
        raise ValueError

    if f.coef[0] == 0 or n == 1:
        return None

    f = poly(np.mod(f.coef, n))

    c = func(f.coef[0], n)
    if c is None:
        return None

    g = poly([c])
    r = int(np.ceil(np.log2(r)))
    x = 2

    for i in range(r):
        g = (2 * g - f * g ** 2).truncate(x)
        g.coef = np.mod(g.coef, n)
        x <<= 1
    return g


def PolyDivModOverQ(a: poly, b: poly) -> (poly, poly):
    if not b.coef.any(): #если полином нулевой
        raise ZeroDivisionError
    a = a.trim()
    b = b.trim()
    n, m = len(a.coef), len(b.coef)  # степени полиномов
    if n < m:
        return poly([0]), a
    else:
        f = poly(b.coef[::-1]).trim() #f:=rev_m(b)
        g = PolyInverseModOverQ(f, n - m + 1) #q inverse of f modulo x^(n-m+1)
        q = (poly(a.coef[::-1]).trim() * g).truncate(n - m + 1)#q:=rev_n(a)g mod x^(n-m+1)
        q = poly(q.coef[::-1]).trim() #q:=rev_(n-m)(q)
        if len(q.coef) < n - m + 1:
            q.coef = np.concatenate([np.zeros(n - len(q)), q.coef])
        r = a - b * q #r:=a-bq
        r = r.trim()
        q = q.trim()
        return q, r


def PolyDivModOverZn(a: poly, b: poly, n: int) -> (poly, poly):
    if n < 1:
        raise ValueError
    if n==1:
        raise ZeroDivisionError
    a = poly(np.mod(a.coef, n))
    b = poly(np.mod(b.coef, n))
    a = a.trim()
    b = b.trim()
    deg_a, deg_b = len(a.coef), len(b.coef)
    if deg_a < deg_b:
        return poly([0]), a
    else:
        f = poly(b.coef[::-1]).trim()
        g = PolyInverseModOverZn(f, deg_a - deg_b + 1, n)
        if g is None:
            raise ZeroDivisionError
        q = (poly(a.coef[::-1]).trim() * g).truncate(deg_a - deg_b + 1)  # q:=rev_n(a)g mod x^(n-m+1)
        q = poly(np.mod(q.coef, n))
        q = poly(q.coef[::-1]).trim()  # q:=rev_(n-m)(q)
        if len(q.coef) < deg_a - deg_b + 1:
            q.coef = np.concatenate([np.zeros(deg_a - deg_b + 1 - len(q)), q.coef])
        bq = poly(np.mod((b * q).coef, n))
        r = a - bq
        r = poly(np.mod(r.coef, n))
        q = poly(np.mod(q.coef, n))
        r = r.trim()
        q = q.trim()
        return q, r

