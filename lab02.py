import math
import numpy as np
from sympy import *
from lab01t2 import func

def IntKaratsuba(x, y):
    if x < 10 or y < 10:
        return x * y
    m = max(len(str(x)), len(str(y)))
    if m % 2 != 0:  # если нечетные
         m += 1
    m_2 = m // 2
    a, b = divmod(x, 10 ** m_2)
    c, d = divmod(y, 10 ** m_2)
    ac = IntKaratsuba(a, c)
    bd = IntKaratsuba(b, d)
    ad_bc = IntKaratsuba((a + b), (c + d)) - ac - bd
    return (ac * (10 ** m)) + bd + (ad_bc * (10 ** m_2))
#print(IntKaratsuba(7331325, 3690280))

def PolyKaratsuba2(a, b):
    if len(a) == 1 and len(b) == 1:
        return a * b

    deg_a, deg_b = len(a), len(b)
    if deg_a < deg_b:
        a, b = b, a
        deg_a, deg_b = deg_b, deg_a

    if deg_a & 1:
        deg_a += 1
        a = np.concatenate([np.zeros(1), a])

    if deg_a != deg_b:
        b = np.concatenate([np.zeros(deg_a - deg_b), b])

    a1 = a[:deg_a // 2]
    b1 = b[:deg_a // 2]
    a2 = a[deg_a // 2:]
    b2 = b[deg_a // 2:]

    result_1 = PolyKaratsuba2(a1, b1) #f
    result_2 = PolyKaratsuba2(a2, b2) #s
    result_3 = PolyKaratsuba2(np.polyadd(a1,a2), np.polyadd(b1,b2)) #m
    result_3 = np.polysub(np.polysub(result_3, result_2), result_1)

    result_1 = np.concatenate([result_1, np.zeros(deg_a)])
    result_3 = np.concatenate([result_3, np.zeros(deg_a//2)])

    # return np.polyadd(result_1, np.polyadd(result_2, result_3))
    res = np.polyadd(result_1, np.polyadd(result_2, result_3))
    #res = np.polymul(res, [1])
    return res

def PolyKaratsuba(a,b):
    res=PolyKaratsuba2(a,b)
    for x in range(len(res)):
        if res[x]:
            break
    return res[x:]

def BinPowMod (a, p, n): #a^p=modn
    if n==0:
        raise ValueError
    if a==0 and p==0:
        return None

    a %= n
    t=1
    if p<0:
        a=func(a,n)
        if a is None:
            return None
        p*=-1
    while p:
     if p & 1:
         t *= a
         t %= n
     a *= a
     a%=n
     p >>= 1
    return t

