#ifndef SMART_FACT_H
#define SMART_FACT_H

#include <gmpxx.h>

inline mpz_class partial_fact(mpz_class a, mpz_class b)
{
    mpz_class res = b;
    for (; a < b; a += 1)
    {
        res *= a;
    }
    return res;
}

inline mpz_class ProdTree(long l, long r)
{
    if (l > r)
        return 1;
    if (l == r)
        return l;
    if (r - l == 1)
        return mpz_class(r * l);
    long m = (l + r) / 2;
    return ProdTree(l, m) * ProdTree(m + 1, r);
}

inline mpz_class FactTree(long n)
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    if (n == 1 || n == 2)
        return n;
    return ProdTree(2, n);
}

#endif // SMART_FACT_H

