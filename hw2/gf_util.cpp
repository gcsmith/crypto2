// File:    gf_util.cpp
// Author:  Garrett Smith
// Created: 04/01/2011

#include <algorithm>
#include <cmath>
#include "hw2_part1.h"

// -----------------------------------------------------------------------------
// Comparison function for sorting Galois pairs by their polynomial.
bool shanks_cmp_gf(gpair a, gpair b)
{
    return a.second < b.second;
}

// -----------------------------------------------------------------------------
// Implements Shanks algorithm for performing discrete logarithms in GF(131^2).
int shanks_gf(int n, const gf2_p &a, const gf2_p &b, int *pa)
{
    int m = ceil(sqrt(n));
    gvec L1(m), L2(m);

    gf2_p am = a.pow(m);
    for (int j = 0; j < m; ++j) {
        gf2_p amj = am.pow(j);
        L1[j] = std::make_pair(j, amj);
    }

    for (int i = 0; i < m; ++i) {
        gf2_p aii = b * a.pow(i).inverse();
        L2[i] = std::make_pair(i, aii);
    }

    std::sort(L1.begin(), L1.end(), shanks_cmp_gf);
    std::sort(L2.begin(), L2.end(), shanks_cmp_gf);

    for (int i = 0, j = 0; (i < m) && (j < m); ) {
        if (L1[i].second == L2[j].second) {
            mpz_class tz(m * L1[i].first + L2[j].first), nz(n), kz;
            mpz_mod(kz.get_mpz_t(), tz.get_mpz_t(), nz.get_mpz_t());
            *pa = kz.get_ui();
            return 0;
        }
        else if (L1[i].second < L2[j].second) ++i;
        else if (L2[j].second < L1[i].second) ++j;
    }

    return -1;
}

// -----------------------------------------------------------------------------
// Compute orders and primitive monic polynomials for GF(p^2)
void compute_GFp2(void)
{
    const int phi = 17160; // (p - 1)(p + 1)

    int result;
    gf2_p base(1, 3), beta(1, 101);
    shanks_gf(phi, base, beta, &result);
    std::cout << "shanks(" << base << ", " << beta << ") = " << result << "\n";
    std::cout << "test: " << base.pow(result) << "\n";

#if 0
    for (int a = 1; a < 100; ++a) {
        for (int b = 1; b < 100; ++b) {
            gf2_p p(a, b);
            int order = p.order(phi);
            cout << "order of " << p << " is " << order;
            cout << ((order == phi) ? " primitive\n" : "\n");
        }
    }
#endif
}

