// File:    zp_util.cpp
// Author:  Garrett Smith
// Created: 04/01/2011

#include <algorithm>
#include <cmath>
#include "hw2_part1.h"

// -----------------------------------------------------------------------------
// Compute Euler's totient function, where p is a prime number raised to k.
int totient_prime(int p, int k)
{
    mpz_class pz(p), a, b, c;
    mpz_pow_ui(a.get_mpz_t(), pz.get_mpz_t(), k);
    mpz_pow_ui(b.get_mpz_t(), pz.get_mpz_t(), k - 1);
    return mpz_class(a - b).get_ui();
}

// -----------------------------------------------------------------------------
// Comparison function for sorting integer pairs by their second element.
bool shanks_cmp(ipair a, ipair b)
{
    return a.second < b.second;
}

// -----------------------------------------------------------------------------
// Implements Shanks algorithm for computing discrete logarithms in Zp.
int shanks_zp(int p, int n, const mpz_class &az, const mpz_class &bz, int *pa)
{
    int m = ceil(sqrt(n));
    mpz_class tz, jz, kz, mz(m), nz(n), pz(p);
    pvec L1(m), L2(m);

    mpz_powm(jz.get_mpz_t(), az.get_mpz_t(), mz.get_mpz_t(), pz.get_mpz_t());
    for (int j = 0; j < m; ++j) {
        mpz_powm_ui(kz.get_mpz_t(), jz.get_mpz_t(), j, pz.get_mpz_t());
        L1[j] = std::make_pair(j, kz.get_ui());
    }

    for (int i = 0; i < m; ++i) {
        mpz_pow_ui(kz.get_mpz_t(), az.get_mpz_t(), i);
        mpz_invert(kz.get_mpz_t(), kz.get_mpz_t(), pz.get_mpz_t());
        mpz_mul(kz.get_mpz_t(), kz.get_mpz_t(), bz.get_mpz_t());
        mpz_mod(kz.get_mpz_t(), kz.get_mpz_t(), pz.get_mpz_t());
        L2[i] = std::make_pair(i, kz.get_ui());
    }

    std::sort(L1.begin(), L1.end(), shanks_cmp);
    std::sort(L2.begin(), L2.end(), shanks_cmp);

    for (int i = 0, j = 0; (i < m) && (j < m); ) {
        if (L1[i].second == L2[j].second) {
            tz = m * L1[i].first + L2[j].first;
            mpz_mod(kz.get_mpz_t(), tz.get_mpz_t(), nz.get_mpz_t());
            *pa = kz.get_ui();
            return 0;
        }
        else if (L1[i].second < L2[j].second) ++i;
        else if (L1[i].second > L2[j].second) ++j;
    }

    return -1;
}

// -----------------------------------------------------------------------------
// Implements Pohlig-Hellman algorithm for computing discrete logarithms.
int pohlig_hellman(int p, int n, int alpha, int beta, int q, int c, ivec &a)
{
    std::vector<mpz_class> bz(c + 1);
    mpz_class az(alpha), cz, dz, lz, z, pz(p);
    int bp = n, i, qj = 1;

    bz[0] = beta;
    a.resize(c);
    mpz_powm_ui(lz.get_mpz_t(), az.get_mpz_t(), n / q, pz.get_mpz_t());
    for (int j = 0; j < c; ++j) {
        bp /= q;
        mpz_powm_ui(dz.get_mpz_t(), bz[j].get_mpz_t(), bp, pz.get_mpz_t());
        if (0 > shanks_zp(p, n, lz, dz, &i)) return -1;
        a[j] = i;
        mpz_pow_ui(cz.get_mpz_t(), az.get_mpz_t(), a[j] * qj);
        mpz_invert(cz.get_mpz_t(), cz.get_mpz_t(), pz.get_mpz_t());
        mpz_mul(cz.get_mpz_t(), cz.get_mpz_t(), bz[j].get_mpz_t());
        bz[j + 1] = cz;
        qj *= q;
    }

    return 0;
}

// -----------------------------------------------------------------------------
// Perform PH algorithm and display both the input parameters and the result.
int pohlig_hellman_print(int p, int n, int alpha, int beta, int q, int c)
{
    ivec an;
    if (0 > pohlig_hellman(p, n, alpha, beta, q, c, an)) {
        std::cout << "PH failed for a=" << alpha << " b=" << beta << "\n";
        return -1;
    }

    int qn = 1, a = 0;
    std::cout << "[ n=" << n << " alpha=" << alpha << " beta=" << beta
              << " q=" << q << " c=" << c << " ]  ->  { ";
    for (int i = 0; i < c; ++i) {
        std::cout << an[i] << " ";
        a += an[i] * qn;
        qn *= q;

    }
    std::cout << "}  ->  " << a << std::endl;
    return 0;
}

// -----------------------------------------------------------------------------
// Calculate the order of element g modulo n.
int calc_order(int g, const mpz_class &n)
{
    mpz_class gz(g), oz;
    for (unsigned int m = 1; m < n.get_ui(); ++m) {
        mpz_powm_ui(oz.get_mpz_t(), gz.get_mpz_t(), m, n.get_mpz_t());
        if (1 == oz) return m;
    }
    return n.get_ui();
}

// -----------------------------------------------------------------------------
// Compute orders, primitives, and powers of elements in Zp^k.
void compute_Zpk(int p, int k, ivec &orders, ivec &prims, ivec &powers)
{
    int phi = totient_prime(p, k);
    
    mpz_class n;
    mpz_pow_ui(n.get_mpz_t(), mpz_class(p).get_mpz_t(), k);

    orders.push_back(1);
    for (int i = 1; i < n; ++i) {
        int o = calc_order(i, n);
        orders.push_back(o);
        if (o == phi) {
            prims.push_back(i);
            // std::cout << (float)i / (float)n.get_ui() << "%\n";
        }
    }

    mpz_class gz(prims[0]), o;
    for (int m = 0; m < n; ++m) {
        mpz_powm_ui(o.get_mpz_t(), gz.get_mpz_t(), m, n.get_mpz_t());
        powers.push_back(o.get_ui());
    }
}

