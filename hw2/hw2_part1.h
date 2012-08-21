// File:    hw2_part1.h
// Author:  Garrett Smith
// Created: 04/01/2011

#ifndef HW2_PART1__H
#define HW2_PART1__H

#include <vector>
#include <map>
#include <gmpxx.h>
#include "gf2_p.h"

typedef std::pair<int, int> ipair;
typedef std::vector<ipair> pvec;
typedef std::vector<int> ivec;

typedef std::pair<int, gf2_p> gpair;
typedef std::vector<gpair> gvec;

// routines for Zp^k

int totient_prime(int p, int k);
int shanks_zp(int p, int n, const mpz_class &az, const mpz_class &bz, int *pa);
int pohlig_hellman(int p, int n, int alpha, int beta, int q, int c, ivec &a);
int pohlig_hellman_print(int p, int n, int alpha, int beta, int q, int c);
void compute_Zpk(int p, int k, ivec &orders, ivec &prims, ivec &powers);

// routines for GF(p^2)

int shanks_gf(int n, const gf2_p &a, const gf2_p &b, int *pa);
void compute_GFp2(void);

#endif // HW2_PART1__H

