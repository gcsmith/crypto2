// File:    gf2_p.h
// Author:  Garrett Smith
// Created: 04/02/2011

#ifndef GF2_P__H
#define GF2_P__H

#include <iostream>
#include <gmpxx.h>

class gf2_p {
public:
    gf2_p() { }
    gf2_p(int _a, int _b) : a(_a), b(_b) { }
    gf2_p(mpz_class _a, mpz_class _b) : a(_a), b(_b) { }
    gf2_p(const gf2_p &p) : a(p.a), b(p.b) { }
    ~gf2_p() { }

    const gf2_p operator+(const gf2_p &rh) const {
        gf2_p r(a + rh.a, b + rh.b);
        mpz_mod_ui(r.a.get_mpz_t(), r.a.get_mpz_t(), 131);
        mpz_mod_ui(r.b.get_mpz_t(), r.b.get_mpz_t(), 131);
        return r;
    }

    const gf2_p operator*(const gf2_p &rh) const {
        gf2_p r(a * rh.b + b * rh.a, b * rh.b - a * rh.a);
        mpz_mod_ui(r.a.get_mpz_t(), r.a.get_mpz_t(), 131);
        mpz_mod_ui(r.b.get_mpz_t(), r.b.get_mpz_t(), 131);
        return r;
    }

    gf2_p &operator+=(const gf2_p &rh) {
        a = (a + rh.a); b = (b + rh.b);
        mpz_mod_ui(a.get_mpz_t(), a.get_mpz_t(), 131);
        mpz_mod_ui(b.get_mpz_t(), b.get_mpz_t(), 131);
        return *this;
    }

    gf2_p &operator*=(const gf2_p &rh) {
        mpz_class a1 = a, b1 = b;
        a = a1 * rh.b + b1 * rh.a; b = b1 * rh.b - a1 * rh.a;
        mpz_mod_ui(a.get_mpz_t(), a.get_mpz_t(), 131);
        mpz_mod_ui(b.get_mpz_t(), b.get_mpz_t(), 131);
        return *this;
    }

    bool operator<(const gf2_p &rh) {
        return (a == rh.a) ? (b < rh.b) : (a < rh.a);
    }

    bool operator==(const gf2_p &rh) {
        return (a == rh.a) && (b == rh.b);
    }

    bool is_one(void) const { return (a == 0) && (b == 1); }
    bool is_monic(void) const { return (a == 1); }

    const gf2_p pow(int n) const {
        if (n == 0) return gf2_p(0, 1);
        gf2_p t(*this);
        for (int i = 1; i < n; ++i) {
            t *= *this;
        }
        return t;
    }

    const gf2_p inverse(void) const {
        for (int a = 0; a < 131; ++a) {
            for (int b = 0; b < 131; ++b) {
                gf2_p t(a, b);
                if ((*this * t).is_one()) return t;
            }
        }
        return gf2_p(0, 0);
    }

    int order(int n) const {
        gf2_p t(*this);
        for (int i = 1; i < n; ++i) {
            if (t.is_one()) return i;
            t *= *this;
        }
        return n;
    }

    mpz_class a, b;
};

// -----------------------------------------------------------------------------
inline std::ostream &operator<<(std::ostream &out, const gf2_p &p)
{
    out << p.a << "x + " << p.b;
    return out;
}

#endif

