#ifndef HW4_UTILITY_H
#define HW4_UTILITY_H

#include <iostream>
#include <vector>
#include <map>

#define Z(x)    (x).get_mpz_t()
#define mZ(x)   mpz_class(x).get_mpz_t()

typedef unsigned int uint_t;
typedef std::pair<uint_t, uint_t> ipair;
typedef std::vector<ipair> ivec;
typedef std::map<int, int> imap;
typedef std::vector<int> NAF;

// represents an elliptic curve x^3 + ax + b mod p
struct curve {
    curve(uint_t _a, uint_t _b, uint_t _p) : a(_a), b(_b), p(_p) {}
    uint_t a;
    uint_t b;
    uint_t p;
};

struct point {
    point() : x(0), y(0), infinity(true) {}
    point(uint_t _x, uint_t _y) : x(_x), y(_y), infinity(false) {}
    point(uint_t _x, uint_t _y, bool inf) : x(_x), y(_y), infinity(inf) {}

    bool operator==(const point &o) {
        return (x == o.x) && (y == o.y) && (infinity == o.infinity);
    }

    uint_t x;
    uint_t y;
    bool infinity;
};

point add(const point &a, const point &b, const curve &ec);
point sub(const point &a, const point &b, const curve &ec);
point pow(const point &a, uint_t exp, const curve &ec);

bool elliptic_add(ipair P, ipair Q, ipair *R, const curve &ec);
ipair point_pow(const ipair &pt, uint_t exp, const curve &ec, bool show);
int calc_order(const point &p, const curve &ec);
void double_and_add(point P, int n, const curve &ec);

NAF int_to_naf(int in);
void print_naf(const NAF &in);
void ext_euclid(int a, int b, int *rr, int *rs, int *rt);

// -----------------------------------------------------------------------------
inline mpz_class mod(const mpz_class &b, const mpz_class &p)
{
    mpz_class out;
    mpz_mod(Z(out), Z(b), Z(p));
    return out;
}

// -----------------------------------------------------------------------------
inline mpz_class inv(const mpz_class &b, const mpz_class &p)
{
    mpz_class out;
    mpz_invert(Z(out), Z(b), Z(p));
    return out;
}

// -----------------------------------------------------------------------------
inline mpz_class pow(const mpz_class &b, const mpz_class &k, const mpz_class &p)
{
    mpz_class out;
    mpz_powm(Z(out), Z(b), Z(k), Z(p));
    return out;
}

// -----------------------------------------------------------------------------
inline mpz_class gcd(const mpz_class &a, const mpz_class &b)
{
    mpz_class out;
    mpz_gcd(Z(out), Z(a), Z(b));
    return out;
}

// -----------------------------------------------------------------------------
inline mpz_class div(const mpz_class &a, const mpz_class &b)
{
    mpz_class out;
    mpz_divexact(Z(out), Z(a), Z(b));
    return out;
}

// -----------------------------------------------------------------------------
inline bool con(const mpz_class &a, const mpz_class &b, const mpz_class &p)
{
    return mpz_congruent_p(Z(a), Z(b), Z(p));
}

// -----------------------------------------------------------------------------
inline std::ostream &operator<<(std::ostream &out, const ipair &p)
{
    out << "(" << p.first << ", " << p.second << ")";
    return out;
}

// -----------------------------------------------------------------------------
inline std::ostream &operator<<(std::ostream &out, const point &p)
{
    if (p.infinity)
        out << "(inf)";
    else
        out << "(" << p.x << ", " << p.y << ")";
    return out;
}

// -----------------------------------------------------------------------------
inline std::ostream &operator<<(std::ostream &out, const NAF &n)
{
    out << "{ ";
    for (NAF::const_reverse_iterator i = n.rbegin(); i != n.rend(); ++i)
        out << *i << " ";
    out << "}";
    return out;
}

#endif // HW4_UTILITY_H

