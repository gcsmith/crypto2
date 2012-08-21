// File:    main.cpp
// Author:  Garrett Smith
// Created: 04/19/2011

#include <iostream>
#include <vector>
#include <map>
#include <gmpxx.h>

typedef unsigned int uint_t;
typedef std::pair<uint_t, uint_t> ipair;
typedef std::vector<ipair> ivec;
typedef std::map<int, int> imap;

// represents an elliptic curve x^3 + ax + b mod p
struct curve {
    curve(uint_t _a, uint_t _b, uint_t _p) : a(_a), b(_b), p(_p) {}

    uint_t a;
    uint_t b;
    uint_t p;
};

using namespace std;

// -----------------------------------------------------------------------------
ostream &operator<<(ostream &out, const ipair &p)
{
    out << "(" << p.first << ", " << p.second << ")";
    return out;
}

// -----------------------------------------------------------------------------
// return true if a is a quadratic residue mod p
bool is_quadratic_residue(mpz_class a, uint_t p)
{
    mpz_class r;
    mpz_pow_ui(r.get_mpz_t(), a.get_mpz_t(), (p - 1) >> 1);
    return 0 != mpz_congruent_ui_p(r.get_mpz_t(), 1, p);
}

// -----------------------------------------------------------------------------
// return true if n is congruent to 3 mod 4
bool is_congruent_3m4(mpz_class n)
{
    return 0 != mpz_congruent_ui_p(n.get_mpz_t(), 3, 4);
}

// -----------------------------------------------------------------------------
// compute the square root of quadratic residue z mod p
bool qr_sqrt(mpz_class z, uint_t p, uint_t *r1, uint_t *r2)
{
    if (!is_congruent_3m4(p)) {
        cerr << p << " is not congruent to 3 mod 4" << endl;
        return false;
    }
    if (!is_quadratic_residue(z, p)) {
        cerr << z << " is not a quadratic residue mod " << p << endl;
        return false;
    }

    mpz_class r, pz(p), x, y;
    mpz_powm_ui(r.get_mpz_t(), z.get_mpz_t(), (p + 1) >> 2, pz.get_mpz_t());

    // compute r1 = +z^(p+1)/4 mod p
    mpz_mod_ui(x.get_mpz_t(), r.get_mpz_t(), p);
    *r1 = x.get_ui();

    // compute r2 = -z^(p+1)/4 mod p
    r = -r;
    mpz_mod_ui(y.get_mpz_t(), r.get_mpz_t(), p);
    *r2 = y.get_ui();

    return true;
}

// -----------------------------------------------------------------------------
// add two points P and Q on the specified elliptic curve mod p
bool elliptic_add(ipair P, ipair Q, ipair *R, const curve &ec)
{
    // check if P + Q = inf
    mpz_class ny2(Q.second); ny2 = -ny2;
    mpz_mod_ui(ny2.get_mpz_t(), ny2.get_mpz_t(), ec.p);
    if ((P.first == Q.first) && (P.second == ny2.get_ui()))
        return true;

    // compute lambda
    mpz_class x1(P.first), y1(P.second), x2(Q.first), y2(Q.second), pz(ec.p), L;
    mpz_class dx, dy;
    if ((P.first != Q.first) || (P.second != Q.second)) {
        mpz_sub(dy.get_mpz_t(), y2.get_mpz_t(), y1.get_mpz_t());
        mpz_sub(dx.get_mpz_t(), x2.get_mpz_t(), x1.get_mpz_t());
        mpz_invert(dx.get_mpz_t(), dx.get_mpz_t(), pz.get_mpz_t());
    }
    else {
        mpz_pow_ui(dx.get_mpz_t(), x1.get_mpz_t(), 2);
        mpz_mul_ui(dx.get_mpz_t(), dx.get_mpz_t(), 3);
        mpz_add_ui(dx.get_mpz_t(), dx.get_mpz_t(), ec.a);
        mpz_mul_ui(dy.get_mpz_t(), y1.get_mpz_t(), 2);
        mpz_invert(dy.get_mpz_t(), dy.get_mpz_t(), pz.get_mpz_t());
    }
    mpz_mul(L.get_mpz_t(), dx.get_mpz_t(), dy.get_mpz_t());
    mpz_mod_ui(L.get_mpz_t(), L.get_mpz_t(), ec.p);

    // compute (x3, y3)
    mpz_class L2, x3, y3;
    mpz_powm_ui(L2.get_mpz_t(), L.get_mpz_t(), 2, pz.get_mpz_t());
    mpz_sub(x3.get_mpz_t(), L2.get_mpz_t(), x1.get_mpz_t());
    mpz_sub(x3.get_mpz_t(), x3.get_mpz_t(), x2.get_mpz_t());
    mpz_mod_ui(x3.get_mpz_t(), x3.get_mpz_t(), ec.p);

    mpz_sub(y3.get_mpz_t(), x1.get_mpz_t(), x3.get_mpz_t());
    mpz_mul(y3.get_mpz_t(), y3.get_mpz_t(), L.get_mpz_t());
    mpz_sub(y3.get_mpz_t(), y3.get_mpz_t(), y1.get_mpz_t());
    mpz_mod_ui(y3.get_mpz_t(), y3.get_mpz_t(), ec.p);

    R->first = x3.get_ui();
    R->second = y3.get_ui();
    return false;
}

// -----------------------------------------------------------------------------
// evaulate the expression y^2 = x^3 + ax + b for a given x
mpz_class evaluate_curve(uint_t x, const curve &ec)
{
    mpz_class zx(x), x3, ax, y2;
    mpz_pow_ui(x3.get_mpz_t(), zx.get_mpz_t(), 3);
    mpz_mul_ui(ax.get_mpz_t(), zx.get_mpz_t(), ec.a);
    mpz_add(y2.get_mpz_t(), x3.get_mpz_t(), ax.get_mpz_t());
    mpz_add_ui(y2.get_mpz_t(), y2.get_mpz_t(), ec.b);
    mpz_mod_ui(y2.get_mpz_t(), y2.get_mpz_t(), ec.p);
    return y2;
}

// -----------------------------------------------------------------------------
// analyze points on the curve y^2 = x^3 + ax + b over Zp
void analyze_curve(const curve &ec)
{
    ivec points;
    for (uint_t x = 0; x < ec.p; ++x) {
        // compute y^2
        mpz_class y2 = evaluate_curve(x, ec);

        // determine if y2 is a quadratic residue
        cout << "x = " << x << ", y^2 = " << y2 << ", qr = ";
        if (is_quadratic_residue(y2, ec.p)) {
            // compute the roots of y^2
            uint_t r1, r2;
            qr_sqrt(y2, ec.p, &r1, &r2);
            points.push_back(make_pair(x, r1));
            points.push_back(make_pair(x, r2));
            cout << "yes, roots = (" << r1 << ", " << r2 << ")\n";
        }
        else {
            cout << "no\n";
        }
    }
    cout << "there are " << points.size() + 1 << " points on E\n";

    // compute the order of each element in E
    imap o;   
    for (ivec::const_iterator i = points.begin(); i != points.end(); ++i) {
        ipair P(*i), Q(*i), R;
        cout << "    (" << i->first << ", " << i->second << ") ";
        uint_t order = 0;
        for (uint_t j = 0; j <= 2 * ec.p; ++j) {
            bool is_inf = elliptic_add(P, Q, &R, ec);
#if 0
            cout << "    " << P << " + " << Q << " = ";
            if (is_inf) cout << "inf\n";
            else cout << R << endl;
#endif
            if (is_inf) {
                order = j + 2;
                break;
            }
            P = R;
        }
        cout << "order = " << order << endl;
        o[order]++;
    }

    cout << "printing the number of elements in each order...\n";
    for (imap::const_iterator i = o.begin(); i != o.end(); ++i) {
        cout << "    " << i->first << " -> " << i->second << endl;
    }
}

// -----------------------------------------------------------------------------
// perform point decompression to compute (x, y) = PD(x, i)
bool point_decompress(ipair *y1, const curve &ec)
{
    mpz_class y2 = evaluate_curve(y1->first, ec);
    if (!is_quadratic_residue(y2, ec.p))
        return false;

    // compute the roots of y^2
    uint_t r1, r2;
    if (!qr_sqrt(y2, ec.p, &r1, &r2))
        return false;

    if ((y1->second & 1) == (r1 & 1))
        y1->second = r1;
    else
        y1->second = r2;
    return true;
}

// -----------------------------------------------------------------------------
ipair point_pow(const ipair &pt, uint_t exp, const curve &ec, bool show)
{
    ipair t, r = pt;
    if (exp < 2) return r;

    for (uint_t i = 2; i <= exp; ++i) {
        if (elliptic_add(pt, r, &t, ec))
            cerr << "point_pow returned infinity...\n";
        if (show)
            cout << i << pt << " = " << pt << " + " << r << " = " << t << endl;
        r = t;
    }
    return r;
}

// -----------------------------------------------------------------------------
bool simplified_ecies(const ipair &y1, uint_t y2, uint_t k, const curve &ec)
{
    ipair y1_d = y1;
    if (!point_decompress(&y1_d, ec))
        return false;

    ipair p0 = point_pow(y1_d, k, ec, false);

    mpz_class dk, x0(p0.first);
    mpz_invert(x0.get_mpz_t(), x0.get_mpz_t(), mpz_class(ec.p).get_mpz_t());
    mpz_mul_ui(dk.get_mpz_t(), x0.get_mpz_t(), y2);
    mpz_mod_ui(dk.get_mpz_t(), dk.get_mpz_t(), ec.p);

    char c = (dk.get_ui() - 1) + 'A';
    cout << "(c, dk, y1_d, p0) = (" << c << ", " << dk << ", "
         << y1_d << ", " << p0 << ")\n";

    return true;
}

// -----------------------------------------------------------------------------
void ex_6_13(void)
{
    cout << "printing results for exercise 6.13 ...\n";

//  analyze_curve(curve(1, 6, 11));
    analyze_curve(curve(1, 28, 71));
//  analyze_curve(curve(2, 7, 31));
}

// -----------------------------------------------------------------------------
void ex_6_17(void)
{
    cout << "printing results for exercise 6.17 ...\n";

    curve ec(2, 7, 31);
    uint_t m = 8;

    ipair P(2, 9);
    ipair Q = point_pow(P, m, ec, true);
    cout << "Q = mP = " << Q << endl;

    simplified_ecies(ipair(18, 1), 21, m, ec);
    simplified_ecies(ipair( 3, 1), 18, m, ec);
    simplified_ecies(ipair(17, 0), 18, m, ec);
    simplified_ecies(ipair(28, 0),  8, m, ec);
}

// -----------------------------------------------------------------------------
// program entry-point
int main(int argc, char *argv[])
{
    ex_6_13();
    ex_6_17();

    return 0;
}

