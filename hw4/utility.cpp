#include <gmpxx.h>
#include <cassert>
#include "utility.h"

using namespace std;

// -----------------------------------------------------------------------------
point add(const point &a, const point &b, const curve &ec)
{
    if (a.infinity)
        return b;
    else if (b.infinity)
        return a;

    ipair P(a.x, a.y), Q(b.x, b.y), out;
    if (elliptic_add(P, Q, &out, ec))
        return point(); // infinity
    else
        return point(out.first, out.second);
}

// -----------------------------------------------------------------------------
point sub(const point &a, const point &b, const curve &ec)
{
    point bn(b.x, -b.y + ec.p, b.infinity);
    return add(a, bn, ec);
}

// -----------------------------------------------------------------------------
point pow(const point &r, uint_t exp, const curve &ec)
{
    point t(r.x, r.y);
    if (exp < 2) return t;
    for (uint_t i = 2; i <= exp; ++i)
        t = add(t, r, ec);
    return t;
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
int calc_order(const point &p, const curve &ec)
{
    int order = 0;
    ipair P(p.x, p.y), Q(p.x, p.y), R;

    for (unsigned int j = 0; j <= 2 * ec.p; ++j) {
        bool is_inf = elliptic_add(P, Q, &R, ec);
        if (is_inf) {
            order = j + 2;
            break;
        }
        P = R;
    }

    return order;
}

// -----------------------------------------------------------------------------
void double_and_add(point P, int n, const curve &ec)
{
    NAF c = int_to_naf(n);
    point Q;

    point check = pow(P, n, ec);
    cout << "expecting " << check << " (long method)\n";

    for (int i = c.size() - 1; i >= 0; --i) {
        Q = add(Q, Q, ec);
        cout << "  [c=" << c[i] << "] 2Q = " << Q << endl;
        if (1 == c[i]) {
            cout << "  [c=1]\tQ = Q + P = " << Q << " + " << P << " = ";
            Q = add(Q, P, ec);
            cout << Q << endl;
        }
        else if (-1 == c[i]) {
            cout << "  [c=-1] Q = Q - P = " << Q << " - " << P << " = ";
            Q = sub(Q, P, ec);
            cout << Q << endl;
        }
    }

    // validate...
    cout << "calculated " << n << P << " = " << Q << endl;
    if (check == Q)
        cout << "success!\n";
    else {
        cout << "failure?\n";
        assert(!"double_and_add failed!");
    }
}

// -----------------------------------------------------------------------------
NAF int_to_naf(int in)
{
    NAF out;
    while (in) {
        out.push_back(in & 1);
        in >>= 1;
    }
    out.push_back(0);

    unsigned int i = 0;
    int beg = -1, end = -1;
    for (;;) {
        if (1 == out[i]) {
            if (-1 == beg) beg = end = i;
            else end = i;
        }
        else {
            if ((-1 != beg) && (end > beg)) {
                out[beg] = -1;
                for (int j = beg + 1; j <= end; ++j) out[j] = 0;
                out[i] = 1;
                beg = end = i;
            }
            else beg = end = -1;
        }
        if (++i >= out.size())
            break;
    }
    return out;
}

// -----------------------------------------------------------------------------
void print_naf(const NAF &in)
{
    for (int i = in.size() - 1; i >= 0; --i)
        cout << in[i] << " ";
    cout << endl;
}

// -----------------------------------------------------------------------------
void ext_euclid(int a, int b, int *rr, int *rs, int *rt)
{
    int a0 = a, b0 = b, t0 = 0, t = 1, s0 = 1, s = 0, q = a0 / b0;
    int r = a0 - q * b0;
    while (r > 0) {
        int temp = t0 - q * t;
        t0 = t;
        t = temp;
        temp = s0 - q * s;
        s0 = s;
        s = temp;
        a0 = b0;
        b0 = r;
        q = a0 / b0;
        r = a0 - q * b0;
    }
    *rr = b0;
    *rs = s;
    *rt = t;
}

