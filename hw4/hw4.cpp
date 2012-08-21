#include <iostream>
#include <cstdlib>
#include <cassert>
#include <gmpxx.h>
#include "utility.h"

using namespace std;

// -----------------------------------------------------------------------------
void ex_6_18(void)
{
    cout << "================== EX 6.18 ==================" << endl;

//  NAF v1 = int_to_naf(3895);
//  print_naf(v1);
//  NAF v2 = int_to_naf(11);
//  print_naf(v2);
//
    cout << "NAF(" << 3895 << ") = " << int_to_naf(3895) << endl;
    cout << "NAF(" << 11 << ") = " << int_to_naf(11) << endl;
    cout << "NAF(" << 87 << ") = " << int_to_naf(87) << endl;

    curve ec(1, 26, 127);
    double_and_add(point(2, 6), 87, ec);
}

// -----------------------------------------------------------------------------
void ex_7_1(void)
{
    cout << "================== EX 7.1 ==================" << endl;

    int p = 31847;
    int alpha = 5;
    int beta = 25703;

    int gamma = 23972;
    int d1 = 31396, x1 = 8990;
    int d2 = 20481, x2 = 31415;

    mpz_class dd = d1 - d2;
    mpz_class dx = x1 - x2;
    mpz_class pm = p - 1;
    mpz_class d, xn, dn, pn, e, k, ak;

    d  = gcd(dd, pm);   // d  = gcd(d1 - d2, p - 1)
    xn = div(dx, d);    // x' = (x1 - x2) / d
    dn = div(dd, d);    // d' = (d1 - d2) / d
    pn = div(pm, d);    // p' = (p - 1) / d
    e  = inv(dn, pn);   // e  = (d')^-1 mod p'

    // compute k ...
    for (unsigned int i = 0; i < d.get_ui(); ++i) {
        k = mod((xn * e) + (pn * i), pm);
        ak = pow(alpha, k, p);
        cout << "testing i=" << i << " k=" << k << " ... ";
        if (con(gamma, ak, p)) {
            cout << "pass\n"; break;
        }
        else cout << "fail\n";
    }

    d = gcd(gamma, pm);
    gmp_printf("gcd(y, p-1) = gcd(%Zd, %Zd) = %Zd\n", mZ(gamma), Z(pm), Z(d));

    mpz_class A = div(gamma, d);
    mpz_class B = div(mod(x1 - (k * d1), pm), d);
    pm = div(pm, d);

    int r, s, t;
    ext_euclid(A.get_ui(), pm.get_ui(), &r, &s, &t);
    mpz_class a;

    // compute a, same method as above ...
    for (unsigned int i = 0; i < d.get_ui(); ++i) {
        a = mod(mod(s, pm) * B + pm * i, pm);
        ak = pow(alpha, a, p);
        if (con(beta, ak, p)) break;
    }

#if 0
    // oops, gcd(gamma, p-1) != 1 ...
    int a1 = compute_a(x1, k, d1, gamma, pm);
    int a2 = compute_a(x2, k, d2, gamma, pm);
#endif

    // double check value of a using brute force approach
    int ac = -1;
    for (int i = 0; i < p; i++) {
        mpz_class test = pow(alpha, i, p);
        if (con(beta, test, p)) {
            ac = i;
            break;
        }
    }
    cout << "k = " << k << "\na = " << a << endl;

    assert(ac >= 0);
    if (a == ac)
        cout << "success!\n";
    else {
        cout << "failure?\n";
        assert(!"ex_7_1 failed!");
    }
}

// -----------------------------------------------------------------------------
void ex_7_7(void)
{
    cout << "================== EX 7.7 ==================" << endl;

    mpz_class q = 101;
    mpz_class p = 7879;
    mpz_class alpha = 170;
    mpz_class a = 75;
    mpz_class beta = 4567;  // example 7.4
    mpz_class k = 49;
    mpz_class sha1 = 52;    // SHA-1(x) = 52

    // calculate the signature (y, d)
    mpz_class y;
    mpz_powm(Z(y), Z(alpha), Z(k), Z(p));
    mpz_mod(Z(y), Z(y), Z(q));

    mpz_class ki, d = sha1;
    mpz_addmul(Z(d), Z(a), Z(y));
    mpz_invert(Z(ki), Z(k), Z(q));
    mpz_mul(Z(d), Z(d), Z(ki));
    mpz_mod(Z(d), Z(d), Z(q));

    gmp_printf("sig(x, k) = (y, d) = (%Zd, %Zd)\n", Z(y), Z(d));

    // perform signature verification
    mpz_class e1 = mod(sha1 * inv(d, q), q);
    mpz_class e2 = mod(y * inv(d, q), q);

    mpz_class t1, t2, v;
    mpz_powm(Z(t1), Z(alpha), Z(e1), Z(p));
    mpz_powm(Z(t2), Z(beta), Z(e2), Z(p));
    mpz_mul(Z(v), Z(t1), Z(t2));
    mpz_mod(Z(v), Z(v), Z(p));
    mpz_mod(Z(v), Z(v), Z(q));

    gmp_printf("(e1, e2) = (%Zd, %Zd)\n", Z(e1), Z(e2));
    gmp_printf("y = %Zd = ((A^e1)(B^e2) mod p) mod q = %Zd\n", Z(y), Z(v));

    // validate ...
    if (y == v)
        cout << "success!\n";
    else {
        cout << "failure?\n";
        assert(!"ex_7_7 failed!");
    }
}

// -----------------------------------------------------------------------------
void compute_ecdsa(const curve &ec, const point &A,
                   mpz_class m, mpz_class sha1, mpz_class k)
{
    mpz_class p = ec.p;
    gmp_printf("computing ECDSA for m=%Zd sha1(x)=%Zd k=%Zd\n",
               Z(m), Z(sha1), Z(k));

    point B = pow(A, m.get_ui(), ec);
    cout << "  B = " << m << A << " = " << B << endl;
    mpz_class q = calc_order(A, ec);
    cout << "  q = order(A) = " << q << endl;

    point uv = pow(A, k.get_ui(), ec);
    mpz_class u = uv.x;
    mpz_class v = uv.y;
    gmp_printf("  (u, v) = (%Zd, %Zd)\n", Z(u), Z(v));

    mpz_class r = mod(u, q);
    mpz_class s = mod(inv(k, q) * (sha1 + m * r), q);
    gmp_printf("  (r, s) = (%Zd, %Zd)\n", Z(r), Z(s));

    cout << "performing signature verification..." << endl;

    mpz_class w = inv(s, q);
    mpz_class i = mod(w * sha1, q);
    mpz_class j = mod(w * r, q);
    gmp_printf("  (w, i, j) = (%Zd, %Zd, %Zd)\n", Z(w), Z(i), Z(j));

    point iA = pow(A, i.get_ui(), ec);
    point jB = pow(B, j.get_ui(), ec);
    uv = add(iA, jB, ec);
    u = uv.x;
    v = uv.y;
    gmp_printf("  (u, v) = (%Zd, %Zd)\n", Z(u), Z(v));

    mpz_class check = mod(u, q);
    gmp_printf("  r = %Zd = (u mod q) = %Zd\n", Z(r), Z(check));

    if (check == r)
        cout << "success!\n";
    else {
        cout << "failure?\n";
        assert(!"compute_ecdsa failed!");
    }
}

// -----------------------------------------------------------------------------
void ex_7_13(void)
{
    cout << "================== EX 7.13 ==================" << endl;

    compute_ecdsa(curve(1, 6, 11), point(2, 7), 7, 4, 3);
    compute_ecdsa(curve(1, 26, 127), point(2, 6), 54, 10, 75);
}

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    ex_6_18();
    ex_7_1();
    ex_7_7();
    ex_7_13();

    return 0;
}

