// File:    hw2_part1.cpp
// Author:  Garrett Smith
// Created: 04/01/2011

#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include <cmath>
#include "hw2_part1.h"

using namespace std;

// -----------------------------------------------------------------------------
// Render the provided data set in an aligned table of the specified width.
void tabulate(ostream &out, const ivec &e, int perline, int width)
{
    int n = 0;
    for (ivec::const_iterator i = e.begin(); i != e.end(); ++i) {
        if (!n) out << "    ";
        out << setw(width) << *i << ((++n == perline) ? "\n" : " ");
        n %= perline;
    }
    if (0 != (n % perline)) out << endl;
}

// -----------------------------------------------------------------------------
// Perform computations for Zp^k and display the results in tabular format.
void test_Zpk(const string &filename, int p, int k)
{
    ofstream fp(filename.c_str());
    int n = totient_prime(p, k);
    fp << "p = " << p << "\nk = " << k << "\nn = " << n << "\n\n";

    ivec orders, prims, powers;
    compute_Zpk(p, k, orders, prims, powers);
    int width = ceil(log10(pow(p, k)));

    fp << "printing order of each element (0 ... " << orders.size()-1 << "):\n";
    tabulate(fp, orders, 20, width);

    fp << "printing primitive elements:\n";
    tabulate(fp, prims, 20, width);

    int mp = prims[0];
    fp << "printing powers of primitive element " << mp << " (" << mp
       << "^0 ... " << mp << "^" << powers.size()-1 << "):\n";
    tabulate(fp, powers, 20, width);

    map<int, int> o;
    for (ivec::const_iterator i = orders.begin(); i != orders.end(); ++i) {
        o[*i]++;
    }

    fp << "printing number of elements in each order ...\n";
    for (map<int, int>::const_iterator i = o.begin(); i != o.end(); ++i) {
        fp << "    order[" << i->first << "] -> " << i->second << endl;
    }

    fp.close();
}

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    int a;
    mpz_class alpha, az, p, b;

    test_Zpk("table_7_1.out", 7, 1);
    test_Zpk("table_13_1.out", 13, 1);
    test_Zpk("table_29_1.out", 29, 1);
    test_Zpk("table_131_1.out", 131, 1);
    // test_Zpk("table_131_2.out", 131, 2);

    // test...
    pohlig_hellman_print(809, 808, 3, 525, 2, 3);
    pohlig_hellman_print(809, 808, 3, 525, 101, 1);

    // test...
    shanks_zp(809, 808, 3, 525, &a);
    cout << "log_3(525) mod 809 = " << a << "\n\n";

    // test the book example: n = p-1 = 28 = (2^7)(7^1), B = 18, alpha = 2
    shanks_zp(29, 28, 2, 18, &a);
    pohlig_hellman_print(29, 28, 2, 18, 2, 2);
    pohlig_hellman_print(29, 28, 2, 18, 7, 1);

    alpha = 2; az = a; p = 29; // after CRT we arrive at a=11
    mpz_powm(b.get_mpz_t(), alpha.get_mpz_t(), az.get_mpz_t(), p.get_mpz_t());
    cout << "a=" << az << " beta=" << b << "\n\n";

    // do the hw problem: n = 17030
    shanks_zp(17161, 17030, 2, 103, &a);
    pohlig_hellman_print(17161, 17030, 2, 103, 2, 1);
    pohlig_hellman_print(17161, 17030, 2, 103, 5, 1);
    pohlig_hellman_print(17161, 17030, 2, 103, 13, 1);
    pohlig_hellman_print(17161, 17030, 2, 103, 131, 1);

    alpha = 2; az = a; p = 17161;
    mpz_powm(b.get_mpz_t(), alpha.get_mpz_t(), az.get_mpz_t(), p.get_mpz_t());
    cout << "a=" << az << " beta=" << b << "\n\n";

    compute_GFp2();

    return 0;
}

