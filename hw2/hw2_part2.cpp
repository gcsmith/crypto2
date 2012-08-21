// File:    hw2_part2.cpp
// Author:  Garrett Smith
// Created: 04/04/2011

#include <string>
#include <iostream>
#include <sstream>
#include <NTL/ZZ_p.h>
#include <NTL/pair_ZZ_pX_long.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace std;

NTL_CLIENT

ZZ_pX pX_from_str(const string &str)
{
    int value, coeff = 0;
    stringstream ss;
    ZZ_pX f;
    ss << str;
    while (ss >> value)
        SetCoeff(f, coeff++, value);
    return f;
}

void print_factors(const string &str, int p)
{
    ZZ pz;
    ZZ_p::init(pz = p);

    vec_pair_ZZ_pX_long factors;
    CanZass(factors, pX_from_str(str));
    
    int num_vecs = factors.length();
    if (num_vecs < 2) {
        cout << "polynomial is irreducible in Z" << p << "[x]\n";
    }
    else {
        for (int i = 0; i < num_vecs; i++) {
            const ZZ_pX &p = factors[i].a;
            int degree = deg(p);
            cout << "(";
            for (int j = degree; j >= 0; --j) {
                ZZ_p x;
                GetCoeff(x, p, j);
                if (0 == x) continue;
                switch (j) {
                case 0: cout << x << ")"; break;
                case 1: cout << x << "x + "; break;
                default: cout << x << "x^" << j << " + "; break;
                }
            }
        }
        cout << endl;
    }
}

int main(int argc, char *argv[])
{
    cout << "Printing factorization of x^5 + x^4 + 1 ... in Z2[x]\n\t";
    print_factors("1 0 0 0 1 1", 2);

    cout << "Printing factorization of x^5 + x^3 + 1 ... in Z2[x]\n\t";
    print_factors("1 0 0 1 0 1", 2);

    cout << "Printing factorization of x^5 + x^4 + x^2 + 1 ... in Z2[x]\n\t";
    print_factors("1 0 1 0 1 1", 2);

    cout << "Printing factorization of x^2 + 1 ... in Z131[x]\n\t";
    print_factors("1 1", 131);

    return 0;
}

