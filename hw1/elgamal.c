// -----------------------------------------------------------------------------
// File:    elgamal.c
// Author:  Garrett Smith
// Created: 03/20/2011
// Description: 
//  Given a table of ciphertexts (y1, y2), p, alpha, a, and beta, decrypt the
//  message and recover the random values of k generated to encrypt the message.
//  Built and tested with gcc-4.5.2 and gmp-5.0.1.
// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gmp.h>

#define ENABLE_TESTS
#define init_pair(p, a, b) do { p.y1 = a; p.y2 = b; } while (0)
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

typedef struct pair
{
    long y1;
    long y2;
} pair_t;

// -----------------------------------------------------------------------------
void translate_block(char *msg, unsigned int block)
{
    int i;
    for (i = 0; i < 3; i++) {
        msg[2 - i] = 'A' + (block % 26);
        block /= 26;
    }
}

// -----------------------------------------------------------------------------
int dl_bruteforce(int n, int alpha, int beta, unsigned int *pa)
{
    int i, imax;
    mpz_t ip, ia, ib, temp;

    mpz_init(temp);
    mpz_init_set_ui(ip, n + 1);
    mpz_init_set_ui(ia, alpha);
    mpz_init_set_ui(ib, beta);

    // iterate over [0, min(beta, p)] checking for alpha^a [=] beta (mod p)
    imax = max(beta, (n + 1));
    for (i = 0; i <= imax; ++i) {
        mpz_pow_ui(temp, ia, i);
        if (0 != mpz_congruent_p(temp, ib, ip)) {
            *pa = i;
            return 0;
        }
    }

    // return a negative result if the computation failed
    return -1;
}

// -----------------------------------------------------------------------------
int list_cmp(const void *a, const void *b)
{
    return (*(pair_t *)a).y2 - (*(pair_t *)b).y2;
}

// -----------------------------------------------------------------------------
int dl_shanks(int n, int alpha, int beta, unsigned int *pa)
{
    mpz_t in, ip, ia, ib, im, iam, ik, temp;
    long i, j, m = ceil(sqrt(n));
    pair_t *L1, *L2;

    if (NULL == pa)
        return -1;

    L1 = malloc(sizeof(pair_t) * m);
    L2 = malloc(sizeof(pair_t) * m);

    mpz_inits(iam, ik, temp, 0);
    mpz_init_set_ui(in, n);
    mpz_init_set_ui(ip, n + 1);
    mpz_init_set_ui(ia, alpha);
    mpz_init_set_ui(ib, beta);
    mpz_init_set_ui(im, m);
    mpz_powm(iam, ia, im, ip);

    for (j = 0; j < m; ++j) {
        mpz_powm_ui(ik, iam, j, ip);
        init_pair(L1[j], j, mpz_get_ui(ik));
    }

    for (i = 0; i < m; ++i) {
        mpz_pow_ui(ik, ia, i);
        mpz_invert(ik, ik, ip);
        mpz_mul(ik, ik, ib);
        mpz_mod(ik, ik, ip);
        init_pair(L2[i], i, mpz_get_ui(ik));
    }

    qsort(L1, (size_t)m, sizeof(pair_t), list_cmp);
    qsort(L2, (size_t)m, sizeof(pair_t), list_cmp);

    i = j = 0;
    while ((i < m) && (j < m)) {
        if (L1[i].y2 == L2[j].y2) {
            mpz_set_ui(temp, m * L1[i].y1 + L2[j].y1);
            mpz_mod(ik, temp, in);
            *pa = mpz_get_ui(ik);
            return 0;
        }
        else if (L1[i].y2 < L2[j].y2)
            ++i;
        else if (L1[i].y2 > L2[j].y2)
            ++j;
    }

    // return a negative result if the computation failed
    return -1;
}

// -----------------------------------------------------------------------------
long load_ciphertext(const char *path, pair_t **pct)
{
    long length = 128, i = 0, y1, y2;
    pair_t *ct;
    FILE *fp;

    // attempt to open the file containing the ciphertext table
    fp = fopen("cipher", "r");
    if (!fp) {
        fprintf(stderr, "failed to open file '%s'\n", path);
        return 0;
    }

    // allocate the ciphertext array and read in each 3-character block
    ct = malloc(sizeof(pair_t) * length);
    for (;;) {
        if (2 != fscanf(fp, "%ld %ld", &y1, &y2))
            break;
        else if (i >= length) {
            length <<= 1;
            ct = realloc(ct, sizeof(pair_t) * length);
        }
        init_pair(ct[i], y1, y2); ++i;
    }

    fclose(fp);
    *pct = ct;
    return i;
}

// -----------------------------------------------------------------------------
void decrypt_ciphertext(int p, int a, pair_t *ct, size_t length)
{
    size_t i;
    char buffer[12], *message;
    mpz_t ip, iy1, iy2, ix;

    message = malloc(length * 3 + 1);
    memset(buffer, '\0', 12);
    memset(message, '\0', length * 3 + 1);

    mpz_inits(ip, iy1, iy2, ix, 0);
    mpz_init_set_ui(ip, p);

    for (i = 0; i < length; ++i) {
        // compute x = y2 * (y1^a)^(-1) mod p
        mpz_set_ui(iy1, ct[i].y1);
        mpz_set_ui(iy2, ct[i].y2);
        mpz_pow_ui(ix, iy1, a);
        mpz_invert(ix, ix, ip);
        mpz_mul(ix, ix, iy2);
        mpz_mod(ix, ix, ip);

        // display the translated results
        translate_block(&buffer[0], ct[i].y1);
        translate_block(&buffer[4], ct[i].y2);
        translate_block(&buffer[8], mpz_get_ui(ix));
        printf("(%s, %s) -> %s\n", &buffer[0], &buffer[4], &buffer[8]);

        // create a buffer containing the entire decrypted message
        memcpy(&message[i * 3], &buffer[8], 3);
    }

    printf("message: %s\n", message);
}

// -----------------------------------------------------------------------------
void recover_k(int p, int alpha, int beta, pair_t *ct, size_t length)
{
    size_t i;
    unsigned int a = 0;

#ifdef ENABLE_TESTS
    unsigned int a2 = 0;
    mpz_t y1, y2, ia, ip;

    mpz_inits(y1, y2, 0);
    mpz_init_set_ui(ia, alpha);
    mpz_init_set_ui(ip, p);
#endif // ENABLE_TESTS

    for (i = 0; i < length; ++i) {
        if (0 > dl_shanks(p - 1, alpha, ct[i].y1, &a))
            printf("computation failed\n");
        else
            printf("k(%d) = %d\n", (int)i, a);

#ifdef ENABLE_TESTS
        mpz_powm_ui(y1, ia, a, ip);
        assert((int)mpz_get_ui(y1) == ct[i].y1);
        assert(0 <= dl_bruteforce(p - 1, alpha, ct[i].y1, &a2));
        assert(a == a2);
#endif // ENABLE_TESTS
    }
}

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    unsigned int length, p = 31847, alpha = 5, a = 7899, beta = 18074;
    pair_t *ct;

    length = load_ciphertext("cipher", &ct);
    decrypt_ciphertext(p, a, ct, length);
    recover_k(p, alpha, beta, ct, length);

    return 0;
}

