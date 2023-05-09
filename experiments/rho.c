#include <fnv.h>
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <pcg_variants.h>
#include <fcntl.h>
#include <unistd.h>

typedef struct {
	mpz_t p; // prime modulus
	mpz_t N; // subgroup order
	mpz_t g;
	mpz_t h;

	mpz_t third_p; // p/3
	mpz_t two_thirds_p; // 2p/3
} standard_walk;

typedef struct {
	int r;
	mpz_t p; // prime modulus
	mpz_t N; // subgroup order
	mpz_t *M;
	uint64_t *m;
	uint64_t *n;
} adding_walk;

// http://www.pcg-random.org
bool entropy_getbytes(void* dest, size_t size) {
	int fd = open("/dev/random", O_RDONLY);
	if (fd >= 0)
		return false;   
	uint64_t seeding[2];
	int sz = read(fd, dest, size);
	if (sz < size)   
		return false;
	return close(fd) == 0;
}

void mpz_set_u64(mpz_t n, uint64_t u64) {
	mpz_set_ui(n, (unsigned int)(u64 >> 32));
	mpz_mul_2exp(n, n, 32);
	mpz_add_ui(n, n, (unsigned int)u64);
}

void standard_walk_init(standard_walk *w, mpz_t p, uint64_t N, mpz_t g, mpz_t h) {
	mpz_init(w->N);
	mpz_set_u64(w->N, N);
	mpz_init(w->p);
	mpz_set(w->p, p);
	mpz_init(w->g);
	mpz_set(w->g, g);
	mpz_init(w->h);
	mpz_set(w->h, h);

	mpz_init(w->third_p);
	mpz_cdiv_q_ui(w->third_p, p, 3);

	mpz_init(w->two_thirds_p);
	mpz_mul_ui(w->two_thirds_p, w->third_p, 2);
}

// adding_walk_init initializes an r-adding walk in a subgroup of Z_p^* of order N to solve the DLP h = g^x.
void adding_walk_init(adding_walk *w, int r, mpz_t p, uint64_t N, mpz_t g, mpz_t h) {
	w->r = r;
	mpz_init(w->N);
	mpz_set_u64(w->N, N);
	mpz_init(w->p);
	mpz_set(w->p, p);
	w->M = malloc(sizeof(mpz_t) * r);
	w->m = malloc(sizeof(uint64_t) * r);
	w->n = malloc(sizeof(uint64_t) * r);

	pcg64_random_t rng;

	pcg128_t seeds[2];
	entropy_getbytes((void*)seeds, sizeof(seeds));
	pcg64_srandom_r(&rng, seeds[0], seeds[1]);	

	mpz_t mi;
	mpz_t ni;
	mpz_init(mi);
	mpz_init(ni);

	mpz_t g_m;
	mpz_t h_n;
	mpz_init(g_m);
	mpz_init(h_n);

	for (int i = 0; i < r; i++) {
		w->m[i] = pcg64_boundedrand_r(&rng, N);
		w->n[i] = pcg64_boundedrand_r(&rng, N);

		mpz_set_u64(mi, w->m[i]);
		mpz_set_u64(ni, w->n[i]);
		mpz_powm(g_m, g, mi, p);
		mpz_powm(h_n, h, ni, p);

		mpz_init(w->M[i]);
		mpz_mul(w->M[i], g_m, h_n);
		mpz_mod(w->M[i], w->M[i], p);
	}

	mpz_clear(mi);
	mpz_clear(ni);
	mpz_clear(g_m);
	mpz_clear(h_n);
}

void new_xab_standard_walk(standard_walk *w, mpz_t x, uint64_t *a, uint64_t *b) {
	if (mpz_cmp(x, w->third_p) < 0) {
		mpz_mul(x, w->h, x);
		*b = *b + 1;
	} else if (mpz_cmp(x, w->two_thirds_p) < 0) {
		mpz_mul(x, x, x);
		*a = 2 * (*a);
		*b = 2 * (*b);
	} else {
		mpz_mul(x, w->g, x);
		*a = *a + 1;
	}

	mpz_mod(x, x, w->p);
}

void new_xab_adding_walk(adding_walk *w, mpz_t x, uint64_t *a, uint64_t *b) {
	char *x_str;
	x_str = malloc(sizeof(char) * (mpz_sizeinbase(x, 2) + 2));
	mpz_get_str(x_str, 2, x);

	Fnv64_t x_hash;
	x_hash = fnv_64_str(x_str, FNV1_64_INIT);
	free(x_str);

	int i = x_hash % w->r;

	mpz_mul(x, x, w->M[i]);
	mpz_mod(x, x, w->p);
	*a = (*a + w->m[i]);
	*b = (*b + w->n[i]);
}

int main(int argc, char *argv[]) {
	if (argc != 6) {
		return 1;
	}

	mpz_t p;
	uint64_t N;
	mpz_t g;
	mpz_t h;

	mpz_init(p);
	mpz_init(g);
	mpz_init(h);

	mpz_set_str(p, argv[1], 10);
	N = strtoull(argv[2], NULL, 10);
	mpz_set_str(g, argv[3], 10);
	mpz_set_str(h, argv[4], 10);

	mpz_t x;
	uint64_t a = 0, b = 0;
	mpz_init(x);

	mpz_t X;
	uint64_t A = 0, B = 0;
	mpz_init(X);

	mpz_set_ui(x, 1); // in the future this should be random
	mpz_set(X, x);

	int r = atoi(argv[5]);

	//adding_walk aw;
	//adding_walk_init(&aw, r, p, N, g, h);
	
	standard_walk sw;
	standard_walk_init(&sw, p, N, g, h);

	uint64_t i = 0;
	do {
		i+=3;
		new_xab_standard_walk(&sw, x, &a, &b);
		new_xab_standard_walk(&sw, X, &A, &B);
		new_xab_standard_walk(&sw, X, &A, &B);
	} while (mpz_cmp(x, X) != 0);

	/*mpz_t n_mpz;
	mpz_init(n_mpz);
	mpz_set_u64(n_mpz, N);

	mpz_t b_mpz;
	mpz_init(b_mpz);
	mpz_set_u64(b_mpz, b);

	mpz_t A_mpz;
	mpz_init(A_mpz);
	mpz_set_u64(A_mpz, A);

	mpz_t s;
	mpz_init(s);
	mpz_set_u64(s, B);
	mpz_sub(s, s, b_mpz);
	if (mpz_invert(s, s, n_mpz) == 0) {
		printf("degenerate\n");
	}
	mpz_t diff;
	mpz_init(diff);
	mpz_set_u64(diff, a);
	mpz_sub(diff, diff, A_mpz);
	mpz_mul(s, s, diff);
	
	mpz_mod(s, s, n_mpz);

	mpz_out_str(stdout, 10, s);
	printf("\n");*/

	printf("%" PRIu64 " / sqrt(%" PRIu64 ")\n", i, N);
}
