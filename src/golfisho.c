/*
 dddddsaa* golfisho.c
 *
 *  Created on: 22/02/2016
 *      Author: ernesto
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>

#define CPX_VACIO &(cpx ) { 0 }

double two_pi;

typedef struct cpx {
	double a;
	double b;
} cpx;

cpx *cpx_init(cpx *caca, double aa, double bb) {
	caca->a = aa;
	caca->b = bb;
	return caca;
}

double cpx_modsq(double a, double b) {
	return a * a + b * b;
}

cpx *cpx_bar(cpx *caca, double a, double b) {
	return cpx_init(caca, a, -b);
}

cpx *cpx_sum(cpx *caca, cpx *a, cpx *b) {
	return cpx_init(caca, a->a + b->a, a->b + b->b);
}

cpx * cpx_mult(cpx *caca, cpx *a, cpx *b) {
	return cpx_init(caca, a->a * b->a - a->b * b->b, a->a * b->b + a->b * b->a);
}

cpx *cpx_div(cpx *caca, cpx *a, cpx *b) {
	cpx *r = cpx_mult(CPX_VACIO, a, cpx_bar(CPX_VACIO, b->a, b->b));
	return cpx_init(caca, r->a / cpx_modsq(b->a, b->b),
			r->b / cpx_modsq(b->a, b->b));
}

cpx *EXP(cpx *caca, double theta) {
	return cpx_init(caca, cos(theta), sin(theta));
}

// in:     input array
// out:    output array
// step:   {SET TO 1} (used internally)
// size:   length of the input/output {MUST BE A POWER OF 2}
// dir:    either plus or minus one (direction of the FFT)
// RESULT: out[k] = \sum_{j=0}^{size - 1} in[j] * exp(dir * 2pi * i * j * k / size)
void FFT(cpx *in, cpx *out, int step, int size, int dir) {
	if (size < 1)
		return;
	if (size == 1) {
		out[0] = in[0];
		return;
	}
	FFT(in, out, step * 2, size / 2, dir);
	FFT(in + step, out + size / 2, step * 2, size / 2, dir);
	for (int i = 0; i < size / 2; i++) {
		cpx even = out[i];
		cpx odd = out[i + (size / 2)];
		cpx_sum(out + i, &even,
				cpx_mult(CPX_VACIO, EXP(CPX_VACIO, dir * two_pi * i / size),
						&odd));
		cpx_sum(out + i + size / 2, &even,
				cpx_mult(CPX_VACIO,
						EXP(CPX_VACIO, dir * two_pi * (i + size / 2) / size),
						&odd));
	}
}

// Usage:
// f[0...N-1] and g[0..N-1] are numbers
// Want to compute the convolution h, defined by
// h[n] = sum of f[k]g[n-k] (k = 0, ..., N-1).
// Here, the index is cyclic; f[-1] = f[N-1], f[-2] = f[N-2], etc.
// Let F[0...N-1] be FFT(f), and similarly, define G and H.
// The convolution theorem says H[n] = F[n]G[n] (element-wise product).
// To compute h[] in O(N log N) time, do the following:
//   1. Compute F and G (pass dir = 1 as the argument).
//   2. Get H by element-wise multiplying F and G.
//   3. Get h by taking the inverse FFT (use dir = -1 as the argument)
//      and *dividing by N*. DO NOT FORGET THIS SCALING FACTOR.

int main(void) {
	printf("If rows come in identical pairs, then everything works.\n");

	cpx a[8] = { { 0 }, { 1, 0 }, { 1, 3 }, { 0, 5 }, { 1, 0 }, { 0 }, { 2, 0 },
			{ 0 } };
	cpx b[8] = { { 1 }, { 0, -2 }, { 0, 1 }, { 3, 0 }, { -1, 0 }, { -3, 0 }, {
			1, 0 }, { -2, 0 } };
	cpx A[8] = { 0 };
	cpx B[8] = { 0 };

	two_pi = 4 * acos(0);

	FFT(a, A, 1, 8, 1);
	FFT(b, B, 1, 8, 1);

	printf("transf a y b con fft\n");
	for (int i = 0; i < 8; i++) {
		printf("%7.2lf%7.2lf", A[i].a, A[i].b);
	}
	printf("\n");
	printf("transf a y b al chingadazo \n");
	for (int i = 0; i < 8; i++) {
		cpx *Ai = &(cpx ) { 0 };
		for (int j = 0; j < 8; j++) {
			cpx_sum(Ai, Ai,
					cpx_mult(CPX_VACIO, a + j,
							EXP(CPX_VACIO, j * i * two_pi / 8)));
		}
		printf("%7.2lf%7.2lf", Ai->a, Ai->b);
	}
	printf("\n");

	cpx AB[8] = { 0 };
	for (int i = 0; i < 8; i++) {
		cpx_mult(AB + i, A + i, B + i);
	}
	cpx aconvb[8] = { 0 };
	FFT(AB, aconvb, 1, 8, -1);
	for (int i = 0; i < 8; i++) {
		cpx_div(aconvb + i, aconvb + i, &(cpx ) { 8, 0 });
	}
	printf("con fft\n");
	for (int i = 0; i < 8; i++) {
		printf("%7.2lf%7.2lf", aconvb[i].a, aconvb[i].b);
	}
	printf("\n");
	printf("al chingadazo \n");
	for (int i = 0; i < 8; i++) {
		cpx *aconvbi = &(cpx ) { 0 };
		for (int j = 0; j < 8; j++) {
			cpx_sum(aconvbi, aconvbi,
					cpx_mult(&(cpx ) { 0 }, a + j, b + ((8 + i - j) % 8)));
		}
		printf("%7.2lf%7.2lf", aconvbi->a, aconvbi->b);
	}
	printf("\n");

	return 0;
}
