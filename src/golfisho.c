/*
 dddddsaa* golfisho.c
 *
 *  Created on: 22/02/2016
 *      Author: ernesto
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stddef.h>
#include <unistd.h>

#define CPX_VACIO &(cpx ) { 0 }

#define TAM_MAX_LINEA 6
#define MAX_NUM 200000

#define CACA_COMUN_ASSERT_DUROTE 0
#define CACA_COMUN_ASSERT_SUAVECITO 1
#define CACA_COMUN_ASSERT_NIMADRES 2

#define CACA_COMUN_TIPO_ASSERT CACA_COMUN_ASSERT_DUROTE

#if CACA_COMUN_TIPO_ASSERT == CACA_COMUN_ASSERT_DUROTE
#define assert_timeout(condition) assert(condition);
#endif
#if CACA_COMUN_TIPO_ASSERT == CACA_COMUN_ASSERT_SUAVECITO
#define assert_timeout(condition) if(!(condition)){printf("fuck\n");sleep(10);}
#endif
#if CACA_COMUN_TIPO_ASSERT == CACA_COMUN_ASSERT_NIMADRES
#define assert_timeout(condition) 0
#endif

/*
 #define caca_log_debug(formato, args...) 0
 */
#define caca_log_debug printf

#define caca_comun_max(x,y) ((x) < (y) ? (y) : (x))
#define caca_comun_min(x,y) ((x) < (y) ? (x) : (y))

typedef unsigned int tipo_dato;

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
void cpx_fft(cpx *in, cpx *out, int step, int size, int dir) {
	if (size < 1)
		return;
	if (size == 1) {
		out[0] = in[0];
		return;
	}
	cpx_fft(in, out, step * 2, size / 2, dir);
	cpx_fft(in + step, out + size / 2, step * 2, size / 2, dir);
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

void cpx_multiplicacion_polinomio(tipo_dato *coeficientes_a,
		tipo_dato *coeficientes_b, tipo_dato *coeficientes_resultado,
		int num_coef_a, int num_coef_b, int num_coef_res) {

	cpx *coef_a_complejos = NULL;
	cpx *coef_b_complejos = NULL;
	cpx *coef_a_complejos_fft = NULL;
	cpx *coef_b_complejos_fft = NULL;
	cpx *coef_res_complejos = NULL;
	cpx *producto_ab = NULL;

	two_pi = 4 * acos(0);

	assert_timeout(num_coef_res >= num_coef_a + num_coef_b - 1);

	coef_a_complejos = calloc(num_coef_a, sizeof(cpx));
	coef_b_complejos = calloc(num_coef_b, sizeof(cpx));
	coef_a_complejos_fft = calloc(num_coef_res, sizeof(cpx));
	coef_b_complejos_fft = calloc(num_coef_res, sizeof(cpx));
	coef_res_complejos = calloc(num_coef_res, sizeof(cpx));
	producto_ab = calloc(num_coef_res, sizeof(cpx));

	for (int i = 0; i < num_coef_a; i++) {
		coef_a_complejos[i].a = coeficientes_a[i];
	}
	for (int i = 0; i < num_coef_b; i++) {
		coef_b_complejos[i].a = coeficientes_b[i];
	}

	cpx_fft(coef_a_complejos, coef_a_complejos_fft, 1, num_coef_res, 1);
	cpx_fft(coef_b_complejos, coef_b_complejos_fft, 1, num_coef_res, 1);

	printf("transf a y b con fft\n");
	for (int i = 0; i < num_coef_res; i++) {
		printf("%7.2lf%7.2lf", coef_a_complejos_fft[i].a,
				coef_a_complejos_fft[i].b);
	}
	printf("\n");
	printf("transf a y b al chingadazo \n");
	for (int i = 0; i < num_coef_res; i++) {
		cpx *Ai = &(cpx ) { 0 };
		for (int j = 0; j < num_coef_res; j++) {
			cpx_sum(Ai, Ai,
					cpx_mult(CPX_VACIO, coef_a_complejos + j,
							EXP(CPX_VACIO, j * i * two_pi / num_coef_res)));
		}
		printf("%7.2lf%7.2lf", Ai->a, Ai->b);
	}
	printf("\n");

	for (int i = 0; i < num_coef_res; i++) {
		cpx_mult(producto_ab + i, coef_a_complejos_fft + i,
				coef_b_complejos_fft + i);
	}
	cpx_fft(producto_ab, coef_res_complejos, 1, num_coef_res, -1);
	for (int i = 0; i < num_coef_res; i++) {
		cpx_div(coef_res_complejos + i, coef_res_complejos + i, &(cpx ) {
						num_coef_res, 0 });
	}
	printf("con fft\n");
	for (int i = 0; i < num_coef_res; i++) {
		printf("%7.2lf%7.2lf", coef_res_complejos[i].a,
				coef_res_complejos[i].b);
	}
	printf("\n");
	printf("al chingadazo \n");
	for (int i = 0; i < num_coef_res; i++) {
		cpx *aconvbi = &(cpx ) { 0 };
		for (int j = 0; j < num_coef_res; j++) {
			cpx_sum(aconvbi, aconvbi,
					cpx_mult(&(cpx ) { 0 }, coef_a_complejos + j,
							coef_b_complejos
									+ ((num_coef_res + i - j) % num_coef_res)));
		}
		printf("%7.2lf%7.2lf", aconvbi->a, aconvbi->b);
	}
	printf("\n");

	for (int i = 0; i < num_coef_res; i++) {
		coeficientes_resultado[i] = coef_res_complejos[i].a;
	}

	free(coef_a_complejos);
	free(coef_b_complejos);
	free(coef_a_complejos_fft);
	free(coef_b_complejos_fft);
	free(coef_res_complejos);
	free(producto_ab);
}

static inline char *caca_arreglo_a_cadena(tipo_dato *arreglo, int tam_arreglo,
		char *buffer) {
	int i;
	char *ap_buffer = NULL;
	int characteres_escritos = 0;

	memset(buffer, 0, 100);
	ap_buffer = buffer;

	for (i = 0; i < tam_arreglo; i++) {
		characteres_escritos += sprintf(ap_buffer + characteres_escritos, "%d",
				*(arreglo + i));
		if (i < tam_arreglo - 1) {
			*(ap_buffer + characteres_escritos++) = ',';
		}
	}
	*(ap_buffer + characteres_escritos) = '\0';
	return ap_buffer;
}

static inline int lee_matrix_long_stdin(tipo_dato *matrix, int *num_filas,
		int *num_columnas, int num_max_filas, int num_max_columnas) {
	int indice_filas = 0;
	int indice_columnas = 0;
	long numero = 0;
	char *siguiente_cadena_numero = NULL;
	char *cadena_numero_actual = NULL;
	char *linea = NULL;

	linea = calloc(TAM_MAX_LINEA, sizeof(char));

	while (indice_filas < num_max_filas && fgets(linea, TAM_MAX_LINEA, stdin)) {
		indice_columnas = 0;
		cadena_numero_actual = linea;
		for (siguiente_cadena_numero = linea;; siguiente_cadena_numero =
				cadena_numero_actual) {
			numero = strtol(siguiente_cadena_numero, &cadena_numero_actual, 10);
			if (cadena_numero_actual == siguiente_cadena_numero) {
				break;
			}
			*(matrix + indice_filas * num_max_columnas + indice_columnas) =
					numero;
			caca_log_debug("en col %d, fil %d, el valor %lu\n", indice_columnas,
					indice_filas, numero);
			indice_columnas++;
			caca_log_debug("las columnas son %d\n", indice_columnas);
		}
		if (num_columnas) {
			num_columnas[indice_filas] = indice_columnas;
		}
		indice_filas++;
		caca_log_debug("las filas son %d, con clos %d\n", indice_filas,
				indice_columnas);
	}

	*num_filas = indice_filas;
	free(linea);
	return 0;
}

int main() {
	tipo_dato a[8] = { 0, 1, 3, 5, 1, 0, 2, 0 };
	tipo_dato b[8] = { 1, 2, 1, 3, 1, 3, 1, 2 };

	tipo_dato c[16] = { 0 };
	char buffer[100] = { '\0' };

	cpx_multiplicacion_polinomio(a, b, c, 8, 8, 16);

	caca_log_debug("el resultado entero %s\n",
			caca_arreglo_a_cadena(c, 16, buffer));

}
