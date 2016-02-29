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

/*#define CPX_VACIO &(cpx ) { 0 }*/

#define TAM_MAX_LINEA 7
#define MAX_NUM 200000
#define MAX_LINEAS MAX_NUM
#define MAX_COEFICIENTES 524288

#define CACA_COMUN_ASSERT_DUROTE 0
#define CACA_COMUN_ASSERT_SUAVECITO 1
#define CACA_COMUN_ASSERT_NIMADRES 2

/*
 #define CACA_COMUN_TIPO_ASSERT CACA_COMUN_ASSERT_SUAVECITO
 */
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

#ifdef CACA_COMUN_LOG
#define caca_log_debug printf
#else
#define caca_log_debug(formato, args...) 0
#endif

#define caca_comun_max(x,y) ((x) < (y) ? (y) : (x))
#define caca_comun_min(x,y) ((x) < (y) ? (x) : (y))

typedef unsigned int tipo_dato;
typedef unsigned long long bitch_vector;

typedef enum BOOLEANOS {
	falso = 0, verdadero
} bool;

double two_pi;

typedef struct cpx {
	double a;
	double b;
} cpx;

static bool caca_comun_checa_bit(bitch_vector *bits, int posicion) {
	bool res = falso;
	int idx_arreglo = 0;
	int idx_registro = 0;

	idx_arreglo = posicion / 64;
	idx_registro = posicion % 64;

	res = !!(bits[idx_arreglo] & (bitch_vector) (1 << idx_registro));

	return res;
}

static void caca_comun_asigna_bit(bitch_vector *bits, int posicion) {
	int idx_arreglo = 0;
	int idx_registro = 0;

	idx_arreglo = posicion / 64;
	idx_registro = posicion % 64;

	bits[idx_arreglo] |= (bitch_vector) (1 << idx_registro);

}

cpx *cpx_init(cpx *caca, double aa, double bb) {
	caca->a = aa;
	caca->b = bb;
	return caca;
}

double cpx_modsq(double a, double b) {
	return a * a + b * b;
}

cpx *cpx_bar(cpx *caca, double a, double b) {
	return cpx_init(caca, a, b * -1);
}

cpx *cpx_sum(cpx *caca, cpx *a, cpx *b) {
	return cpx_init(caca, a->a + b->a, a->b + b->b);
}

cpx * cpx_mult(cpx *caca, cpx *a, cpx *b) {
	return cpx_init(caca, a->a * b->a - a->b * b->b, a->a * b->b + a->b * b->a);
}

cpx *cpx_div(cpx *caca, cpx *a, cpx *b) {
	cpx tmp1 = { 0 };
	cpx tmp2 = { 0 };
	cpx *r = cpx_mult(&tmp2, a, cpx_bar(&tmp1, b->a, b->b));
	return cpx_init(caca, r->a / cpx_modsq(b->a, b->b),
			r->b / cpx_modsq(b->a, b->b));
}

cpx *EXP(cpx *caca, double theta) {
	return cpx_init(caca, cos(theta), sin(theta));
}

void cpx_fft(cpx *in, cpx *out, int step, int size, int dir) {
	int i = 0;
	cpx tmp1 = { 0 };
	cpx tmp2 = { 0 };
	cpx tmp3 = { 0 };
	cpx tmp4 = { 0 };
	if (size < 1)
		return;
	if (size == 1) {
		out[0] = in[0];
		return;
	}
	cpx_fft(in, out, step * 2, size / 2, dir);
	cpx_fft(in + step, out + size / 2, step * 2, size / 2, dir);
	for (i = 0; i < size / 2; i++) {
		cpx even = out[i];
		cpx odd = out[i + (size / 2)];
		cpx_sum(out + i, &even,
				cpx_mult(&tmp1, EXP(&tmp2, dir * two_pi * i / size), &odd));
		cpx_sum(out + i + (size / 2), &even,
				cpx_mult(&tmp3,
						EXP(&tmp4, dir * two_pi * (i + size / 2) / size),
						&odd));
	}
}

void cpx_multiplicacion_polinomio(tipo_dato *coeficientes_a,
		tipo_dato *coeficientes_b, tipo_dato *coeficientes_resultado,
		int num_coef_a, int num_coef_b, int num_coef_res) {

	int i = 0;
	int j = 0;
	cpx tmp1 = { 0 };
	cpx tmp2 = { 0 };
	cpx tmp3 = { 0 };
	cpx tmp4 = { 0 };

	char *num_buf = NULL;

	cpx *coef_a_complejos = NULL;
	cpx *coef_b_complejos = NULL;
	cpx *coef_a_complejos_fft = NULL;
	cpx *coef_b_complejos_fft = NULL;
	cpx *coef_res_complejos = NULL;
	cpx *producto_ab = NULL;

	two_pi = 4 * acos(0);

	coef_a_complejos = calloc(num_coef_a, sizeof(cpx));
	coef_b_complejos = calloc(num_coef_b, sizeof(cpx));
	coef_a_complejos_fft = calloc(num_coef_res, sizeof(cpx));
	coef_b_complejos_fft = calloc(num_coef_res, sizeof(cpx));
	coef_res_complejos = calloc(num_coef_res, sizeof(cpx));
	producto_ab = calloc(num_coef_res, sizeof(cpx));

	for (i = 0; i < num_coef_a; i++) {
		coef_a_complejos[i].a = (double) coeficientes_a[i];
	}
	for (i = 0; i < num_coef_b; i++) {
		coef_b_complejos[i].a = (double) coeficientes_b[i];
	}
	for (i = 0; i < num_coef_a; i++) {
		caca_log_debug("mierda comple %f %f\n", coef_a_complejos[i].a,
				coef_b_complejos[i].b);
	}

	cpx_fft(coef_a_complejos, coef_a_complejos_fft, 1, num_coef_res, 1);
	cpx_fft(coef_b_complejos, coef_b_complejos_fft, 1, num_coef_res, 1);

	caca_log_debug("transf a y b con fft\n");
	caca_log_debug("coef a\n");
	for (i = 0; i < num_coef_res; i++) {
		caca_log_debug("%7.2lf%7.2lf", coef_a_complejos_fft[i].a,
				coef_a_complejos_fft[i].b);
	}
	caca_log_debug("\n");
	caca_log_debug("coef b\n");
	for (i = 0; i < num_coef_res; i++) {
		caca_log_debug("%7.2lf%7.2lf", coef_b_complejos_fft[i].a,
				coef_b_complejos_fft[i].b);
	}
	caca_log_debug("\n");
#ifndef ONLINE_JUDGE
	caca_log_debug("transf a y b al chingadazo \n");
	for (i = 0; i < num_coef_res; i++) {
		cpx Ai_e = { 0 };
		cpx *Ai = &Ai_e;
		for (j = 0; j < num_coef_a; j++) {
			cpx_sum(Ai, Ai,
					cpx_mult(&tmp1, coef_a_complejos + j,
							EXP(&tmp2, j * i * two_pi / num_coef_res)));
		}
		caca_log_debug("%7.2lf%7.2lf", Ai->a, Ai->b);
	}
	caca_log_debug("\n");
#endif

	for (i = 0; i < num_coef_res; i++) {
		cpx_mult(producto_ab + i, coef_a_complejos_fft + i,
				coef_b_complejos_fft + i);
	}

	caca_log_debug("producto ab\n");
	for (i = 0; i < num_coef_res; i++) {
		caca_log_debug("%7.2lf%7.2lf", producto_ab[i].a, producto_ab[i].b);
	}
	caca_log_debug("\n");

	cpx_fft(producto_ab, coef_res_complejos, 1, num_coef_res, -1);

	for (i = 0; i < num_coef_res; i++) {
		cpx tmp_a = { 0 };
		tmp_a.a = num_coef_res;
		cpx_div(coef_res_complejos + i, coef_res_complejos + i, &tmp_a);
	}
	caca_log_debug("con fft\n");
	for (i = 0; i < num_coef_res; i++) {
		caca_log_debug("%7.2lf%7.2lf", coef_res_complejos[i].a,
				coef_res_complejos[i].b);
	}
	caca_log_debug("\n");
#ifndef ONLINE_JUDGE
	caca_log_debug("al chingadazo \n");
	for (i = 0; i < num_coef_res; i++) {
		cpx aconvbi_e = { 0 };
		cpx *aconvbi = &aconvbi_e;
		for (j = 0; j < num_coef_a; j++) {
			cpx_sum(aconvbi, aconvbi,
					cpx_mult(&tmp3, coef_a_complejos + j,
							coef_b_complejos
									+ ((num_coef_a + i - j) % num_coef_res)));
		}
		caca_log_debug("%7.2lf%7.2lf", aconvbi->a, aconvbi->b);
	}
	caca_log_debug("\n");
#endif

	num_buf = calloc(1000, sizeof(char));

	for (i = 0; i < num_coef_res; i++) {
		sprintf(num_buf, "%6.0lf", coef_res_complejos[i].a);
		caca_log_debug("rep numerica de %7.10lf %s\n", coef_res_complejos[i].a,
				num_buf);
		coeficientes_resultado[i] = atoi(num_buf);
		/*		coeficientes_resultado[i] = trunc(coef_res_complejos[i].a);*/
		caca_log_debug("%7.10lf convertido a %u en %u\n",
				coef_res_complejos[i].a, coeficientes_resultado[i], i);
	}
	free(coef_a_complejos);
	free(coef_b_complejos);
	free(coef_a_complejos_fft);
	free(coef_b_complejos_fft);
	free(coef_res_complejos);
	free(producto_ab);
	free(num_buf);
}

static char *caca_arreglo_a_cadena(tipo_dato *arreglo, int tam_arreglo,
		char *buffer) {
	int i;
	char *ap_buffer = NULL;
	int characteres_escritos = 0;
#ifdef ONLINE_JUDGE
	return NULL;
#endif

	memset(buffer, 0, 100);
	ap_buffer = buffer;

	for (i = 0; i < tam_arreglo; i++) {
		characteres_escritos += sprintf(ap_buffer + characteres_escritos, "%2d",
				*(arreglo + i));
		if (i < tam_arreglo - 1) {
			*(ap_buffer + characteres_escritos++) = ',';
		}
	}
	*(ap_buffer + characteres_escritos) = '\0';
	return ap_buffer;
}

static char *caca_bitch_vector_a_cadena(bitch_vector *sektor, int tam_arreglo,
		char *buffer) {
	int i;
	char *ap_buffer = NULL;
	int characteres_escritos = 0;
#ifdef ONLINE_JUDGE
	return NULL;
#endif

	memset(buffer, 0, 100);
	ap_buffer = buffer;

	for (i = 0; i < tam_arreglo; i++) {
		if (caca_comun_checa_bit(sektor, i)) {
			characteres_escritos +=
					sprintf(ap_buffer + characteres_escritos, "%2d",
							i);
			if (i < tam_arreglo - 1) {
				*(ap_buffer + characteres_escritos++) = ',';
			}
		}
	}
	*(ap_buffer + characteres_escritos) = '\0';
	return ap_buffer;
}

void caca_comun_strreplace(char s[], char chr, char repl_chr) {
	int i = 0;
	while (s[i] != '\0') {
		if (s[i] == chr) {
			s[i] = repl_chr;
		}
		i++;
	}
}

static int lee_matrix_long_stdin(tipo_dato *matrix, int *num_filas,
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
		caca_comun_strreplace(linea, '\n', '\0');
		if (!strlen(linea)) {
			caca_log_debug("weird, linea vacia\n");
			continue;
		}
		for (siguiente_cadena_numero = linea;; siguiente_cadena_numero =
				cadena_numero_actual) {
			caca_log_debug("el numero raw %s\n", linea);
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

static void golfisho_main() {
	int i = 0;
	int j = 0;
	int k = 0;
	int num_filas = 0;
	int tam_coeficientes_perrisha = 0;
	int tam_coeficientes_resultado_redondeado = 0;
	int tam_coeficientes_resultado = 0;
	int conteo_oyos_factibles = 0;
	bitch_vector *distancias_oyos_presentes = NULL;
	bitch_vector *distancias_oyos_factibles = NULL;
	tipo_dato max_distancia_perrisha = 0;
	tipo_dato max_distancia_oio = 0;
	tipo_dato numero_distancias_perisha = 0;
	tipo_dato numero_distancias_oyos = 0;
	tipo_dato numero_distancias_oyos_unicos = 0;
	tipo_dato dist_act = 0;
	tipo_dato *matrix = NULL;
	tipo_dato *distancias_perrilla = NULL;
	tipo_dato *distancias_oyos = NULL;
	tipo_dato *distancias_oyos_unicas = NULL;
	tipo_dato *a = NULL;
	tipo_dato *b = NULL;
	tipo_dato *c = NULL;
	tipo_dato *c1 = NULL;
	char buffer[100] = { '\0' };
	char buffer1[100] = { '\0' };

	matrix = calloc(MAX_LINEAS * 2 + 2, sizeof(tipo_dato));
	assert_timeout(matrix);

	lee_matrix_long_stdin(matrix, &num_filas, NULL, MAX_LINEAS * 2 + 2, 1);

	numero_distancias_perisha = *matrix;
	distancias_perrilla = matrix + 1;
	numero_distancias_oyos = *(matrix + 1 + numero_distancias_perisha);
	distancias_oyos = matrix + 1 + numero_distancias_perisha + 1;

	caca_log_debug("el numero de distancias perrilla %d el de oyos %d\n",
			numero_distancias_perisha, numero_distancias_oyos);
	caca_log_debug("las dstancias perrilla %s de oyos %s\n",
			caca_arreglo_a_cadena(distancias_perrilla, numero_distancias_perisha, buffer),
			caca_arreglo_a_cadena(distancias_oyos, numero_distancias_oyos, buffer1));

	assert_timeout(numero_distancias_perisha);
	assert_timeout(numero_distancias_perisha<=MAX_NUM);
	assert_timeout(numero_distancias_oyos);
	assert_timeout(numero_distancias_oyos<=MAX_NUM);

	distancias_oyos_factibles = calloc(
			(MAX_NUM + 1) / (sizeof(bitch_vector) * 8) + 1,
			sizeof(bitch_vector));
	assert_timeout(distancias_oyos_factibles);

	distancias_oyos_presentes = calloc(
			(MAX_NUM / (sizeof(bitch_vector) * 8)) + 1, sizeof(bitch_vector));
	assert_timeout(distancias_oyos_presentes);

	distancias_oyos_unicas = calloc(MAX_NUM, sizeof(tipo_dato));
	assert_timeout(distancias_oyos_unicas);

	for (i = 0; i < numero_distancias_oyos; i++) {
		dist_act = distancias_oyos[i];
		assert_timeout(dist_act);
		if (!caca_comun_checa_bit(distancias_oyos_presentes, dist_act)) {
			distancias_oyos_unicas[numero_distancias_oyos_unicos++] = dist_act;
			caca_comun_asigna_bit(distancias_oyos_presentes, dist_act);
		}
		if (dist_act > max_distancia_oio) {
			max_distancia_oio = dist_act;
		}
	}

	caca_log_debug("distancias unicas %s\n",
			caca_arreglo_a_cadena(distancias_oyos_unicas,numero_distancias_oyos_unicos,buffer));

	/*
	 printf("el numero de distancias perrilla %hu el de oyos %hu\n",
	 numero_distancias_perisha, numero_distancias_oyos);
	 */

	for (i = 0; i < numero_distancias_perisha; i++) {
		dist_act = distancias_perrilla[i];
		assert_timeout(dist_act);
		if (dist_act > max_distancia_perrisha) {
			max_distancia_perrisha = dist_act;
		}
	}

	caca_log_debug("maxima distancia perrilla %u\n", max_distancia_perrisha);

	assert_timeout(max_distancia_perrisha<=MAX_NUM);
	assert_timeout(max_distancia_oio<=MAX_NUM);

	tam_coeficientes_perrisha = max_distancia_perrisha + 1;
	tam_coeficientes_resultado = tam_coeficientes_perrisha * 2;

	tam_coeficientes_resultado_redondeado = 1;
	while (tam_coeficientes_resultado_redondeado
			<= tam_coeficientes_resultado - 1) {
		tam_coeficientes_resultado_redondeado <<= 1;
	}
	caca_log_debug("tamano coef %u y redond %u\n", tam_coeficientes_resultado,
			tam_coeficientes_resultado_redondeado);

	assert_timeout(tam_coeficientes_resultado_redondeado<=MAX_COEFICIENTES);

	/*
	 printf("tamano coef %u y redond %u\n", tam_coeficientes_resultado,
	 tam_coeficientes_resultado_redondeado);
	 */

	a = calloc(tam_coeficientes_resultado_redondeado, sizeof(tipo_dato));
	assert_timeout(a);

	b = calloc(tam_coeficientes_resultado_redondeado, sizeof(tipo_dato));
	assert_timeout(b);

	for (i = 0; i < numero_distancias_perisha; i++) {
		tipo_dato distancia_act = 0;
		distancia_act = distancias_perrilla[i];
		caca_log_debug("prendiendo la caca %u \n", distancia_act);
		*(a + distancia_act) = 1;
		*(b + distancia_act) = distancia_act;
	}

	caca_log_debug("los coeficientes son %s\n",
			caca_arreglo_a_cadena(a, tam_coeficientes_resultado_redondeado, buffer));

#ifdef GOLFISHO_PUTEADO

	c = calloc(tam_coeficientes_resultado_redondeado, sizeof(tipo_dato));
	assert_timeout(c);
	c1 = calloc(tam_coeficientes_resultado_redondeado, sizeof(tipo_dato));
	assert_timeout(c1);

	cpx_multiplicacion_polinomio(a, b, c, tam_coeficientes_resultado_redondeado,
			tam_coeficientes_resultado_redondeado,
			tam_coeficientes_resultado_redondeado);

	caca_log_debug("el resultado entero %s\n",
			caca_arreglo_a_cadena(c, tam_coeficientes_resultado, buffer));

	cpx_multiplicacion_polinomio(a, a, c1,
			tam_coeficientes_resultado_redondeado,
			tam_coeficientes_resultado_redondeado,
			tam_coeficientes_resultado_redondeado);

	caca_log_debug("el resultado conteo %s\n",
			caca_arreglo_a_cadena(c1, tam_coeficientes_resultado, buffer));

	for (i = 0; i < tam_coeficientes_resultado_redondeado; i++) {
		tipo_dato dista_act = 0;
		tipo_dato conteo_act = 0;
		dista_act = c[i];
		conteo_act = c1[i];
		if (dista_act && i < MAX_NUM + 1) {
			if (conteo_act && dista_act) {
				caca_log_debug("por fft activando %u\n", i);
			}
			caca_comun_asigna_bit(distancias_oyos_factibles, i);
		}
	}
	caca_log_debug("las distancias factilbes de fft %s\n",
			caca_bitch_vector_a_cadena(distancias_oyos_factibles, (MAX_NUM + 1) / (sizeof(bitch_vector) * 8) + 1, buffer));

	for (i = 0; i < numero_distancias_perisha; i++) {
		tipo_dato dista_act = 0;
		dista_act = distancias_perrilla[i];
		caca_log_debug("por inputs activando %u\n", dista_act);
		caca_comun_asigna_bit(distancias_oyos_factibles, dista_act);
		if (dista_act * 2 < MAX_NUM + 1) {
			caca_log_debug("por dobles de inputs activando %u\n",
					dista_act * 2);
			caca_comun_asigna_bit(distancias_oyos_factibles, dista_act * 2);
		}
	}

	caca_log_debug("las distancias factilbes %s\n",
			caca_bitch_vector_a_cadena(distancias_oyos_factibles, (MAX_NUM + 1) / (sizeof(bitch_vector) * 8) + 1, buffer));

	for (i = 0; i < numero_distancias_oyos_unicos; i++) {
		tipo_dato distan_act = 0;
		distan_act = distancias_oyos_unicas[i];
		if (caca_comun_checa_bit(distancias_oyos_factibles, distan_act)) {
			caca_log_debug("el # %u se puede formar\n", distan_act);
			conteo_oyos_factibles++;
		}
	}

	caca_log_debug("el num d oios fact %u\n", conteo_oyos_factibles);
	printf("%u\n", conteo_oyos_factibles);

	free(c);
	free(c1);
	free(distancias_oyos_presentes);
#else
	for (i = 0; i < tam_coeficientes_resultado_redondeado; i++) {
		int acumulador_distancia = 0;
		for (j = 0; j <= i; j++) {
			int numero_arriba = 0;
			int numero_abajo = 0;

			numero_arriba = a[j];
			numero_abajo = a[i - j];

			caca_log_debug(
					"con tu carota de chiste arriba %u (%u) abajo %u (%u)\n", j,
					numero_arriba, i-j, numero_abajo);

			acumulador_distancia += numero_arriba * numero_abajo;
			caca_log_debug("acumulador distan de %u asta aora %u\n", i,
					acumulador_distancia);

		}
		if (acumulador_distancia) {
			caca_log_debug("distancia %u es factivle\n", i);
			caca_comun_asigna_bit(distancias_oyos_factibles, i);
		}
	}

	for (i = 0; i < numero_distancias_perisha; i++) {
		tipo_dato dista_act = 0;
		dista_act = distancias_perrilla[i];
		caca_log_debug("por inputs activando %u\n", dista_act);
		caca_comun_asigna_bit(distancias_oyos_factibles, dista_act);
		if (dista_act * 2 < MAX_NUM + 1) {
			caca_log_debug("por dobles de inputs activando %u\n",
					dista_act * 2);
			caca_comun_asigna_bit(distancias_oyos_factibles, dista_act * 2);
		}
	}

	caca_log_debug("las distancias factilbes %s\n",
			caca_bitch_vector_a_cadena(distancias_oyos_factibles, (MAX_NUM + 1) / (sizeof(bitch_vector) * 8) + 1, buffer));

	for (i = 0; i < numero_distancias_oyos_unicos; i++) {
		tipo_dato distan_act = 0;
		distan_act = distancias_oyos_unicas[i];
		if (caca_comun_checa_bit(distancias_oyos_factibles, distan_act)) {
			caca_log_debug("el # %u se puede formar\n", distan_act);
			conteo_oyos_factibles++;
		}
	}

	caca_log_debug("el num d oios fact %u\n", conteo_oyos_factibles);
	printf("%u\n", conteo_oyos_factibles);

#endif

	free(b);
	free(a);
	free(distancias_oyos_unicas);
	free(distancias_oyos_factibles);
	free(matrix);
}

int main(int argc, char *argv[]) {
	golfisho_main();
	return 0;
}
