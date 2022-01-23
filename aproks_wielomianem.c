#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdlib.h>
#include <math.h>

#define POLYN_DEG 4

double fi(points_t *pts, double *a, int point){
    double value = a[0];
    int i;
    for(i=1; i<POLYN_DEG+1; i++)
        value+=a[i]*pow(pts->x[point], i);
    return value;
}

double dfi(points_t *pts, double *a, int point){
    double power, value = 0.0;
    int i;
    for(i=1; i<POLYN_DEG+1; i++){
        if(pts->x[point]<0.001 && pts->x[point]>-0.001) power = 0.0; /* obejście sytuacji 0 do potęgi 0, gdzie pow zwraca 1, a potrzebne nam 0 */
        else power = pow(pts->x[point], i-1);
        value+=a[i]*i*power;
    }
    return value;
}

double d2fi(points_t *pts, double *a, int point){
    double power, value = 0.0;
    int i;
    for(i=2; i<POLYN_DEG+1; i++){
        if(pts->x[point]<0.001 && pts->x[point]>-0.001) power = 0.0; /* obejście sytuacji 0 do potęgi 0, gdzie pow zwraca 1, a potrzebne nam 0 */
        else power = pow(pts->x[point], i-2);
        value+=a[i]*i*(i-1)*power;
    }
    return value;
}

double d3fi(points_t *pts, double *a, int point){
    double power, value = 0.0;
    int i;
    for(i=3; i<POLYN_DEG+1; i++){
        if(pts->x[point]<0.001 && pts->x[point]>-0.001) power = 0.0; /* obejście sytuacji 0 do potęgi 0, gdzie pow zwraca 1, a potrzebne nam 0 */
        else power = pow(pts->x[point], i-3);
        value+=a[i]*i*(i-1)*(i-2)*power;
    }
    return value;
}

void make_spl(points_t * pts, spline_t * spl){

    double *s = calloc(((POLYN_DEG * 2) + 1), sizeof *s);
    double *t = calloc((POLYN_DEG + 1), sizeof *t);
    int i, n;
    matrix_t *stMat = make_matrix(POLYN_DEG+1, POLYN_DEG+2); /* macierz współczynników s i t reprezentująca układ równań do wyznaczenia współczynników a0..am */

    /* obliczanie współczynników s0 .. s2m sumując wartości argumentów wszystkich punktów funkcji podniesionych do potęgi indeksu s */
    for(i=0; i<(POLYN_DEG * 2)+1; i++){
        for(n=0; n<pts->n; n++){
            s[i]+=pow(pts->x[n], (double)i);
        }
    }
    /* obliczanie współczynników t0 .. tm sumując wartości argumentów wszystkich punktów funkcji podniesionych do potęgi indeksu t, pomnożonych przez wartość funkcji dla tych argumentów */
    for(i=0; i<POLYN_DEG+1; i++){
        for(n=0; n<pts->n; n++){
            t[i]+=pow(pts->x[n], (double)i) * pts->y[n];
        }
    }
    /* konstrukcja macierzy do układu równań - wstawianie s */
    for(i=0; i<POLYN_DEG+1; i++){
        for(n=0; n<POLYN_DEG+1; n++){
            put_entry_matrix(stMat, i, n, s[i+n]);
        }
    }
    /* wstawianie t do ostatniej kolumny macierzy */
    for(i=0; i<POLYN_DEG+1; i++)
        put_entry_matrix(stMat, i, stMat->cn-1, t[i]);

    if(piv_ge_solver(stMat)!=0){
        spl->n = 0;
        return;
    }

    double *a = malloc((POLYN_DEG + 1) * sizeof *a); /* tablica współczynników a, dla wygody */
    for(i=0; i<POLYN_DEG+1; i++){
        a[i] = get_entry_matrix(stMat, i, stMat->cn-1);
    }
    alloc_spl(spl, pts->n-1);
    for(i=0; i<spl->n; i++){
        spl->x[i] = pts->x[i];
        spl->f[i] = fi(pts, a, i);
        spl->f1[i] = dfi(pts, a, i);
        spl->f2[i] = d2fi(pts, a, i);
        spl->f3[i] = d3fi(pts, a, i);
    }

    free_matrix(stMat);
    free(s);
    free(t);
    free(a);
}