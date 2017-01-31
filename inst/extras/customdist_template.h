#ifndef CDIST_H
#define CDIST_H

#include <gsl/gsl_randist.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
void customdist(double* fitness, struct FitnessParameters *fit_params, gsl_rng* rng);
#ifdef __cplusplus
}
#endif
#endif
