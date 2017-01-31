#ifndef CCLONE_H
#define CCLONE_H
#include "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/siapopr/include/constantGlobalStructs.h"
#include "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/siapopr/include/constantRVFunctions.h"

#include <gsl/gsl_randist.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
void customclone(struct clone *new_clone, struct clone *parent_clone,
struct FitnessParameters* fit_params, struct MutationParameters* mut_params,
struct PunctuationParameters* punct_params,
struct EpistaticParameters* epi_params, int* number_mutations, gsl_rng* rng,
void (*ConstantGenerateFitness)(double *, struct FitnessParameters*, gsl_rng*));
#ifdef __cplusplus
}
#endif
#endif
