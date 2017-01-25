/*
 * =====================================================================================
 *
 *       Filename:  rvfunctions.h
 *
 *    Description: Header for functions for generating random variables
 *
 *        Version:  1.0
 *        Created:  08/24/2016 16:50:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thomas McDonald (), mcdonald@jimmy.harvard.edu
 *   Organization:  DFCI
 *
 * =====================================================================================
 */

#ifndef __TIMEDEPRVFUNCTIONS_TD_H_INCLUDED__
#define __TIMEDEPRVFUNCTIONS_TD_H_INCLUDED__

#include <gsl/gsl_randist.h>
#include <cmath>
#include "timedepGlobalStructs.h"

double tdoubleexp(struct FitnessParameters fit_params, gsl_rng* rng);
double tnormal(struct FitnessParameters fit_params, gsl_rng* rng);
double tuniform(struct FitnessParameters fit_params, gsl_rng* rng);

double TDGenerateMutationProb(MutationParameters mut_params, gsl_rng* rng);
int TDGeneratePunctuation(PunctuationParameters punct_params, gsl_rng* rng);

#endif // __TIMEDEPRVFUNCTIONS_TD_H_INCLUDED__
