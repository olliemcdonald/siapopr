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


extern GlobalParameters gptime;

double TDGenerateFitness(FitnessParameters fit_params);
double TDGenerateMutationProb(MutationParameters mut_params);
int TDGeneratePunctuation(PunctuationParameters punct_params);

#endif // __TIMEDEPRVFUNCTIONS_TD_H_INCLUDED__
