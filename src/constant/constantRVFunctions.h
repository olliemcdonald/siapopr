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
#ifndef __CONSTANTRVFUNCTIONS_H_INCLUDED__
#define __CONSTANTRVFUNCTIONS_H_INCLUDED__

#include <gsl/gsl_randist.h>
#include <cmath>
#include "constantGlobalStructs.h"


extern GlobalParameters gpcons;

double ConstantGenerateFitness(FitnessParameters fit_params);
double ConstantGenerateMutationProb(MutationParameters mut_params);
int ConstantGeneratePunctuation(PunctuationParameters punct_params);

#endif // __CONSTANTRVFUNCTIONS_H_INCLUDED__
