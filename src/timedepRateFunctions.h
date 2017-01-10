/*
 * =====================================================================================
 *
 *       Filename:  rvfunctions.h
 *
 *    Description: header for functions for generating random variables
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

#ifndef __TIMEDEPRATEFUNCTIONS_TD_H_INCLUDED__
#define __TIMEDEPRATEFUNCTIONS_TD_H_INCLUDED__

#include "timedepGlobalStructs.h"
#include <gsl/gsl_math.h>


class RateFunctions
{
public:
  static double constant(double t, void *params);
  static double linear(double t, void *params);
  static double logistic(double t, void *params);
  static double Gompertz(double t, void *params);
  //double piecewise(double t, void *params);
  static double custom(double t, void *params);
};


double MaximizeRate(gsl_function rate_function, double start_time, double end_time, int bins = 1000);

#endif // __TIMEDEPRATEFUNCTIONS_TD_H_INCLUDED__
