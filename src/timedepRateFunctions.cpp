/*
 * =====================================================================================
 *
 *       Filename:  ratefunctions.cpp
 *
 *    Description: Contains different time-dependent functions that may be
 *                 implemented to adjust rates over time.
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

#include "timedepRateFunctions.h"


double RateFunctions::constant(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;
  double additional_rate = params->additional_rate;
  double a = params->coefs[0];
  double y = additional_rate + a;

  return y;
}

double RateFunctions::linear(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;

  double additional_rate = params->additional_rate;
  double a = params->coefs[0];
  double b = params->coefs[1];
  double min = params->coefs[2];

  double y = (a + t * b);
  y = additional_rate + GSL_MAX(y, min);
  return y;
}

double RateFunctions::logistic(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;

  double additional_rate = params->additional_rate;
  double A = params->coefs[0];
  double K = params->coefs[1];
  double B = params->coefs[2];
  double nu = params->coefs[3];
  double Q = params->coefs[4];
  double M = params->coefs[5];

  double y = additional_rate +
    (A + (K - A) / pow(1 + Q * exp(-B * (t-M)), 1 / nu) );
  return y;
}


double RateFunctions::Gompertz(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;

  double additional_rate = params->additional_rate;
  double asymptote = params->coefs[0];
  double alpha = params->coefs[1];
  double beta = params->coefs[2];

  double y = additional_rate +
    (asymptote + beta * exp(- alpha * t));
  return y;
}


double RateFunctions::custom(double t, void *p)
{
  struct TimeDependentParameters *params = (struct TimeDependentParameters *)p;
  double additional_rate = params->additional_rate;

  double a = params->coefs[0];
  double b = params->coefs[1];
  double c = params->coefs[2];
  // currently a cosine function just for fun
  double y = additional_rate + a + b * cos(c * t);
  return y;
}



double MaximizeRate(gsl_function rate_function, double start_time, double end_time, int bins)
{
  double max = GSL_FN_EVAL(&(rate_function), start_time);
  double test_max;
  double delta_t = (end_time - start_time) / bins;

  if(delta_t > 0.1)
  {
    delta_t = 0.1;
    bins = ceil((end_time - start_time) / delta_t);
  }

  for(int step = 1; step <= bins; ++step)
  {
    test_max = GSL_FN_EVAL(&(rate_function), start_time + delta_t * step);
    max = fmax(max, test_max);
  }
  return max;
}
