/*
 * =====================================================================================
 *
 *       Filename:  rvfunctions.cpp
 *
 *    Description: Contains functions for generating random variables
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

#include "constantRVFunctions.h"

// Functions for generating distributions
double ConstantGenerateFitness(FitnessParameters fit_params, gsl_rng* rng)
{
  double fitness;
  double z = gsl_ran_flat(rng, 0, 1);


  if( (z > fit_params.pass_prob) &&
    ((fit_params.beta_fitness == 0) || z <= fit_params.pass_prob + (1 - fit_params.pass_prob) / 2) )
  {
    fitness = gsl_ran_exponential(rng, 1 / fit_params.alpha_fitness);
    return fitness;

  }
  else if( (fit_params.alpha_fitness == 0) || (z > fit_params.pass_prob + (1 - fit_params.pass_prob) / 2) )
  {
    fitness = -1 * gsl_ran_exponential(rng, 1 / fit_params.beta_fitness);
    return fitness;

  }
  else
  {
    // If both values are 0 or z < fit_params.pass_prob then return 0.
    fitness = 0.0;
    return fitness;
  }
}

/*
// For creating a custom distribution for fitness
double ConstantGenerateFitness(FitnessParameters fit_params)
{
  double fitness;
  double z = gsl_ran_flat(rng, 0, 1);

  if( z < fit_params.pass_prob )
  {
    return 0;
  }
  else
  {
    fitness = "favorite distribution";
    return fitness;
  }
}
*/


// Function for generating the mutation probability from a beta distribution
double ConstantGenerateMutationProb(MutationParameters mut_params, gsl_rng* rng)
{
  // might want to customize later
  double mut_prob;

  // allow for passenger probability
  double z = gsl_ran_flat(rng, 0, 1);

  if(z < mut_params.pass_prob)
  {
    mut_prob = 0;
  }
  else
  {
    mut_prob = gsl_ran_beta(rng, mut_params.alpha_mutation, mut_params.beta_mutation);
  }

  return mut_prob;
}

// Generate zero-truncated poisson random variable
int ConstantGeneratePunctuation(PunctuationParameters punct_params, gsl_rng* rng)
{
  // Generates a zero-truncated Poisson
  double lambda = punct_params.poisson_param;
  double rand_pois = gsl_ran_flat(rng, 0, 1);
  rand_pois = -log(1 - rand_pois * (1 - exp(-lambda)));
  rand_pois = lambda - rand_pois;

  int zero_trun_pois = gsl_ran_poisson(rng, rand_pois) + 1;

  return zero_trun_pois;
}
