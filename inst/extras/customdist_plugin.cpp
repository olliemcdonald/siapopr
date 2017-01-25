/* To create a shared library for use as a plugin,
    alter any lines and assign your desired additional fitness to
    (*fitness). Keep the function name as "customdist"
*/

// do not touch this line
void customdist(double* fitness, struct FitnessParameters *fit_params, gsl_rng* rng)
{

  double z = gsl_ran_flat(rng, 0, 1);

  if( (z > fit_params->pass_prob) )
  {
    (*fitness) = gsl_ran_flat(rng, fit_params->alpha_fitness, fit_params->beta_fitness);
//    (*fitness) = 1;
  }
}
