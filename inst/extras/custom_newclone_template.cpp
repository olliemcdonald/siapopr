// Do not touch these first lines
void customclone(struct clone *new_clone, struct clone *parent_clone,
  struct FitnessParameters* fit_params_ptr, struct MutationParameters* mut_params_ptr,
  struct PunctuationParameters* punct_params_ptr,
  struct EpistaticParameters* epi_params_ptr, int number_mutations*, gsl_rng* rng,
  void (*ConstantGenerateFitness)(double *, struct FitnessParameters*, gsl_rng*))
{
  // additional code

  // To sample from a fitness distribution and assign the new fitness to
  // additional_rate (not required)
  double additional_rate = 0;
  (*ConstantGenerateFitness)(&additional_rate, fit_params_ptr, rng);

  //additional code

}
