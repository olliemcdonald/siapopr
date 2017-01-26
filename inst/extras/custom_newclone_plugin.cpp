void customclone(struct clone *new_clone, struct clone *parent_clone,
  struct FitnessParameters* fit_params_ptr, struct MutationParameters* mut_params_ptr,
  struct PunctuationParameters* punct_params_ptr,
  struct EpistaticParameters* epi_params_ptr, gsl_rng* rng,
  void (*ConstantGenerateFitness)(double *, struct FitnessParameters*, gsl_rng*))
{
  bool did_count_driver = false;
  // generation of additive rate to the fitness
    double additional_rate = 0;
    (*ConstantGenerateFitness)(&additional_rate, fit_params_ptr, rng);

    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      did_count_driver = true;
      new_clone->is_driver = true;
    }
    new_clone->birth_rate = fmax(0, additional_rate + parent_clone->birth_rate);

}
