// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// siapop
int siapop(double tot_life, int max_pop, double start_time, int ancestors, int ancestor_clones, int num_sims, bool allow_extinction, int num_samples, int sample_size, double detection_threshold, double observation_frequency, SEXP observation_times, double birth_rate, double death_rate, double mutation_prob, SEXP fitness_distribution, SEXP custom_distribution_file, double alpha_fitness, double beta_fitness, double pass_prob, SEXP upper_fitness, SEXP lower_fitness, double alpha_mutation, double beta_mutation, bool trace_ancestry, bool count_alleles, double punctuated_prob, double poisson_param, double punctuated_multiplier, double punctuated_advantageous_prob, double epistatic_mutation_thresh, double epistatic_multiplier, SEXP seed, SEXP custom_model_file, SEXP input_file, SEXP output_dir, SEXP ancestor_file);
RcppExport SEXP _siapopr_siapop(SEXP tot_lifeSEXP, SEXP max_popSEXP, SEXP start_timeSEXP, SEXP ancestorsSEXP, SEXP ancestor_clonesSEXP, SEXP num_simsSEXP, SEXP allow_extinctionSEXP, SEXP num_samplesSEXP, SEXP sample_sizeSEXP, SEXP detection_thresholdSEXP, SEXP observation_frequencySEXP, SEXP observation_timesSEXP, SEXP birth_rateSEXP, SEXP death_rateSEXP, SEXP mutation_probSEXP, SEXP fitness_distributionSEXP, SEXP custom_distribution_fileSEXP, SEXP alpha_fitnessSEXP, SEXP beta_fitnessSEXP, SEXP pass_probSEXP, SEXP upper_fitnessSEXP, SEXP lower_fitnessSEXP, SEXP alpha_mutationSEXP, SEXP beta_mutationSEXP, SEXP trace_ancestrySEXP, SEXP count_allelesSEXP, SEXP punctuated_probSEXP, SEXP poisson_paramSEXP, SEXP punctuated_multiplierSEXP, SEXP punctuated_advantageous_probSEXP, SEXP epistatic_mutation_threshSEXP, SEXP epistatic_multiplierSEXP, SEXP seedSEXP, SEXP custom_model_fileSEXP, SEXP input_fileSEXP, SEXP output_dirSEXP, SEXP ancestor_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tot_life(tot_lifeSEXP);
    Rcpp::traits::input_parameter< int >::type max_pop(max_popSEXP);
    Rcpp::traits::input_parameter< double >::type start_time(start_timeSEXP);
    Rcpp::traits::input_parameter< int >::type ancestors(ancestorsSEXP);
    Rcpp::traits::input_parameter< int >::type ancestor_clones(ancestor_clonesSEXP);
    Rcpp::traits::input_parameter< int >::type num_sims(num_simsSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_extinction(allow_extinctionSEXP);
    Rcpp::traits::input_parameter< int >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type detection_threshold(detection_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type observation_frequency(observation_frequencySEXP);
    Rcpp::traits::input_parameter< SEXP >::type observation_times(observation_timesSEXP);
    Rcpp::traits::input_parameter< double >::type birth_rate(birth_rateSEXP);
    Rcpp::traits::input_parameter< double >::type death_rate(death_rateSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_prob(mutation_probSEXP);
    Rcpp::traits::input_parameter< SEXP >::type fitness_distribution(fitness_distributionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type custom_distribution_file(custom_distribution_fileSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_fitness(alpha_fitnessSEXP);
    Rcpp::traits::input_parameter< double >::type beta_fitness(beta_fitnessSEXP);
    Rcpp::traits::input_parameter< double >::type pass_prob(pass_probSEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_fitness(upper_fitnessSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_fitness(lower_fitnessSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_mutation(alpha_mutationSEXP);
    Rcpp::traits::input_parameter< double >::type beta_mutation(beta_mutationSEXP);
    Rcpp::traits::input_parameter< bool >::type trace_ancestry(trace_ancestrySEXP);
    Rcpp::traits::input_parameter< bool >::type count_alleles(count_allelesSEXP);
    Rcpp::traits::input_parameter< double >::type punctuated_prob(punctuated_probSEXP);
    Rcpp::traits::input_parameter< double >::type poisson_param(poisson_paramSEXP);
    Rcpp::traits::input_parameter< double >::type punctuated_multiplier(punctuated_multiplierSEXP);
    Rcpp::traits::input_parameter< double >::type punctuated_advantageous_prob(punctuated_advantageous_probSEXP);
    Rcpp::traits::input_parameter< double >::type epistatic_mutation_thresh(epistatic_mutation_threshSEXP);
    Rcpp::traits::input_parameter< double >::type epistatic_multiplier(epistatic_multiplierSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type custom_model_file(custom_model_fileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type output_dir(output_dirSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ancestor_file(ancestor_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(siapop(tot_life, max_pop, start_time, ancestors, ancestor_clones, num_sims, allow_extinction, num_samples, sample_size, detection_threshold, observation_frequency, observation_times, birth_rate, death_rate, mutation_prob, fitness_distribution, custom_distribution_file, alpha_fitness, beta_fitness, pass_prob, upper_fitness, lower_fitness, alpha_mutation, beta_mutation, trace_ancestry, count_alleles, punctuated_prob, poisson_param, punctuated_multiplier, punctuated_advantageous_prob, epistatic_mutation_thresh, epistatic_multiplier, seed, custom_model_file, input_file, output_dir, ancestor_file));
    return rcpp_result_gen;
END_RCPP
}
// siapopNoMut
int siapopNoMut(double tot_life, int ancestors, int ancestor_clones, int num_sims, int num_samples, int sample_size, bool allow_extinction, double detection_threshold, double birth_rate, double death_rate, SEXP seed, SEXP input_file, SEXP output_dir, SEXP ancestor_file);
RcppExport SEXP _siapopr_siapopNoMut(SEXP tot_lifeSEXP, SEXP ancestorsSEXP, SEXP ancestor_clonesSEXP, SEXP num_simsSEXP, SEXP num_samplesSEXP, SEXP sample_sizeSEXP, SEXP allow_extinctionSEXP, SEXP detection_thresholdSEXP, SEXP birth_rateSEXP, SEXP death_rateSEXP, SEXP seedSEXP, SEXP input_fileSEXP, SEXP output_dirSEXP, SEXP ancestor_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tot_life(tot_lifeSEXP);
    Rcpp::traits::input_parameter< int >::type ancestors(ancestorsSEXP);
    Rcpp::traits::input_parameter< int >::type ancestor_clones(ancestor_clonesSEXP);
    Rcpp::traits::input_parameter< int >::type num_sims(num_simsSEXP);
    Rcpp::traits::input_parameter< int >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_extinction(allow_extinctionSEXP);
    Rcpp::traits::input_parameter< double >::type detection_threshold(detection_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type birth_rate(birth_rateSEXP);
    Rcpp::traits::input_parameter< double >::type death_rate(death_rateSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type output_dir(output_dirSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ancestor_file(ancestor_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(siapopNoMut(tot_life, ancestors, ancestor_clones, num_sims, num_samples, sample_size, allow_extinction, detection_threshold, birth_rate, death_rate, seed, input_file, output_dir, ancestor_file));
    return rcpp_result_gen;
END_RCPP
}
// siapopTD
int siapopTD(double tot_life, int max_pop, double start_time, int ancestors, int ancestor_clones, int num_sims, bool allow_extinction, bool is_custom_model, int num_samples, int sample_size, double detection_threshold, double observation_frequency, SEXP observation_times, int birth_function, Rcpp::NumericVector birth_coefs, int death_function, Rcpp::NumericVector death_coefs, double mutation_prob, SEXP fitness_distribution, SEXP custom_distribution_file, double alpha_fitness, double beta_fitness, double pass_prob, SEXP upper_fitness, SEXP lower_fitness, double alpha_mutation, double beta_mutation, bool trace_ancestry, bool count_alleles, double punctuated_prob, double poisson_param, double punctuated_multiplier, double punctuated_advantageous_prob, double epistatic_mutation_thresh, double epistatic_multiplier, SEXP seed, SEXP input_file, SEXP output_dir, SEXP ancestor_file);
RcppExport SEXP _siapopr_siapopTD(SEXP tot_lifeSEXP, SEXP max_popSEXP, SEXP start_timeSEXP, SEXP ancestorsSEXP, SEXP ancestor_clonesSEXP, SEXP num_simsSEXP, SEXP allow_extinctionSEXP, SEXP is_custom_modelSEXP, SEXP num_samplesSEXP, SEXP sample_sizeSEXP, SEXP detection_thresholdSEXP, SEXP observation_frequencySEXP, SEXP observation_timesSEXP, SEXP birth_functionSEXP, SEXP birth_coefsSEXP, SEXP death_functionSEXP, SEXP death_coefsSEXP, SEXP mutation_probSEXP, SEXP fitness_distributionSEXP, SEXP custom_distribution_fileSEXP, SEXP alpha_fitnessSEXP, SEXP beta_fitnessSEXP, SEXP pass_probSEXP, SEXP upper_fitnessSEXP, SEXP lower_fitnessSEXP, SEXP alpha_mutationSEXP, SEXP beta_mutationSEXP, SEXP trace_ancestrySEXP, SEXP count_allelesSEXP, SEXP punctuated_probSEXP, SEXP poisson_paramSEXP, SEXP punctuated_multiplierSEXP, SEXP punctuated_advantageous_probSEXP, SEXP epistatic_mutation_threshSEXP, SEXP epistatic_multiplierSEXP, SEXP seedSEXP, SEXP input_fileSEXP, SEXP output_dirSEXP, SEXP ancestor_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tot_life(tot_lifeSEXP);
    Rcpp::traits::input_parameter< int >::type max_pop(max_popSEXP);
    Rcpp::traits::input_parameter< double >::type start_time(start_timeSEXP);
    Rcpp::traits::input_parameter< int >::type ancestors(ancestorsSEXP);
    Rcpp::traits::input_parameter< int >::type ancestor_clones(ancestor_clonesSEXP);
    Rcpp::traits::input_parameter< int >::type num_sims(num_simsSEXP);
    Rcpp::traits::input_parameter< bool >::type allow_extinction(allow_extinctionSEXP);
    Rcpp::traits::input_parameter< bool >::type is_custom_model(is_custom_modelSEXP);
    Rcpp::traits::input_parameter< int >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type detection_threshold(detection_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type observation_frequency(observation_frequencySEXP);
    Rcpp::traits::input_parameter< SEXP >::type observation_times(observation_timesSEXP);
    Rcpp::traits::input_parameter< int >::type birth_function(birth_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type birth_coefs(birth_coefsSEXP);
    Rcpp::traits::input_parameter< int >::type death_function(death_functionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type death_coefs(death_coefsSEXP);
    Rcpp::traits::input_parameter< double >::type mutation_prob(mutation_probSEXP);
    Rcpp::traits::input_parameter< SEXP >::type fitness_distribution(fitness_distributionSEXP);
    Rcpp::traits::input_parameter< SEXP >::type custom_distribution_file(custom_distribution_fileSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_fitness(alpha_fitnessSEXP);
    Rcpp::traits::input_parameter< double >::type beta_fitness(beta_fitnessSEXP);
    Rcpp::traits::input_parameter< double >::type pass_prob(pass_probSEXP);
    Rcpp::traits::input_parameter< SEXP >::type upper_fitness(upper_fitnessSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lower_fitness(lower_fitnessSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_mutation(alpha_mutationSEXP);
    Rcpp::traits::input_parameter< double >::type beta_mutation(beta_mutationSEXP);
    Rcpp::traits::input_parameter< bool >::type trace_ancestry(trace_ancestrySEXP);
    Rcpp::traits::input_parameter< bool >::type count_alleles(count_allelesSEXP);
    Rcpp::traits::input_parameter< double >::type punctuated_prob(punctuated_probSEXP);
    Rcpp::traits::input_parameter< double >::type poisson_param(poisson_paramSEXP);
    Rcpp::traits::input_parameter< double >::type punctuated_multiplier(punctuated_multiplierSEXP);
    Rcpp::traits::input_parameter< double >::type punctuated_advantageous_prob(punctuated_advantageous_probSEXP);
    Rcpp::traits::input_parameter< double >::type epistatic_mutation_thresh(epistatic_mutation_threshSEXP);
    Rcpp::traits::input_parameter< double >::type epistatic_multiplier(epistatic_multiplierSEXP);
    Rcpp::traits::input_parameter< SEXP >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type output_dir(output_dirSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ancestor_file(ancestor_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(siapopTD(tot_life, max_pop, start_time, ancestors, ancestor_clones, num_sims, allow_extinction, is_custom_model, num_samples, sample_size, detection_threshold, observation_frequency, observation_times, birth_function, birth_coefs, death_function, death_coefs, mutation_prob, fitness_distribution, custom_distribution_file, alpha_fitness, beta_fitness, pass_prob, upper_fitness, lower_fitness, alpha_mutation, beta_mutation, trace_ancestry, count_alleles, punctuated_prob, poisson_param, punctuated_multiplier, punctuated_advantageous_prob, epistatic_mutation_thresh, epistatic_multiplier, seed, input_file, output_dir, ancestor_file));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_siapopr_siapop", (DL_FUNC) &_siapopr_siapop, 37},
    {"_siapopr_siapopNoMut", (DL_FUNC) &_siapopr_siapopNoMut, 14},
    {"_siapopr_siapopTD", (DL_FUNC) &_siapopr_siapopTD, 39},
    {NULL, NULL, 0}
};

RcppExport void R_init_siapopr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
