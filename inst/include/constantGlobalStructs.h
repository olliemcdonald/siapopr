/*
 * =====================================================================================
 *
 *       Filename:  globalstructs.h
 *
 *    Description:  Header for structures that are accessed by multiple classes
 *                  throughout the simulation
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
#ifndef __CONSTANTGLOBALSTRUCTS_H_INCLUDED__
#define __CONSTANTGLOBALSTRUCTS_H_INCLUDED__

#include <string>

struct GlobalParameters
{
  double tot_life;
  long unsigned int max_pop;
  int ancestors;
  int ancestor_clones;
  int num_sims;
  int num_samples;
  int sample_size;
  double detection_threshold;
  double observation_frequency;
  double start_time;
  bool allow_extinction;
  bool is_custom_model;
  double birth_rate;
  double death_rate;
  double mutation_prob;
  bool trace_ancestry;
  bool count_alleles;
  double seed;
};


// Clone Structure and Linked List Class containing list of clones
struct clone
{
  std::string clone_id;
  int subclone_count;
  long unsigned int cell_count;
  int allele_count;
  int mut_count;
  int driver_count;
  bool is_driver;
  double birth_rate;
  double death_rate;
  double mut_prob;
  double clone_time;

  // linking other nodes
  struct clone *parent;
  struct clone *nextnode;
  struct clone *prevnode;
};

struct FitnessParameters
{
  std::string fitness_distribution;
  bool is_randfitness;
  double alpha_fitness;
  double beta_fitness;
  double pass_prob;
  double upper_fitness;
  double lower_fitness;
};

struct MutationParameters
{
  bool is_mutator;
  double alpha_mutation;
  double beta_mutation;
  double pass_prob;
};

struct PunctuationParameters
{
    bool is_punctuated;
    double punctuated_prob;
    double decay_rate;
    double poisson_param;
    double punctuated_multiplier;
    double punctuated_advantageous_prob;
};

struct EpistaticParameters
{
    bool is_epistasis;
    int epistatic_mutation_thresh;
    double epistatic_multiplier;
};

#endif // __CONSTANTGLOBALSTRUCTS_H_INCLUDED__
