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
#ifndef __SIMPLEGLOBALSTRUCTS_H_INCLUDED__
#define __SIMPLEGLOBALSTRUCTS_H_INCLUDED__

#include <string>

struct GlobalParameters
{
  double tot_life;
  int ancestors;
  int ancestor_clones;
  int num_sims;
  int num_samples;
  int sample_size;
  double detection_threshold;
  bool allow_extinction;
  double birth_rate;
  double death_rate;
  gsl_rng* rng;
  double seed;
};


// Clone Structure and Linked List Class containing list of clones
struct clone
{
  std::string clone_id;
  int cell_count;
  double birth_rate;
  double death_rate;

  // linking other nodes
  struct clone *nextnode;
  struct clone *prevnode;
};

#endif // __SIMPLEGLOBALSTRUCTS_H_INCLUDED__
