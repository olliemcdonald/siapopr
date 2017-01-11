/*
 * =====================================================================================
 *
 *       Filename:  clonelist.h
 *
 *    Description: Header for class SimpleCloneList which creates a linked list
 *                 of clone structures along with the methods that modify,
 *                 add, delete, and output clones.
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

#ifndef __SIMPLECLONELIST_H_INCLUDED__
#define __SIMPLECLONELIST_H_INCLUDED__

// dependencies
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include "simpleGlobalStructs.h"


extern GlobalParameters gpsimp;

class SimpleCloneList
{
private:
  struct clone *root;
  struct clone *currnode;

public:
  double tot_rate;
  int num_clones;
  int tot_cell_count;

  SimpleCloneList() { init(); };
  void init();


  void InsertAncestor(struct clone* ancestor);

  void DeleteNode();
  // Output Functions
  void Traverse(std::ofstream &F, int sim_number);
  void SampleAndTraverse(std::ofstream &F, int run, int sample_size, int nsamples);
  void DeleteList();
};

#endif // __SIMPLECLONELIST_H_INCLUDED__
