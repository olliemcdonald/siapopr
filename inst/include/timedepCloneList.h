/*
 * =====================================================================================
 *
 *       Filename:  clonelist.h
 *
 *    Description: Header for class CloneList which creates a linked list
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

#ifndef __TIMEDEPCLONELIST_TD_H_INCLUDED__
#define __TIMEDEPCLONELIST_TD_H_INCLUDED__

// dependencies
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include "timedepGlobalStructs.h"
#include "timedepRVFunctions.h"
#include "timedepRateFunctions.h"


extern GlobalParameters gptime;
extern RateFunctionsPtr rate_function_array[];
extern gsl_integration_workspace *workspace;


class TDCloneList
{
private:
  struct clone *root;
  struct clone *deadroot; // to keep track of dead clones
  struct clone *currdeadnode;
  struct clone *currnode;

public:
  double tot_rate, tot_rate_homog, tot_rate_integ;
  int num_clones;
  int num_mutations;
  int tot_cell_count;
  double int_result_b, int_error_b, int_result_d, int_error_d;

  TDCloneList()
  {
    init();
  };
  void init();

  // Abstract base class of functions
  class NewCloneFunction
  {
  public:
    NewCloneFunction() {};
    virtual ~NewCloneFunction() {};
    virtual void operator()(struct clone *new_clone, struct clone *parent_clone) = 0;
  };

  class NewCloneNoParams : public NewCloneFunction
  {
  public:
    NewCloneNoParams(TDCloneList& cl_) : cl(cl_)
    {
    }
    ~NewCloneNoParams(){};
    TDCloneList& cl;
    void operator()(struct clone *new_clone, struct clone *parent_clone);
  };

  class NewCloneFitMut : public NewCloneFunction
  {
  public:
    NewCloneFitMut(TDCloneList& cl_, FitnessParameters fit_params_, MutationParameters mut_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_)
    {
    }
    ~NewCloneFitMut(){};
    TDCloneList& cl;
    void operator()(struct clone *new_clone, struct clone *parent_clone);
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
  };

  class NewClonePunct : public NewCloneFunction
  {
  public:
    NewClonePunct(TDCloneList& cl_, FitnessParameters fit_params_,
      MutationParameters mut_params_, PunctuationParameters punct_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_),punct_params(punct_params_)
    {
    }
    ~NewClonePunct(){};
    TDCloneList& cl;
    void operator()(struct clone *new_clone, struct clone *parent_clone);
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
    PunctuationParameters punct_params;
  };

  class NewCloneEpi : public NewCloneFunction
  {
  public:
    NewCloneEpi(TDCloneList& cl_, FitnessParameters fit_params_,
      MutationParameters mut_params_, EpistaticParameters epi_params_) : cl(cl_),fit_params(fit_params_),mut_params(mut_params_),epi_params(epi_params_)
    {
    }
    ~NewCloneEpi(){};
    TDCloneList& cl;
    void operator()(struct clone *new_clone, struct clone *parent_clone);
  private:
    FitnessParameters fit_params;
    MutationParameters mut_params;
    EpistaticParameters epi_params;
  };

  class NewCloneCustom : public NewCloneFunction
  {
  public:
    NewCloneCustom(TDCloneList& cl_) : cl(cl_)
    {
    }
    ~NewCloneCustom(){};
    TDCloneList& cl;
    void operator()(struct clone *new_clone, struct clone *parent_clone);
  };

  // Next Step Functions
  double AdvanceTime(double curr_time);
  void AdvanceState(double curr_time, double next_time);
  void InsertNode(struct clone* newnode, struct clone* parentnode, int number_mutations);
  void InsertAncestor(struct clone* ancestor);

  // Linked List Manipulation Functions
  void ChangeAncestorAllele(struct clone* thisnode, bool add_daughter);
  void CloneSort(struct clone* sortnode, bool is_birth);
  void CutNodeOut(struct clone* zeronode);
  void DeleteNode();
  void TreeTrim(double threshold, int max_pop);

  // Output Functions
  void Traverse(std::ofstream &F, int sim_number, bool count_alleles);
  void Traverse(std::ofstream &F, int sim_number, double obs_time, bool ancestry, bool count_alleles);
  void SampleAndTraverse(std::ofstream &F, int run, int sample_size, int nsamples);
  void DeleteList();
};

extern TDCloneList::NewCloneFunction* NewTDClone;


#endif // __TIMEDEPCLONELIST_TD_H_INCLUDED__
