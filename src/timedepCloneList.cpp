/*
 * =====================================================================================
 *
 *       Filename:  clonelist.cpp
 *
 *    Description:  Contains the class TDCloneList which creates a linked list
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
#include "timedepCloneList.h"

// Constructor
void TDCloneList::init()
{
    root = NULL;
    currnode = NULL;
    deadroot = NULL;
    currdeadnode = NULL;
    num_clones = 0; // number of types
    num_mutations = 0; // number of mutations (types and deadtypes)
    tot_rate = 0; // total rate (birth + death)
    tot_rate_homog = 0;
    tot_rate_integ = 0;
    tot_cell_count = 0; // total number of cells in all types
}

// InsertNode attaches the newly created node (in NewClone) to the end of the
// linked list. Number of mutations was included for the punctuated scenario
// to add to the ID.
void TDCloneList::InsertNode(struct clone* newnode, struct clone* parentnode, int number_mutations)
{
  num_clones++;

  if( root == NULL )
  {
    newnode->nextnode = NULL;
    newnode->parent = NULL;
    newnode->prevnode = NULL;

    root = newnode;
    currnode = newnode;
  }
  else
  {
    // Declare extra pieces to append to id
    std::string add_id;
    if(number_mutations == 1)
    {
      add_id = add_id + '>' + std::to_string(num_clones);
    }
    else
    {
      for(int mutation_counter = 1; mutation_counter <= number_mutations; mutation_counter ++)
      {
        add_id = add_id + '>' + std::to_string(num_clones) + '.' + std::to_string(mutation_counter);
      }
    }
    // append new id to list of ancestors
    newnode->clone_id = newnode->clone_id + add_id;
    newnode->nextnode = NULL;
    newnode->parent = parentnode;

    // move to end of list to attach newnode to the end
    while (currnode->nextnode != NULL)
    {
      currnode = currnode->nextnode;
    }

    newnode->prevnode = currnode;
    currnode->nextnode = newnode;
  }
}


/*
  Insert the ancestors into the linked list. If ancestors are defined
  in a separate file it parses the file and inserts. Otherwise it
  creates nodes based on the number of clones and individuals defined
  in input file.
*/
void TDCloneList::InsertAncestor(struct clone* ancestor)
{
  if( ancestor->clone_id.empty() )
  {
    num_clones++;
    ancestor->clone_id = std::to_string(num_clones) + ".a";
  } else
  {
    ancestor->clone_id = ancestor->clone_id + ".a";
  }

  ancestor->nextnode = NULL;
  ancestor->parent = NULL;

  if( root == NULL )
  {
    ancestor->prevnode = NULL;
    root = ancestor;
    currnode = ancestor;
  }
  else
  {
    while (currnode->nextnode != NULL)
    {
      currnode = currnode->nextnode;
    }

    ancestor->prevnode = currnode;
    currnode->nextnode = ancestor;
  }
}


/*
  ChangeAncestorAllele occurs after determining which event occurs. It travels
  through parent pointer sin the list to add or subtract a clone's
  numallele value for each ancestor. This generates a count for the number
  of individuals a particular allele is represented by (allele is the last
  term of unique_id). Implemented if count_alleles is TRUE as an input.
  add_daughter argument is based on whether a birth or death occurred.
*/
void TDCloneList::ChangeAncestorAllele(struct clone* thisnode, bool add_daughter)
{
  if(add_daughter) // if added a daughter add to each allele_count
  {
    while(thisnode != NULL)
    {
      // add 1 to all allele_count for each ancestor
      thisnode->allele_count++;
      thisnode = thisnode->parent;
    }
  }
  else // else subtract from each allele_count
  {
    while(thisnode != NULL)
    {
      // add 1 to all allele_count for each ancestor
      thisnode->allele_count--;
      thisnode = thisnode->parent;
    }
  }
}


/*
  Method to sort the current node based on the number of cells in it
  relative to its prevnode and nextnode. If death, move toward nextnode
  if birth, move to prevnode direction second argument is true if birth.
  Method should help speed up later runs when subclones have larger fitness
  and start to dominate the process.
*/
void TDCloneList::CloneSort(struct clone* sortnode, bool is_birth)
{
  if(is_birth && sortnode->prevnode != NULL)
  {
    if( sortnode->cell_count * (sortnode->birth_rate + sortnode->death_rate) >=
      sortnode->prevnode->cell_count * (sortnode->prevnode->birth_rate + sortnode->prevnode->death_rate) )
    {
      struct clone* newprevnode;
      struct clone* newnextnode;

      newnextnode = sortnode->prevnode;
      newprevnode = newnextnode->prevnode;
      // attach prevnode and nextnode
      if(sortnode->nextnode != NULL)
      {
        sortnode->nextnode->prevnode = sortnode->prevnode;
      }
      sortnode->prevnode->nextnode = sortnode->nextnode;
      // insert clone between newprevnode and newnextnode
      sortnode->prevnode = newprevnode;
      sortnode->nextnode = newnextnode;

      if(newprevnode != NULL)
      {
        newprevnode->nextnode = sortnode;
      } else
      {
        // reroot to the current node
        root = sortnode;
      }
      newnextnode->prevnode = sortnode;
    }
  } else if(!is_birth && sortnode->nextnode != NULL)
  {
    // move newprevnode to the right
    if( sortnode->cell_count * (sortnode->birth_rate + sortnode->death_rate) <=
      sortnode->nextnode->cell_count  * (sortnode->nextnode->birth_rate + sortnode->nextnode->death_rate) )
    {
      struct clone* newprevnode;
      struct clone* newnextnode;

      newprevnode = sortnode->nextnode;
      newnextnode = newprevnode->nextnode;
      // attach prevnode and nextnode
      if(sortnode->prevnode != NULL)
      {
        sortnode->prevnode->nextnode = sortnode->nextnode;
      } else
      {
        root = sortnode->nextnode;
      }
      sortnode->nextnode->prevnode = sortnode->prevnode;

      // insert clone between newprevnode and newnextnode
      sortnode->prevnode = newprevnode;
      sortnode->nextnode = newnextnode;
      if(newnextnode != NULL)
      {
        newnextnode->prevnode = sortnode;
      }
      newprevnode->nextnode = sortnode;
    }
  }
}

/*
  AdvanceTime determines the next time in the process. It does so using adaptive
  thinning, which generates exponential random variables with a rate greater
  than the value of rate function at time t and accepts/reject the r.v. with
  according to a probability based on the value of the rate
*/
double TDCloneList::AdvanceTime(double curr_time, gsl_rng* rng)
{
  struct clone *pnode;
  double rand_next_time;

  while(true) // loop until a r.v. is generated that is accepted
  {
    Rcpp::checkUserInterrupt();

    tot_rate = 0;
    rand_next_time = gsl_ran_exponential(rng, 1 / tot_rate_homog);
    pnode = root;

    while(pnode) // add the rates of all nodes at the next time to get total rate
    {
      //    time dependent rate to test
      tot_rate = tot_rate + (GSL_FN_EVAL(&(pnode->B), curr_time + rand_next_time) +
        GSL_FN_EVAL(&(pnode->D), curr_time + rand_next_time)) * pnode->cell_count;
      pnode = pnode->nextnode;
    }

    // if total rate is less than the total rate assuming constant rates, accept
    // the next time.
    double u_thin = gsl_ran_flat(rng, 0, 1);
    double beta_ratio = tot_rate / tot_rate_homog;
    //std::cout << beta_ratio << "\n";

    if(u_thin <= beta_ratio)
    {
      break;
    }
  }

  //std::cout << tot_rate << "\t" << curr_time + rand_next_time << "\t" << tot_cell_count << "\n";

  return rand_next_time;
}

/*
  Advance state determines the next event that occurs in the process and which
  clone it occurs in (birth w/o mutation, birth w/ mutation, death) and updates
  the process accordingly
*/
void TDCloneList::AdvanceState(double curr_time, double next_time, gsl_rng* rng)
{
  double summand = 0;
  tot_rate_integ = 0;
  bool flag = false;

  struct clone *pnode;
  pnode = root;

  // Collect total rates over time between last and next event
  while(pnode)
  {
    // integrals for determining next clone
    gsl_integration_qags(&(pnode->B), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(int_result_b), &(int_error_b));
    pnode->birth_rate = int_result_b;
    gptime.tot_error = gptime.tot_error + int_error_b;

    gsl_integration_qags(&(pnode->D), curr_time, curr_time + next_time,
      0, 1e-12, 1000, workspace, &(int_result_d), &(int_error_d));
    pnode->death_rate = int_result_d;
    gptime.tot_error = gptime.tot_error + int_error_d;


    tot_rate_integ = tot_rate_integ + (pnode->birth_rate + pnode->death_rate) * pnode->cell_count;
    pnode = pnode->nextnode;
  }

  double rand_next_event = gsl_ran_flat(rng, 0, tot_rate_integ);
  double rand_mut_occur;
  // Put pnode back to the root
  pnode = root;

  // cycle through nodes to determine the next event that occurs
  while ( (pnode) && !(flag) )
  {
    // Condition for new birth
    if ( rand_next_event <= summand + (pnode->cell_count) * (pnode->birth_rate) )
    {
      rand_mut_occur = gsl_ran_flat(rng, 0, 1);
      // Condition to determine if mutation occurs in new daughter
      if (rand_mut_occur <= pnode->mut_prob)
      {
        // Creation of new clone
        struct clone *new_mut_node;
        new_mut_node = new struct clone;
        new_mut_node->clone_id = pnode->clone_id;
        new_mut_node->driver_count = pnode->driver_count;
        new_mut_node->subclone_count = 0;
        new_mut_node->cell_count = 1;
        new_mut_node->allele_count = 0;
        new_mut_node->clone_time = curr_time + next_time;
        new_mut_node->mut_count = pnode->mut_count + 1;
        new_mut_node->is_driver = false;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;

        new_mut_node->mut_prob = pnode->mut_prob;
        new_mut_node->birth_params = pnode->birth_params;
        new_mut_node->death_params = pnode->death_params;
        (new_mut_node->B).function = rate_function_array[new_mut_node->birth_params.type];
        (new_mut_node->B).params = &(new_mut_node->birth_params);
        (new_mut_node->D).function = rate_function_array[new_mut_node->death_params.type];
        (new_mut_node->D).params = &(new_mut_node->death_params);

        // Function class to add any modifications to rates and insert clone
        (*NewTDClone)(new_mut_node, pnode);

        // Reset birth and death rates (recalc at beginning of next AdvanceState)
        new_mut_node->birth_rate = 0;
        new_mut_node->death_rate = 0;

        // Get max values of function to add to homogeneous rate
        new_mut_node->birth_params.homogeneous_rate = MaximizeRate((new_mut_node->B), curr_time, gptime.tot_life, 1000);
        //(new_mut_node->B).params = &(new_mut_node->birth_params);
        new_mut_node->death_params.homogeneous_rate = MaximizeRate((new_mut_node->D), curr_time, gptime.tot_life, 1000);
        //(new_mut_node->D).params = &(new_mut_node->death_params);

        tot_rate_homog = tot_rate_homog + new_mut_node->birth_params.homogeneous_rate +
          new_mut_node->death_params.homogeneous_rate;

        tot_cell_count++;
        num_clones++;
        num_mutations++;


        pnode = new_mut_node;
        // add allele count to all ancestors
        if(gptime.count_alleles)
        {
          ChangeAncestorAllele(pnode, true);
        }

      }
      else // no mutation, increment clone by 1
      {
        pnode->cell_count++;
        tot_cell_count++;
        // Add to homogeneous process that undergoes thinning
        tot_rate_homog = tot_rate_homog + pnode->birth_params.homogeneous_rate +
          pnode->death_params.homogeneous_rate;
      }

      if(gptime.count_alleles) ChangeAncestorAllele(pnode, true);
      // Reprioritize clone based on size by moving to left
      CloneSort(pnode, true);

      flag = true;
      break;

    }
    // Condition if a death occurs
    if(rand_next_event <= summand + (pnode->cell_count) * ((pnode->birth_rate) + (pnode->death_rate)) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      tot_cell_count = tot_cell_count - 1;
      // decrease allele count in ancestors
      if(gptime.count_alleles) ChangeAncestorAllele(pnode, false);
      tot_rate_homog = tot_rate_homog - pnode->birth_params.homogeneous_rate -
        pnode->death_params.homogeneous_rate;

      /* Need to figure out
      // Clean up of clones with zero to speed up later runs
      if(pnode->cell_count == 0)
      {
        // remove pnode by attaching it to deadlist
        CutNodeOut(pnode);
        flag = true;
        break;
      }
      */

      // sort by moving to the right until fits
      CloneSort(pnode, false);
      flag = true;
      break;
    }

    summand = summand + (pnode->cell_count) * ((pnode->birth_rate + pnode->death_rate));
    pnode = pnode->nextnode;

  }

  if (!flag)
  {
    std::cout << "error: step not completed" << "\n";
    exit(0);
  }
}

/* *********************************************************************
  Beginning of function classes for modifying and inserting a new clone. The
  classes all require InsertNode(new_clone, mutation_count) at the end, but
  rates can be adjusted in the new clones to account for changes in rates or
  other fitness parameters
  *********************************************************************
*/

/*
  NewCloneNoParams is a simple birth death mutation process with only passenger
  mutations.
*/
void TDCloneList::NewCloneNoParams::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  NewCloneFitMut is a simple birth death mutation process where the new clone
  can have a different birth rate from a fitness distribution or a different
  mutation probability based on the parameters imported.
*/
void TDCloneList::NewCloneFitMut::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  // updating time-dependent parameters
  new_clone->birth_params = parent_clone->birth_params;
  new_clone->death_params = parent_clone->death_params;

  bool did_count_driver = false;
  // generation of additive rate to the fitness
  if(fit_params.is_randfitness)
  {
    double additional_rate = (*TDGenerateFitness)(fit_params, rng);

    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      did_count_driver = true;
      new_clone->is_driver = true;
    }
    new_clone->birth_params.homogeneous_rate = fmax(0, additional_rate + new_clone->birth_params.homogeneous_rate);
    new_clone->birth_params.additional_rate = fmax(0, additional_rate + new_clone->birth_params.additional_rate);
  }

  if(mut_params.is_mutator)
  {
    double additional_mut_prob = TDGenerateMutationProb(mut_params, rng);
    if(additional_mut_prob > 0)
    {
      if(!did_count_driver) new_clone->driver_count++;
      new_clone->is_driver = true;
    }
    new_clone->mut_prob = fmin(1, parent_clone->mut_prob + additional_mut_prob);
  }
  else
  {
    new_clone->mut_prob = parent_clone->mut_prob;
  }

  (new_clone->B).function = rate_function_array[new_clone->birth_params.type];
  (new_clone->D).function = rate_function_array[new_clone->death_params.type];
  (new_clone->B).params = &(new_clone->birth_params);
  (new_clone->D).params = &(new_clone->death_params);

  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  NewClonePunct allows for punctuated equilibrium, where a burst of mutations
  can occur which are denoted in the ID. A burst can decrease or the fitness the
  fitness with some probability parameter
*/
void TDCloneList::NewClonePunct::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  // generation of punctuated number of mutations
  int number_mutations = 1;

  double rand_punct = gsl_ran_flat(rng, 0, 1);
  double rand_advantage = 0;
  if(rand_punct < punct_params.punctuated_prob)
  {
    number_mutations = TDGeneratePunctuation(punct_params, rng);
    rand_advantage = gsl_ran_flat(rng, 0, 1);
  }

  // updating time-dependent parameters
  new_clone->birth_params = parent_clone->birth_params;
  new_clone->death_params = parent_clone->death_params;

  bool did_count_driver = false;
  // generation of additive rate to the fitness
  if(fit_params.is_randfitness)
  {
    double additional_rate = (*TDGenerateFitness)(fit_params, rng);
    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      did_count_driver = true;
      new_clone->is_driver = true;
      if(rand_punct < punct_params.punctuated_prob)
      {
        additional_rate = additional_rate * punct_params.punctuated_multiplier;
        // the additional rate can go to the birth or the death rate
        if( rand_advantage < punct_params.punctuated_advantageous_prob )
        {
          new_clone->birth_params.homogeneous_rate = fmax(0, additional_rate + new_clone->birth_params.homogeneous_rate);
          new_clone->birth_params.additional_rate = fmax(0, additional_rate + new_clone->birth_params.additional_rate);
        }
        else
        {
          new_clone->death_params.homogeneous_rate = fmax(0, additional_rate + new_clone->death_params.homogeneous_rate);
          new_clone->death_params.additional_rate = fmax(0, additional_rate + new_clone->death_params.additional_rate);
        }
      }
    }
  }

  if(mut_params.is_mutator)
  {
    double additional_mut_prob = TDGenerateMutationProb(mut_params, rng);
    if(additional_mut_prob > 0)
    {
      if(!did_count_driver) new_clone->driver_count++;
      new_clone->is_driver = true;
    }
    new_clone->mut_prob = fmin(1, parent_clone->mut_prob + additional_mut_prob);
  }
  else
  {
    new_clone->mut_prob = parent_clone->mut_prob;
  }

  (new_clone->B).function = rate_function_array[new_clone->birth_params.type];
  (new_clone->B).params = &(new_clone->birth_params);
  (new_clone->D).function = rate_function_array[new_clone->death_params.type];
  (new_clone->D).params = &(new_clone->death_params);
  // end of update rates

  cl.InsertNode(new_clone, parent_clone, number_mutations);
}

/*
  NewCloneEpi changes the fitness after k mutations have exists in the popoulation
*/
void TDCloneList::NewCloneEpi::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  // updating time-dependent parameters
  new_clone->birth_params = parent_clone->birth_params;
  new_clone->death_params = parent_clone->death_params;

  bool did_count_driver = false;
  // generation of additive rate to the fitness
  if(fit_params.is_randfitness)
  {
    double additional_rate = (*TDGenerateFitness)(fit_params, rng);
    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      did_count_driver = true;
      new_clone->is_driver = true;
      if(new_clone->mut_count == epi_params.epistatic_mutation_thresh)
      {
        additional_rate = additional_rate * epi_params.epistatic_multiplier;
      }
      new_clone->birth_params.homogeneous_rate = fmax(0, additional_rate + new_clone->birth_params.homogeneous_rate);
    }
  }

  if(mut_params.is_mutator)
  {
    double additional_mut_prob = TDGenerateMutationProb(mut_params, rng);
    if(additional_mut_prob > 0)
    {
      if(!did_count_driver) new_clone->driver_count++;
      new_clone->is_driver = true;
    }
    new_clone->mut_prob = fmin(1, parent_clone->mut_prob + additional_mut_prob);
  }
  else
  {
    new_clone->mut_prob = parent_clone->mut_prob;
  }

  (new_clone->B).function = rate_function_array[new_clone->birth_params.type];
  (new_clone->B).params = &(new_clone->birth_params);
  (new_clone->D).function = rate_function_array[new_clone->death_params.type];
  (new_clone->D).params = &(new_clone->death_params);
  // end of update rates

  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  For inserting a new custom function
*/
void TDCloneList::NewCloneCustom::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  Final output to file after running through a simulation
*/
void TDCloneList::Traverse(std::ofstream &F, int sim_number, bool count_alleles = false)
{
  struct clone *pnode;
  if(count_alleles)
  {
    for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->allele_count << "\t" <<
           pnode->birth_params.homogeneous_rate << "\t" <<
           pnode->death_params.homogeneous_rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }

    for (pnode = deadroot; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->allele_count << "\t" <<
           pnode->birth_params.homogeneous_rate << "\t" <<
           pnode->death_params.homogeneous_rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }
  }
  else
  {
    for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->birth_params.homogeneous_rate << "\t" <<
           pnode->death_params.homogeneous_rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }

    for (pnode = deadroot; pnode != NULL; pnode = pnode->nextnode)
    {
      F << sim_number << "\t" <<
           pnode->clone_id << "\t" <<
           pnode->cell_count << "\t" <<
           pnode->birth_params.homogeneous_rate << "\t" <<
           pnode->death_params.homogeneous_rate << "\t" <<
           pnode->mut_prob << "\t" <<
           pnode->clone_time << "\t" <<
           pnode->subclone_count << "\t" <<
           pnode->mut_count << "\t" <<
           pnode->driver_count << "\t" <<
           pnode->is_driver << "\n";
    }
  }
}

/*
  For traversing and outputting at a designated observation time not the end
*/
void TDCloneList::Traverse(std::ofstream &F, int sim_number, double obs_time, bool ancestry = false, bool count_alleles = false)
{
  struct clone *pnode;
  if(ancestry)
  {
    if(count_alleles)
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" << pnode->allele_count << "\t" <<
             GSL_FN_EVAL(&(pnode->B), obs_time) - GSL_FN_EVAL(&(pnode->D), obs_time) << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << GSL_FN_EVAL(&(pnode->parent->B), obs_time) - GSL_FN_EVAL(&(pnode->parent->D), obs_time) << "\t" <<
               pnode->parent->clone_time << "\n";
        }
      }
    }
    else
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" <<
             GSL_FN_EVAL(&(pnode->B), obs_time) - GSL_FN_EVAL(&(pnode->D), obs_time) << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << GSL_FN_EVAL(&(pnode->parent->B), obs_time) - GSL_FN_EVAL(&(pnode->parent->D), obs_time) << "\t" <<
               pnode->parent->clone_time << "\n";
        }
      }
    }
  }
  else
  {
    if(count_alleles)
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" << pnode->allele_count << "\t" <<
             GSL_FN_EVAL(&(pnode->B), obs_time) - GSL_FN_EVAL(&(pnode->D), obs_time) << "\t" <<
             pnode->clone_time << "\n";
      }
    }
    else
    {
      for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
      {
        F << sim_number << "\t" <<
             obs_time << "\t" <<
             pnode->clone_id << "\t" <<
             pnode->cell_count << "\t" <<
             GSL_FN_EVAL(&(pnode->B), obs_time) - GSL_FN_EVAL(&(pnode->D), obs_time) << "\t" <<
             pnode->clone_time << "\n";
      }
    }
  }
}

/*
  For sampling individuals from the population multiple times.
*/
void TDCloneList::SampleAndTraverse(std::ofstream &F, int run, int sample_size, int nsamples, gsl_rng* rng)
{
  // loop through to repeat with all samples
  for(int sample_counter = 1; sample_counter <= nsamples; sample_counter++)
  {
    int samples_to_place = sample_size;
    int cells_left = tot_cell_count;

    // assign cells from samples_to_place to nodes using binomial
    // probabilities (same as multinomial r.v. but don't need to
    // know all beforehand
    struct clone *pnode = root;
    while(samples_to_place > 0)
    {
      // calculate prob and simulate number to sample from this clone
      double prob = (double)pnode->cell_count / (double)cells_left;
      int samples_placed = gsl_ran_binomial(rng, prob, samples_to_place);

      // only write if sampled any
      if(samples_placed > 0)
      {
        F << run << "\t" << sample_counter << "\t" <<
          pnode->clone_id << "\t" << samples_placed << "\n";
      }

      samples_to_place = samples_to_place - samples_placed;
      cells_left = cells_left - pnode->cell_count;

      pnode = pnode->nextnode;
    }
  }
}

/*
  Run through all clones and check that they have a large enough popsize -
  remove those that don't by putting their count into their parent pop
  then deleting the node
*/
void TDCloneList::TreeTrim(double threshold, int max_pop)
{
  //Move to last clone
  while(currnode->nextnode != NULL)
  {
    currnode = currnode->nextnode;
  }

  double cell_cutoff = threshold * max_pop;
  // Starting with current node and working back until at root
  while( currnode->prevnode != NULL )
  {
    if( currnode->cell_count < cell_cutoff )
    {
      // add current clones numbers to that of parent
      currnode->parent->cell_count = currnode->parent->cell_count + currnode->cell_count;
      // remove the node and update pointers
      DeleteNode();
    }
    else
    {
      currnode = currnode->prevnode;
    }

  }
}

/*
  When a node hits 0 count, move to dead side of list
  if at end of list, no need to link the next node (NULL value), Otherwise
  cut out value
*/
void TDCloneList::CutNodeOut(struct clone* zeronode)
{
  if(zeronode->nextnode != NULL)
  {
    zeronode->nextnode->prevnode = zeronode->prevnode;
  }
  /* if at beginning of list, need to reroot the list to the next node
     otherwise, cut out node*/
  if(zeronode->prevnode != NULL)
  {
    zeronode->prevnode->nextnode = zeronode->nextnode;
  }
  else
  {
    root = zeronode->nextnode;
  }

  /* if dead size not rooted, root with first zero node,
     otherwise add to end of list*/
  if( deadroot == NULL)
  {
    deadroot = zeronode;
    currdeadnode = zeronode;
  }
  else
  {
    zeronode->prevnode = currdeadnode;
    currdeadnode->nextnode = zeronode;
    currdeadnode = zeronode;
    zeronode->nextnode = NULL;
  }
}

// This function used for TreeTrim to remove nodes that are too small
void TDCloneList::DeleteNode()
{
  struct clone *tmpcurrnode;
  tmpcurrnode = currnode;

  if(currnode->nextnode == NULL)
  {
    currnode = currnode->prevnode;
    currnode->nextnode = NULL;
  }
  else
  {
    struct clone *tmpnextnode;
    tmpnextnode = currnode->nextnode;
    currnode = currnode->prevnode; // set previous node as current
    currnode->nextnode = tmpnextnode; // set the nextnode of the current node to skip over tmpnode
    tmpnextnode->prevnode = currnode;
  }

  delete tmpcurrnode;
  num_clones = num_clones - 1; // increase numtypes (have 1 type in pop)
}

// Delete all members of the linked list
void TDCloneList::DeleteList()
{
  struct clone *pnode;
  while (root != NULL)
  {
    pnode = root;
    root = root->nextnode;
    delete pnode;
  }
}
