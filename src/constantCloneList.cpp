/*
 * =====================================================================================
 *
 *       Filename:  clonelist.cpp
 *
 *    Description: Contains the class ConstantCloneList which creates a linked list
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
#include "constantCloneList.h"
#include <Rcpp.h>

// Constructor
void ConstantCloneList::init()
{
    root = NULL; // node points to first ancestor clone
    deadroot = NULL; // node points to first extinct clone (to keep track)
    currnode = NULL; // for moving around the list and keeping track of location
    currdeadnode = NULL;
    num_clones = 0; // number of types
    num_mutations = 0; // number of mutations (types and deadtypes)
    tot_rate = 0; // total rate (birth + death)
    tot_cell_count = 0; // total number of cells in all types
}

// InsertNode attaches the newly created node (in NewConstantClone) to the end of the
// linked list. Number of mutations was included for the punctuated scenario
// to add to the ID.
void ConstantCloneList::InsertNode(struct clone* newnode, struct clone* parentnode, int number_mutations)
{
  if( root == NULL ) // if linked list not yet rooted (no clones created)
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
    currnode = newnode;
  }
}

/*
  Insert the ancestors into the linked list. If ancestors are defined
  in a separate file it parses the file and inserts. Otherwise it
  creates nodes based on the number of clones and individuals defined
  in input file.
*/
void ConstantCloneList::InsertAncestor(struct clone* ancestor)
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
    ancestor->nextnode = NULL;
    ancestor->parent = NULL;
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
    currnode = ancestor;
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
void ConstantCloneList::ChangeAncestorAllele(struct clone* thisnode, bool add_daughter)
{
  if(add_daughter) // if added a daughter then add to each allele_count
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
void ConstantCloneList::CloneSort(struct clone* sortnode, bool is_birth)
{
  if(is_birth && sortnode->prevnode != NULL) // if a birth occurred and not root node
  {
    // if the rates are greater in node than the previous node, switch them
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
  } else if(!is_birth && sortnode->nextnode != NULL) // if death and not final node
  {
    // if the rates are greater in node than the previous node, switch them
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
  AdvanceTime determines the next time in the process.
*/
double ConstantCloneList::AdvanceTime(double curr_time, gsl_rng* rng)
{
  double rand_next_time;
  rand_next_time = gsl_ran_exponential(rng, 1 / tot_rate);
  return rand_next_time;
}

/*
  Advance state determines the next event that occurs in the process and which
  clone it occurs in (birth w/o mutation, birth w/ mutation, death) and updates
  the process accordingly
*/
void ConstantCloneList::AdvanceState(double curr_time, double next_time, gsl_rng* rng)
{
  double summand = 0;
  double rand_next_event = gsl_ran_flat(rng, 0, tot_rate);
  double rand_mut_occur;
  bool flag = false;

  struct clone *pnode;
  pnode = root;


  // cycle through nodes to determine the next event that occurs
  while ( (pnode) && !(flag) )
  {
    // Condition for new birth
    if ( rand_next_event <= summand + (pnode->cell_count * pnode->birth_rate) )
    {
      rand_mut_occur = gsl_ran_flat(rng, 0, 1);
      bool is_mutation = false;
      // Check which mutation rate to use based on population size
      if(tot_cell_count < gpcons.max_pop_mutation)
      {
        is_mutation = (rand_mut_occur <= pnode->mut_prob);
      }
      else
      {
        is_mutation = (rand_mut_occur <= gpcons.max_pop_mut_rate;
      }
      // Condition to determine if mutation occurs in new daughter
      if ( is_mutation )
      {
        // Creation of new clone
        struct clone *new_mut_node;
        new_mut_node = new struct clone;
        new_mut_node->clone_id = pnode->clone_id;
        new_mut_node->driver_count = pnode->driver_count;
        new_mut_node->subclone_count = 0;
        new_mut_node->cell_count = 1;
        new_mut_node->allele_count = 1;
        new_mut_node->clone_time = curr_time + next_time;
        new_mut_node->mut_count = pnode->mut_count + 1;
        new_mut_node->is_driver = false;
        // Update parent subclones
        pnode->subclone_count = pnode->subclone_count + 1;
        new_mut_node->mut_prob = pnode->mut_prob;
        new_mut_node->birth_rate = pnode->birth_rate;
        new_mut_node->death_rate = pnode->death_rate;

        // Function class to add any modifications to rates and insert clone
        (*NewConstantClone)(new_mut_node, pnode);

        tot_rate = tot_rate + new_mut_node->birth_rate + new_mut_node->death_rate;
        tot_cell_count++;
        num_clones++;
        num_mutations++;

        currnode = new_mut_node;
        // add allele count to all ancestors
        if(gpcons.count_alleles)
        {
          ChangeAncestorAllele(pnode, true);
        }
      }
      else // no mutation, increment clone by 1
      {
        pnode->cell_count++;
        tot_cell_count++;
        tot_rate = tot_rate + (pnode->birth_rate) + (pnode->death_rate);
        currnode = pnode;
        if(gpcons.count_alleles) ChangeAncestorAllele(pnode, true);
        // Reprioritize clone based on size by moving to left
        CloneSort(pnode, true);
      }

      flag = true;
      break;
    }
    // Condition if a death occurs
    if(rand_next_event <= summand + (pnode->cell_count * (pnode->birth_rate + pnode->death_rate) ) )
    {
      // Death occurs
      pnode->cell_count = pnode->cell_count - 1;
      tot_cell_count = tot_cell_count - 1;
      tot_rate = tot_rate - pnode->birth_rate - pnode->death_rate;
      currnode = pnode;
      // decrease allele count in ancestors
      if(gpcons.count_alleles) ChangeAncestorAllele(pnode, false);

      /* Couldn't figure out issue
      // Clean up of clones with zero to speed up later runs
      if(pnode->cell_count == 0)
      {
        struct clone *zeronode;
        zeronode = pnode;
        // remove pnode by attaching it to deadlist
        CutNodeOut(zeronode);
      }
      */

      // sort by moving to the right until fits
      CloneSort(pnode, false);

      flag = true;
      break;

    }

    summand = summand + (pnode->cell_count) * (pnode->birth_rate + pnode->death_rate);
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
void ConstantCloneList::NewCloneNoParams::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  // Do Nothing - nothing happens in this case
  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  NewCloneFitMut is a simple birth death mutation process where the new clone
  can have a different birth rate from a fitness distribution or a different
  mutation probability based on the parameters imported.
*/
void ConstantCloneList::NewCloneFitMut::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  // updating rate parameters
  bool did_count_driver = false;
  // generation of additive rate to the fitness
  if(fit_params.is_randfitness)
  {
    double additional_rate = 0;
    (*ConstantGenerateFitness)(&additional_rate, &fit_params, rng);

    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      did_count_driver = true;
      new_clone->is_driver = true;
    }
    new_clone->birth_rate = fmax(0, additional_rate + parent_clone->birth_rate);
  }
  else
  {
    new_clone->birth_rate = parent_clone->birth_rate;
  }

  if(mut_params.is_mutator)
  {
    double additional_mut_prob = ConstantGenerateMutationProb(mut_params, rng);
    if(additional_mut_prob > 0)
    {
      if(!did_count_driver) new_clone->driver_count++;
      new_clone->is_driver = true;
    }
    new_clone->mut_prob = fmin(1, parent_clone->mut_prob + additional_mut_prob);
  }

  // Insert new clone
  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  NewClonePunct allows for punctuated equilibrium, where a burst of mutations
  can occur which are denoted in the ID. A burst can decrease or the fitness the
  fitness with some probability parameter
*/
void ConstantCloneList::NewClonePunct::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  int number_mutations = 1;
  // generation of punctuated number of mutations
  double rand_punct = gsl_ran_flat(rng, 0, 1);
  double rand_advantage = 0;
  if(rand_punct < punct_params.punctuated_prob)
  {
    number_mutations = ConstantGeneratePunctuation(punct_params, rng);
    rand_advantage = gsl_ran_flat(rng, 0, 1);
  }

  bool did_count_driver = false;
  // generation of additive rate to the fitness
  if(fit_params.is_randfitness)
  {
    double additional_rate = 0;
    (*ConstantGenerateFitness)(&additional_rate, &fit_params, rng);
    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      new_clone->is_driver = true;
      did_count_driver = true;
      new_clone->birth_rate = fmax(0, additional_rate + new_clone->birth_rate);

      if(rand_punct < punct_params.punctuated_prob)
      {
        additional_rate = additional_rate * punct_params.punctuated_multiplier;
        // the additional rate can go to the birth or the death rate
        if( rand_advantage < punct_params.punctuated_advantageous_prob)
        {
          new_clone->birth_rate = fmax(0, additional_rate + new_clone->birth_rate);
        }
        else
        {
          new_clone->death_rate = fmax(0, additional_rate + new_clone->death_rate);
        }
      }
      else
      {
        new_clone->birth_rate = fmax(0, additional_rate + new_clone->birth_rate);
      }
    }
  }

  if(mut_params.is_mutator)
  {
    double additional_mut_prob = ConstantGenerateMutationProb(mut_params, rng);
    if(additional_mut_prob > 0)
    {
      if(!did_count_driver) new_clone->driver_count++;
      new_clone->is_driver = true;
    }
    new_clone->mut_prob = fmin(1, new_clone->mut_prob + additional_mut_prob);
  }

  // Insert new clone
  cl.InsertNode(new_clone, parent_clone, number_mutations);
}

/*
  NewCloneEpi changes the fitness after k mutations have exists in the popoulation
*/
void ConstantCloneList::NewCloneEpi::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  new_clone->death_rate = parent_clone->death_rate;
  new_clone->birth_rate = parent_clone->birth_rate;

  bool did_count_driver = false;
  // generation of additive rate to the fitness
  if(fit_params.is_randfitness)
  {
    double additional_rate = 0;
    (*ConstantGenerateFitness)(&additional_rate, &fit_params, rng);
    if (additional_rate > 0)
    {
      new_clone->driver_count++;
      did_count_driver = true;
      new_clone->is_driver = true;
      if(new_clone->mut_count == epi_params.epistatic_mutation_thresh)
      {
        additional_rate = additional_rate * epi_params.epistatic_multiplier;
      }
    }
    new_clone->birth_rate = fmax(0, additional_rate + parent_clone->birth_rate);
  }

  if(mut_params.is_mutator)
  {
    double additional_mut_prob = ConstantGenerateMutationProb(mut_params, rng);
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
  // Insert new clone
  cl.InsertNode(new_clone, parent_clone, 1);
}

/*
  For inserting a new custom function
*/
void ConstantCloneList::NewCloneCustom::operator()(struct clone *new_clone, struct clone *parent_clone)
{
  // Insert custom code here for how to update a new clone - need to put all params
  (*CreateNewCustomClone)(new_clone, parent_clone, &fit_params, &mut_params,
    &punct_params, &epi_params, &num_new_muts, rng, ConstantGenerateFitness);

  // End with this piece - Insert new clone
  cl.InsertNode(new_clone, parent_clone, num_new_muts);
}

/*
  Final output to file after running through a simulation
*/
void ConstantCloneList::Traverse(std::ofstream &F, int sim_number, bool count_alleles = false)
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
           pnode->birth_rate << "\t" <<
           pnode->death_rate << "\t" <<
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
void ConstantCloneList::Traverse(std::ofstream &F, int sim_number, double obs_time, bool ancestry = false, bool count_alleles = false)
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << pnode->parent->birth_rate - pnode->parent->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
             pnode->clone_time << "\t";

        if(pnode->parent == NULL)
        {
          F << "NA" << "\t" << "NA" << "\n";
        }
        else
        {
          F << pnode->parent->birth_rate - pnode->parent->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
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
             pnode->birth_rate - pnode->death_rate << "\t" <<
             pnode->clone_time << "\n";
      }
    }
  }
}

/*
  For sampling individuals from the population multiple times.
*/
void ConstantCloneList::SampleAndTraverse(std::ofstream &F, int sim_number, int sample_size, int nsamples, gsl_rng* rng)
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
        F << sim_number << "\t" << sample_counter << "\t" <<
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
void ConstantCloneList::TreeTrim(double threshold, int max_pop)
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
      if(currnode->parent != NULL)
      {
        currnode->parent->cell_count = currnode->parent->cell_count + currnode->cell_count;
      }
      // remove the node and update pointers
      DeleteNode();
    }
    else
    {
      currnode = currnode->prevnode;
    }

  }
}

/****  CURRENTLY NOT USED ****
  When a node hits 0 count, move to dead side of list
  if at end of list, no need to link the next node (NULL value), Otherwise
  cut out value
*/
void ConstantCloneList::CutNodeOut(struct clone* zeronode)
{
  if(zeronode->prevnode == NULL) // if root
  {
    if(zeronode->nextnode == NULL)
    {
      root = NULL;
    }
    else
    {
      root = zeronode->nextnode;
      root->prevnode = NULL;
    }
  }
  else
  {
    if(zeronode->nextnode != NULL)
    {
      zeronode->nextnode->prevnode = zeronode->prevnode;
      zeronode->prevnode->nextnode = zeronode->nextnode;
    }
    else
    {
      zeronode->prevnode->nextnode = NULL;
    }
  }

  // if dead not rooted, root with first zero node,
  //   otherwise add to end of list
  if( deadroot == NULL)
  {
    deadroot = zeronode;
    zeronode->prevnode = NULL;
  }
  else
  {
    zeronode->prevnode = currdeadnode;
    currdeadnode->nextnode = zeronode;
  }
  currdeadnode = zeronode;
  zeronode->nextnode = NULL;
}

// This function used for TreeTrim to remove nodes that are too small
void ConstantCloneList::DeleteNode()
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
void ConstantCloneList::DeleteList()
{
  struct clone *pnode;
  while (root != NULL)
  {
    pnode = root;
    root = root->nextnode;
    delete pnode;
  }
}
