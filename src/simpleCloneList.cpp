/*
 * =====================================================================================
 *
 *       Filename:  clonelist.cpp
 *
 *    Description: Contains the class SimpleCloneList which creates a linked list
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

#include "simpleCloneList.h"

// Constructor
void SimpleCloneList::init()
{
    root = NULL; // node points to first ancestor clone
    currnode = NULL; // for moving around the list and keeping track of location
    num_clones = 0; // number of types
    tot_rate = 0; // total rate (birth + death)
    tot_cell_count = 0; // total number of cells in all types
}

/*
  Insert the ancestors into the linked list. If ancestors are defined
  in a separate file it parses the file and inserts. Otherwise it
  creates nodes based on the number of clones and individuals defined
  in input file.
*/
void SimpleCloneList::InsertAncestor(struct clone* ancestor)
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

  if( root == NULL )
  {
    ancestor->prevnode = NULL;
    ancestor->nextnode = NULL;
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
  Final output to file after running through a simulation
*/
void SimpleCloneList::Traverse(std::ofstream &F, int sim_number)
{
  struct clone *pnode;
  for (pnode = root; pnode != NULL; pnode = pnode->nextnode)
  {
    F << sim_number << "\t" <<
         pnode->clone_id << "\t" <<
         pnode->cell_count << "\t" <<
         pnode->birth_rate << "\t" <<
         pnode->death_rate << "\n";
  }
}

/*
  For sampling individuals from the population multiple times.
*/
void SimpleCloneList::SampleAndTraverse(std::ofstream &F, int sim_number, int sample_size, int nsamples, gsl_rng* rng)
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

// This function used for TreeTrim to remove nodes that are too small
void SimpleCloneList::DeleteNode()
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
void SimpleCloneList::DeleteList()
{
  struct clone *pnode;
  while (root != NULL)
  {
    pnode = root;
    root = root->nextnode;
    delete pnode;
  }
}
