/*
 * =============================================================================
 *
 *       Filename:  SIApop.cpp
 *
 *    Description: Birth-Death-Mutation process simulation for infinite-allele
 *                 model with random fitness contributions using Gillespie
 *                 Algorithm. Imports data, runs SSA and outputs to designated
 *                 folder.
 *
 *        Version:  1.0
 *        Created:  08/24/2016 16:50:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thomas McDonald (), mcdonald@jimmy.harvard.edu
 *   Organization:  DFCI
 *
 * =============================================================================
 */
 #include <RcppGSL.h>
 #include <Rcpp.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <vector>
#include <cstring>
#include <string>
#include <cmath>
#include <gsl/gsl_randist.h>

#include "simpleGlobalStructs.h"
#include "simpleCloneList.h"
#include "simpleParameterList.h"

// structure contains all global parameters used in multiple source files
GlobalParameters gpsimp;

//' SIApop for non-mutating processes
//'
//' @param input input character vector of input file
//' @param output_dir input character vector of output location
//' @param ancestor_file input character vector of ancestor file
//' @export
// [[Rcpp::export]]
int siapopSimple(Rcpp::Nullable<Rcpp::CharacterVector> input = R_NilValue,
                 Rcpp::Nullable<Rcpp::CharacterVector> output_dir = R_NilValue,
                 Rcpp::Nullable<Rcpp::CharacterVector> ancestor_file = R_NilValue)
{

  //  declaring random number generator and setting seed
  gpsimp.rng = gsl_rng_alloc(gsl_rng_mt19937);
  gpsimp.seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  /*
   VARIABLE INPUT AND CONVERSION
  */
  // parsing arguments in command line by searching for argument options
  // default output is current directory
  const char *output_folder;
  if( output_dir.isNull() )
  {
    output_folder = "./";
  }
  else
  {
    std::vector<std::string> output_dir_ = Rcpp::as<std::vector <std::string> > (output_dir);
    output_folder = output_dir_[0].c_str();
  }


  // declare and open output stream for simulation statistics
  char fn[100];
  std::ofstream sim_stats;
  sprintf(fn,"%s/sim_stats.txt", output_folder);
  sim_stats.open(fn);
  sim_stats.setf(std::ios::fixed);
  sim_stats.precision(8);


  // declare and initialize parameter list for simulation
  SimpleParameterList params;
  params.init();
  int count_extinct = 0;


  // parsing through the input file and converting/adding to parameter list
  if( input.isNotNull() )
  {
    std::vector<std::string> input_ = Rcpp::as<std::vector <std::string> > (input);
    const char* input_params = input_[0].c_str();

    std::string s = input_params;
    std::ifstream infile(input_params);
    if(infile.is_open())
    {
        while(getline(infile, s))
        {
            // allow comments
            if(s[0] == '#') continue;
            params.SplitAndFill(s);
        }
    }
  }

  for(std::map<std::string, std::string>::iterator it=params.begin(); it!=params.end(); ++it)
  {
    sim_stats << it->first << ", " << it->second << "\n";
  }

  // convert all parameters imported from file into respective values in gpsimp
  params.convert("tot_life", gpsimp.tot_life);
  params.convert("max_pop", gpsimp.max_pop);
  params.convert("ancestors", gpsimp.ancestors);
  params.convert("ancestor_clones", gpsimp.ancestor_clones);
  params.convert("num_sims", gpsimp.num_sims);
  params.convert("num_samples", gpsimp.num_samples);
  params.convert("sample_size", gpsimp.sample_size);
  params.convert("detection_threshold", gpsimp.detection_threshold);
  params.convert("allow_extinction", gpsimp.allow_extinction);
  params.convert("birth_rate", gpsimp.birth_rate);
  params.convert("death_rate", gpsimp.death_rate);

  /*
    END OF VARIABLE INPUT AND CONVERSION
  */


  // declare and open other output streams for time and end of sim clone list
  std::ofstream clonedata;
  sprintf(fn, "%s/clonedata.txt", output_folder);
  clonedata.open(fn);
  clonedata.setf(std::ios::fixed);
  clonedata.precision(12);
  clonedata << "run\tunique_id\tnumcells\tbirthrate\tdeathrate\n";


  // Open output stream for sampling data
  std::ofstream sample_data;
  if(gpsimp.sample_size > 0 & gpsimp.num_samples > 0)
  {
    sprintf(fn, "%s/sampledata.txt", output_folder);
    sample_data.open(fn);
    sample_data.setf(std::ios::fixed);
    sample_data << "run\tsample_number\tunique_id\tnumber_obs\n";
  }

  // set RNG seed
  gsl_rng_set(gpsimp.rng, gpsimp.seed);


  // Beginning of simulation that has "sim" number of runs
  for (int sim = 1; sim <= gpsimp.num_sims; sim++)
  {
    Rcpp::checkUserInterrupt();
    // Define CloneList population and initialize variables;
    SimpleCloneList population;
    population.init();
    population.tot_cell_count = 0;

    if( ancestor_file.isNull() ) // if no ancestor file exists
    {
      // go through all ancestors and initialize clones for each
      for(int ance__clone_count = 1; ance__clone_count <= gpsimp.ancestor_clones; ance__clone_count++)
      {
        // Create Ancestor Node
        struct clone* ancestor;
        ancestor = new struct clone;

        ancestor->cell_count = 0;
        ancestor->birth_rate = gpsimp.birth_rate;
        ancestor->death_rate = gpsimp.death_rate;

        double meangrowth = exp((ancestor->birth_rate - ancestor->death_rate) * gpsimp.tot_life);
        double alpha = (ancestor->death_rate * meangrowth - ancestor->death_rate) /
                       (ancestor->birth_rate * meangrowth - ancestor->death_rate);
        double beta = (ancestor->birth_rate * meangrowth - ancestor->birth_rate) /
                       (ancestor->birth_rate * meangrowth - ancestor->death_rate);

        // Update ancestors to time
        for(int i = 0; i < gpsimp.ancestors; ++i)
        {
          int not_zero = gsl_ran_bernoulli(gpsimp.rng, 1 - alpha);
          if(not_zero == 1)
          {
            int new_cells = gsl_ran_geometric(gpsimp.rng, 1 - beta);
            ancestor->cell_count = ancestor->cell_count + new_cells;
            population.tot_cell_count = population.tot_cell_count + new_cells;
          }
        }

        population.InsertAncestor(ancestor);
      }
    }
    else // ancestor file exists to read from
    {
      Rcpp::Rcout << "Reading ancestor file...";

      // import first line of file to a vector of keys
      std::vector<std::string> ancestor_keys;
      std::vector<std::string>::iterator it;


      std::vector<std::string> ancestor_file_ = Rcpp::as<std::vector <std::string> > (ancestor_file);
      const char *ancestors = ancestor_file_[0].c_str();
      std::string a = ancestors;
      std::ifstream ancfile(ancestors);

      if(ancfile.is_open())
      {
        getline(ancfile, a);
        std::istringstream ss( a );
        std::string s2;
        while(ss >> s2)
        {
          ancestor_keys.push_back(s2);
        }

        // Loop through all lines of the file
        while(getline(ancfile, a))
        {
          // create a map between column names (from vector) and values
          std::map<std::string, std::string> ancestor_map;
          std::istringstream ss( a );
          std::string s2;

          for(it = ancestor_keys.begin(); it < ancestor_keys.end(); it++)
          {
            ss >> s2;
            ancestor_map[*it];
            ancestor_map[*it] = s2;
          }
            // Move map values to structures
            struct clone* ancestor;
            ancestor = new struct clone;

            ancestor->clone_id = !ancestor_map["unique_id"].empty() ? ancestor_map["unique_id"] : "";
            ancestor->cell_count = 0;
            int num_ancestors = !ancestor_map["numcells"].empty() ? stoi(ancestor_map["numcells"]) : 1;
            ancestor->birth_rate = !ancestor_map["birthrate"].empty() ? stof(ancestor_map["birthrate"]) : gpsimp.birth_rate;
            ancestor->death_rate = !ancestor_map["deathrate"].empty() ? stof(ancestor_map["deathrate"]) : gpsimp.death_rate;


            double meangrowth = exp((ancestor->birth_rate - ancestor->death_rate) * gpsimp.tot_life);
            double alpha = (ancestor->death_rate * meangrowth - ancestor->death_rate) /
                           (ancestor->birth_rate * meangrowth - ancestor->death_rate);
            double beta = (ancestor->birth_rate * meangrowth - ancestor->birth_rate) /
                           (ancestor->birth_rate * meangrowth - ancestor->death_rate);

            // Update ancestors to time
            for(int i = 0; i < num_ancestors; ++i)
            {
              int not_zero = gsl_ran_bernoulli(gpsimp.rng, 1 - alpha);
              if(not_zero == 1)
              {
                int new_cells = gsl_ran_geometric(gpsimp.rng, 1 - beta);
                ancestor->cell_count = ancestor->cell_count + new_cells;
                population.tot_cell_count = population.tot_cell_count + new_cells;
              }
            }

            population.InsertAncestor(ancestor);

            ancestor_map.clear();
        }
      }
    }

    // nonextinction checker - if set to false and goes extinction, restart that sim
    if( (population.tot_cell_count == 0) && (gpsimp.allow_extinction == false) )
    {
      population.DeleteList();
      count_extinct++;
      gpsimp.num_sims++; // increase number of sims
      Rcpp::Rcout << "Population went extinct. Restarting.\n";
      continue;
    }

    // Sampling from population
    if(gpsimp.sample_size > 0 & gpsimp.num_samples > 0)
    {
      population.SampleAndTraverse(sample_data, sim, gpsimp.sample_size, gpsimp.num_samples);
    }
    // Output of end state with clone info
    Rcpp::Rcout << "Traversing and outputting run " << sim << "\n";
    population.Traverse(clonedata, sim);
    Rcpp::Rcout << "Traversal Done\n";

    population.DeleteList();

  }

  //avg_sim_endtime = avg_sim_endtime * (double)gpsimp.num_sims / (double)count_detect;

  gsl_rng_free(gpsimp.rng);

  sim_stats << "count_extinct, " << count_extinct  << "\n";

  clonedata.close();
  sample_data.close();
  sim_stats.close();

  return 0;
}
