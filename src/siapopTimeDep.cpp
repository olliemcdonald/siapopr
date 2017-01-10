/*
 * =============================================================================
 *
 *       Filename:  SIApop.cpp
 *
 *    Description: Time-Dependent Birth-Death-Mutation process simulation for
 *                 infinite-allele model with random fitness contributions using
 *                 Gillespie Algorithm. Imports data, runs SSA and outputs to
 *                 designated folder.
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
#include <cmath>
#include <string>
#include <gsl/gsl_randist.h>

#include "timedepGlobalStructs.h"
#include "timedepCloneList.h"
#include "timedepParameterList.h"
#include "timedepRateFunctions.h"

// structure contains all global parameters used in multiple source files
GlobalParameters gptime;
// Function array for different time-dependent function types
RateFunctionsPtr rate_function_array[] = {&RateFunctions::constant,
  &RateFunctions::linear, &RateFunctions::logistic,
  &RateFunctions::Gompertz, &RateFunctions::custom};
gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
// Pointer to Function class which will point to an instance of one based on parameters
TDCloneList::NewCloneFunction* NewTDClone;

//’ SIApop for time-dependent processes.
//’
//’ @param input input character vector of input file
//’ @param output_dir input character vector of output location
//’ @param ancestor_file input character vector of ancestor file
// [[Rcpp::export]]
int siapopTimeDep(Rcpp::Nullable<Rcpp::CharacterVector> input = R_NilValue,
                  Rcpp::Nullable<Rcpp::CharacterVector> output_dir = R_NilValue,
                  Rcpp::Nullable<Rcpp::CharacterVector> ancestor_file = R_NilValue)
{

  // track total error from integration - FOR TESTING
  gptime.tot_error = 0;

  //  declaring random number generator and setting seed
  gptime.rng = gsl_rng_alloc(gsl_rng_mt19937);
  gptime.seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

  /*
    VARIABLE INPUT AND CONVERSION
  */
  // parsing arguments in command line by searching for argument options
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
  TDParameterList params;
  params.init();

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

  // convert all parameters imported from file into respective values in gptime
  params.convert("tot_life", gptime.tot_life);
  params.convert("max_pop", gptime.max_pop);
  params.convert("start_time", gptime.start_time);
  params.convert("ancestors", gptime.ancestors);
  params.convert("ancestor_clones", gptime.ancestor_clones);
  params.convert("num_sims", gptime.num_sims);
  params.convert("num_samples", gptime.num_samples);
  params.convert("sample_size", gptime.sample_size);
  params.convert("detection_threshold", gptime.detection_threshold);
  params.convert("observation_frequency", gptime.observation_frequency);
  if (gptime.observation_frequency == 0) gptime.observation_frequency = gptime.tot_life;

  std::vector<double> observation_times;
  params.ParseVector("observation_times", observation_times);
  if (observation_times.size() == 0)
  {
    int num_obs =  ceil(gptime.tot_life / gptime.observation_frequency);
    for (int i = 1; i <= num_obs; i++)
    {
      observation_times.push_back(i * gptime.observation_frequency);
    }
  }
  if(observation_times.back() != gptime.tot_life) observation_times.push_back(gptime.tot_life);


  params.convert("allow_extinction", gptime.allow_extinction);
  params.convert("trace_ancestry", gptime.trace_ancestry);
  params.convert("count_alleles", gptime.count_alleles);
  params.convert("mutation_prob", gptime.mutation_prob);

  FitnessParameters fit_params;
  params.convert("alpha_fitness", fit_params.alpha_fitness);
  params.convert("beta_fitness", fit_params.beta_fitness);
  params.convert("pass_prob", fit_params.pass_prob);
  params.convert("upper_fitness", fit_params.upper_fitness);
  params.convert("lower_fitness", fit_params.lower_fitness);
  fit_params.is_randfitness = false;
  if( ((fit_params.alpha_fitness > 0 && fit_params.beta_fitness > 0) ||
      (fit_params.upper_fitness != fit_params.lower_fitness)) &&
      (fit_params.pass_prob < 1) )
  {
    fit_params.is_randfitness = true;
  }

  MutationParameters mut_params;
  params.convert("alpha_mutation", mut_params.alpha_mutation);
  params.convert("beta_mutation", mut_params.beta_mutation);
  params.convert("pass_prob", mut_params.pass_prob);
  mut_params.is_mutator = false;
  if( (mut_params.alpha_mutation > 0 && mut_params.beta_mutation > 0) &&
      (mut_params.pass_prob < 1) )
  {
    mut_params.is_mutator = true;
  }

  PunctuationParameters punct_params;
  params.convert("punctuated_prob", punct_params.punctuated_prob);
  params.convert("poisson_param", punct_params.poisson_param);
  params.convert("punctuated_multiplier", punct_params.punctuated_multiplier);
  params.convert("punctuated_advantageous_prob", punct_params.punctuated_advantageous_prob);
  punct_params.is_punctuated = punct_params.punctuated_prob > 0;

  EpistaticParameters epi_params;
  params.convert("epistatic_mutation_thresh", epi_params.epistatic_mutation_thresh);
  params.convert("epistatic_multiplier", epi_params.epistatic_multiplier);
  epi_params.is_epistasis = epi_params.epistatic_mutation_thresh > 0;
  if (epi_params.epistatic_multiplier == 1.0)
  {
    epi_params.is_epistasis = false;
  }


  TimeDependentParameters td_birth_params;
  params.convert("birth_function", td_birth_params.type);
  params.ParseVector("td_birth_params", td_birth_params.coefs);
  td_birth_params.homogeneous_rate = 1;
  td_birth_params.additional_rate = 0;

  TimeDependentParameters td_death_params;
  params.convert("death_function", td_death_params.type);
  params.ParseVector("td_death_params", td_death_params.coefs);
  td_death_params.homogeneous_rate = 1;
  td_death_params.additional_rate = 0;

  /*
    END OF VARIABLE INPUT AND CONVERSION
  */

  // declare and open other output streams for time and end of sim clone list
  std::ofstream clonedata;
  sprintf(fn, "%s/clonedata.txt", output_folder);
  clonedata.open(fn);
  clonedata.setf(std::ios::fixed);
  clonedata.precision(12);
  if(gptime.count_alleles)
  {
    clonedata << "run\tunique_id\tnumcells\tallelefreq\tbirthrate\tdeathrate\tmutprob\t"
      "initialtime\tsubclone_count\tnum_mut\tnum_drivers\tis_driver" << "\n";
  }
  else
  {
    clonedata << "run\tunique_id\tnumcells\tbirthrate\tdeathrate\tmutprob\t"
      "initialtime\tsubclone_count\tnum_mut\tnum_drivers\tis_driver" << "\n";
  }

  std::ofstream timedata;
  sprintf(fn, "%s/timedata.txt", output_folder);
  timedata.open(fn);
  timedata.setf(std::ios::fixed);
  timedata.precision(12);
  if(gptime.trace_ancestry)
  {
    if(gptime.count_alleles)
    {
      timedata << "run\ttime\tunique_id\tnumcells\tallelefreq\tgrowth_rate\tinitialtime\t"
        "parent_growth_rate\tparent_initialtime" << "\n";
    }
    else
    {
      timedata << "run\ttime\tunique_id\tnumcells\tgrowth_rate\tinitialtime\t"
        "parent_growth_rate\tparent_initialtime" << "\n";
    }
  }
  else
  {
    if(gptime.count_alleles)
    {
      timedata << "run\ttime\tunique_id\tnumcells\tallelefreq\tgrowth_rate\tinitialtime" << "\n";
    }
    else
    {
      timedata << "run\ttime\tunique_id\tnumcells\tgrowth_rate\tinitialtime" << "\n";
    }
  }

  // Open output stream for sampling data
  std::ofstream sample_data;
  if(gptime.sample_size > 0 & gptime.num_samples > 0)
  {
    sprintf(fn, "%s/sampledata.txt", output_folder);
    sample_data.open(fn);
    sample_data.setf(std::ios::fixed);
    sample_data << "run\tsample_number\tunique_id\tnumber_observed\n";
  }


  // set RNG seed
  gsl_rng_set(gptime.rng, gptime.seed);

  // simulation variables
  double avg_sim_endtime = 0;
  int count_detect = 0;
  double current_time;
  int curr_observation;
  double rand_next_time;
  int count_extinct = 0;

  // Beginning of simulation that has "sim" number of runs.
  for (int sim = 1; sim <= gptime.num_sims; sim++)
  {
    // initialize time to zero
    current_time = gptime.start_time;

    // initialize current obseration output to 0
    curr_observation = 0;

    // Define TDCloneList population and initialize variables;
    TDCloneList population;
    population.init();

    // Determine Advance function class to use based on the parameters
    if (gptime.is_custom_model)
    {
      NewTDClone = new TDCloneList::NewCloneCustom(population);
    }
    else if( punct_params.is_punctuated )
    {
      NewTDClone = new TDCloneList::NewClonePunct(population, fit_params, mut_params, punct_params);
    }
    else if( fit_params.is_randfitness || mut_params.is_mutator )
    {
      if ( epi_params.is_epistasis )
      {
        NewTDClone = new TDCloneList::NewCloneEpi(population, fit_params, mut_params, epi_params);
      }
      else
      {
        NewTDClone = new TDCloneList::NewCloneFitMut(population, fit_params, mut_params);
      }
    }
    else
    {
      NewTDClone = new TDCloneList::NewCloneNoParams(population);
    }


    if( ancestor_file.isNull() ) // if no ancestor file exists
    {
      /*
        total rate for time-homogeneous population should be the max for the
        birth and death functions times their associated birth rates times
        the number of ancestors
      */
      gsl_function B_max;
      gsl_function D_max;
      B_max.function = rate_function_array[td_birth_params.type];
      D_max.function = rate_function_array[td_death_params.type];
      B_max.params = &td_birth_params;
      D_max.params = &td_death_params;

      // crude global max finder just to give us a homogeneous rate value
      td_birth_params.homogeneous_rate = MaximizeRate(B_max, gptime.start_time, gptime.tot_life, 1000);
      td_death_params.homogeneous_rate = MaximizeRate(D_max, gptime.start_time, gptime.tot_life, 1000);
      Rcpp::Rcout << "default max birth rate:\t" << td_birth_params.homogeneous_rate << "\t" <<
                   "default max death rate:\t" << td_death_params.homogeneous_rate << "\n";


      population.tot_rate_homog = (td_birth_params.homogeneous_rate +
        td_death_params.homogeneous_rate) * gptime.ancestors * gptime.ancestor_clones;
      population.tot_cell_count = gptime.ancestors * gptime.ancestor_clones;

      // go through all ancestors and initialize clones for each
      for(int ance_clone_count = 1; ance_clone_count <= gptime.ancestor_clones; ance_clone_count++)
      {
        // Create Ancestor Node
        struct clone* ancestor;
        ancestor = new struct clone;

        ancestor->cell_count = gptime.ancestors;
        ancestor->allele_count = gptime.ancestors;
        ancestor->mut_prob = gptime.mutation_prob;
        ancestor->clone_time = current_time;
        ancestor->subclone_count = 0;
        ancestor->mut_count = 0;
        ancestor->driver_count = 0;
        ancestor->is_driver = false;

        // input default types
        ancestor->birth_params = td_birth_params;
        ancestor->death_params = td_death_params;
        (ancestor->B).function = rate_function_array[td_birth_params.type];
        (ancestor->D).function = rate_function_array[td_death_params.type];
        (ancestor->B).params = &(ancestor->birth_params);
        (ancestor->D).params = &(ancestor->death_params);

        // initializing the rate accumulation values
        ancestor->birth_rate = 0;
        ancestor->death_rate = 0;

        population.InsertAncestor(ancestor);
      }
    }
    else // ancestor file exists to read from
    {
      Rcpp::Rcout << "Reading ancestor file...\n";

      population.tot_rate_homog = 0;
      population.tot_rate = 0;
      population.tot_cell_count = 0;

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
            ancestor->cell_count = !ancestor_map["numcells"].empty() ? stoi(ancestor_map["numcells"]) : 0;
            ancestor->allele_count = ancestor->cell_count;
            ancestor->mut_prob = !ancestor_map["mutprob"].empty() ? stof(ancestor_map["mutprob"]) : gptime.mutation_prob;
            ancestor->clone_time = !ancestor_map["initialtime"].empty() ? stof(ancestor_map["initialtime"]) : current_time;
            ancestor->subclone_count = !ancestor_map["subclone_count"].empty() ? stoi(ancestor_map["subclone_count"]) : 0;
            ancestor->mut_count = !ancestor_map["num_mut"].empty() ? stoi(ancestor_map["num_mut"]) : 0;
            ancestor->driver_count = !ancestor_map["num_drivers"].empty() ? stoi(ancestor_map["num_drivers"]) : 0;
            ancestor->is_driver = !ancestor_map["is_driver"].empty() ? stoi(ancestor_map["is_driver"]) : false;

            ancestor->birth_params.type = !ancestor_map["birth_function"].empty() ? stoi(ancestor_map["birth_function"]) : td_birth_params.type;
            ancestor->death_params.type = !ancestor_map["death_function"].empty() ? stoi(ancestor_map["death_function"]) : td_death_params.type;

            if(!ancestor_map["bf_params"].empty()) // If the ancestor file contains a parameter list, parse and add to function
            {
              std::string bf_input = ancestor_map["bf_params"];
              std::istringstream ss(bf_input);
              std::string token;
              double value;
              while(getline(ss, token, ','))
              {
                value = stof(token);
                ancestor->birth_params.coefs.push_back(value);
              }
              // Add params to birth function
              (ancestor->B).params = &(ancestor->birth_params);

            }
            else // Use defaults if not present
            {
              ancestor->birth_params = td_birth_params;
              (ancestor->B).params = &(ancestor->birth_params);
            }
            (ancestor->B).function = rate_function_array[ancestor->birth_params.type];


            if(!ancestor_map["df_params"].empty())
            {
              std::string df_input = ancestor_map["df_params"];
              std::istringstream ss(df_input);
              std::string token;
              double value;
              while(getline(ss, token, ','))
              {
                value = stof(token);
                ancestor->death_params.coefs.push_back(value);
              }
              (ancestor->D).params = &(ancestor->death_params);
            }
            else
            {
              ancestor->death_params = td_death_params;
              (ancestor->D).params = &(ancestor->death_params);
            }
            (ancestor->D).function = rate_function_array[ancestor->death_params.type];

            // crude global max finder just to give us a homogeneous rate value
            ancestor->birth_params.homogeneous_rate = MaximizeRate((ancestor->B), gptime.start_time, gptime.tot_life, 1000);
            //(ancestor->B).params = &(ancestor->birth_params);
            ancestor->death_params.homogeneous_rate = MaximizeRate((ancestor->D), gptime.start_time, gptime.tot_life, 1000);
            //(ancestor->D).params = &(ancestor->death_params);

            Rcpp::Rcout << "ancestor " << ancestor->clone_id << " max birth rate:\t" << ancestor->birth_params.homogeneous_rate << "\t" <<
                         "ancestor " << ancestor->clone_id << " max death rate:\t" << ancestor->death_params.homogeneous_rate << "\n";


            ancestor->birth_rate = ancestor->cell_count * ancestor->birth_params.homogeneous_rate;
            ancestor->death_rate = ancestor->cell_count * ancestor->death_params.homogeneous_rate;

            population.InsertAncestor(ancestor);

            population.tot_rate_homog = population.tot_rate_homog +
                (ancestor->birth_params.homogeneous_rate + ancestor->death_params.homogeneous_rate) *
                ancestor->cell_count;
            population.tot_cell_count = population.tot_cell_count + ancestor->cell_count;

            ancestor_map.clear();
        }
      }
      Rcpp::Rcout << "Ancestor File Read...";
    }

    Rcpp::Rcout << "Output Ancestor Population...";
    population.Traverse(timedata, sim, current_time, gptime.trace_ancestry, gptime.count_alleles);
    Rcpp::Rcout << "Initial Traverse Done\n";

    // Begin single simulation with while loop that exists when hit max time, max pop, or extinction
    while ( (population.tot_cell_count < gptime.max_pop) &&
            (population.tot_cell_count > 0) &&
            (current_time < gptime.tot_life) )
    {
      Rcpp::checkUserInterrupt();
      //Rcpp::Rcout << current_time << "\n";
      // Get next event time by advancing with adaptive thinning
      rand_next_time = population.AdvanceTime(current_time);

      // Advance Simulation State (choose next event)
      population.AdvanceState(current_time, rand_next_time);

      // update current_time
      current_time = current_time + rand_next_time;

      // Method to output data at designated observation times
      while(current_time > observation_times[curr_observation])
      {
        if( (current_time < observation_times[curr_observation + 1]) ||
            (observation_times.size() == curr_observation + 1) )
        {
          population.Traverse(timedata, sim, observation_times[curr_observation], gptime.trace_ancestry, gptime.count_alleles);
          if((observation_times.size() == curr_observation + 1))
          {
            break;
          }
        }
        curr_observation++;
      }

      /*
      if((population.tot_cell_count % 50000) == 0)
      {
        Rcpp::Rcout << "Time: " << current_time << "\t size: " << population.tot_cell_count << "\n";
      }
      //*/
    }

    // nonextinction checker - if set to false and goes extinction, restart that sim
    if( (population.tot_cell_count == 0) && (gptime.allow_extinction == false) )
    {
      population.DeleteList();
      count_extinct++;
      gptime.num_sims++; // increase number of sims
      Rcpp::Rcout << "Population went extinct. Restarting.\n";
      continue;
    }

    // If simulation made it to the gptime.max_pop in the alotted time
    if( (population.tot_cell_count >= gptime.max_pop) &&
        (current_time < gptime.tot_life) )
    {
      count_detect = count_detect + 1;
      avg_sim_endtime = avg_sim_endtime + (current_time) / (double)gptime.num_sims;
    }

    // Final Timed Output
    population.Traverse(timedata, sim, current_time, gptime.trace_ancestry, gptime.count_alleles);
    // Sampling from population
    if(gptime.sample_size > 0 & gptime.num_samples > 0)
    {
      population.SampleAndTraverse(sample_data, sim, gptime.sample_size, gptime.num_samples);
    }
    // Trim tree if threshold is higher. Otherwise, Traverse
    population.TreeTrim(gptime.detection_threshold, gptime.max_pop);
    // Output of end state with clone info
    Rcpp::Rcout << "Traversing and outputting run " << sim << "\n";
    population.Traverse(clonedata, sim, gptime.count_alleles);
    Rcpp::Rcout << "Traversal Done\n";

    population.DeleteList();

  }

  //avg_sim_endtime = avg_sim_endtime * (double)num_sims / (double)count_detect;

  gsl_rng_free(gptime.rng);
  delete NewTDClone;

  sim_stats << "avg_sim_endtime, " << avg_sim_endtime << "\n" <<
               "count_detect, " << count_detect << "\n" <<
               "count_extinct, " << count_extinct  << "\n" <<
               "tot_integration_error, " << gptime.tot_error;

  clonedata.close();
  timedata.close();
  sample_data.close();
  sim_stats.close();

  return 0;
}
