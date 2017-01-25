#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <vector>
#include <cstring>
#include <string>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <dlfcn.h>


#include "constantGlobalStructs.h"
#include "constantCloneList.h"
#include "constantParameterList.h"

// structure contains all global parameters used in multiple source files
GlobalParameters gpcons;
gsl_rng* constant_rng;
// Function class ptr defined in main() but used in clonelist.cpp
ConstantCloneList::NewCloneFunction* NewConstantClone;
void (*ConstantGenerateFitness)(double*, struct FitnessParameters*, gsl_rng*);
void *lib_handle;


//' siapopConstant
//'
//' SIApop for time-homogeneous populations. Simulates a
//' Birth-Death-Mutation process simulation for infinite-allele
//' model with random fitness contributions using the Gillespie Algorithm and
//' outputs data to the location of \code{output_dir}.
//'
//' The infinite-allele birth death process assumes anytime a mutation that
//' occurs, it is unique and not seen in the population. \code{siapopConstant}
//' can generates an infinite-allele birth-death time homogeneous process along
//' with additional scenarios that affect new clones.
//'
//' A fitness distribution is provided so that the birth rate of new clones is
//' additive with the additional fitness coming from this distribution. The
//' mutation distribution has a similar effect on the mutation rate of a new
//' clone. The punctuated parameters assume a new clone may arise that contains
//' multiple new alleles instead of a single new allele. These clones have
//' different fitnesses as well and a Poisson number of new alleles. The
//' epistatic parameters refer to a simple epistatic model where a change in
//' fitness occurs in a clone when it reaches \code{epistatic_mutation_thresh}
//' total mutations.
//'
//' The process can be run with defaults that result in a simple model without
//' mutations. The defaults for most parameters that would affect fitness of
//' new clones are set such that they have no effect on the process.
//'
//' Simulations are output as text files and input can be in the form of text
//' files or a comma-delimeted input file. Currently input for the ancestors
//' requires a text file.
//'
//' @param tot_life total lifetime to run a simulation for
//' @param max_pop maximum population to stop simulation
//' @param start_time time index to begin each simulation
//' @param ancestors number of ancestors in a clone to initialize simulation
//' @param ancestor_clones number of ancestor clones each containing
//'   \code{ancestors} individuals to initialize simulation with
//' @param num_sims number of simulations to run
//' @param allow_extinction if TRUE then each simulation restarts when
//'   extinction occurs. The run counter is incremented and the data is still
//'   recorded in \emph{timedata.txt}
//' @param is_custom_model (not implemented yet) allows user to provide a shared
//'   object file to advance the clone
//' @param num_samples number of single cell samples to take from each
//'   simulation
//' @param sample_size size of each sample of single cells
//' @param detection_threshold minimum threshold to report clones. If a clone is
//'   below minimum its number is added to its parent count.
//' @param observation_frequency the frequency of time to output the current
//'   population to \emph{timedata.txt}
//' @param observation_times a vector of specific time points to output the
//'   current population to \emph{timedata.txt}
//' @param birth_rate ancestor birth rate
//' @param death_rate ancestor birth rate
//' @param mutation_prob ancestor mutation probability, probability that a
//' daughter is a new mutant allele
//' @param distribution_function one of "doubleexp", "normal", or "uniform"
//' @param alpha_fitness fitness distribution (right-side) rate parameter. When
//'   a new clone arises, the fitness of the new clone is a double exponential
//'   with the positive side having rate \code{alpha}
//' @param beta_fitness fitness distribution (left-side) rate parameter. When a
//'   new clone arises, the fitness of the new clone is a double exponential
//'   with the negative side having rate \code{beta}
//' @param pass_prob probability of no change in fitness
//' @param upper_fitness upper bound to fitness distribution
//' @param lower_fitness lower bound to fitness distribution
//' @param alpha_mutation mutation distribution alpha parameter. A new clone
//'   can have a new mutation rate coming from a distribution with parameters
//'   \code{alpha} and \code{beta}.
//' @param beta_mutation mutation distribution beta parameter. A new clone can
//'   have a new mutation rate coming from a distribution with parameters
//'   \code{alpha} and \code{beta}
//' @param trace_ancestry if \code{TRUE} then \emph{clonedata.txt} reports the
//'   parent clone id and time
//' @param count_alleles if \code{TRUE} then reports the number of cells with
//'   allele equal to the final id of the clone. For example, a single ancestor
//'   should have an allele frequency equal to the sum of the current living
//'   population.
//' @param punctuated_prob probability that a new clone has a burst of multiple
//'   mutations
//' @param poisson_param parameter of a poisson distribution used to determine
//'   the number of mutations in a new clone given a burst occurs
//' @param punctuated_multiplier = given a burst occurs in the new clone, the
//'   fitness is multiplied by this amount to affect the fitness distribution.
//' @param punctuated_advantageous_prob the probability that a punctuated clone
//'   has a positive fitness
//' @param epistatic_mutation_thresh = number of mutations before an epistatic
//'   event occurs affecting the fitness
//' @param epistatic_multiplier = given an epistatic event occurs, the fitness
//'   is multiplied by this amount
//' @param input_file input character vector of input file
//' @param output_dir input character vector of output location
//' @param ancestor_file input character vector of ancestor file
//' @examples
//' \dontrun{
//' # Use default values
//' siapopConstant()
//' siapopConstant(outputdir = "./")
//' siapopConstant(ancestor_file = "./ancestors.txt")
//' siapopConstant(tot_life = 10, max_pop = 1000, birth_rate = 1.1,
//'                death_rate = 0.99, mutation_prob = 0.01,
//'                allow_extinction = FALSE, num_sims = 1, num_samples = 1,
//'                alpha_fitness = 100, beta_fitness = 100,
//'                sample_size = 100, observation_times = c(1, 5, 10))
//' }
//' @export
// [[Rcpp::export]]
int siapopConstant(double tot_life = 40000.0,
                   int max_pop = 10000,
                   double start_time = 0.0,
                   int ancestors = 1,
                   int ancestor_clones = 1,
                   int num_sims = 1,
                   bool allow_extinction = true,
                   bool is_custom_model = false,
                   int num_samples = 0,
                   int sample_size = 0,
                   double detection_threshold = 0.0,
                   double observation_frequency = 0.0,
                   SEXP observation_times = R_NilValue,
                   double birth_rate = 1.5,
                   double death_rate = 1.0,
                   double mutation_prob = 0.0,
                   SEXP fitness_distribution = R_NilValue,
                   double alpha_fitness = 0.0,
                   double beta_fitness = 0.0,
                   double pass_prob = 1.0,
                   double upper_fitness = 0.0,
                   double lower_fitness = 0.0,
                   double alpha_mutation = 0.0,
                   double beta_mutation = 0.0,
                   bool trace_ancestry = true,
                   bool count_alleles = true,
                   double punctuated_prob = 0.0,
                   double poisson_param = 1.0,
                   double punctuated_multiplier = 1.0,
                   double punctuated_advantageous_prob = 1.0,
                   double epistatic_mutation_thresh = 1.0,
                   double epistatic_multiplier = 1.0,
                   SEXP seed = R_NilValue,
                   SEXP input_file = R_NilValue,
                   SEXP output_dir = R_NilValue,
                   SEXP ancestor_file = R_NilValue)
{
  // for now, hardcode is_custom_model to false;
  is_custom_model = false;

  //  declaring random number generator and setting seed
  constant_rng = gsl_rng_alloc(gsl_rng_mt19937);
  if( Rf_isNull(seed) )
  {
    gpcons.seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  }
  else
  {
    gpcons.seed = Rf_asInteger(seed);
  }
  /*
    VARIABLE INPUT AND CONVERSION
  */
  // parsing arguments in command line by searching for argument options
  // default output is current directory
  const char *output_folder;
  if( Rf_isNull(output_dir) )
  {
    output_folder = "./";
  }
  else
  {
    output_folder = CHAR(Rf_asChar(output_dir));
  }


  // declare and open output stream for simulation statistics
  char fn[100];
  std::ofstream sim_stats;
  sprintf(fn,"%s/sim_stats.txt", output_folder);
  sim_stats.open(fn);
  sim_stats.setf(std::ios::fixed);
  sim_stats.precision(8);

  FitnessParameters fit_params;
  MutationParameters mut_params;
  PunctuationParameters punct_params;
  EpistaticParameters epi_params;
  std::vector<double> observation_times_;

  // if the input_file path exists, use this to import variables
  // parsing through the input file and converting/adding to parameter list
  if( !Rf_isNull(input_file) )
  {
    // declare and initialize parameter list for simulation
    ConstantParameterList params;
    params.init();

    const char* input_params = CHAR(Rf_asChar(input_file));

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

    for(std::map<std::string, std::string>::iterator it=params.begin(); it!=params.end(); ++it)
    {
      sim_stats << it->first << ", " << it->second << "\n";
    }

    // convert all parameters imported from file into respective values in gpcons
    params.convert("tot_life", gpcons.tot_life);
    params.convert("max_pop", gpcons.max_pop);
    params.convert("start_time", gpcons.start_time);
    params.convert("ancestors", gpcons.ancestors);
    params.convert("ancestor_clones", gpcons.ancestor_clones);
    params.convert("num_sims", gpcons.num_sims);
    params.convert("num_samples", gpcons.num_samples);
    params.convert("sample_size", gpcons.sample_size);
    params.convert("detection_threshold", gpcons.detection_threshold);
    params.convert("observation_frequency", gpcons.observation_frequency);
    if (gpcons.observation_frequency == 0) gpcons.observation_frequency = gpcons.tot_life;

    params.ParseVector("observation_times", observation_times_);
    if (observation_times_.size() == 0)
    {
      int num_obs =  ceil(gpcons.tot_life / gpcons.observation_frequency);
      for (int i = 1; i <= num_obs; i++)
      {
        observation_times_.push_back(i * gpcons.observation_frequency);
      }
    }
    if(observation_times_.back() != gpcons.tot_life) observation_times_.push_back(gpcons.tot_life);


    params.convert("allow_extinction", gpcons.allow_extinction);
    params.convert("trace_ancestry", gpcons.trace_ancestry);
    params.convert("count_alleles", gpcons.count_alleles);
    params.convert("birth_rate", gpcons.birth_rate);
    params.convert("death_rate", gpcons.death_rate);
    params.convert("mutation_prob", gpcons.mutation_prob);
    params.convert("is_custom_model", gpcons.is_custom_model);

    std::string fitness_distribution;
    params.convert("fitness_distribution", fit_params.fitness_distribution);
    params.convert("alpha_fitness", fit_params.alpha_fitness);
    params.convert("beta_fitness", fit_params.beta_fitness);
    params.convert("pass_prob", fit_params.pass_prob);
    params.convert("upper_fitness", fit_params.upper_fitness);
    params.convert("lower_fitness", fit_params.lower_fitness);
    fit_params.is_randfitness = false;
    if( !fit_params.fitness_distribution.empty() )
    {
      fit_params.is_randfitness = true;
      if( fit_params.fitness_distribution.compare("doubleexp") == 0 )
      {
        Rcpp::Rcout << "Double Exponential Fitness Distribution\n";
        ConstantGenerateFitness = &cdoubleexp;
      }
      else if( fit_params.fitness_distribution.compare("normal") == 0 )
      {
        Rcpp::Rcout << "Normal Fitness Distribution\n";
        ConstantGenerateFitness = &cnormal;
      }
      else if( fit_params.fitness_distribution.compare("uniform") == 0 )
      {
        Rcpp::Rcout << "Uniform Fitness Distribution\n";
        ConstantGenerateFitness = &cuniform;
      }
      else
      {
        Rcpp::stop("Not a valid fitness distribution");
      }
    }

    params.convert("alpha_mutation", mut_params.alpha_mutation);
    params.convert("beta_mutation", mut_params.beta_mutation);
    params.convert("pass_prob", mut_params.pass_prob);
    mut_params.is_mutator = false;
    if( (mut_params.alpha_mutation > 0 && mut_params.beta_mutation > 0) &&
        (mut_params.pass_prob < 1) )
    {
      mut_params.is_mutator = true;
    }

    params.convert("punctuated_prob", punct_params.punctuated_prob);
    params.convert("poisson_param", punct_params.poisson_param);
    params.convert("punctuated_multiplier", punct_params.punctuated_multiplier);
    params.convert("punctuated_advantageous_prob", punct_params.punctuated_advantageous_prob);
    punct_params.is_punctuated = punct_params.punctuated_prob > 0;


    params.convert("epistatic_mutation_thresh", epi_params.epistatic_mutation_thresh);
    params.convert("epistatic_multiplier", epi_params.epistatic_multiplier);
    epi_params.is_epistasis = epi_params.epistatic_mutation_thresh > 0;
    if (epi_params.epistatic_multiplier == 1.0)
    {
      epi_params.is_epistasis = false;
    }
    /*
      END OF VARIABLE INPUT AND CONVERSION FROM input_file
    */
  }
  else
  {
    // convert all parameters imported from file into respective values in gpcons
    gpcons.tot_life = tot_life;
    gpcons.max_pop = max_pop;
    gpcons.start_time = start_time;
    gpcons.ancestors = ancestors;
    gpcons.ancestor_clones = ancestor_clones;
    gpcons.num_sims = num_sims;
    gpcons.num_samples = num_samples;
    gpcons.sample_size = sample_size;
    gpcons.detection_threshold = detection_threshold;
    gpcons.observation_frequency = observation_frequency;
    if ( observation_frequency == 0)
    {
      gpcons.observation_frequency = gpcons.tot_life;
    }

    if ( Rf_isNull(observation_times) )
    {
      int num_obs =  ceil(gpcons.tot_life / gpcons.observation_frequency);
      for (int i = 1; i <= num_obs; i++)
      {
        observation_times_.push_back(i * gpcons.observation_frequency);
      }
    }
    else
    {
      int obs_time_length = Rf_length(observation_times);
      observation_times_.resize(obs_time_length);

      observation_times = Rf_coerceVector(observation_times, REALSXP);
      for(int iter = 0; iter < obs_time_length; ++iter)
      {
        observation_times_[iter] = REAL(observation_times)[iter];
      }
    }
    if(observation_times_.back() != gpcons.tot_life)
    {
      observation_times_.push_back(gpcons.tot_life);
    }

    gpcons.allow_extinction = allow_extinction;
    gpcons.trace_ancestry = trace_ancestry;
    gpcons.count_alleles = count_alleles;
    gpcons.birth_rate = birth_rate;
    gpcons.death_rate = death_rate;
    gpcons.mutation_prob = mutation_prob;
    gpcons.is_custom_model = is_custom_model;

    fit_params.alpha_fitness = alpha_fitness;
    fit_params.beta_fitness = beta_fitness;
    fit_params.pass_prob = pass_prob;
    fit_params.upper_fitness = upper_fitness;
    fit_params.lower_fitness = lower_fitness;
    fit_params.is_randfitness = false;
    if( !Rf_isNull(fitness_distribution) )
    {
      fit_params.is_randfitness = true;
      std::string fitness_distribution_ = CHAR(STRING_ELT(fitness_distribution, 0));
      fit_params.fitness_distribution = fitness_distribution_;
      if( fit_params.fitness_distribution.compare("doubleexp") == 0 )
      {
        Rcpp::Rcout << "Double Exponential Fitness Distribution\n";
        ConstantGenerateFitness = &cdoubleexp;
      }
      else if( fit_params.fitness_distribution.compare("normal") == 0 )
      {
        Rcpp::Rcout << "Normal Fitness Distribution\n";
        ConstantGenerateFitness = &cnormal;
      }
      else if( fit_params.fitness_distribution.compare("uniform") == 0 )
      {
        Rcpp::Rcout << "Uniform Fitness Distribution\n";
        ConstantGenerateFitness = &cuniform;
      }
      else if( fit_params.fitness_distribution.compare("custom") == 0 )
      {
        lib_handle = dlopen("/Users/mcdonald/Desktop/testplugin/plugin.so", RTLD_LAZY);
        ConstantGenerateFitness = (void (*)(double *, struct FitnessParameters*, gsl_rng*))dlsym(lib_handle, "custom");
      }
      else
      {
        Rcpp::stop("Not a valid fitness distribution");
      }
    }

    mut_params.alpha_mutation = alpha_mutation;
    mut_params.beta_mutation = beta_mutation;
    mut_params.pass_prob = pass_prob;
    mut_params.is_mutator = false;
    if( (mut_params.alpha_mutation > 0 && mut_params.beta_mutation > 0) &&
        (mut_params.pass_prob < 1) )
    {
      mut_params.is_mutator = true;
    }

    punct_params.punctuated_prob = punctuated_prob;
    punct_params.poisson_param = poisson_param;
    punct_params.punctuated_multiplier = punctuated_multiplier;
    punct_params.punctuated_advantageous_prob = punctuated_advantageous_prob;
    punct_params.is_punctuated = punct_params.punctuated_prob > 0;

    epi_params.epistatic_mutation_thresh = epistatic_mutation_thresh;
    epi_params.epistatic_multiplier = epistatic_multiplier;
    epi_params.is_epistasis = epi_params.epistatic_mutation_thresh > 0;
    if (epi_params.epistatic_multiplier == 1.0)
    {
      epi_params.is_epistasis = false;
    }
    /*
      END OF VARIABLE INPUT AND CONVERSION
    */
  }

  // declare and open other output streams for time and end of sim clone list
  std::ofstream clonedata;
  sprintf(fn, "%s/clonedata.txt", output_folder);
  clonedata.open(fn);
  clonedata.setf(std::ios::fixed);
  clonedata.precision(12);
  if(gpcons.count_alleles)
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
  if(gpcons.trace_ancestry)
  {
    if(gpcons.count_alleles)
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
    if(gpcons.count_alleles)
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
  if( (gpcons.sample_size > 0) & (gpcons.num_samples > 0) )
  {
    sprintf(fn, "%s/sampledata.txt", output_folder);
    sample_data.open(fn);
    sample_data.setf(std::ios::fixed);
    sample_data << "run\tsample_number\tunique_id\tnumber_obs\n";
  }

  // set RNG seed
  gsl_rng_set(constant_rng, gpcons.seed);

  // simulation variables
  double avg_sim_endtime = 0;
  int count_detect = 0;
  double current_time;
  unsigned int curr_observation;
  double rand_next_time;
  int count_extinct = 0;

  // Beginning of simulation that has "sim" number of runs
  for (int sim = 1; sim <= gpcons.num_sims; sim++)
  {
    // initialize time to zero
    current_time = gpcons.start_time;

    // initialize current obseration output to 0
    curr_observation = 0;

    // Define ConstantCloneList population and initialize variables;
    ConstantCloneList population;
    population.init();

    // Determine Advance function class to use based on the parameters
    if (gpcons.is_custom_model)
    {
      Rcpp::Rcout << "Custom model\n";
      NewConstantClone = new ConstantCloneList::NewCloneCustom(population, constant_rng);
    }
    else if( punct_params.is_punctuated )
    {
      Rcpp::Rcout << "Punctuated model\n";
      NewConstantClone = new ConstantCloneList::NewClonePunct(population, fit_params, mut_params, punct_params, constant_rng);
    }
    else if( fit_params.is_randfitness || mut_params.is_mutator )
    {
      if ( epi_params.is_epistasis )
      {
        Rcpp::Rcout << "Epistatic model\n";
        NewConstantClone = new ConstantCloneList::NewCloneEpi(population, fit_params, mut_params, epi_params, constant_rng);
      }
      else
      {
        Rcpp::Rcout << "Fitness model\n";
        NewConstantClone = new ConstantCloneList::NewCloneFitMut(population, fit_params, mut_params, constant_rng);
      }
    }
    else
    {
      Rcpp::Rcout << "No parameters model\n";
      NewConstantClone = new ConstantCloneList::NewCloneNoParams(population);
    }

    if( Rf_isNull(ancestor_file) ) // if no ancestor file exists
    {
      // total rate for SSA is equal to number of individuals alive times
      // respective rates
      population.tot_rate = (gpcons.birth_rate + gpcons.death_rate) * gpcons.ancestors * gpcons.ancestor_clones;
      population.tot_cell_count = gpcons.ancestors * gpcons.ancestor_clones;

      // go through all ancestors and initialize clones for each
      for(int ance__clone_count = 1; ance__clone_count <= gpcons.ancestor_clones; ance__clone_count++)
      {
        // Create Ancestor Node
        struct clone* ancestor;
        ancestor = new struct clone;

        ancestor->cell_count = gpcons.ancestors;
        ancestor->allele_count = gpcons.ancestors;
        ancestor->birth_rate = gpcons.birth_rate;
        ancestor->death_rate = gpcons.death_rate;
        ancestor->mut_prob = gpcons.mutation_prob;
        ancestor->clone_time = current_time;
        ancestor->subclone_count = 0;
        ancestor->mut_count = 0;
        ancestor->driver_count = 0;
        ancestor->is_driver = false;

        population.InsertAncestor(ancestor);
      }
    }
    else // ancestor file exists to read from
    {
      Rcpp::Rcout << "Reading ancestor file...";

      population.tot_rate = 0;
      population.tot_cell_count = 0;

      // import first line of file to a vector of keys
      std::vector<std::string> ancestor_keys;
      std::vector<std::string>::iterator it;

      const char *ancestors = CHAR(Rf_asChar(ancestor_file));
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
            ancestor->birth_rate = !ancestor_map["birthrate"].empty() ? stof(ancestor_map["birthrate"]) : gpcons.birth_rate;
            ancestor->death_rate = !ancestor_map["deathrate"].empty() ? stof(ancestor_map["deathrate"]) : gpcons.death_rate;
            ancestor->mut_prob = !ancestor_map["mutprob"].empty() ? stof(ancestor_map["mutprob"]) : gpcons.mutation_prob;
            ancestor->clone_time = !ancestor_map["initialtime"].empty() ? stof(ancestor_map["initialtime"]) : current_time;
            ancestor->subclone_count = !ancestor_map["subclone_count"].empty() ? stoi(ancestor_map["subclone_count"]) : 0;
            ancestor->mut_count = !ancestor_map["num_mut"].empty() ? stoi(ancestor_map["num_mut"]) : 0;
            ancestor->driver_count = !ancestor_map["num_drivers"].empty() ? stoi(ancestor_map["num_drivers"]) : 0;
            ancestor->is_driver = !ancestor_map["is_driver"].empty() ? stoi(ancestor_map["is_driver"]) : false;

            population.InsertAncestor(ancestor);

            population.tot_rate = population.tot_rate + (ancestor->birth_rate + ancestor->death_rate) * ancestor->cell_count;
            population.tot_cell_count = population.tot_cell_count + ancestor->cell_count;

            ancestor_map.clear();
        }
      }
      Rcpp::Rcout << "Ancestor File Read...";
    }

    Rcpp::Rcout << "Output Ancestor Population...";
    population.Traverse(timedata, sim, current_time, gpcons.trace_ancestry, gpcons.count_alleles);
    Rcpp::Rcout << "Ancestor Output Written...\n";

    // Begin single simulation with while loop that exists when hit max time, max pop, or extinction
    while ( (population.tot_cell_count < gpcons.max_pop) &&
            (population.tot_cell_count > 0) &&
            (current_time < gpcons.tot_life) )
    {
      Rcpp::checkUserInterrupt();
      // Advance Simulation Time (choose next event time)
      rand_next_time = population.AdvanceTime(current_time, constant_rng);

      // Advance Simulation State (choose next event)
      population.AdvanceState(current_time, rand_next_time, constant_rng);

      // update current_time
      current_time = current_time + rand_next_time;

      // Method to output data at designated observation times
      while(current_time > observation_times_[curr_observation])
      {
        if( (current_time < observation_times_[curr_observation + 1]) ||
            (observation_times_.size() == curr_observation + 1) )
        {
          population.Traverse(timedata, sim, observation_times_[curr_observation], gpcons.trace_ancestry, gpcons.count_alleles);
          if((observation_times_.size() == curr_observation + 1) )
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
    if( (population.tot_cell_count == 0) && (gpcons.allow_extinction == false) )
    {
      population.DeleteList();
      count_extinct++;
      gpcons.num_sims++; // increase number of sims
      Rcpp::Rcout << "Population went extinct. Restarting.\n";
      continue;
    }

    // If simulation made it to the max_pop in the alotted time
    if( (population.tot_cell_count >= gpcons.max_pop) &&
        (current_time < gpcons.tot_life) )
    {
      count_detect = count_detect + 1;
      avg_sim_endtime = avg_sim_endtime + (current_time) / (double)gpcons.num_sims;
    }

    // Final Timed Output
    population.Traverse(timedata, sim, current_time, gpcons.trace_ancestry, gpcons.count_alleles);
    // Sampling from population
    if( (gpcons.sample_size > 0) & (gpcons.num_samples > 0) )
    {
      population.SampleAndTraverse(sample_data, sim, gpcons.sample_size, gpcons.num_samples, constant_rng);
    }
    // Trim tree if threshold is higher. Otherwise, Traverse
    population.TreeTrim(gpcons.detection_threshold, gpcons.max_pop);
    // Output of end state with clone info
    Rcpp::Rcout << "Traversing and outputting run " << sim << "\n";
    population.Traverse(clonedata, sim, gpcons.count_alleles);
    Rcpp::Rcout << "Traversal Done\n";

    population.DeleteList();

  }

  //avg_sim_endtime = avg_sim_endtime * (double)gpcons.num_sims / (double)count_detect;

  gsl_rng_free(constant_rng);
  delete NewConstantClone;

  sim_stats << "avg_sim_endtime, " << avg_sim_endtime << "\n" <<
               "count_detect, " << count_detect << "\n" <<
               "count_extinct, " << count_extinct  << "\n";

  clonedata.close();
  timedata.close();
  sample_data.close();
  sim_stats.close();
  dlclose(lib_handle);



  return 0;
}
