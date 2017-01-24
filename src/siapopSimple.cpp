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

#include "simpleGlobalStructs.h"
#include "simpleCloneList.h"
#include "simpleParameterList.h"

// structure contains all global parameters used in multiple source files
GlobalParameters gpsimp;
gsl_rng* simple_rng;


//' siapopSimple
//'
//' SIApop for non-mutating processes. Runs an exact process by simulating
//' binomial and negative binomial random values at each time step. The process
//' is only ran for a total length of time since the simulation is exact during
//' that time instead of running until a specific population size.
//'
//' Simple birth-death processes can be simulated exactly by generating a
//' binomial random variable to determine the number of ancestors of a specific
//' clone that give rise to descendant, then a negative binomial random
//' variable to determine the total number of descendant that arise out of those
//' replicating ancestors. This can be done for each ancestor clone with
//' various rates. Mutation is not permitted in this scheme.
//'
//' Simulations are output as text files and input can be in the form of text
//' files or a comma-delimeted input file. Currently input for the ancestors
//' requires a text file.
//'
//' @param tot_life total lifetime to run a simulation for
//' @param ancestors number of ancestors in a clone to initialize simulation
//' @param ancestor_clones number of ancestor clones each containing
//'   \code{ancestors} individuals to initialize simulation with
//' @param num_sims number of simulations to run
//' @param allow_extinction if TRUE then each simulation restarts when
//'   extinction occurs. The run counter is incremented and the data is still
//'     recorded in \emph{timedata.txt}
//' @param num_samples number of single cell samples to take from each
//'   simulation
//' @param sample_size size of each sample of single cells
//' @param detection_threshold minimum threshold to report clones. If a clone is
//'   below minimum its number is added to its parent count.
//' @param input input character vector of input file
//' @param output_dir input character vector of output location
//' @param ancestor_file input character vector of ancestor file
//' @param birth_rate ancestor birth rate
//' @param death_rate ancestor birth rate
//' @examples
//' \dontrun{
//' # Use default values
//' siapopSimple()
//' siapopSimple(input_file = "./input.txt, "output_dir = "./",
//'              ancestor_file = "ancestors.txt")
//' siapopSimple(tot_life = 5, ancestors = 10, ancestor_clones = 5,
//'              num_sims = 1, allow_extinction = FALSE,
//'              detection_threshold = 0, birth_rate = 1.1, death_rate = 1)
//' }
//' @export
// [[Rcpp::export]]
int siapopSimple(double tot_life = 40000.0,
                 int ancestors = 1,
                 int ancestor_clones = 1,
                 int num_sims = 1,
                 int num_samples = 0,
                 int sample_size = 0,
                 bool allow_extinction = true,
                 double detection_threshold = 0.0,
                 double birth_rate = 1.5,
                 double death_rate = 1.0,
                 SEXP seed = R_NilValue,
                 SEXP input_file = R_NilValue,
                 SEXP output_dir = R_NilValue,
                 SEXP ancestor_file = R_NilValue)
{

  //  declaring random number generator and setting seed
  simple_rng = gsl_rng_alloc(gsl_rng_mt19937);
  if( Rf_isNull(seed) )
  {
    gpsimp.seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  }
  else
  {
    gpsimp.seed = Rf_asInteger(seed);
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



  int count_extinct = 0;


  // if the input_file path exists, use this to import variables
  // parsing through the input file and converting/adding to parameter list
  if( !Rf_isNull(input_file) )
  {
    // declare and initialize parameter list for simulation
    SimpleParameterList params;
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

    // convert all parameters imported from file into respective values in gpsimp
    params.convert("tot_life", gpsimp.tot_life);
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
  }
  else
  {
    // convert all parameters imported from file into respective values in gpsimp
    gpsimp.tot_life = tot_life;
    gpsimp.ancestors = ancestors;
    gpsimp.ancestor_clones = ancestor_clones;
    gpsimp.num_sims = num_sims;
    gpsimp.num_samples = num_samples;
    gpsimp.sample_size = sample_size;
    gpsimp.detection_threshold = detection_threshold;
    gpsimp.allow_extinction = allow_extinction;
    gpsimp.birth_rate = birth_rate;
    gpsimp.death_rate = death_rate;
  }

  // declare and open other output streams for time and end of sim clone list
  std::ofstream clonedata;
  sprintf(fn, "%s/clonedata.txt", output_folder);
  clonedata.open(fn);
  clonedata.setf(std::ios::fixed);
  clonedata.precision(12);
  clonedata << "run\tunique_id\tnumcells\tbirthrate\tdeathrate\n";


  // Open output stream for sampling data
  std::ofstream sample_data;
  if( (gpsimp.sample_size > 0) & (gpsimp.num_samples > 0) )
  {
    sprintf(fn, "%s/sampledata.txt", output_folder);
    sample_data.open(fn);
    sample_data.setf(std::ios::fixed);
    sample_data << "run\tsample_number\tunique_id\tnumber_obs\n";
  }

  // set RNG seed
  gsl_rng_set(simple_rng, gpsimp.seed);


  // Beginning of simulation that has "sim" number of runs
  for (int sim = 1; sim <= gpsimp.num_sims; sim++)
  {
    Rcpp::checkUserInterrupt();
    // Define CloneList population and initialize variables;
    SimpleCloneList population;
    population.init();
    population.tot_cell_count = 0;

    if( Rf_isNull(ancestor_file) ) // if no ancestor file exists
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
          int not_zero = gsl_ran_bernoulli(simple_rng, 1 - alpha);
          if(not_zero == 1)
          {
            int new_cells = gsl_ran_geometric(simple_rng, 1 - beta);
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
              int not_zero = gsl_ran_bernoulli(simple_rng, 1 - alpha);
              if(not_zero == 1)
              {
                int new_cells = gsl_ran_geometric(simple_rng, 1 - beta);
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
    if( (gpsimp.sample_size > 0) & (gpsimp.num_samples > 0) )
    {
      population.SampleAndTraverse(sample_data, sim, gpsimp.sample_size, gpsimp.num_samples, simple_rng);
    }
    // Output of end state with clone info
    Rcpp::Rcout << "Traversing and outputting run " << sim << "\n";
    population.Traverse(clonedata, sim);
    Rcpp::Rcout << "Traversal Done\n";

    population.DeleteList();

  }

  //avg_sim_endtime = avg_sim_endtime * (double)gpsimp.num_sims / (double)count_detect;

  gsl_rng_free(simple_rng);

  sim_stats << "count_extinct, " << count_extinct  << "\n";

  clonedata.close();
  sample_data.close();
  sim_stats.close();

  return 0;
}
