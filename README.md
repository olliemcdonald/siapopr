---
title: "SIApopr: A computational method to simulate evolutionary branching trees for analysis of tumor clonal evolution"
author: "Thomas McDonald and Franziska Michor"
output: pdf_document
# output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# Siapopr

Siapopr is an R package that wraps the C++ functions SIApop. These functions
simulate birth-death-mutation processes with mutations having random fitnesses
to simulate clonal evolution.

## Dependencies

* [GNU Scientific Library](https://www.gnu.org/software/gsl/)
    + (OSX) `brew install gsl` with Homebrew or from
  [here](http://ftpmirror.gnu.org/gsl/).
    + (Windows) download and extract
  the file [local###.zip](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/)
  and create an environmental variable LIB_GSL to add the directory (see notes
  about Windows installation below for more details).
    + (Linux) install libgsl0-dev and gsl-bin.
* [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (*Windows only*)
* [devtools](https://github.com/hadley/devtools)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
    + Issues arise with installation of recursive dependencies using
    install_github. Installing this R package first solves this issue.

### Important Notes about Windows installation
Rtools contains the necessary resources to compile C++ files when installing
packages in R. GSL is also required which can be downloaded using msys2 which is in the Rtools directory. To install GSL, open the msys2 application and type/run the following commands:

`pacman -Sy mingw-w64-x86_64-gsl`

`pacman -Sy mingw-w64-i686-gsl`

Verify that the GSL libraries are installed by navigating to __[Your Rtools Directory]>mingw64>include>gsl__ and verifying this folder exists.

Add these locations to your PATH and LIB_GSL environmental variables in R by adding the following lines to your __.Renviron__ file. If you do not have a .Renviron file, create one in your home directory and add the following lines:
`PATH=“${RTOOLS40_HOME}\usr\bin;${PATH}”`

`LIB_GSL=“${RTOOLS40_HOME}\mingw64”`

## Recommended R packages
The following R packages are required for certain functions.

* [ggmuller](https://github.com/robjohnnoble/ggmuller)
    + `devtools::install_github("robjohnnoble/ggmuller")`
* [fishplot](https://github.com/chrisamiller/fishplot)
    + `devtools::install_github("chrisamiller/fishplot")`
* ape
* igraph
* phangorn

# Installation

To install in R, type:
```{r, eval = F}
install.packages("devtools")
devtools::install_git("git://github.com/olliemcdonald/siapopr.git")
install.packages("dplyr")
```

Installing the library should compile all necessary functions so that SIApopr
can be run as an R function.

# Uses
SIApopr (Simulating Infinite-Allele populations in R) is an R package that
uses C++ to simulate time homogeneous and inhomogeneous birth-death-mutation
processes under a very flexible set of assumptions.

The software simulates clonal evolution with the emergence of driver and
passenger mutations under the infinite-allele assumption. The software uses an
application of the direct Gillespie Stochastic Simulation Algorithm expanded to
a large number of cell types and scenarios with the intention of allowing
users to easily modify existing models or create their own. Conversion functions
to other visualization packages in R are included to show results of
individual simulations.

A branching process is a stochastic process used to model the growth and
composition of reproducing populations. Assumptions made in branching processes are
individuals live for a random amount of time before splitting into a random number
of individuals (both random variables are dictated by their given distribution
functions). Individuals of the same type are independent and identically
distributed. These processes are useful for modeling cell growth and evolution,
as in a tumor. Mutations may occur that lead to different types of individuals
with different probability laws dictating those individuals' birth and death
rates.

# Using SIApopr in R

SIApopr contains three functions to simulate birth-death processes based on the
desired model. To run a birth-death process without mutation where birth and
death rates are constant, `siapopNoMut` simulates the number of individuals
with respect to each ancestor clone at the exact time given by the parameter.
`siapopNoMut` generates the exact number of descendants for each clone by
generating the number of ancestors that have nonextinct lines as a binomial
random variable, then generating the number of descendants of those lines
as a negative binomial random variable. Since mutation is not allowed,
`siapopNoMut` does not simulate evolution, and only allows growth of
preexisting clones.

If a model with mutation is desired and rates are constant, `siapop`
simulates a birth-death-mutation process under the infinite-allele assumption,
that any mutation leads to a new allele that has not yet been observed in the
population. This is distinct from a traditional multitype branching process
where types are accessible from each other more than once. Within this general
model contains the ability to include fitness effects, where a new mutation
has a birth rate equal to the sum of its mother's and a random variable from
a double exponential distribution. Details about the distribution are given
below. A neutral model can also be created, where mutations are passengers only
and have no effect on fitness of new clones.

The final function allows for time-dependent birth and death rates, named
`siapopTD`. For now, we included preselected time-dependent kernels
for the rate functions and the user can input the selected parameters.
Parameters are input as a numeric vector. More details are given below.

To run, type in the R console:
```{r, eval = F}
siapop(...)
```
and the simulation proceeds. Results from the simulation are output as a file
instead of being loaded into the environment, and the user imports the data
after the simulation using the import functions provided. Default values
output data to the current working directory and import from that directory as
well.

# Model Parameters

| Variable Name       | Variable Type | Description |
| ------------------- | ------------- | ------------------------------------- |
| tot_life            | numeric        | total lifetime of branching process |
| max_pop             | numeric        | maximum population to end process at|
| start_time          | numeric        | starting time (typically 0)|
| ancestors           | int           | number of individuals per ancestor clone (if ancestor file not provided)|
| ancestor_clones     | int           | number of initial types |
| num_sims            | int           | total number of simulations of the same process|
| allow_extinction    | TRUE/FALSE    | 1 if allow a simulation to go extinct |
| detection_threshold | numeric        | minimum proportion of the total population such that a clone is output|
| num_samples         | int           | number of samples to take|
| sample_size         | int           | size of each sample|
| birth_rate          | numeric > 0    |  starting birth rate |
| death_rate          | numeric > 0    | starting death rate |
| mutation_prob       | numeric [0, 1] | default mutation probability for new clone|
| trace_ancestry      | TRUE/FALSE    | Track info on parent of each clone |
| count_alleles       | TRUE/FALSE    | adds/subtracts and individual to allele_count of individual and all ancestors |
| custom_model_file    | string    | A custom model shared object file |

### FITNESS DISTRIBUTION PARAMETERS

The fitness distribution provided is a double exponential with an atom
at 0 indicated by the probability of a passenger mutation, `pass_prob`. The
rate for the right side of the distribution is `alpha` and the left side is
`beta`. Bounds are provided as `upper_fitness` and `lower_fitness` to bound
the probability distribution to this range. These parameters are not required.
The default model assumes passenger mutations only.

| Variable Name | Variable Type | Description |
| ------------- | ------------- | ----------------------------------------------- |
| distribution_function | doubleexp, normal, uniform, custom | the type of distribution to use (see below)
| custom_distribution | string | if using a custom distribution, the filepath for the .so file used (see below)
| alpha_fitness | numeric > 0             | parameter 1 for the distribution |
| beta_fitness  | numeric > 0             | parameter 2 for the distribution |
| pass_prob     | numeric [0,1]           | probability that additional fitness of new mutant is 0|
| upper_fitness | numeric                 | upper bound to fitness distribution|
| lower_fitness | numeric <= upper_fitness | lower bound to fitness distribution|

#### Fitness Distribution Functions
There is an option to specify different distributions using the parameter `distribution_function`. For
a double exponential, `alpha_fitness` refers to the rate parameter with respect to an exponential distribution to the right of 0. `beta_fitness` is the rate parameter for the left of 0.

If `distribution_function = "normal"`, fitnesses are generated from a N(`alpha_fitness`, `beta_fitness`)
distribution. If `distribution_function = "normal"`, fitnesses are generated from a U(`alpha_fitness`, `beta_fitness`) distribution.

The user has the option to specify a custom distribution. A distribution should
be written as a '.cpp' file and include the following:
```{r engine='Rcpp', eval = FALSE}
void customdist(double* fitness, struct FitnessParameters *fit_params, gsl_rng* rng)
{
  (*fitness) = //insert fitness function;
}
```
An example for the uniform distribution with a passenger probability is:
```{r engine='Rcpp', eval = FALSE}
void customdist(double* fitness, struct FitnessParameters *fit_params, gsl_rng* rng)
{
  double z = gsl_ran_flat(rng, 0, 1);
  if( (z > fit_params->pass_prob) )
  {
    (*fitness) = gsl_ran_flat(rng, fit_params->alpha_fitness, fit_params->beta_fitness);
  }
}
```
The file should be saved as a ".cpp" file and then use the R function
`compile_custom_fitness(cppfile)` which creates the necessary headers, compiles, and
builds a shared library out of the distribution function file so that SIApopr can
load the shared object.

### MUTATION DISTRIBUTION PARAMETERS

If these parameters are defined, any new clone has a mutation probability coming
from a beta distribution with parameters `alpha_mutation` and `beta_mutation`.

| Variable Name  | Variable Type | Description |
| -------------- | ------------- | ----------------------------------------------- |
| alpha_mutation | numeric > 0    | alpha parameter for Beta distribution for additional mutation probability in new mutant clone|
| beta_mutation  | numeric > 0    | beta parameter for Beta distribution for additional mutation probability in new mutant clone|

### PUNCTUATED EVOLUTION PARAMETERS

If the parameters are defined, a punctuated evolution model is used, which assumes
with probability `punctuated_prob` that a new mutation event has a burst of mutations,
or multiple alleles arise at once. When this occurs, the number of new
alleles is Poisson with parameter `poisson_param`. Punctuated models should
provide a disadvantage to the new clone, so these clones have a very high death
rate unless an advantageous clone arises with probability defined by the parameter
`punctuated_advantageous_prob`. Finally, given an advantageous clone arises, its
new birth rate is multiplied by a factor of `punctuated_fitness_multiplier`.

| Variable Name                 | Variable Type | Description |
| ----------------------------- | ------------- | ---------------------------- |
| punctuated_prob               | numeric [0,1]  | probability of mutation burst |
| poisson_param                 | numeric > 0    | rate parameter for zero-truncated Poisson distribution number of mutations in burst |
| punctuated_fitness_multiplier | numeric        | amount to multiply additional fitness by |
| punctuated_advantageous_prob  | numeric [0,1]  | probability that burst affects birth rate instead of death rate |


### EPISTATIC PARAMETERS

The epistatic model only assumes the birth rate is multiplied
by a factor of `epistatic_multiplier` after the $k^{th}$ mutation arises. The
variable `epistatic_mutation_threshold` is $k$.

| Variable Name               | Variable Type | Description |
| --------------------------- | ------------- | ------------------------------ |
|epistatic_mutation_threshold | int > 0       | number of mutation required before burst in fitness due to epistasis|
|epistatic_multiplier         | numeric        | amount to multiply fitness contribution in new clone by due to epistasis occurring|


### TIME-DEPENDENT PARAMETERS (siapopTD only)

The following parameters define the birth and death rate functions and
associated parameters. More information is given in the next section.

| Variable Name  | Variable Type                       | Description |
| -------------- | ----------------------------- | --------------------------- |
|birth_function  | 0, 1, 2, 3, 4                       | see below |
|death_function  | 0, 1, 2, 3, 4                       | see below |
|td_birth_params | numeric vector | see below |
|td_death_params | numeric vector | see below |

# Ancestor File

The ancestor file is a tab-delimited file with the same structure as the output.
The first line contains variable names and each line contains information for a
single clone to serve as an ancestor population. The only requirement for this
file is a column containing the number of cells. If not provided, the program
will look at the arguments ancestors and ancestor_clones described above to run
with nonunique clones. If those are not provided a default of a single ancestor
individual is used. The following table describes the possible variables for the
ancestor file. Some parameters below are included since they are present in the
output and store information when continuing a previous simulation.

### PARAMETERS

| Variable Name | Variable Type | Description |
| ------------  | ------------- | -------------------------------------- |
| unique_id     | string        | id for each ancestor |
| numcells      | int           | the number of cells for the ancestor |
| mutprob       | numeric [0,1]  | the probability of initiating a new clone given a birth occurs |

### TIME-DEPENDENT PARAMETERS

| Variable Name  | Variable Type                       | Description |
| -------------- | ----------------------------------- | ------------------------------- |
|birth_function  | 0, 1, 2, 3, 4                       | see below |
|death_function  | 0, 1, 2, 3, 4                       | see below |
|bf_params       | vector of doubles (space-delimited) | see below |
|df_params       | vector of doubles (space-delimited) | see below |


## Time-Dependent Rate Functions

The functions all are predefined and allow the user to provide a list of
parameters. The parameters should be listed the same way regardless of which
function is used in the form of “x1,x2,x3,x4,x5…” in the tab-delimeted file
(see example). An extra function is included called "custom" to allow the user
to define a unique function, but recompiling the program is required after.
The curves are parameterized as follows:

| Function Parameter  | Function Name     | Mathematical Description |
| ------------------- | ----------------- | ------------------------ |
| 0                   | constant          | ![constant](README/README-image-1.png) |
| 1                   | linear            | ![linear](README/README-image-2.png)  |
| 2                   | logistic          | ![logistic](README/README-image-3.png)  |
| 3                   | Gompertz growth   | ![Gompertz](README/README-image-4.png)  |
| 4                   | Custom            | Include own parameters  |
