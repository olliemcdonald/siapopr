
<!-- README.md is generated from README.Rmd. Please edit that file -->
siapopr
=======

siapopr is an R package that wraps the C++ functions SIApop. These functions simulate birth-death-mutation processes with mutations having random fitnesses to simulate clonal evolution.

Dependencies
------------

-   [GNU Scientific Library](https://www.gnu.org/software/gsl/)
    -   (OSX) `brew install gsl` with Homebrew or from [here](http://ftpmirror.gnu.org/gsl/).
    -   (Windows) download and extract the file [local\#\#\#.zip](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/) and create an environmental variable LIB\_GSL to add the directory (see notes about Windows installation below for more details).
    -   (Linux) install libgsl0-dev and gsl-bin.
-   [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (*Windows only*)
-   [devtools](https://github.com/hadley/devtools)
-   [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
    -   Issues arise with installation of recursive dependencies using install\_github. Installing this R package first solves this issue.

### Important Notes about Windows installation

Rtools contains the necessary resources to compile C++ files when installing packages in R. GSL is also required which can be downloaded from [here](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/local323.zip). After downloading, unzip to your R directory. Depending on if your computer is 32 or 64-bit, move the library files from **local\#\#\#/lib/i386** (32-bit) or **local\#\#\#/lib/x64** (64-bit) to **local\#\#\#/lib**.

To set the environmental variable LIB\_GSL on a Windows 7 computer, go to "Advanced system settings" in *Control Panel &gt; System and Security &gt; System* and click *Environmental Variables*. Create a new system variable with

-   Variable Name: **LIB\_GSL**
-   Variable Value: **"C:/path/to/local323"** (include quotes)

Recommended R packages
----------------------

The following R packages are required for certain functions.

-   [ggmuller](https://github.com/robjohnnoble/ggmuller)
    -   `devtools::install_github("robjohnnoble/ggmuller")`
-   [fishplot](https://github.com/chrisamiller/fishplot)
    -   `devtools::install_github("chrisamiller/fishplot")`
-   ape
-   igraph
-   phangorn

Uses
====

SIApopr (Simulating Infinite-Allele populations in R) is an R package that uses C++ to simulate homogeneous and inhomogeneous birth-death-mutation processes under a very flexible set of assumptions.

The software simulates clonal evolution with the emergence of driver and passenger mutations under the infinite-allele assumption. The software is a application of the direct Gillespie Stochastic Simulation Algorithm expanded to a large number of cell types and scenarios, with the intention of allowing users to easily modify existing models or create their own. Visualization functions in R are included to show results of individual simulations.

A branching process is a stochastic process used to model the growth and composition of reproducing populations. Assumptions made in branching processes are individuals live for a random amount of time before splitting into a random number of individuals (both dictated by distribution functions). Individuals of the same type are independent and identically distributed. These processes are useful for modeling cell growth and evolution, as in a tumor. Mutations may occur that lead to different types of individuals with different probability laws dictating those individuals' birth and death rates.

Installation
============

To install in R, type:

``` r
devtools::install_git("https://github.com/olliemcdonald/siapopr")
install.packages("dplyr")
```

Installing the library should compile all necessary functions so that SIApopr can be run as an R function.

Using SIApop in R
=================

SIApopr contains three functions to simulate birth-death processes based on the desired model. To run a birth-death process without mutation where birth and death rates are constant, `siapopSimple` simulates the number of individuals with respect to each ancestor clone at the exact time given by the parameter. `siapopSimple` generates the exact number of descendants for each clone by generating the number of ancestors that have nonextinct lines as a binomial random variable, then generating the number of descendants of those lines as a negative binomial random variable. Since mutation is not allowed, `siapopConstant` does not simulate evolution, and only allows growth of preexisting clones.

If a model with mutation is desired and rates are constant, `siapopConstant` simulates a birth-death-mutation process under the infinite-allele assumption, that any mutation leads to a new allele that has not yet been observed in the population. This is distinct from a traditional multitype branching process where types are accessible from each other more than once. Within this general model contains the ability to include fitness effects, where a new mutation has a birth rate equal to the sum of its mother's and a random variable from a double exponential distribution. Details about the distribution are given below. A neutral model can also be created, where mutations are passengers only and have no effect on fitness of new clones.

The final function allows for time-dependent birth and death rates, named `siapopTimeDep`. For now, we included preselected time-dependent kernels for the rate functions and the user can input the selected parameters. Parameters are input as a numeric vector. More details are given below.

To run, type in the R console:

``` r
siapopConstant(...)
```

and the simulation proceeds. Results from the simulation are output as a file instead of being loaded into the environment, and the user imports the data after the simulation using the import functions provided. Default values output data to the current working directory and import from that directory as well.

Model Parameters
================

<table style="width:64%;">
<colgroup>
<col width="27%" />
<col width="19%" />
<col width="16%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable Name</th>
<th>Variable Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>tot_life</td>
<td>numeric</td>
<td>total lifetime of branching process</td>
</tr>
<tr class="even">
<td>max_pop</td>
<td>numeric</td>
<td>maximum population to end process at</td>
</tr>
<tr class="odd">
<td>start_time</td>
<td>numeric</td>
<td>starting time (typically 0)</td>
</tr>
<tr class="even">
<td>ancestors</td>
<td>int</td>
<td>number of individuals per ancestor clone (if ancestor file not provided)</td>
</tr>
<tr class="odd">
<td>ancestor_clones</td>
<td>int</td>
<td>number of initial types</td>
</tr>
<tr class="even">
<td>num_sims</td>
<td>int</td>
<td>total number of simulations of the same process</td>
</tr>
<tr class="odd">
<td>allow_extinction</td>
<td>TRUE/FALSE</td>
<td>1 if allow a simulation to go extinct</td>
</tr>
<tr class="even">
<td>detection_threshold</td>
<td>numeric</td>
<td>minimum proportion of the total population such that a clone is output</td>
</tr>
<tr class="odd">
<td>num_samples</td>
<td>int</td>
<td>number of samples to take</td>
</tr>
<tr class="even">
<td>sample_size</td>
<td>int</td>
<td>size of each sample</td>
</tr>
<tr class="odd">
<td>birth_rate</td>
<td>numeric &gt; 0</td>
<td>starting birth rate</td>
</tr>
<tr class="even">
<td>death_rate</td>
<td>numeric &gt; 0</td>
<td>starting death rate</td>
</tr>
<tr class="odd">
<td>mutation_prob</td>
<td>numeric [0, 1]</td>
<td>default mutation probability for new clone</td>
</tr>
<tr class="even">
<td>trace_ancestry</td>
<td>TRUE/FALSE</td>
<td>Track info on parent of each clone</td>
</tr>
<tr class="odd">
<td>count_alleles</td>
<td>TRUE/FALSE</td>
<td>adds/subtracts and individual to allele_count of individual and all ancestors</td>
</tr>
</tbody>
</table>

### FITNESS DISTRIBUTION PARAMETERS

The fitness distribution provided is a double exponential with an atom at 0 indicated by the probability of a passenger mutation, `pass_prob`. The rate for the right side of the distribution is `alpha` and the left side is `beta`. Bounds are provided as `upper_fitness` and `lower_fitness` to bound the probability distribution to this range. These parameters are not required. The default model assumes passenger mutations only.

<table style="width:56%;">
<colgroup>
<col width="19%" />
<col width="19%" />
<col width="16%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable Name</th>
<th>Variable Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>alpha_fitness</td>
<td>numeric &gt; 0</td>
<td>exponential distribution parameter for positive side of fitness distribution</td>
</tr>
<tr class="even">
<td>beta_fitness</td>
<td>numeric &gt; 0</td>
<td>exponential distribution parameter for negative side of fitness distribution</td>
</tr>
<tr class="odd">
<td>pass_prob</td>
<td>numeric [0,1]</td>
<td>probability that additional fitness of new mutant is 0</td>
</tr>
<tr class="even">
<td>upper_fitness</td>
<td>numeric</td>
<td>upper bound to fitness distribution</td>
</tr>
<tr class="odd">
<td>lower_fitness</td>
<td>numeric &lt;= upper_fitness</td>
<td>lower bound to fitness distribution</td>
</tr>
</tbody>
</table>

### MUTATION DISTRIBUTION PARAMETERS

If these parameters are defined, any new clone has a mutation probability coming from a beta distribution with parameters `alpha_mutation` and `beta_mutation`.

<table style="width:57%;">
<colgroup>
<col width="20%" />
<col width="19%" />
<col width="16%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable Name</th>
<th>Variable Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>alpha_mutation</td>
<td>numeric &gt; 0</td>
<td>alpha parameter for Beta distribution for additional mutation probability in new mutant clone</td>
</tr>
<tr class="even">
<td>beta_mutation</td>
<td>numeric &gt; 0</td>
<td>beta parameter for Beta distribution for additional mutation probability in new mutant clone</td>
</tr>
</tbody>
</table>

### PUNCTUATED EVOLUTION PARAMETERS

If the parameters are defined, a punctuated evolution model is used, which assumes with probability `punctuated_prob` that a new mutation event has a burst of mutations, or multiple alleles arise at once. When this occurs, the number of new alleles is Poisson with parameter `poisson_param`. Punctuated models should provide a disadvantage to the new clone, so these clones have a very high death rate unless an advantageous clone arises with probability defined by the parameter `punctuated_advantageous_prob`. Finally, given an advantageous clone arises, its new birth rate is multiplied by a factor of `punctuated_fitness_multiplier`.

<table style="width:78%;">
<colgroup>
<col width="41%" />
<col width="19%" />
<col width="16%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable Name</th>
<th>Variable Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>punctuated_prob</td>
<td>numeric [0,1]</td>
<td>probability of mutation burst</td>
</tr>
<tr class="even">
<td>poisson_param</td>
<td>numeric &gt; 0</td>
<td>rate parameter for zero-truncated Poisson distribution number of mutations in burst</td>
</tr>
<tr class="odd">
<td>punctuated_fitness_multiplier</td>
<td>numeric</td>
<td>amount to multiply additional fitness by</td>
</tr>
<tr class="even">
<td>punctuated_advantageous_prob</td>
<td>numeric [0,1]</td>
<td>probability that burst affects birth rate instead of death rate</td>
</tr>
</tbody>
</table>

### EPISTATIC PARAMETERS

The epistatic model is a basic one that assumes the birth rate is multiplied by a factor of `epistatic_multiplier` after the *k*<sup>*t**h*</sup> mutation arises. The variable `epistatic_mutation_threshold` is *k*.

<table style="width:75%;">
<colgroup>
<col width="38%" />
<col width="19%" />
<col width="16%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable Name</th>
<th>Variable Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>epistatic_mutation_threshold</td>
<td>int &gt; 0</td>
<td>number of mutation required before burst in fitness due to epistasis</td>
</tr>
<tr class="even">
<td>epistatic_multiplier</td>
<td>numeric</td>
<td>amount to multiply fitness contribution in new clone by due to epistasis occurring</td>
</tr>
</tbody>
</table>

### TIME-DEPENDENT PARAMETERS (siapopTimeDep only)

The following parameters define the birth and death rate functions and associated parameters. More information is given in the next section.

| Variable Name     | Variable Type  | Description |
|-------------------|----------------|-------------|
| birth\_function   | 0, 1, 2, 3, 4  | see below   |
| death\_function   | 0, 1, 2, 3, 4  | see below   |
| td\_birth\_params | numeric vector | see below   |
| td\_death\_params | numeric vector | see below   |

Ancestor File
=============

The ancestor file is a tab-delimited file with the same structure as the output. The first line contains variable names and each line contains information for a single clone to serve as an ancestor population. The only requirement for this file is a column containing the number of cells. If not provided, the program will look at the arguments ancestors and ancestor\_clones described above to run with nonunique clones. If those are not provided a default of a single ancestor individual is used. The following table describes the possible variables for the ancestor file. Some parameters below are included since they are present in the output and store information when continuing a previous simulation.

### PARAMETERS

<table style="width:54%;">
<colgroup>
<col width="18%" />
<col width="19%" />
<col width="16%" />
</colgroup>
<thead>
<tr class="header">
<th>Variable Name</th>
<th>Variable Type</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>unique_id</td>
<td>string</td>
<td>id for each ancestor</td>
</tr>
<tr class="even">
<td>numcells</td>
<td>int</td>
<td>the number of cells for the ancestor</td>
</tr>
<tr class="odd">
<td>mutprob</td>
<td>numeric [0,1]</td>
<td>the probability of initiating a new clone given a birth occurs</td>
</tr>
</tbody>
</table>

### TIME-DEPENDENT PARAMETERS

| Variable Name   | Variable Type                       | Description |
|-----------------|-------------------------------------|-------------|
| birth\_function | 0, 1, 2, 3, 4                       | see below   |
| death\_function | 0, 1, 2, 3, 4                       | see below   |
| bf\_params      | vector of doubles (space-delimited) | see below   |
| df\_params      | vector of doubles (space-delimited) | see below   |

Time-Dependent Rate Functions
-----------------------------

The functions all are predefined and allow the user to provide a list of parameters. The parameters should be listed the same way regardless of which function is used in the form of “x1,x2,x3,x4,x5…” in the tab-delimeted file (see example). An extra function is included called "custom" to allow the user to define a unique function, but recompiling the program is required after. The curves are parameterized as follows:

| Function Parameter | Function Name   | Mathematical Description               |
|--------------------|-----------------|----------------------------------------|
| 0                  | constant        | ![constant](README/README-image-1.png) |
| 1                  | linear          | ![linear](README/README-image-2.png)   |
| 2                  | logistic        | ![logistic](README/README-image-3.png) |
| 3                  | Gompertz growth | ![Gompertz](README/README-image-4.png) |
| 4                  | Custom          | Include own parameters                 |
