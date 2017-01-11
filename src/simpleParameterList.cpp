/*
 * =====================================================================================
 *
 *       Filename:  parameterlist.cpp
 *
 *    Description: Contains the class that imports the input file, parses
 *                 and converts to parameters for the model.
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

#include "simpleParameterList.h"

/*
  ParseVector turns space-delimited text into a double vector
*/
void SimpleParameterList::ParseVector(const std::string& s, std::vector<double>& v)
{
    std::istringstream s2((*this)[s]);
    double tmp;

    while(s2 >> tmp) v.push_back(tmp);
}

/*
  SplitAndFill splits each line into the parameter name and its value and stores
  it into a map.
*/
void SimpleParameterList::SplitAndFill(const std::string& s)
{
    std::string::size_type pos1 = s.find_first_of(" ,=\t");
    std::string::size_type pos2 = s.find_first_not_of(" ,=\t", pos1+1);
    if (pos2 == std::string::npos) return;

    std::string k = s.substr(0, pos1);
    std::string val =  s.substr(pos2);

    if (!val.compare("NA")) return;

    (*this)[k] = val;
}



void SimpleParameterList::init()
{
    //default values
    insert(std::make_pair("tot_life", "40000"));
    insert(std::make_pair("ancestors", "1"));
    insert(std::make_pair("ancestor_clones", "1"));
    insert(std::make_pair("num_sims", "1"));
    insert(std::make_pair("allow_extinction", "1"));
    insert(std::make_pair("num_samples", "0"));
    insert(std::make_pair("sample_size", "0"));
    insert(std::make_pair("birth_rate", "1.5"));
    insert(std::make_pair("death_rate", "1"));
}
