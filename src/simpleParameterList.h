/*
 * =====================================================================================
 *
 *       Filename:  parameterlist.h
 *
 *    Description: Header for class that imports the input file, parses
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
#ifndef __SIMPLEPARAMETERLIST_H_INCLUDED__
#define __SIMPLEPARAMETERLIST_H_INCLUDED__

// dependencies
#include <map>
#include <string>
#include <vector>
#include <sstream>


class SimpleParameterList : public std::map<std::string, std::string>
{
public:
    SimpleParameterList() { init(); };
    void init();
    template <class T> bool convert(const std::string s, T& result);
    void ParseVector(const std::string& s, std::vector<double>& v);
    void SplitAndFill(const std::string& s);
};

/*
  SimpleParameterList is a function template to convert the SimpleParameterList map values
  to their appropriate types
*/
template <class T> bool SimpleParameterList::convert(const std::string s, T& result)
{
    std::string val = (*this)[s];
    std::stringstream ss(val);
    ss >> result;
    return val.empty();
}

#endif // __SIMPLEPARAMETERLIST_H_INCLUDED__
