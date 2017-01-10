#include <Rcpp.h>
#include "siapopConstant.h"
#include "siapopSimple.h"

// [[Rcpp::export]]
int siapop_constant(Rcpp::CharacterVector input)
{
  siapopConstant(input);
  return 0;
}

// [[Rcpp::export]]
int siapop_simple(Rcpp::CharacterVector input)
{
  siapopSimple(input);
  return 0;
}
