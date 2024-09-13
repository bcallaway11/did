#include <Rcpp.h>
using namespace Rcpp;

// This code contains a series of sanity checks for preprocessing.
// These checks are used to ensure that the input data is consistent
// Date: 2024-09-12

// [[Rcpp::export]]
void checkTreatmentUniqueness(DataFrame dta, String idname, String treatName) {
  CharacterVector id = dta[idname];
  NumericVector cohort = dta[treatName];

  std::map<String, std::set<double>> cohortMap;

  // Building a map where each key is an idname and the value is a set of cohort.names
  for (int i = 0; i < id.size(); ++i) {
    cohortMap[id[i]].insert(cohort[i]);
  }

  // Check if any id has more than one unique cohort.name
  for (const auto &pair : cohortMap) {
    if (pair.second.size() > 1) {
      stop("The value of dname (treatment variable) must be the same across all periods for each particular unit");
    }
  }
  // If we reach here, all dname are unique by idname
}

// [[Rcpp::export]]
void checkWeightsUniqueness(DataFrame dta, String idname) {
  CharacterVector id = dta[idname];
  NumericVector weights = dta["weights"];

  std::map<String, std::set<double>> weightsMap;

  // Building a map where each key is an idname and the value is a set of weights
  for (int i = 0; i < id.size(); ++i) {
    weightsMap[id[i]].insert(weights[i]);
  }

  // Check if any id has more than one unique weight
  for (const auto &pair : weightsMap) {
    if (pair.second.size() > 1) {
      stop("The value of weights must be the same across all periods for each particular unit");
    }
  }
  // If we reach here, all weights are unique by idname
}
