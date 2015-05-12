#include <Rcpp.h>
using namespace Rcpp;

void nested_loop_operation_haplotype_probabilities_sum(double* res[], int counters[], int start_alleles[], int stop_alleles[], int level, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int loci = y.ncol();
  
  if (level == loci) {
    int clusters = y.nrow();

    /*
    Rcout << "(";
    for (int l = 0; l < loci; l++) {
      Rcout << counters[l] << ", ";
    }
    Rcout << ")" << std::endl;
    */
    
    double hprob = 0.0;
        
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      double component_prob = tau(c);
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        component_prob *= pow(pcl, abs(counters[l] - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      hprob += component_prob;
    }
    
    //Rcout << (*res)[0] << std::endl;
    
    (*res)[0] += hprob;
    (*res)[1] += hprob * hprob;
  } else {
    for (int allele = start_alleles[level]; allele <= stop_alleles[level]; allele++) {
      counters[level] = allele;
      nested_loop_operation_haplotype_probabilities_sum(res, counters, start_alleles, stop_alleles, level + 1, y, p, tau);
    }
  }
}

// [[Rcpp::export]]
NumericVector rcpp_calculate_haplotype_probabilities_sum(IntegerMatrix allele_range, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int loci = allele_range.ncol();

  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }

  if (allele_range.nrow() != 2) {
    throw std::range_error("Exactly two rows in allele_range required");
  }
  
  double* res = new double[2];
  int* counters = new int[loci];
  int* start_alleles = new int[loci];
  int* stop_alleles = new int[loci];
  
  res[0] = 0.0;
  res[1] = 0.0;
  
  for (int l = 0; l < loci; l++) {
    start_alleles[l] = allele_range(0, l);
    stop_alleles[l] = allele_range(1, l);
  }
  
  nested_loop_operation_haplotype_probabilities_sum(&res, counters, start_alleles, stop_alleles, 0, y, p, tau);
  
  // Need to convert to vector, or else the return is not OK. Weird.
  std::vector<double> ans(2);
  ans[0] = res[0];
  ans[1] = res[1];
  
  delete[] res;
  delete[] counters;
  delete[] start_alleles;
  delete[] stop_alleles;
  
  NumericVector ret = NumericVector::create(
    _["sum_prob"] = ans[0],
    _["sum_sq_prob"] = ans[1]);
  
  return ret;
}

/****************************/

void nested_loop_operation_match_quantities(NumericMatrix res, int counters[], int start_alleles[], int stop_alleles[], int level, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int loci = y.ncol();
  
  if (level == loci) {
    int clusters = y.nrow();
    
    NumericVector happrobs(clusters);
    
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      double component_prob = 1.0;
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        component_prob *= pow(pcl, abs(counters[l] - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      happrobs(c) = component_prob;
    }
    
    for (int c1 = 0; c1 < (clusters - 1); c1++) {
      for (int c2 = c1+1; c2 < clusters; c2++) {
        res(c1, c2) += happrobs(c1) * happrobs(c2);
      }
    }
  } else {
    for (int allele = start_alleles[level]; allele <= stop_alleles[level]; allele++) {
      counters[level] = allele;
      nested_loop_operation_match_quantities(res, counters, start_alleles, stop_alleles, level + 1, y, p, tau);
    }
  }
}

// [[Rcpp::export]]
NumericVector rcpp_match_quantities(IntegerMatrix allele_range, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int clusters = y.nrow();
  int loci = allele_range.ncol();

  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }

  if (allele_range.nrow() != 2) {
    throw std::range_error("Exactly two rows in allele_range required");
  }
  
  NumericMatrix res(clusters, clusters);
  int* counters = new int[loci];
  int* start_alleles = new int[loci];
  int* stop_alleles = new int[loci];
  
  for (int l = 0; l < loci; l++) {
    start_alleles[l] = allele_range(0, l);
    stop_alleles[l] = allele_range(1, l);
  }
  
  nested_loop_operation_match_quantities(res, counters, start_alleles, stop_alleles, 0, y, p, tau);
  
  delete[] counters;
  delete[] start_alleles;
  delete[] stop_alleles;
  
  return res;
}

