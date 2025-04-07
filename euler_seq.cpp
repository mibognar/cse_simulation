#include <Rcpp.h>
#include <algorithm>
#include <cstdlib>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector generate_sequence(int n_trials) {
  int n_each = n_trials / 4;
  int a = n_each, b = n_each, c = n_each, d = n_each;

  std::vector<int> circuit;
  std::vector<int> stack;

  int current = 1;
  stack.push_back(current);

  while (!stack.empty()) {
    if ((current == 1 && (a > 0 || b > 0)) || (current == -1 && (c > 0 || d > 0))) {
      stack.push_back(current);
      if (current == 1) {
        int total = a + b;
        int r = rand() % total;
        current = (r < a) ? (--a, 1) : (--b, -1);
      } else {
        int total = c + d;
        int r = rand() % total;
        current = (r < c) ? (--c, -1) : (--d, 1);
      }
    } else {
      circuit.push_back(current);
      stack.pop_back();
      if (!stack.empty()) current = stack.back();
    }
  }

  std::reverse(circuit.begin(), circuit.end());

  IntegerVector is_congruent(n_trials);
  for (int i = 0; i < n_trials; i++) {
    is_congruent[i] = circuit[i + 1];
  }
  return is_congruent;
}
