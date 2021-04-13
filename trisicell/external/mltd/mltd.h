#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <string>
#include <queue>
#include <map>

typedef struct {
  int distance;
  int similarity;
  double normalized_similarity;
} MLTDResult;
MLTDResult calc_mltd(const char *tree1, const char *tree2);