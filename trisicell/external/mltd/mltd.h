#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

typedef struct {
  int distance;
  int similarity;
  double normalized_similarity;
} MLTDResult;
MLTDResult calc_mltd(const char *tree1, const char *tree2);
