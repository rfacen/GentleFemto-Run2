#include "GetCorrelationsNano.C"
#include "GetCorrelationsNanoBBar.C"
#include "GetCorrelationsNanoBenchmark.C"
#include "METoSEReweightingBBar.C"
#include "TSystem.h"
int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  GetCorrelationsNanoBenchmark(filename, prefix, addon, 0.24, 0.34);
  return 1;
}
