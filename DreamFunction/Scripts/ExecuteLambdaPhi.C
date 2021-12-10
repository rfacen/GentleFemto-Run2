#include "GetCorrelationsLambdaPhi.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  GetCorrelationsLambdaPhi(filename, prefix, addon, 0.8, 1.);
  return 1;
}
