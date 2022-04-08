#include "GetCorrelationsLKCommonAncestors.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  //GetCorrelationsLambdaKaonRossana(filename, prefix, addon, 0.24, 0.34);
  GetCorrelationsLKCommonAncestors(filename, prefix, addon, 0.24, 0.34);
  
  return 1;
}