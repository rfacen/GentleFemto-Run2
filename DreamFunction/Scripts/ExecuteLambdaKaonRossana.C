#include "GetCorrelationsLambdaKaonRossana.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  //GetCorrelationsLambdaKaonRossana(filename, prefix, addon, 0.24, 0.34);
  GetCorrelationsLambdaKaonRossana(filename, prefix, addon, 0.50, 0.80);
  
  return 1;
}