
#include <fstream>
#include <iostream>
#include <string>

void readFASTA(std::string FileName, std::vector<char> &Seq) {
  std::ifstream inf(FileName);
  char c = inf.get();
  while (inf.good()) {
    if (c=='>') do { c = inf.get(); } while(c!='\r' && c!='\n' && inf.good());
    if (c==';') do { c = inf.get(); } while(c!='\r' && c!='\n' && inf.good());
    if (std::isalpha(c)) Seq.push_back(std::toupper(c));
    c = inf.get();
  }
  inf.close();
}

