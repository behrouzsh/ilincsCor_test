#include "ilincs_corx.h"

//===================================================================
std::vector<std::vector<std::string> > read_tsv(std::string fname) {
    ifstream ifs(fname);
    std::vector<std::vector<std::string> > output;
    // std::cerr << "Inside read_tsv(" << fname << ")\n";
    if (ifs.fail()) {
        cerr << "error opening tsv file " << fname << std::endl;
        return output;
    }
    string line;
    while (getline(ifs, line)) {
       if ( line.find ("\r\n") != string::npos )
       {
          line.erase ( line.find ("\r\n"), 2 );
       }
       if ( line.find ("\r") != string::npos )
       {
          line.erase ( line.find ("\r"), 1 );
       }
       if (line[0]!=0) {
          stringstream ss(line);
          vector<string> item;
          string tmp;
          while(getline(ss, tmp, '\t')) {
             item.push_back(tmp);
          }
          output.push_back(item);
       }
    }
    return output;
}
