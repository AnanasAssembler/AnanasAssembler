#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "math/Spearman.h"


class FPKM
{
public:
  FPKM() {}
  
  void AddTissue(const string & s) {
    svec<double> tmp;
    //char n[1024];
    //strcpy(n, s.c_str());
    //n[strlen(n)-5] = 0;
    m_data.push_back(tmp);
    m_name.push_back(s);
  }

  void Add(int i, double val) {
    m_data[i].push_back(val);
  }

  int Num() const {return m_data.isize();}
  const svec<double> & Get(int i) {
    return m_data[i];
  }
  const string & Name(int i) const {
    return m_name[i];
  }

private:
  svec < svec < double > > m_data;
  svec<string> m_name;
};

void Read(FPKM & out, const string & fileName) 
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  parser.ParseLine();
  int i;
  svec<int> index;
  for (i=0; i<parser.GetItemCount(); i++) {
    if (strstr(parser.AsString(i).c_str(), "FPKM") != NULL) {
      out.AddTissue(parser.AsString(i));
      index.push_back(i);
    }
  }

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    for (i=0; i<index.isize(); i++) {
      out.Add(i, parser.AsFloat(index[i]));
    }
  }
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  //commandArg<string> file2Cmmd("-i2","input file 2");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing Spearman's rho.");
  P.registerArg(fileCmmd);
  //P.registerArg(file2Cmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  //string fileName2 = P.GetStringValueFor(file2Cmmd);
  
  FPKM a;
  Read(a, fileName);
  // Read(b, fileName2);
 
  printf("%d\n", a.Num());

  int i, j, k;
  for (i=0; i<a.Num(); i++) {
    const svec<double> & one = a.Get(i);
    char tmp[256];
    strcpy(tmp, a.Name(i).c_str());
    tmp[strlen(tmp)-5] = 0;
    printf("%s", tmp);
    for (int y=strlen(a.Name(i).c_str()); y<14; y++)
      printf(" ");
    for (j=0; j<a.Num(); j++) {
      const svec<double> & two = a.Get(j);
      svec<double> x, y;
      for (k=0; k<one.isize(); k++) {
	if (one[k] > 0.1 && two[k] > 0.1) {
	  x.push_back(one[k]);
	  y.push_back(two[k]);
	}
      }
      SpearmansRho r;
      double rho = r.Compute(x, y);
      if (j > 0)
	printf(" ");
      double dd = 1. - rho;
      if (dd > 0.99)
	dd = 0.99;

      printf("%1.3f", dd);
      //cout << a.Name(i) << " " << b.Name(j) << " ";
      //cout << "Rho " << rho << " Significance=" << r.Significance() << " data points: " << x.isize() <<  endl;
    
    }
    printf("\n");
  }

  
 
  return 0;
}
