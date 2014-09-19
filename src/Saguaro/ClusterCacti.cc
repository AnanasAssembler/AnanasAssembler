#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Saguaro/Cactus.h"

#define FORCE_DEBUG

class DistMatrix
{
public:
  DistMatrix(int n) {
    m_dist.resize(n);
    for (int i=0; i<n; i++)
      m_dist[i].resize(n, 0.);
  }
  double Get(int i, int j) const {
    const svec<double> & d = m_dist[i];
    return d[j];
  }
  void Set(int i, int j, double v) {
    svec<double> & d = m_dist[i];
    d[j] = v;
  }

private:
  svec< svec<double> > m_dist;
};

void Iterate(svec<Cactus> & cacti)
{
  int i, j;
  DistMatrix dm(cacti.isize());
  
  double min = 99999.;
  int minI = -1;
  int minJ = -1;
  for (i=0; i<cacti.isize(); i++) {
    for (j=i+1; j<cacti.isize(); j++) {
      double d = cacti[i].Distance(cacti[j]);
      dm.Set(i, j, d);
      if (d < min) {
	min = d;
	minJ = j;
	minI = i;
      }
    }
  }
  cout << "Merging " << cacti[minI].Name() << " weight " << cacti[minI].Weight() << " with " << cacti[minJ].Name() << " weight " << cacti[minJ].Weight() << " dist=" << min << endl;

  cacti[minI].Merge(cacti[minJ]);
  cacti[minJ] = cacti[cacti.isize()-1];
  cacti.resize(cacti.isize()-1);
  
  //svec<Cactus> out;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-o","output file");
  commandArg<string> cCmmd("-i","cacti file");
  commandArg<int> nCmmd("-n","number of final cacti");
  commandLineParser P(argc,argv);
  P.SetDescription("Clusters cacti by similarity.");
  P.registerArg(cCmmd);
  P.registerArg(fileCmmd);
  P.registerArg(nCmmd);
  
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string cName = P.GetStringValueFor(cCmmd);
  int n = P.GetIntValueFor(nCmmd);
  
  int i, j;
  svec<Cactus> cacti;
  LoadCacti(cacti, cName);

  while(cacti.isize() > n) {
    cout << "Iterating: " << cacti.isize() << endl;
    Iterate(cacti);
  }

  cout << "Saving cacti." << endl;
  for (i=0; i<cacti.isize(); i++)
    cacti[i].SetName(i);

  SaveCacti(cacti, fileName);

  return 0;
}
