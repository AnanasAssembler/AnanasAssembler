#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

void Read(svec<string> & l, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    l.push_back(parser.AsString(0));
  }
  Sort(l);
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","satsuma input file");
  commandArg<string> listCmmd("-l","list of sequences");
  commandLineParser P(argc,argv);
  P.SetDescription("Clusters sequences by sequence similarity.");
  P.registerArg(fileCmmd);
  P.registerArg(listCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string listName = P.GetStringValueFor(listCmmd);
  
  svec<string> all;
  
  Read(all, listName);

  svec<int> cluster;
  cluster.resize(all.isize(), -1);
  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int k = 0;
  int i;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & n1 = parser.AsString(3);
    const string & n2 = parser.AsString(0);
    if (n1 == n2 )
      continue;

    double ident = parser.AsFloat(6);
    int len = parser.AsInt(2) - parser.AsInt(1);

    if (len < 200 ||  ident < 0.52)
      continue;
  
    long long i1 = BinSearch(all, n1);
    long long i2 = BinSearch(all, n2);
    
    int c1 = cluster[i1];
    int c2 = cluster[i2];
    if (c1 == -1) {
      if (c2 != -1) {
	cluster[i1] = c2;	
	//cout << "Case 1 " << all[i1] << " " << all[i2] << " " << c2 << endl;
      } else {
	cluster[i1] = k;
	cluster[i2] = k;
	//cout << "New " << all[i1] << " " << all[i2] << " " << k << endl;
 	k++;
      }
      continue;
    }
    if (c2 == -1) {
      cluster[i2] = c1;
      //cout << "Case 2 " << all[i1] << " " << all[i2] << " " << c1 << endl;
      continue;
    }
    // Update all indices
    if (c1 != c2) {
      //cout << "Update " << all[i1] << " " << all[i2] << " " << c2 << " -> " << c1 << endl;
      for (i=0; i<all.isize(); i++) {
	if (cluster[i] == c2)
	  cluster[i] = c1;
      }
    }
  }
  /*
  cout << "Print clusters: " << endl;
  for (int j=0; j<all.isize(); j++) {
    cout << all[j] << " " << cluster[j] << endl;
    }*/

  for (i=0; i<k; i++) {
    cout << "Cluster: " << i << endl;
    int n = 0;
    for (int j=0; j<all.isize(); j++) {
      if (cluster[j] == i) {
	cout << all[j] << endl;
	n++;
      }
    }
    cout << "Members: " << n << endl << endl;
  }

  return 0;
}
