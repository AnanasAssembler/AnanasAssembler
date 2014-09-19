#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "math/Spearman.h"

/*void LoadTissues(svec<string> & t) {
  FlatFileParser parser;
  parser.Open("tissues_short");
  parser.ParseLine();
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    t.push_back(parser.Line());
  }
  }*/

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<bool> ogCmmd("-outgroup","add dummy outgroup", false);
  //commandArg<string> file2Cmmd("-i2","input file 2", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing Spearman's rho.");
  P.registerArg(fileCmmd);
  P.registerArg(ogCmmd);
  //P.registerArg(file2Cmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  bool bOut = P.GetBoolValueFor(ogCmmd);
  //string fileName2 = P.GetStringValueFor(file2Cmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<string> t;
  //LoadTissues(t);

  int i, j;
  //for (i=0; i<t.isize(); i++)
  //  cout << t[i] << "\t";
  //cout << endl;
  //cout << t.isize() << endl;
  svec< svec< double > > v;
  //v.resize(t.isize());
  //int k = 0;
  int off = 1;

  svec<string> output;

  //if (bOut)
  //t.push_back("out");

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.GetItemCount() == 1) {
      t.push_back(parser.AsString(0));
      parser.ParseLine();
      svec<double> tmp;
      v.push_back(tmp);
      //cout << "Added " << t[t.isize()-1] << endl;
      continue;
    }
    svec<double> & s = v[t.isize()-1];
    for (i=1; i<parser.GetItemCount(); i++) {
      s.push_back(parser.AsFloat(i));
      //if (i == 12+off)
      //cout << "get " << parser.AsFloat(i) << " " << i << endl;
    }
    //k++;
  }

  if (bOut) {
    char one[4096];
    string line;
    //cout << t[i];
    sprintf(one, "outgroup");
    line += one;
    for (j=strlen("outgroup"); j<11; j++) {
      sprintf(one, " ");
      line += one;
    }
    for (j=0; j<t.isize()+1; j++) {
      if (j == 0)
	sprintf(one, "0.000");
      else
	sprintf(one, "1.000");
      line += one;
      if (j<t.isize()) {
	sprintf(one, " ");
       line += one;
      }
    }
    sprintf(one, "\n");
    line += one;
    output.push_back(line);
  }




  //cout << "Loaded." << endl;
  for (i=0; i<t.isize(); i++) {
    char one[4096];
    string line;
    //cout << t[i];
    sprintf(one, "%s", t[i].c_str());
    line += one;
    for (j=strlen(t[i].c_str()); j<11; j++) {
      sprintf(one, " ");
      line += one;
    }

    if (bOut) {
      sprintf(one, "1.000");
      line += one;     
      sprintf(one, " ");
      line += one;      
    }

    for (j=0; j<t.isize(); j++) {
      //svec<double> x, y;
      svec<double> x = v[i];
      svec<double> y = v[j];
      //for (int l = 0; l<k; l++) {
	//if (a[l] > 1.0 && b[l] > 1.0) {
	//x.push_back(a[l]);
	//y.push_back(b[l]);
	  //if (i==12 && j == 10) {
	  //  cout << "push " << a[l] << " " << b[l] << endl;
	  //}
	  //}
      //}

      SpearmansRho r;
      double rho = r.Compute(x, y);
      double dist = 1 - rho;
      //double dist = rho;
      //cout << "\t" << dist;
      //cout << endl << t[i] << " " << t[j] << " points: " << x.isize() << " " << i << " " << j << endl;
      //if (t[i] == "lung-d" && t[j] == "liver-d") {
      //for (int z=0; z<x.isize(); z++)
      //	  cout << x[z] <<  " " << y[z] << endl;
      //}


      sprintf(one, "%1.3f", dist);
      line += one;
      if (j<t.isize()-1) {
	sprintf(one, " ");
       line += one;
      }
     //cout << "Rho=" << rho << " Significance=" << r.Significance() << endl;
      
    }
    sprintf(one, "\n");
    line += one;
    output.push_back(line);
  }
  
  cout << output.isize() << endl;
  for (i=0; i<output.isize(); i++)
    cout << output[i];

  return 0;
}
