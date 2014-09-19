#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

void Load (svec<double> & v, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  double max = 0.;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    v.push_back(parser.AsFloat(0));
    if (parser.AsFloat(0) > max)
      max = parser.AsFloat(0);
  }
  int i;
  double scale = 0.98;
  for (i=0; i<v.isize(); i++) {
    v[i] /= max;
    v[i] *= scale;    
  }
  
}

void Solve(double & x, double & y, svec<double> & m, svec<double> & pure)
{
  int i;
  double a1 = 0.;
  double a2 = 0.;
  double a3 = 0.;
  double b1 = 0.;
  double b2 = 0.;
  double b3 = 0.;

  double s1 = 0;
  double s2 = 0.;

  for (int i=0; i<m.isize(); i++) {
    double p = pure[i];
    double q = 1-p;
    a1 += p*p;
    a2 += p*q;
    a3 += m[i]*p;
    b1 += p*q;
    b2 += q*q;
    b3 += m[i]*q;

    s1 += p*p - m[i]*p;
    s2 += p*q - m[i]*q;
  }
  //cout << "s1=" << s1 << " s2=" << s2 << endl;

  //cout << " a1=" << a1 << " a2=" << a2 << " a3=" << a3;
  //cout << " b1=" << b1 << " b2=" << b2 << " b3=" << b3;
 
  //cout << endl;

  b2 /= b1;
  b3 /= b1;
  b1 /= b1;

  a2 /= a1;
  a3 /= a1;
  a1 /= a1;

  double r = a3 - b3;
  double s = a2 - b2;
  y = r/s;

  // y = (b3*a1 - a3) / (b2*a1 - a2);
  x = (a3 - y*a2) / a1;
  // x = 1.;
  //y = 0.;

  //cout << "x = " << x << " y = " << y << endl;
  // for (i=0; i<m.isize(); i++) {
  //  cout << "   " << pure[i] * x << " " << (1-pure[i]) * y << " sum= ";
  // cout << pure[i] * x + (1-pure[i]) * y << " obs=" << m[i] << endl;
  //}
    
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file (FPKM)");
  commandArg<string> pureCmmd("-p","purity file");
  commandArg<bool> bCmmd("-s","second population", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Normalizes FPKM values numerically by cell type.");
  P.registerArg(fileCmmd);
  P.registerArg(pureCmmd);
  P.registerArg(bCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  svec<double> purity;
  Load(purity, P.GetStringValueFor(pureCmmd));
  bool b = P.GetBoolValueFor(bCmmd);


  //svec<double> test;
  //double x1, y1;
 
  //Solve(x1, y1, purity, purity);
  //return 0;

  FlatFileParser parser;
  
  parser.Open(fileName);
  parser.ParseLine();
  cout << parser.AsString(0) << "\t" << parser.AsString(1);
  int i;
  for (i=2; i<parser.GetItemCount(); i++) {
    cout << "\t" << parser.AsString(i) << "\t" << parser.AsString(i) << "_norm";
  }
  cout << "\t" << "endocrine\tnon_endocrine" << endl;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    svec<double> n;
    for (i=2; i<parser.GetItemCount(); i++)
      n.push_back(parser.AsFloat(i));

    cout << parser.AsString(0) << "\t" << parser.AsString(1);
    double x, y;
    Solve(x, y, n, purity);
    for (i=0; i<n.isize(); i++) {
      cout << "\t" << parser.AsString(2+i);
      cout << "\t" << y * (1. - purity[i]) + x * purity[i];      
    }
    cout << "\t" << x << "\t" << y;
    cout << endl;
  }
  return 0;
}
