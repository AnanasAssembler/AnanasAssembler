#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>

double CDF(double x, double m, double s)
{
  
  double r = 0.5 * (1. + erf((x-m)/s/1.414213562));
  return r;
}

double Sigma(double p, int N)
{
  return sqrt(p * (1. - p) * (double)N);
}

double GetExpect(double ident,
		 int length,
		 double expect_ident,
		 int tries)
{
  double p_match = expect_ident;
    
  double scale = 4.0; // Adjust standard deviation for normal distribution 
  int scaledLength = (double)length/scale;
  if(length<2) { scaledLength = 1; } // Added to prevent 0 sigma for length<2 TODO - review
  double s = expect_ident * scale * Sigma(p_match, scaledLength);

  //cout << "Sigma=" << s << endl; 
      
  double m = p_match * (double)length;
  double x = (double)length * ident;
    
  //cout << "len=" << length << " ident=" << ident << " s=" << s << endl;

  double cdf = CDF(m, x, s); // Cumulative distribution function

  double targetSize = tries;

  double expect = cdf * targetSize; // Adjust for sequence lengths
  
  double p_val = -expm1(-expect); // Poisson distribution

  cout << "s=" << s << " cdf=" << cdf << " size=" << targetSize << " expect=" << expect << " p=" << p_val << endl; 
  cout << "expect=" << expect << " p-val=" << cdf << endl;
  return p_val;
}


double FPRate(double len, double ident, int qLen, int tLen) 
{

  //double ident_mod = (len * ident - 5)/len;
  //if (ident_mod < 0.)
  //ident_mod = 0.;
  double tries = qLen * tLen;
  tries *= 2.;
  if (tries < 0)
    tries = 1;

  double dd = GetExpect(ident, len, 0.43, tries);
  cout << ident << " " << len << " " << dd << " " << tries << endl;
  return dd;
}

int main( int argc, char** argv )
{

 
  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  
  int i, j;
  
  /*for (int i=100; i<4000; i+=100) {
    double r = FPRate(100, 0.3, i, 4000);  
    cout << i << " " << r << endl;
  }  
  return 0;*/


  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  double len = 0.;
  double ident = 0.;

  svec<double> scores;
  scores.resize(201, 0.);
  double scale = (double)scores.isize()-1.;
  double n = 0.;
  int tLen = 0;
  int qLen = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "Target" && parser.AsString(1) == "sequence") {
      tLen = parser.AsInt(3);
    }
    if (parser.AsString(0) == "Query" && parser.AsString(1) == "sequence") {
      qLen = parser.AsInt(3);
    }
    if (parser.AsString(0) == "Target" && parser.AsString(1) == "aligned") {
      len = parser.AsFloat(3);
    }
    if (parser.AsString(0) == "Identity") {
      ident = parser.AsFloat(2);
      double fp = FPRate(len, ident, qLen, tLen);    
      //cout << "FP: " << fp << endl;
      int index = (int)(fp * scale);
      scores[index] += 1.;
      n += 1.;
    }

  }
  double sum = 0.;

  for (i=0; i<scores.isize(); i++) {
    sum += scores[i]/n;
    cout << ((double)i)/scale << "\t" << sum << endl;
  }

  return 0;
}
