#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "math/Spearman.h"

void Read(svec<string> & names, const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    names.push_back(parser.AsString(0));
  }
  Sort(names);
}


double Tau(svec<double> & d) 
{
  double n = d.isize();
  double m = 0.;
  int i;
  for (i=0; i<d.isize(); i++) {
    if (d[i] > m)
      m = d[i];
  }
  if (m < 1.0)
    return -1;
  double t = 0.;
  for (i=0; i<d.isize(); i++) {
    t += (1. - d[i]/m);
  }
  //cout << n << " " << t << endl;
  //for (i=0; i<d.isize(); i++)
  //cout << d[i] << "\t";
  //cout << "Tau: " << t /(n-1) <<endl;
  return t /(n-1);
  
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","FPKM file");
  commandArg<string> listCmmd("-l","list file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints FPKM stats.");
  P.registerArg(fileCmmd);
  P.registerArg(listCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string listName = P.GetStringValueFor(listCmmd);
  
  svec<string> names;
  Read(names, listName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<double> polyA_all, dsn_all;
  int i;
  int s = 4;
  
  double polyA_spec = 0.;
  double dsn_spec = 0.;
  double nn = 0;

  parser.ParseLine();

  cout << "Read " << names.isize() << " entries." << endl;

  svec<string> label;
  for (i=0; i<parser.GetItemCount(); i++) {
    label.push_back(parser.AsString(i));
  }

  svec<double> a_brainA, a_brainA2, a_brainD, a_kidneyA, a_kidneyA2, a_kidneyD;
  
  double ex_brain = 0.;
  double ex_brain_div = 0.;

  double polyA_count = 0.;
  double dsn_count = 0.;
  int polyA_only = 0;
  int dsn_only = 0;
  int both = 0;

  double tau_sum = 0.;
  double tau_div = 0.;
  double tau_1 = 0.;
  double tau_98 = 0.;
  double tau_95 = 0.;
  double tau_60 = 0.;

  svec<int> superspec;
  

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const string & n = parser.AsString(0);
    int index = BinSearch(names, n);
    if (index < 0)
      continue;

    svec<double> polyA, dsn;
    svec<double> for_tau;
   //cout << parser.GetItemCount() << endl;

    double brainA = 0.;
    double brainA2 = 0; 
    double brainD = 0.;
    double kidneyA= 0.;
    double kidneyA2 = 0.; 
    double kidneyD = 0.;
    
    for (i=0; i<parser.GetItemCount(); i++) {
      if (strstr(label[i].c_str(), "polyA_FPKM") != NULL) {
	//cout << label[i] << " " << parser.AsFloat(i) << endl;
	polyA.push_back(parser.AsFloat(i));
	if (label[i] != "brain2.polyA_FPKM" && label[i] != "kidney2.polyA_FPKM") {
	  for_tau.push_back(parser.AsFloat(i));
	}
      }
      if (strstr(label[i].c_str(), "dsn_FPKM") != NULL) {
	dsn.push_back(parser.AsFloat(i));
      }

      if (label[i] == "brain2.polyA_FPKM") 
	brainA2 = parser.AsFloat(i);
      if (label[i] == "brain.polyA_FPKM")
      //if (label[i] == "liver.polyA_FPKM")
	brainA = parser.AsFloat(i);
 	
      if (label[i] == "brain.dsn_FPKM")
	brainD = parser.AsFloat(i);

      if (label[i] == "kidney2.polyA_FPKM")
	kidneyA2 = parser.AsFloat(i);
      if (label[i] == "kidney.polyA_FPKM")
	kidneyA = parser.AsFloat(i);
      if (label[i] == "kidney.dsn_FPKM")
	kidneyD = parser.AsFloat(i);

    }
    if (kidneyA > 1. || kidneyA2 > 1. /*|| kidneyD > 1.*/) {
      a_kidneyA.push_back(kidneyA);
      a_kidneyA2.push_back(kidneyA2);
      a_kidneyD.push_back(kidneyD);
      //cout << "Kidney " << kidneyA << " " << kidneyA2 << " " << kidneyD << endl;
    } 
    if (brainA > 1. || brainA2 > 1. /*|| brainD > 1.*/) {      
      a_brainA.push_back(brainA);
      a_brainA2.push_back(brainA2);
      a_brainD.push_back(brainD);
    }
    if (brainA > 1. || brainA2 > 1.) {
      ex_brain_div += 1.;
      if (brainA * brainA2 == 0.)
	ex_brain += 1.;
    }

    //cout << polyA.isize() << " " << dsn.isize() << endl;
    nn += 1.;

    double sum = 0.;
    double hi = 0.;
    for (i=0; i<polyA.isize(); i++) {
      sum += polyA[i];
      //cout << polyA[i] << " ";
      if (polyA[i] > hi)
	hi = polyA[i];
    }
    //cout << endl;
    //sum -= hi;
    sum /= (double)polyA.isize();
    if (hi > 5000000) {
      sum = 1.;
      hi = 1.;
    }

    bool polyA_b = false;
    bool dsn_b = false;
    bool both_b = false;
    //cout << sum << endl;
    if (hi > 0.1) {
      //polyA_spec += 1.-sum/hi/(1.-1/(double)polyA.isize());
      //double val = 1.-sum/hi/(1.-1/(double)polyA.isize());
      double val = sum/hi;
      polyA_spec += val;
      //if (val < 0)
      //cout << "Value: " << val << endl;
      polyA_count += 1.;
      polyA_b = true;
    }
    //sum = log(1.+sum);
    polyA_all.push_back(sum);

    //cout << sum << endl;

    hi = 0.;
    sum = 0.;
    for (i=0; i<dsn.isize(); i++) {
      sum += dsn[i];
      //cout << dsn[i] << " ";
      if (dsn[i] > hi)
	hi = dsn[i];
    }
    //cout << endl;
    sum /= dsn.isize();
    if (hi > 5000000) {
      sum = 1.;
      hi = 1.;
    }
    if (hi > 0.1) {
      dsn_spec += 1.-sum/hi/(1.-1/(double)dsn.isize());
      //cout << 1.-sum/hi/(1.-1/(double)dsn.isize()) << endl;
      dsn_count += 1.;
      dsn_b = true;
    }

    if (polyA_b && dsn_b) {
      both ++;
    } else {
      if (polyA_b) 
	polyA_only++;
      if (dsn_b) 
	dsn_only++;
    }
    //sum = log(1. + sum);
    dsn_all.push_back(sum);
    double t = Tau(for_tau);

    
    superspec.resize(for_tau.isize(), 0);
    if (t > 0.999) {
      double max_ex = 0.;
      int max_ex_i = -1;
      for (int zz=0; zz<for_tau.isize(); zz++) {
	if (for_tau[zz] > max_ex) {
	  max_ex = for_tau[zz];
	  max_ex_i = zz;
	}
      }
      if (max_ex_i >= 0)
	superspec[max_ex_i]++;
    }


    if (t >= 0.) {
      tau_sum += Tau(for_tau);
      tau_div += 1.;
      if (t > 0.9999)
	tau_1 += 1.;
      if (t > 0.98)
	tau_98 += 1.;
      if (t > 0.95)
	tau_95 += 1.;

      if (t < 0.6)
	tau_60 += 1.;

    }
   
  }

  
  cout << "Specificity: polyA " << polyA_spec/polyA_count << " DSN " << dsn_spec/dsn_count << endl;
  cout << "Tau avg. (polyA): " << tau_sum / tau_div << endl;
  cout << "Tau 100  (polyA): " << tau_1 / tau_div << endl;
  cout << "Tau  98  (polyA): " << tau_98 / tau_div << endl;
  cout << "Tau  95  (polyA): " << tau_95 / tau_div << endl;
  cout << "Tau <60  (polyA): " << tau_60 / tau_div << endl;
  cout << endl;
  cout << "tau_div: " << tau_div << endl;
  Sort(polyA_all);
  Sort(dsn_all);
  cout << nn/2 << endl;
  cout << polyA_all.isize() << endl;
  cout << dsn_all.isize() << endl;

  cout << "Superspecific: " << endl;
  for (i=0; i<superspec.isize(); i++)
    cout << superspec[i] << endl;
  cout << endl;

  int percent = (75*nn)/100;
  cout << "75%ile FPKM: polyA " << polyA_all[percent] << " DSN " << dsn_all[percent] << endl;
  double ssA = 0.;
  double ssD = 0.;

  double nnn = 0.;
  for (i=0; i<polyA_all.isize(); i++) {
    //if (polyA_all[i] > 1.0 || dsn_all[i] > 1.0) {
    ssA += polyA_all[i];
    ssD += dsn_all[i];
    nnn +=1. ;
    //}
  }

  cout << "PERCENT polyA-only %: " << 100.*((double)polyA_only)/(double)names.isize() << " ";
  cout << "dsn-only %: " << 100.*((double)dsn_only)/(double)names.isize() << " ";
  cout << "both %: " << 100.*((double)both)/(double)names.isize() << endl;
  cout << "FULL: " << polyA_only << " " << dsn_only << " " << both << endl;

  cout << "Average FPKM: polyA " << ssA/nnn << " DSN " << ssD/(double)nnn << endl;
  cout << "Correlations:" << endl;
  SpearmansRho r;
  double rho = r.Compute(a_brainA, a_brainA2);
  cout << "Brain polyA vs. polyA2: Rho=" << rho << " Significance=" << r.Significance() << endl;

  rho = r.Compute(a_brainA, a_brainD);
  cout << "Brain polyA vs. DSN: Rho=" << rho << " Significance=" << r.Significance() << endl;

  rho = r.Compute(a_kidneyA, a_kidneyA2);
  cout << "Kidney polyA vs. polyA2: Rho=" << rho << " Significance=" << r.Significance() << endl;

  rho = r.Compute(a_kidneyA, a_kidneyD);
  cout << "Kidney polyA vs. DSN: Rho=" << rho << " Significance=" << r.Significance() << endl;

  cout << "Null ratio: " << ex_brain / ex_brain_div << endl;

  FILE * x = fopen("x", "w");
  FILE * y = fopen("y", "w");

  for (i=0; i<a_brainA.isize(); i++) {
    fprintf(x, "%f\n", log(a_brainA[i]+1.));
    fprintf(y, "%f\n", log(a_brainA2[i]+1.));
  }

  fclose(x);
  fclose(y);
  return 0;
}
