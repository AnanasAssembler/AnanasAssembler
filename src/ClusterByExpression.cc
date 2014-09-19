#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "math/SquareCorr.h"
#include "src/DNAVector.h"
#include <math.h>


int sub = 0;

class Locus
{
public: 
  Locus() {
    m_exp.resize(12, 0.);
  }

  bool Read(FlatFileParser & f) {
    static int count = 0;

    //cout << "Count: " << count << endl;
    count++;
    if (f.GetItemCount() != 19-sub) {
      //cout << "Not 19." << endl;
      return false;
    }

    m_chr = f.AsString(2-sub);
    m_start = f.AsInt(3-sub);
    m_stop = f.AsInt(4-sub);
    m_ori = f.AsInt(5-sub);

    //if (m_stop - m_start < 200)
    //return false;

    bool b = false;
    m_name = f.AsString(1-sub);   
    for (int i=6-sub; i<f.GetItemCount()-1; i++) {
      m_exp[i-(6-sub)] = f.AsFloat(i);
      if (m_exp[i-(6-sub)] > 1.)
	b = true;
      if (m_exp[i-(6-sub)] > 200000) {
	//cout << "Too high." << endl;
	return false;
      }
    }  
    //if (!b)
    //cout << "Failed." << endl;
    return b;
  }

  const svec<double> & GetExp() {return m_exp;}

  double Dist(const Locus & l) const {
    SquareCorrelation corr;
    return corr.Compute(m_exp, l.m_exp);
    double d = 0.;
    for (int i=0; i<m_exp.isize(); i++) {
      //if (i == 7)
      //continue;
      d += (m_exp[i] - l.m_exp[i]) * (m_exp[i] - l.m_exp[i]);
    }
    return d;
  }
 
  void Sub(double v) {
    for (int i=0; i<m_exp.isize(); i++)
      m_exp[i] += v;
  }
  void Add(const Locus & v) {
    double val = 0.;
    for (int i=0; i<m_exp.isize(); i++) {
      if (v.m_exp[i] > val)
	val = v.m_exp[i];
    }
    val = sqrt(val);
    for (int i=0; i<m_exp.isize(); i++) {
      if (v.m_exp[i] < 0)
	continue;
      m_exp[i] += v.m_exp[i] /*/val*/;
    }
  }

  double Specific() const {
    double avg = 0.;
    double hi = 0.;
    for (int i=0; i<m_exp.isize(); i++) {
      avg += m_exp[i];
      if (m_exp[i] > hi)
	hi = m_exp[i];
    }
    avg /= (double)m_exp.isize();
    return hi /avg;
 
  }

  void Div(double v) {
    if (v < 1.)
      return;
    for (int i=0; i<m_exp.isize(); i++)
      m_exp[i] /= v;
  }
  void Print() const {
    cout << m_name;
    for (int i=0; i<m_exp.isize(); i++) {
      cout << "\t" << m_exp[i];
    }
    cout << endl;
  }
  const string & Name() const {return m_name;}
  const string & Chr() const {return m_chr;}
  int Start() const {return m_start;}
  int Stop() const {return m_stop;}
  int Ori() const {return m_ori;}
  double Highest() const {
    double h = 0.;
    for (int i=0; i<m_exp.isize(); i++) {
      if (m_exp[i] > h)
	h = m_exp[i];
    }
    return h;
  }

private:
  svec<double> m_exp;
  string m_name;
  string m_chr;
  int m_start;
  int m_stop;
  int m_ori;

};


class Means
{
public: 
  Means() {
    Locus l;   
    m_points.push_back(l);
  }

  void Print(const svec<Locus> & data, const vecDNAVector & dna) {
    int i, j;
  
    svec<int> num, ind;
    int highest = Assign(num, ind, data);
    
    for (i=0; i<m_points.isize(); i++) {
      cout << "Center_" << i << " ";
      m_points[i].Print();
      cout << endl;
      
      FILE * p = NULL;
      if (dna.isize() > 0) {
	char name[256];
	sprintf(name, "cluster_%d.fasta", i);
	p = fopen(name, "w");
      }
      
      for (j=0; j<data.isize(); j++) {
	if (ind[j] == i) {
	  data[j].Print();
	  if (dna.isize() > 0) {
	    int c = PrintDNA(p, data[j], dna);
	    //if (c >= 0) 
	    //cout << "CpG stats: " << data[j].Specific() << " " << c << endl;
	  }
	}
      }
      if (p != NULL)
	fclose(p);
    }
  }

  int PrintDNA(FILE * p, const Locus & l, const vecDNAVector & dna) const {
    
    fprintf(p, ">%s %f\n", l.Name().c_str(), l.Highest());

    int index = dna.NameIndex(l.Chr());
    //cout << "Index: " << index << " name: " << l.Chr() << endl;
    const DNAVector & d = dna[index];
    int left = 1000;
    int right = 500;
    int start = l.Start() - left;
    int stop = l.Start() + right;
    if (l.Ori() == -1) {
      stop = l.Stop() + left;
      start = l.Stop() - right;
    }

    if (start < 0 || stop >= d.isize())
      return -1;
    DNAVector sub;
    sub.SetToSubOf(d, start, stop-start);
    if (l.Ori() == -1) {
      sub.ReverseComplement();
    }    
    int cpg = 0;
    for (int i=1; i<sub.isize(); i++) {
      if (sub[i-1] == 'C' && sub[i] == 'G')
	cpg++;
    }
    for (int i=0; i<sub.isize(); i++)
      fprintf(p, "%c", sub[i]);
    fprintf(p, "\n");
    return cpg;
  }

  void Split(const svec<Locus> & data) {
    int i, j;
    svec<int> num, ind;
    int highest = Assign(num, ind, data);
    cout << "Highest: " << highest << endl;
    Locus n = m_points[highest];
    n.Sub(0.1);
    m_points[highest].Sub(-0.1);
    m_points.push_back(n);
    svec<int> a;
    Assign(num, a, data);
    svec<double> counts;
    counts.resize(m_points.isize(), 0.);
    for (i=0; i<m_points.isize(); i++) {
      for (j=0; j<data.isize(); j++) {
	if (a[j] == i) {
	  counts[i] += 1.;
	  m_points[i].Add(data[j]);
	}
      }
    }
    for (int i=0; i<m_points.isize(); i++) {
      m_points[i].Div(counts[i]);
      cout << "Centroid: " << i << " counts: " << counts[i] << endl;
    }
    
  }
  
  int Assign(svec<int> & num, svec<int> & ind, const svec<Locus> & data) {
    int i, j;
    ind.resize(data.isize(), -1);
    num.resize(m_points.isize());
    int hi = 0;
    int highest = -1;
    for (i=0; i<data.isize(); i++) {
      double min = data[i].Dist(m_points[0]);
      int index = 0;
      for (j=0; j<m_points.isize(); j++) {
	double d = data[i].Dist(m_points[j]);
	if (d < min) {
	  min = d;
	  index = j;
	}
      }
      //cout << index << " " << min << endl;
      //data[i].Print();
      ind[i] = index;
      num[index]++;
      if (num[index] > hi) {
	hi = num[index];
	highest = index;
      }
    }
    return highest;
  }
 
private:
  svec<Locus> m_points;
};

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> fastaCmmd("-f","fasta file (prints promoters)", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the file parser.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);

  vecDNAVector dna;
  if (fastaName != "")
    dna.Read(fastaName);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  parser.ParseLine();

  svec<Locus> data;
  svec<Locus> means;

  while (parser.ParseLine()) {
    Locus l;
    if (l.Read(parser))
      data.push_back(l);
  }

  cout << "Data points: " << data.isize() << endl;

  Means m;
  for (int i=0; i<100; i++) {
    cout << "Split " << i << endl;
    m.Split(data);
  }
  m.Print(data, dna);

  return 0;
}
