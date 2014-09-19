#include "src/Smorgas/ProtNameCluster.h"
#include "src/Cola/Cola.h"
#include "base/StreamParser.h"


class Pair
{
public:
  Pair(const string & s, int n){
    m_s = s;
    m_n = n;
  }

  const string & S() const {return m_s;}
  int N() const {return m_n;}

  bool operator < (const Pair & p) const {
    return m_n < p.m_n;
  }

private:
  int m_n;
  string m_s;
};

void ProteinNameClusterer::Add(const string & name)
{
  string clean;
  int i;
  int n = 0;
  for (i=0; i<(int)name.size(); i++) {
    if (name[i] == '|') {
      n++;
      if (n == 4)
	break;
    }
  }
  i++;
  if (name[i] == '_')
    i++;
  for (; i<(int)name.size(); i++) {
    if (name[i] == '[')
      break;
    clean += (char)tolower(name[i]);
  }

  m_names.push_back(clean);
}

void ProteinNameClusterer::Cluster()
{

  int i, j;
  for (i=0; i<m_names.isize(); i++) {
    for (j=i+1; j<m_names.isize(); j++) {
      string result;
      if (AlignTwo(result, m_names[i], m_names[j])) {
	AddReplace(result);
      }
    }
  }

  svec<Pair> pair;
  for (i=0; i<m_result.isize(); i++) {
    for (j=0; j<m_names.isize(); j++) {
      string rrr;
      if (AlignTwo(rrr, m_result[i], m_names[j])) {
	m_counts[i]++;
      }
    }
    pair.push_back(Pair(m_result[i], m_counts[i]));
  }
  Sort(pair);
  for (i=0; i<pair.isize(); i++) {
    m_result[i] = pair[i].S();
    m_counts[i] = pair[i].N();
  }

  for (i=0; i<m_result.isize(); i++) {
    string & tmp = m_result[i];
    for (j=0; j<(int)tmp.size(); j++) {
      if (tmp[j] == '_')
	tmp[j] = ' ';
    }
  }
}

void ProteinNameClusterer::AddReplace(const string & s)
{
  int i;
  //cout << "Enter AddReplace w/ " << s << endl;
  for (i=0; i<m_result.isize(); i++) {
    string r;
    if (AlignTwo(r, m_result[i], s)) {
      if (s.size() < m_result[i].size())
	m_result[i] = s;
      //m_counts[i]++;
      //cout << "Found: " << s << endl;
      return;
    }
  }
  //cout << "Adding: " << s << endl; 
  m_result.push_back(s);
  m_counts.push_back(0);
}

bool ProteinNameClusterer::AlignTwo(string & result, const string & a, const string & b)
{
  if (a == b) {
    result = a;
    return true;
  }
  //cout << "Trying " << a << " " << b << endl;
  Cola cola1 = Cola();
  result = "";
  DNAVector target, query;
  target.SetFromBases(a);
  query.SetFromBases(b);

  cola1.createAlignment(target, query, AlignerParams(-1, AlignerType(1)));
  Alignment cAlgn = cola1.getAlignment();
  double maxP = 0.1;
  double minIdent = 0.5;
  
  int i;

  if (cAlgn.calcPVal()<=maxP && cAlgn.getIdentityScore()>=minIdent) {
    //cout << a << " vs " << b << endl;
    stringstream tmp;
    cAlgn.print(0,1,tmp,100);
    //cAlgn.print(0,1,cout,100);
    //cout << tmp.str() << endl;
    StreamParser p;
    p.Set(tmp);
    while (p.ParseLine()) {
      if (p.GetItemCount() == 0)
	continue;
      if (p.AsString(0) == "Query:") {
	const string & s = p.AsString(2);
	result = s;
      }
    }

    return true;
  } else {
    //cout<<"No Alignment at given significance threshold"<<endl;
    return false;
  }

}
