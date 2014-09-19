#ifndef CACTUS_H_
#define CACTUS_H_

#include "base/SVector.h"
#include <string>
#include <stdio.h>
#include "base/StringUtil.h"

class FlatFileParser;

class Cactus
{
 public:
  Cactus() {
    m_weight = 1.;
  }

  bool Read(FlatFileParser & f);

  const string & Name() const {return m_name;}
  void Print();
  void Write(FILE * p);

  double Distance(const Cactus & c);
  void Merge(const Cactus & c);
  double Weight() const {return m_weight;}

  double Get(int i, int j) const {
    const svec<double> & d = m_data[i];
    return d[j];
  } 
  void Set(int i, int j, double v) {
    svec<double> & d = m_data[i];
    d[j] = v;
  } 

  int Size() const {return m_label.isize();}

  void SetSize(int n) {
    m_data.resize(n);
    for (int i=0; i<n; i++)
      m_data[i].resize(n, 0.);
  }
 

  void SetName(const string & n) {
    m_name = n;
  }
  void SetName(int n) {
    m_name = "cactus" + Stringify(n);
  }
 private:
  string m_name;
  svec<string> m_label;
  svec< svec<double> > m_data;
  double m_weight;
};

void LoadCacti(svec<Cactus> & all, const string & name);
void SaveCacti(svec<Cactus> & all, const string & name);

//void Merge(Cactus & out, svec<Cactus> & all, svec<int> & indices);

#endif 

