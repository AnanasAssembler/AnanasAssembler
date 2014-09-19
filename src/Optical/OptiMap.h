#ifndef OPTIMAP_H
#define OPTIMAP_H

#include "base/SVector.h"
#include "src/CrossCorr.h"

class OpticalSeq
{
 public:
  OpticalSeq() {}
  void SetName(const string & n) {
    m_name = n;
  }
  const string & Name() const {return m_name;}

  int isize() const {return m_data.isize();}
  int operator [] (int i) const {return m_data[i];}
  void push_back(int n) {
    m_data.push_back(n);
  }
  
  void Reverse() {
    throw;
    int i; 
    svec<int> tmp;
    tmp.resize(m_data.isize(), 0);
    int k = m_data.isize()-1;
    int last = m_data[k];
    for (i=0; i<m_data.isize(); i++) {
      tmp[i] = m_data[k];
      k--;
    }
    m_data = tmp;
 }

 private:
  svec<int> m_data;
  string m_name;
};

class OpticalMap
{
 public:
  OpticalMap() {}

  int isize() const {return m_data.isize();}
  const OpticalSeq & operator[] (int i) {return m_data[i];}
  void push_back(const OpticalSeq & s) {m_data.push_back(s);}

  double Average();
 private:
  svec<OpticalSeq> m_data;
};

//=====================================================
class OptiAligner
{
 public:
  OptiAligner() {
    m_size = 2048;
    //m_scale = 0.0025;
    m_decay = 0.1;
  }
  
  double Align(const OpticalSeq & a, const OpticalSeq & b, bool bPrint = false);
  double AlignAuto(const OpticalSeq & a, const OpticalSeq & b);
 private:
  double AlignInt(const OpticalSeq & a, const OpticalSeq & b, double scale, double thresh);
  void Fit(svec<float> & o, const OpticalSeq & a, double scale);

  int m_size;
  CrossCorrelation m_xc;
  // double m_scale;
  double m_decay;
};


#endif //OPTIMAP_H
