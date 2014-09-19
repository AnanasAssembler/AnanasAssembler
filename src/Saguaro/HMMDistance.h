#ifndef HMMDISTANCE_H_
#define HMMDISTANCE_H_

#include "base/SVector.h"
#include <string>
#include <math.h>
#include <iostream>

class CMReadFileStream;
class CMWriteFileStream;

class HMMFeature
{
 public:
  HMMFeature() {
    m_n = 0;
    m_weight = 1.;
  }

  
  double Weight() const {return m_weight;}

  const string & Name() const {return m_chr;}
  int GetPosition() const {return m_pos;}
  int size() const {return m_data.isize();}
  int isize() const {return m_data.isize();}
  const double & operator[] (int i) const {return m_data[i];}

  void resize(int size, double v=0.) {
    m_data.resize(size, v);
  }

  double & operator[] (int i) {return m_data[i];}
  void SetName(const string & n) {m_chr = n;}
  void SetPosition(int pos) {m_pos = pos;}
  
  void Read(CMReadFileStream & f, int ver);
  void Write(CMWriteFileStream & f);

  void Reset() {m_n = 0;}

//  // Average over incoming features
//  void Update(const HMMFeature & f, HMMFeature & weights) {
//    int i;
//    //double weight = 1./((double)m_n + 1);
//    //cout << "Update weight=" << weight << endl;
//    for (i=0; i<m_data.isize(); i++) {
//      if (f[i] < -0.5)
//	continue;     
//      
//      double w = 1./(weights[i] + 1.);
//      m_data[i] = f[i] * w + m_data[i] * (1 - w);
//      weights[i] += 1.;
//
//      //if (m_data[i] < 0.) {
//      //cout << "ERROR: i=" << i << " m_data=" << m_data[i];
//      //cout << " w=" << weight << endl; 
//      //}
//    }
//    m_n++;
//  } 



  void Update(const HMMFeature & f, HMMFeature & same, HMMFeature & diff) {
    
    int i;
    
    for (i=0; i<m_data.isize(); i++) {
      if (f[i] < 0.)
	continue;     
      
      /*
      if (f[i] == 0) {
	same[i] += 1.;
      }

      if (f[i] == 1) {
        diff[i] += 1.;
	}*/

      diff[i] += f[i];
      same[i] += (1. - f[i]);

      if (same[i] == 0) {
        diff[i] == 0 ? m_data[i]=0 : m_data[i]=log(diff[i]);
      } else {
        m_data[i] = -log(same[i]/(same[i]+diff[i]));
      }
     
//     cout << i << " same = " << same[i] << " diff = " << diff[i] << " dist = " << m_data[i];
//     if (m_data[i] > 1) {cout << " * ";}
//     cout << endl;

    }
    m_n++;
  } 


private:
  svec<double> m_data;
  double m_weight;
  int m_pos;
  string m_chr;
  int m_n;
};



class CompressedFeature
{
 public:
  CompressedFeature() {}

  void SetSize(int n);
  void Set(const char * p);
  void SetFloat(const svec<double> & d);

  void MakeFeature(HMMFeature & f) const;
  void SimpleFeature(svec<double> & f);

  void Read(CMReadFileStream & f, int ver);
  void Write(CMWriteFileStream & f);

  //svec<char> & Data() {return m_data;} 
  int isize() const {return m_data.isize() + m_floats.isize();}
 private:
  svec<char> m_data;
  svec<double> m_floats;
};






class HMMFeatureVector
{
 public:
  HMMFeatureVector() {}
  ~HMMFeatureVector() {}


  int size() const {return m_data.isize();}
  int isize() const {return m_data.isize();}
  void SimpleFeature(svec<double> & f, int i) {
    m_comp[i].SimpleFeature(f);
  }

  const HMMFeature & operator[] (int i) const {
    if (m_data[i].isize() > 0 || m_comp.isize() == 0 ) {
      //cout << "Warning (const)" << endl;
      return m_data[i];
    } else {
      //cout << "OK (const)" << endl;
      m_comp[i].MakeFeature( ((HMMFeatureVector*)this)->m_temp);
      ( ((HMMFeatureVector*)this)->m_temp).SetName(m_data[i].Name());
      ( ((HMMFeatureVector*)this)->m_temp).SetPosition(m_data[i].GetPosition());
      return m_temp;
    }
  }
  
  void resize(int size) {
    m_data.resize(size);
    //m_comp.resize(size);
  }
  
  void push_back(const HMMFeature & f) {
    m_data.push_back(f);
  }

  HMMFeature & operator[] (int i) {
    if (m_data[i].isize() > 0 || m_comp.isize() == 0 || m_comp[i].isize() == 0) {
      //cout << "OK (1)" << endl;
      return m_data[i];
    } else {
      //cout << "Warning (1)" << endl;
      //cout << m_data[i].isize() << " " << m_comp.isize() << " " << m_comp[i].isize() << endl;
      m_comp[i].MakeFeature(m_temp);
      m_temp.SetName(m_data[i].Name());
      m_temp.SetPosition(m_data[i].GetPosition());
     return m_temp;
    }    
  }
 
  void SetCompressed(int i, const char * p, int n, const string & chr, int pos) {
    if (i >= m_comp.isize())
      m_comp.resize(m_data.isize());
    m_comp[i].SetSize(n);
    m_comp[i].Set(p);
    m_data[i].SetName(chr);
    m_data[i].SetPosition(pos);
  }

  void SetCompressedFloat(int i, const svec<double> & data, int n, const string & chr, int pos) {
    if (i >= m_comp.isize())
      m_comp.resize(m_data.isize());
    m_comp[i].SetSize(n);
    m_comp[i].SetFloat(data);
    m_data[i].SetName(chr);
    m_data[i].SetPosition(pos);
  }

  void Read(const string & name);

  // Reads more data but keeps what's in there
  void MergeRead(const string & name);
  void Write(const string & name);

  void SetName(int i, const string & name) {
    if (i>=m_names.isize())
      m_names.resize(i+1);
    m_names[i] = name;
  }

  int GetNameCount() const {return m_names.isize();}
  const string & GetName(int i) const {return m_names[i];}


 private:
  svec<HMMFeature> m_data;
  svec<CompressedFeature> m_comp;
  HMMFeature m_temp;
  svec<string> m_names;
};





class HMMDistance
{
 public:
  HMMDistance() {}
  virtual ~HMMDistance() {}

  virtual double Distance(const HMMFeature & v1, const HMMFeature & v2) = 0;

};

class HMMEuclidianDistance : public HMMDistance
{
 public:
  HMMEuclidianDistance() {}
  virtual ~HMMEuclidianDistance() {}

  virtual double Distance(const HMMFeature & v1, const HMMFeature & v2);

};

class HMMTreeDistance : public HMMDistance
{
 public:
  HMMTreeDistance();
  virtual ~HMMTreeDistance() {}

  virtual double Distance(const HMMFeature & v1, const HMMFeature & v2);
 private:
  svec<float> m_cache;
  int m_n;
  int m_shift;
  int Index(double d1, double d2);
};

class NameIndex
{
 public:
  NameIndex() {m_i = -1;}
  NameIndex(const string & s, int i) {
    m_s = s;
    m_i = i;
  }

  int Index() const {return m_i;}
  void SetIndex(int i) {m_i = i;}
  void SetName(const string & s) {m_s = s;}
  const string & String() const {return m_s;}

  bool operator < (const NameIndex & d) const {
    return m_s < d.m_s;
  }
  
 private:
  string m_s;
  int m_i;
};

//#define USE_CACHE

class HMMTrees
{
 public:
  HMMTrees();
  ~HMMTrees();


  void Reset() {
    for (int i=0; i<m_model.isize(); i++) {
      m_model[i].Reset();
    }
  }

  double GetScore(const string & name, const HMMFeature & f) {
    int i = Index(name);

#ifndef USE_CACHE
    return m_pDist->Distance(m_model[i], f);
#else
    if (m_cache.isize() == 0)
      return m_pDist->Distance(m_model[i], f);

    if (m_cache[i] < 0.) {
      m_cache[i] = m_pDist->Distance(m_model[i], f);
    }
    return m_cache[i];
#endif
  }

  void ClearCache() {
    if (m_cache.isize() == 0)
      m_cache.resize(m_model.isize());
    for (int i=0; i<m_cache.isize(); i++)
      m_cache[i] = -1;
  }


  int GetNumTrees() const {return m_names.isize();}
  const string & GetTreeName(int i) const {return m_names[i].String();}

//  void Update(const string & name, const HMMFeature & f) {
//    int i = Index(name);
//    m_model[i].Update(f, m_weights);
//  }

  void Update(const string & name, const HMMFeature & f) {
    int i = Index(name);
    m_model[i].Update(f, m_same[i], m_diff[i]);
  }

  void Update(int i, const HMMFeature & f) {
    m_model[i].Update(f, m_same[i], m_diff[i]);
  }

  void PrettyPrintOne(const string & name);

  void PrettyPrint();

  void PrettyPrint(const string & file);

  // Read and add from ASCII file
  void AddRead(const string & file);


  void Add(const HMMFeature & f, const string & name) {
    m_model.push_back(f);
    m_names.push_back(NameIndex(name, m_names.isize()));
    HMMFeature w;
    w.resize(f.isize(), 0.);
    m_weights.push_back(w);
    m_same.push_back(w); 
    m_diff.push_back(w); 
    m_bSorted = false;
  }


  void ResetWeights() {
    int i, j;
    for (i=0; i<m_weights.isize(); i++) {
      HMMFeature & f = m_weights[i];
      for (j=0; j<f.isize(); j++)
	f[j] = 0.;
    }
  }

  void ResetSameDiff() {
    int i, j;
    for (i=0; i<m_same.isize(); i++) {
      HMMFeature & f = m_same[i];
      for (j=0; j<f.isize(); j++)
	f[j] = 0.;
    }
    for (i=0; i<m_diff.isize(); i++) {
      HMMFeature & f = m_diff[i];
      for (j=0; j<f.isize(); j++)
	f[j] = 0.;
    }
  }
  
  void SetMatrixName(int i, const string & name) {
    m_model.SetName(i, name);
  }

  void Sort();

 private:
  int Index(const string & name) {
    if (!m_bSorted)
      Sort();

    /*
    int realIndex = -1;
    for (int i=0; i<m_names.isize(); i++) {
      if (m_names[i].String() == name)
        return i;
    }
    return -1;
    */

    int index = BinSearch(m_searchNames, NameIndex(name, -1));

    //if (realIndex != m_searchNames[index].Index())
    //cout << "ERROR!!" << endl;

    return m_searchNames[index].Index();
  }
		

  bool m_bSorted;
  svec<NameIndex> m_names;
  svec<NameIndex> m_searchNames;
  HMMFeatureVector m_model;
  HMMFeatureVector m_weights;
  HMMFeatureVector m_same;
  HMMFeatureVector m_diff;

  HMMDistance * m_pDist;

  svec<double> m_cache;

};






#endif //HMMDISTANCE_H_





