#ifndef PROTNAMECLUSTER_H
#define PROTNAMECLUSTER_H

#include <string>
#include "base/SVector.h"

class ProteinNameClusterer
{
 public:
  ProteinNameClusterer() {}
  void Clear() {
    m_names.clear();
    m_result.clear();
    m_counts.clear();
  }

  void Add(const string & name);

  void Cluster();

  int Num() const {return m_result.isize();}
  int Count(int i) const {return m_counts[i];}
  const string & Class(int i) const {return m_result[i];}

 private:
  bool AlignTwo(string & result, const string & a, const string & b);
  void AddReplace(const string & s);

  svec<string> m_names;
  svec<string> m_result;
  svec<int> m_counts;

  

};

#endif

