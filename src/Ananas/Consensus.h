#ifndef CONSENSUS_H
#define CONSENSUS_H

#include <string>
#include "src/DNAVector.h"
#include "base/SVector.h"

class Hypothesis;
class ConsensOverlapUnit;
class Contig;
class Assembled;


class Consensus
{
public:
  Consensus() {m_data.resize(4);} // 4 characters : ACGT

  void resize(int n) {
    for (int i=0; i<m_data.isize(); i++) {
      m_data[i].resize(n, 0);
    }
  }

  void Add(int i, char c) {
    int l = NucIndex(c);
    if (l < m_data.isize() && l >= 0) {
      (m_data[l])[i] += 1;
    }
  }

  int isize() const {return m_data[0].isize();}
  char operator [] (int pos) const {
    int i;
    int max = 0;
    int index = -1;
    for (i=0; i<m_data.isize(); i++) {
      int m = (m_data[i])[pos];
      if (m > max) {
	index = i;
	max = m;
      }	
    }
    if (index == -1)
      return 'N';
    return NucLetter(index);
  }
 
private:
  svec< svec<int> > m_data;
};


class ConsensusBuilder
{
 public:
  ConsensusBuilder() {}

  void Build(DNAVector & d, const Contig& cont, const ConsensOverlapUnit & reads);

};



class LayoutSink
{
 public:
  LayoutSink();
  ~LayoutSink();
 
  void SetLayoutFile(const string & layoutFile);
  void SetPrefix(int i) {
    m_prefix = i;
  }
 
  void Dump(const Hypothesis & hyp, const ConsensOverlapUnit & reads, bool bPrev = false, int minLen = 0);

  void fastaFromAssembly(const string& fastaFile, const Assembled& asmb, const ConsensOverlapUnit & COUnit);

  const string & LastContig() const {return m_lastContigName;}

 private:
  bool IsSame(const DNAVector & d);
  void End();

  FILE * m_pLayout;
  int m_counter;
  int m_minor;
  map<int, bool> m_currReads;
  map<int, bool> m_currPairedReads;
  string m_lastContigName;
  string m_lastScaffName;
  DNAVector m_lastCons;
  double m_minIdent;
  int m_prefix;
};


#endif //CONSENSUS_H
