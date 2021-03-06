#ifndef CONSENSUS_H
#define CONSENSUS_H

#include <string>
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/SVector.h"

class Hypothesis;
class ConsensOverlapUnit;
class Contig;
class Assembled;


class Consensus
{
public:
  Consensus() {
    m_data.resize(4);
  } // 4 characters : ACGT

  void resize(int n) {
    for (int i=0; i<m_data.isize(); i++) {
      m_data[i].resize(n, 0);
    }
    m_first.resize(n, -1);
  }

  void Add(int i, char c) {
    
    int l = NucIndex(c);
    if (l < m_data.isize() && l >= 0) {
      (m_data[l])[i] += 1;
    }
    if (m_first[i] == -1)
      m_first[i] = c;
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
    if (index == -1 || max==0)
      return 'N';
    return NucLetter(index);
  }
  char GetFirst(int i) {return m_first[i];}
  
private:
  svec< svec<int> > m_data;
  svec<char> m_first;
};


class ConsensusBuilder
{
 public:
  ConsensusBuilder() {}

  void Build(DNAVector & d, const Contig& cont, const ConsensOverlapUnit & reads);
  void BuildWithGaps(DNAVector & d, const Contig& cont, const ConsensOverlapUnit & reads);

};



class LayoutSink
{
 public:
  LayoutSink();
  ~LayoutSink();
 
  void SetLayoutFile(const string & layoutFile);

  void SetPrefix(string pref) {
    m_prefix = pref;
  }

  void SetIndex(int i) { 
    m_index = i;
  }
 
  void Dump(const Hypothesis & hyp, const ConsensOverlapUnit & reads, bool bPrev = false, int minLen = 0);

  void fastaFromAssembly(const string& fastaFile, Assembled& asmb, const ConsensOverlapUnit & COUnit, int minLen, bool bUseGaps, bool onlyTop);

  const string & LastContig() const {return m_lastContigName;}

 private:
  bool IsSame(const DNAVector & d);
  void End();

  FILE * m_pLayout;
  int m_counter;
  int m_minor;
  map<int, bool> m_currReads;        ///Keep count of uique raw reads
  map<int, bool> m_currPairedReads;  ///Keep count of unique raw pairs
  string m_lastContigName;
  string m_lastScaffName;
  DNAVector m_lastCons;
  double m_minIdent;
  int m_index;
  string m_prefix;
};


#endif //CONSENSUS_H
