#ifndef CONTSCAFF_H_
#define CONTSCAFF_H_

#include "ryggrad/src/general/DNAVector.h"


class ReadPlacement
{
 public:
  ReadPlacement() {
    m_start = -1;
    m_stop = -1;
    m_ori = 0;
    m_read = -1;
    m_paired = -1;
    m_pairOrient = 0;
  }
 
  ReadPlacement(int read, int ori, int start, int stop, int paired, int pairOrient, const string & name) {
    m_start = start;
    m_stop = stop;
    m_ori = ori;
    m_read = read;
    m_paired = paired;
    m_pairOrient = pairOrient;
    m_name = name;
  }

  
  int Start() const  {return m_start;}
  int Stop()  const  {return m_stop;}
  int Ori()   const  {return m_ori;}
  int Read()  const  {return m_read;}
  int Pair() const { return m_paired;} 
  int PairOrient() const { return m_pairOrient;} 
  const string & Name()  const {return m_name;}

  bool operator < (const ReadPlacement & r) const {
    return m_start < r.m_start;
  }

  void Reverse(int len) {
    int tmp = m_start;
    m_start = len - m_stop;
    m_stop = len - tmp;
    m_ori = -m_ori;
  }

 private:
  int m_start;
  int m_stop;
  int m_ori;
  int m_read;
  int m_paired;
  int m_pairOrient;
  string m_name;
};


class Contig
{
 public:
  Contig() {
    m_ori = 1;
    m_numOfReads = 0;
    m_numOfPaired = 0;
    m_discard = false;
  }
  ~Contig() {}
 
  int isize() const {return m_placements.isize();}
  const ReadPlacement & operator [] (int i) const {return m_placements[i];}
  ReadPlacement & operator [] (int i) {return m_placements[i];}
  void push_back(const ReadPlacement & s) {m_placements.push_back(s);} 

  const string & Name() const {return m_name;}
  void SetName(const string & n) {m_name = n;}

  int NumReads() const { return m_numOfReads; }
  void SetNumReads(int readCnt) { m_numOfReads = readCnt; }

  int NumPairs() const { return m_numOfPaired/2; }
  float Score() const  { return ((float) m_numOfPaired)/m_numOfReads; }
  void SetNumPairs(int pairCnt) { m_numOfPaired = pairCnt*2; }
  

  bool Discard() const { return m_discard; }
  void SetDiscard(bool dis) { m_discard = dis; }

  void clear() {
    m_placements.clear();
    m_name = "";
    m_numOfReads = 0;
    m_numOfPaired = 0;
    m_discard = false;
  }

  void Reverse() {
    int len = 0;
    int i;
    for (i=0; i<m_placements.isize(); i++) {
      if (m_placements[i].Stop() > len)
	len = m_placements[i].Stop();
    }
    for (i=0; i<m_placements.isize(); i++) {
      m_placements[i].Reverse(len);
    }
    Sort(m_placements);
    m_ori = -m_ori;
  }

  int Ori() const {return m_ori;}
  void SetOri(int o) {
    m_ori = o;
  }

  bool operator < (const Contig & c) const {
    return m_placements.isize() < c.m_placements.isize();
  }

  int Highest() const {
    int h = 0;
    for (int i=0; i<m_placements.isize(); i++) {
      if (m_placements[i].Stop() > h)
	h = m_placements[i].Stop();
    }
    return h;
  }

 private:

  svec<ReadPlacement> m_placements;
  string m_name;
  int    m_ori;
  int    m_numOfReads;
  int    m_numOfPaired;
  bool   m_discard;   /// Flag for allowing to discard/ignore the contig in different contexts
};

class Scaffold
{
 public:
  Scaffold() {}
  ~Scaffold() {}

  int isize() const {return m_contigs.isize();}
  const Contig & operator [] (int i) const {return m_contigs[i];}
  Contig & operator [] (int i) {return m_contigs[i];}
  
  void resize(int n) {
    m_contigs.resize(n);
    m_offset.resize(n);
  }


  void push_back(const Contig & s, int off, int ori) {
    m_offset.push_back(off);
    m_contigs.push_back(s);
    //if (ori == -1) {
    //m_contigs[m_contigs.isize()-1].Reverse();
    //}
  }
  
  void push_back(const Scaffold & s, int off, int ori) {
    if (ori == 1) {
      for (int i=0; i<s.isize(); i++) {
	push_back(s[i], off + s.Offset(i), ori);
      } 
    } else {
      throw;
    }
  }

  void clear() {
    m_contigs.clear();
    m_offset.clear();
  }

  const string & Name() const {return m_name;}
  void SetName(const string & n) {m_name = n;}

  int Offset(int i) const {return m_offset[i];}
  void SetOffset(int i, int o) {m_offset[i] = o;}

  int Length() const {
    int n = 0;
    for (int i=0; i<m_contigs.isize(); i++) {
      n += m_contigs[i].Highest();
    }
    return n;
  }

  int NumUniqReads() const          { return m_numUniqReads;    }
  void SetNumUniqReads(int readCnt) { m_numUniqReads = readCnt; }

  int NumUniqPairs() const          { return m_numUniqPairs;     }
  void SetNumUniqPairs(int pairCnt) { m_numUniqPairs = pairCnt;  }

  int NumReads() const {
    int n = 0;
    for (int i=0; i<m_contigs.isize(); i++) {
      n += m_contigs[i].NumReads();
    }
    return n;
  }

  bool operator < (const Scaffold & s) const {
    return NumReads() > s.NumReads();
  }

  void Reverse() {
    for (int i=0; i<m_contigs.isize(); i++)
      m_contigs[i].Reverse();
  }


 private:

  svec<Contig> m_contigs;
  svec<int> m_offset;
  string m_name;
  int  m_numUniqReads;
  int  m_numUniqPairs;
};


class Assembled
{
 public:
  Assembled() {}

  int isize() const {return m_scaffolds.isize();}
  const Scaffold & operator [] (int i) const {return m_scaffolds[i];}
  Scaffold & operator [] (int i) {return m_scaffolds[i];}
  void push_back(const Scaffold & s) {m_scaffolds.push_back(s);} 

  void Sort() {
    ::Sort(m_scaffolds);
  }

 private:
  svec<Scaffold> m_scaffolds;
};






#endif //CONTSCAFF_H_

