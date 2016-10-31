#ifndef SEQUENCE_GRAPH
#define SEQUENCE_GRAPH

#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/aligns/KmerDynProg.h"

class SeqCoords
{
 public:
  SeqCoords() {
    m_first = -1;
    m_last = -1;
    m_seq = -1;
  }
  SeqCoords(int seq, int first, int last) {
    Set(seq, first, last);
  }

  void Set(int seq, int first, int last) {
    m_first = first;
    m_last = last;
    m_seq = seq;
  }

  int First() const {return m_first;}
  int Last() const {return m_last;}
  int Seq() const {return m_seq;}

 private:
  int m_seq;
  int m_first;
  int m_last;
};

class SeqNode
{
 public:
  SeqNode() {}

  void AddContrib(const SeqCoords & s) {
    m_contrib.push_back(s);
  }
  void push_back(const SeqCoords & s) {
    m_contrib.push_back(s);
  }
  
  int GetNumContrib() const {return m_contrib.isize();}
  const SeqCoords & Contrib(int i) const {return m_contrib[i];}
  int NumOut() const {return m_out.isize();}
  int Out(int i) const {return m_out[i];}
  int NumIn() const {return m_in.isize();}
  int In(int i) const {return m_in[i];}

  void AddOut(int i) {m_out.push_back(i);}
  void AddIn(int i) {m_in.push_back(i);}

 private:
 
  svec<int> m_out;
  svec<int> m_in;
  svec<SeqCoords> m_contrib;
  
};


class SeqGraph
{
 public:
  SeqGraph() {}
  
  void AddCoords(int i, int seq = -1, int from = -1, int to = -1) {
    cout << "Adding " << seq << " " << from << " - " << to << endl;
    m_nodes[i].AddContrib(SeqCoords(seq, from, to));
  }

  int AddSeqNode() {
    m_nodes.push_back(SeqNode());
    return m_nodes.isize()-1;
  }

  void Connect(int a, int b) {
    if (a < 0 || b < 0)
      return;
    m_nodes[a].AddOut(b);
    m_nodes[b].AddIn(a);
  }
  

  SeqNode & operator[] (int i) {return m_nodes[i];}
  const SeqNode & operator[] (int i) const {return m_nodes[i];}
  int isize() const {return m_nodes.isize();}
  void clear() {m_nodes.clear();}
  void push_back(const SeqNode & s) {
    m_nodes.push_back(s);
  }
  void pop() {
    m_nodes.resize(m_nodes.isize()-1);
  }

 private:
  svec<SeqNode> m_nodes;

};


class GraphEnumerate
{
 public:
  GraphEnumerate() {}

  void Enumerate(svec<SeqGraph> & linear, const SeqGraph & in);

 private:
  bool Search(svec<SeqGraph> & linear, const SeqGraph & in, int index);
  
  SeqGraph m_work;
};

class SequenceBuild
{
 public:
  SequenceBuild() {}

  void Build(vecDNAVector & out, const vecDNAVector & in, const svec<SeqGraph> & linear);
  void Build(DNAVector & out, const vecDNAVector & in, const SeqGraph & linear);
};

class GraphConstructor
{
 public:
  GraphConstructor(const vecDNAVector * pDna);

  void Construct(SeqGraph & out);

  void Sequence(DNAVector & out, const SeqGraph & in);

 private:
  //KmerSuperAligner m_sa;
  const vecDNAVector * m_pDna;
};

#endif //SEQUENCE_GRAPH
