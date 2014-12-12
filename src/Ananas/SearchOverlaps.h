#ifndef SEARCHOVERLAPS_H
#define SEARCHOVERLAPS_H

#include <string>
#include "src/DNAVector.h"
#include "base/SVector.h"
#include "src/Ananas/Consensus.h"
#include "src/Ananas/CachedVec.h"
#include "src/Ananas/ConsensOverlapUnit.h"

class GlobalUsageHandler;

class SearchNode
{
 public:
  SearchNode() {
    m_read = -1;
    m_back = -1;
    m_ex = false;
    m_ori = 0;
    m_lap = 0;
    m_pos = 0;
    m_counter = 0;
    m_nodesSoFar = 0;
  }

  SearchNode(int i, int b = -1, int ori = 1, int lap = 0, int pos = 0, int counter = 0, int nodesSoFar=0) {
    m_read = i;
    m_back = b;
    m_ex = false;
    m_ori = ori;
    m_lap = lap;
    m_pos = pos;
    m_counter = counter;
    m_nodesSoFar = nodesSoFar;
  }
  
  void SetRead(int i) {m_read = i;}
  void SetOri(int o) {m_ori = o;}
  void SetPos(int i) {m_pos = i;}
  void SetNodeCount(int n) {m_nodesSoFar = n; }

  int Overlap() const {return m_lap;}
  int Offset() const {return m_lap;}
  int Ori()  const {return m_ori;}
  int Read() const {return m_read;}
  int Back() const {return m_back;}
  int Pos() const {return m_pos;}
  bool Ext() const {return m_ex;}
  void SetExt() {m_ex = true;}
  int Counter() const {return m_counter;}
  int NodeCount() const {return m_nodesSoFar;}

  void IncCounter() {m_counter++;}

 protected:
  int m_read;
  int m_back;
  bool m_ex;
  int m_ori;
  int m_lap;
  int m_pos;
  int m_counter;
  int m_nodesSoFar;
};

class SearchStack
{
public:
  SearchStack() {
    m_size = 0;
    m_len = 0;
    m_pairs = 0;
  }
  
  void Resize(int size) {
    m_stack.resize(size);
  }

  void Reset() {
    m_size = 0;
  }

  void Push(const SearchNode & n) {
    if (m_stack.isize() <= m_size)
      m_stack.resize(m_size + 4096);
    m_stack[m_size] = n;
    m_size++;
  }

  int Curr() const {return m_size-1;}

  bool Pop(SearchNode & n) {
    if (m_size == 0)
      return false;
    m_size--;
    n = m_stack[m_size];    
    return true;
  }
  bool Pop() {
    if (m_size == 0)
      return false;
    m_size--;
    return true;
  }

  SearchNode & Top() {return m_stack[m_size-1];}
  bool Empty() const {return (m_size == 0);    }


  int Size() const {return m_size;}

  void Minimal(SearchStack & s) {
    if (m_size == 0)
      return;
    int totSize = Top().NodeCount();
    s.Resize(totSize);
    int next  = m_size - 1;
    int count = totSize - 1;
    while (count>=0) {
      s.m_stack[count] = m_stack[next];
      next = m_stack[next].Back();
      count--;
    }
    s.m_size = totSize;
  }

  void SetPairs(int i) {m_pairs = i;}
  int Pairs() const {return m_pairs;}
 
  
  void Trim(const ConsensOverlapUnit & COUnit, int max) {
    int len = 0;
    int i;
    for (i=m_size-1; i>=0; i--) {
      len += m_stack[i].Overlap();
      if (len > max) {
	i++;
	break;
      }
    }
    int j = 0;
    if (i > 0) {
      for (; i<m_size; i++) {
	m_stack[j] = m_stack[i];
	j++;
      }
      m_size = j;
      //m_stack.resize(j);
    }
    //cout << "Old length: " << m_len << endl;
    int l = Length(COUnit);
    //cout << "New length: " << m_len << endl;
  }

  int Length(const ConsensOverlapUnit & COUnit) {
    int len = 0;
    const ConsensReads & consReads = COUnit.getConsReads();
    //cout << "Items in stack: " << m_size << endl;
    for (int i=m_size-1; i>=0; i--) {
      len += m_stack[i].Overlap();
    }
    
    if (m_size > 0)
      len += consReads.getSize(m_stack[0].Read());

    m_len = len;
    return len;
  }

  int GetNodes(svec<int> & nodes, const ConsensOverlapUnit & COUnit) const {
    int len = 0;
    for (int i=m_size-1; i>=0; i--) {
      int r = m_stack[i].Read();
      nodes.push_back(r);
      len += m_stack[i].Overlap();
    }
    //m_len = len;
    return len;
  }

  int Print(const ConsensOverlapUnit & COUnit) const {
    int len = 0;
    DNAVector consensus;
    int i;
    int prevLen = 0;
    for (i=m_size-1; i>=0; i--) {
      const ConsensReads & consReads = COUnit.getConsReads();
      int r = m_stack[i].Read();
      len += m_stack[i].Overlap();
      DNAVector d = consReads[r];
      if (m_stack[i].Ori() == -1)
	d.ReverseComplement();
 
      DNAVector sub;
      int snippet = prevLen - m_stack[i].Overlap();
      sub.SetToSubOf(d, snippet, d.isize()-snippet);
      consensus += sub;

      prevLen = d.isize();

    }
    return len;
  }
  const svec<SearchNode> & Raw() const {return m_stack;}

  bool operator < (const SearchStack & s) const {
    if (m_pairs != s.m_pairs)
      return m_pairs < s.m_pairs;
    return m_len < s.m_len;
  }

  int Length() const {return m_len;}


  void Reverse() {
    svec<SearchNode> tmp;
    tmp.resize(m_size);
    int k = 0;
    for (int i=m_size-1; i>=0; i--) {
      tmp[k] = m_stack[i];
      k++;
    }
    m_stack = tmp;
  }


private:
  int m_size;
  svec<SearchNode> m_stack;
  int m_len;
  int m_pairs;
  int m_maxKeep;
};

//=====================================================================
class UsageTracker
{
 public:
  UsageTracker() {
    m_div         = 10;
    m_maxPerRead  = 10;
    m_usedCounter = 0;
  }

  void Resize(int totReads, int totAvailReads) {
    m_full.resize(totReads);
    m_usedList.resize(totAvailReads);
  }

  void Clear() {
    for(int i=0; i<m_usedCounter ; i++) {
      m_full[m_usedList[i]].clear();
      m_usedCounter = 0;
    }
  }

  bool IsUsed(int read, int pos) {
    unsigned int groupPos = pos/m_div;
    if(m_full[read].size()==m_maxPerRead){ //If maximum number of positions has been registered mark as used
      return true;
    } else {
      if(m_full[read].find(groupPos)==m_full[read].end()) {
        return false;
      } else {
        return true;
      }
    }
  }

  void SetUsed(int read, int pos) {
    unsigned int posCnt = m_full[read].size();
    if(posCnt<m_maxPerRead) {
      if(posCnt==0) { //First time a read is used (update used counter)
         m_usedList[m_usedCounter] = read;
         m_usedCounter++; 
      } 
      unsigned int groupPos = pos/m_div;
      m_full[read][groupPos] = true;
    }
  }

 private:
  vector< map<int, bool> > m_full; ///List of positions per used read (indexed over all reads) 
  vector<int> m_usedList;          ///List of read indexes that have been used and need to be cleared in the next round
  int m_usedCounter;               ///Counts the number of reads so far used (at least once) 
  int m_div;                       ///Divide the position values into blocks 
  unsigned int m_maxPerRead;       ///Maximum number of positions registered per read
};
//=====================================================================
class HypothesisNode : public SearchNode
{
public:
  HypothesisNode() {
    m_start = m_stop = -1;
    m_paired = -1;
  }
  HypothesisNode(const SearchNode & s) {
    m_start = m_stop = -1;
    m_paired = -1;
    *this = s;
  }

  void Print() const {
    //cout << "Read " << m_read << " ori " << m_ori << " " << m_start << " - " << m_stop << endl;
  }

  HypothesisNode & operator = (const SearchNode & s) {
    m_start = m_stop = -1;
    m_lap = s.Overlap();
    m_ori = s.Ori();
    m_read = s.Read();
    m_back = s.Back();
    m_ex = s.Ext();     
    m_pos = s.Pos();     
    m_counter = s.Counter();
    return *this;
  }
  
  int Start() const {return m_start;}
  int Stop() const {return m_stop;}
  
  int Pair() const {return m_paired;}
  void SetPaired(int i) {m_paired = i;}

  void SetCoords(int start, int stop, bool reverse = false) {
    m_start = start;
    m_stop = stop;
    if (reverse)
      m_ori *= -1;
  }

  bool operator < (const HypothesisNode & n) const {
    if (m_start != n.m_start)
      return m_start < n.m_start;
    return m_read < n.m_read;
  }
  bool operator == (const HypothesisNode & n) const {
    if (m_read == n.m_read && m_start == n.m_start && m_ori == n.m_ori)
      return true;
    return false;
  }
  bool CloseEnough(const HypothesisNode & n) const {
    if (m_read == n.Read() && m_ori == n.Ori()) {
      int x = m_start - n.Start();
      if (x <= 3 && x >= -3)
	return true;
    }
    return false;
  }


private:
  int m_start;
  int m_stop;
  int m_paired;
};


class Hypothesis
{
public:
  Hypothesis() {
    m_pairs = 0;
  }

  bool operator < (const Hypothesis & h) const {
    if (m_pairs != h.m_pairs)
      return m_pairs < h.m_pairs;
    return Length() < h.Length();
  }

  void SetPairs(int i) {m_pairs = i;}
  int Pairs() const {return m_pairs;}
  int Length() const { 
    if(Size()==0) { return 0; } 
    else          { return m_main[Size()-1].Stop(); } 
  } 
  
  void Reserve(int size) { 
    m_main.reserve(size);
  } 

  void Add(HypothesisNode & n) {
    m_main.push_back(n);
  }

  void Sort() {
    sort(m_main.begin(), m_main.end());
  }

  void PrettyPrint(const ConsensOverlapUnit & COUnit) const {
    int i;

    //=================================================
    //return;
    //=================================================

    const ConsensReads & consReads = COUnit.getConsReads();
    for (i=0; i<m_main.isize(); i++) {
      int r = m_main[i].Read();
      cout << r << "\t" << m_main[i].Ori() << "\t" << m_main[i].Start() << " - " << m_main[i].Stop() << "\t";
      cout << consReads[r].Name();
      cout << endl;
    }
    for (i=0; i<m_main.isize(); i++) {
      int r = m_main[i].Read();
      DNAVector d = consReads[r];
      if (m_main[i].Ori() == -1)
	d.ReverseComplement();
      for (int y=0; y<m_main[i].Start(); y++)
	cout << " ";
      for (int x=0; x<d.isize(); x++)
	cout << d[x];
      cout << endl;      
    }    
  }

  int GetNodes(svec<int> & nodes, const ConsensOverlapUnit & COUnit) const {
    int len = 0;
    int numOfNodes = Size();
    nodes.resize(numOfNodes);
    for (int i=numOfNodes-1; i>=0; i--) {
      int r = m_main[i].Read();
      nodes[i] = r;
      len += m_main[i].Overlap();
    }
    return len;
  }
  bool ContainsSubset(const Hypothesis & h);

  void Reverse(int totalLen) {
    int i;
    for (i=0; i<m_main.isize(); i++) {
      int start = m_main[i].Start();
      int stop = m_main[i].Stop();
      int startR = totalLen - stop; 
      int stopR = totalLen - start;
      m_main[i].SetCoords(startR, stopR, true);
    }

    Sort();
  }

  void RemoveUnpaired() {
    int n = 0;
    int i;
    for (i=0; i<m_main.isize(); i++) {
      if (m_main[i].Pair() != -1)
	n++;
    }
    if (n == 0)
      return;

    Hypothesis tmp = *this;
    m_main.clear();
    for (i=0; i<tmp.Size(); i++) {
      if (tmp[i].Pair() != -1)
	Add(tmp[i]);
    }
  }


  void TrimLeft(int from) {
    if (from == 0)
      return;
    int i;
    //int off = 0;
    int lowest = -1;
    for (i=0; i<m_main.isize(); i++) {
      int start = m_main[i].Start();
      if (start >= from) {
	if (lowest == -1) {
	  lowest = start;
	} else {
	  if (start < lowest)
	    lowest = start;
	}
      } 
    }
    //cout << "Left trimming hypothesis, offset: " << lowest << endl;
    Hypothesis tmp = *this;
    m_main.clear();
    //cout << "Before: " << tmp.isize() << endl;
    for (i=0; i<tmp.Size(); i++) {
      int start = tmp[i].Start();
      //cout << "start=" << start << " from=" << from << endl;
      if (start >= from) {
	//cout << "Add." << endl;
	tmp[i].SetCoords(tmp[i].Start()-lowest, tmp[i].Stop()-lowest);
	Add(tmp[i]);
      }
    }
    //cout << "After: " << m_main.isize() << endl;
  }
  void TrimRight(int to) {
    if (to < 0)
      return;
    int i;
    //int off = 0;
    int highest = -1;
    for (i=0; i<m_main.isize(); i++) {
      int start = m_main[i].Start();
      if (start < to) {
	highest = i;
      } else {
	break;
      }
    }
    m_main.resize(highest+1);
 
    //cout << "After: " << m_main.isize() << endl;
  }

  int Size() const {return m_main.isize();}
  const HypothesisNode & operator[] (int i) const {return m_main[i];}
  HypothesisNode & operator[] (int i) {return m_main[i];}
  void clear() {m_main.clear();}

private:
  svec<HypothesisNode> m_main;
  int m_pairs;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Search
{
 
public:
  Search() {
    m_exhaust = false;
    m_discount = 30;
    m_pGlob = NULL;
    m_report = 0;
    m_lastNoPairs = -1;
    m_pairDir = -1;
    m_override = false;
    m_maxResults = 500; //TODO depends on the strategy taken in structure chosen for m_results
    m_minAltKeep = 200;
  }

  void SetGlobalUsage(GlobalUsageHandler * p) {
    m_pGlob = p;
  }

  void SetMinAltKeep(int i) {
    m_minAltKeep = i;
  }

  void SetDir(const string & d) {
    if (d == "fr") {
      m_pairDir = -1;
      return;
    }
    if (d == "ff") {      
      m_pairDir = 1;
      return;
    }
    if (d == "na") {      
      m_pairDir = 0;
      return;
    }
    throw;
  }

  bool DoSearchAll(const ConsensOverlapUnit & COUnit, int numAvailableReads, int startWithRead = 0);

  void SetOutput(const string & layout) {
    m_sink.SetLayoutFile(layout);
  }
  
  void SetExhaustive(bool b) {
    m_exhaust = b;
  }

  void SetUsedAll(const ConsensOverlapUnit & COUnit) {
    m_globalUsed.clear();
    m_globalUsed.resize(COUnit.GetNumReads(), 1);
  }
  void SetUsed(int i, bool b) {
    if (b == false)
      m_globalUsed.Set(i, 0);
    else
      m_globalUsed.Set(i, 1);
       
  }

  void SetIndex(int i) {
    m_sink.SetPrefix(i);
  }

protected:
  int DoSearch(const ConsensOverlapUnit & COUnit, int index, bool rc = false);

  int CountPairs(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint = false);
  int CountPairs_fullStat(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint = false);
  int CountUnPairs(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint = false);
  void SetPairs(Hypothesis & hyp, const ConsensOverlapUnit & COUnit);

  bool IsUsedGlobal(const SearchNode & n) const {    
    return m_globalUsed[n.Read()];
  }

  bool HasExtensions(const ConsensOverlapUnit & COUnit, int id) const;

  bool IsUsed(const SearchNode & n) const {
    if (m_globalUsed[n.Read()] > 0)
      return true;
    return (m_used[n.Read()] > 0);
  }
  void SetUsed(const SearchNode & n, bool bSet = true)  {
    if (bSet)
      m_used.Set(n.Read(), m_used[n.Read()]+1);
    else
      m_used.Set(n.Read(), 0);      
  }

  void SetGlobalUsed(const SearchNode & n, bool bSet = true)  {
   if (bSet)
     m_globalUsed.Set(n.Read(), m_globalUsed[n.Read()]+1);
    else
      m_globalUsed.Set(n.Read(), 0);      
  }
  
  void MakeHypothesis(Hypothesis & hyp, const SearchStack & r, const ConsensOverlapUnit & COUnit, bool bAddIn = false);

  void Commit(const Hypothesis & hyp);

  void SelectTopN(const ConsensOverlapUnit & COUnit, bool rc);
  int SelectLeftest(const ConsensOverlapUnit & COUnit,  bool rc);
  int Evaluate(SearchStack & stack, int diffNodeCount, const ConsensOverlapUnit & COUnit);
  
  bool IsNew(const SearchStack & test, const ConsensOverlapUnit & COUnit);
  
private:
  //void WeedOut();

  svec<SearchStack> m_results;
  int m_maxResults;
  VecInt m_used;
  VecIntInc m_globalUsed;

  VecInt m_present;
  VecInt m_localUsed;
  svec<int> m_ids;
  VecInt m_cov_seq;
  VecInt m_cov_pair;

  LayoutSink m_sink;
  bool m_exhaust;
  int m_discount;
  Hypothesis m_workHyp;
  
  UsageTracker m_usage;
  GlobalUsageHandler * m_pGlob;
  int m_report;
  int m_lastNoPairs;
  int m_pairDir;

  bool m_override;
  int m_minAltKeep;

};


#endif //SEARCHOVERLAPS_H
