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
  }

  SearchNode(int i, int b = -1, int ori = 1, int lap = 0, int pos = 0, int counter = 0) {
    m_read = i;
    m_back = b;
    m_ex = false;
    m_ori = ori;
    m_lap = lap;
    m_pos = pos;
    m_counter = counter;
  }
  
  void SetRead(int i) {m_read = i;}
  void SetOri(int o) {m_ori = o;}
  void SetPos(int i) {m_pos = i;}

  int Overlap() const {return m_lap;}
  int Offset() const {return m_lap;}
  int Ori()  const {return m_ori;}
  int Read() const {return m_read;}
  int Back() const {return m_back;}
  int Pos() const {return m_pos;}
  bool Ext() const {return m_ex;}
  void SetExt() {m_ex = true;}
  int Counter() const {return m_counter;}

  void IncCounter() {m_counter++;}

 protected:
  int m_read;
  int m_back;
  bool m_ex;
  int m_ori;
  int m_lap;
  int m_pos;
  int m_counter;
};

class SearchStack
{
public:
  SearchStack() {
    m_size = 0;
    m_len = 0;
    m_pairs = 0;
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
  bool Empty() const {return (m_size == 0);}


  int Size() const {return m_size;}

  void Minimal(SearchStack & s) {
    if (m_size == 0)
      return;
    int i = m_size - 1;
    while (i >= 0) {
      s.Push(m_stack[i]);
      i = m_stack[i].Back();
    }
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


class UsageItem
{
 public:
  UsageItem() {
    m_read = -1;
    m_pos = 0;
  }
  UsageItem(int read, int pos) {
    m_read = read;
    m_pos = pos/10;
  }

  bool operator < (const UsageItem & u) const {
    if (m_read != u.m_read)
      return m_read < u.m_read;
    return m_pos < u.m_pos;
  }

  int Read() const {return m_read;}
  int Pos()  const {return m_pos;}

 private:
  int m_read;
  int m_pos;
};

//=====================================================================
class UsageTracker
{
 public:
  UsageTracker() {
    m_div = 10;
    m_max = 100000; // Limit space 
  }
  UsageTracker(int n) {
    m_div = 10;
    Resize(n);
  }
  
  void Resize(int n) {
    m_cache.resize(n, -1);
  }
  void Clear() {
    int n = m_cache.isize();
    m_cache.clear();
    m_cache.resize(n, -1);
    m_full.clear();
  }

  bool IsUsed(int read, int pos) {
    pos /= m_div;
    int p = m_cache[read];
    if (p < 0)
      return false;
    // Limit the buffer
    if (p == pos || m_full.isize() > m_max)
      return true;
    if (!m_bSorted) {
      m_bSorted = true;
      m_full.Sort();
    }
    //cout << "Need bin search..." << endl;
    //int index = BinSearch(m_full, UsageItem(read, pos));
    int index = m_full.BinSearch(UsageItem(read, pos));
    if (index < 0) {
      return false;
    } else {
      return true;
    }
  }

  void SetUsed(int read, int pos) {
    if (IsUsed(read, pos))
      return;
    pos /= m_div;
    //if (m_full.isize() > m_max) {
    //  pos = 0;
    //}
    m_full.push_back(UsageItem(read, pos));
    m_cache.Set(read, pos);
    m_bSorted = false;
    //cout << "USED Size: " << m_full.isize() << endl;
  }

 private:
  bool m_bSorted;
  svec_buff<UsageItem> m_full;
  VecInt m_cache;
  int m_div;
  int m_max;
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
  Hypothesis() {}

  void Add(HypothesisNode & n) {
    m_main.push_back(n);
  }

  void Sort() {
    //::Sort(m_main);
    m_main.Sort();
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
    for (i=0; i<tmp.isize(); i++) {
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
    for (i=0; i<tmp.isize(); i++) {
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

  int isize() const {return m_main.isize();}
  const HypothesisNode & operator[] (int i) const {return m_main[i];}
  HypothesisNode & operator[] (int i) {return m_main[i];}
  void clear() {m_main.clear();}

private:
  //svec<HypothesisNode> m_main;
  svec_buff<HypothesisNode> m_main;
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
    m_maxResults = 50;
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

  bool DoSearchAll(const ConsensOverlapUnit & COUnit, int startWithRead = 0);

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
  int CountUnPairs(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint = false);
  void SetPairs(Hypothesis & hyp, const ConsensOverlapUnit & COUnit);

  bool IsUsedGlobal(const SearchNode & n) const {    
    //cout << m_globalUsed.isize() << endl;
    return m_globalUsed[n.Read()];
  }

  bool HasExtensions(const ConsensOverlapUnit & COUnit, int id) const;

  bool IsUsed(const SearchNode & n) const {
    //return false;
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
  int Evaluate(SearchStack & stack, const ConsensOverlapUnit & COUnit);
  
  bool IsNew(const SearchStack & test, const ConsensOverlapUnit & COUnit);
  
private:
  //void WeedOut();

  svec_buff<SearchStack> m_results;
  int m_maxResults;
  VecInt m_used;
  //VecInt m_usedFW;
  //VecInt m_usedRC;
  VecIntInc m_globalUsed;

  VecInt m_present;
  VecInt m_localUsed;
  svec_buff<int> m_ids;
  VecInt m_cov_seq;
  VecInt m_cov_pair;
  //VecInt m_cov_seq_strict;
  //VecInt m_cov_pair_strict;

  //svec<int> m_used;
  //svec<int> m_globalUsed;
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
