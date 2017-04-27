#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>
#include <ctime>
#include "ryggrad/src/base/FileParser.h"
#include "SearchOverlaps.h"
#include "GlobUsage.h"

bool Hypothesis::ContainsSubset(const Hypothesis & h)
{
    int i, j;
    int l1 = 0;
    int l2 = 0;

    if (h.Size() == 0)
        return true;

    // Stupid heuristics!!!!!
    int max = 3;
    int k = 0;
    for (int j=h.Size()-1; j>=0; j--) {
        const HypothesisNode & n = h[j];
        for (i=0; i<Size(); i++) {
            if (n.CloseEnough(m_main[i]))
                return true;
        }
        k++;
        if (k >= max)
            break;
    }
    return false;
}

bool Hypothesis::IsNew(const Hypothesis & h, const ConsensOverlapUnit & COUnit) 
{

  svec<int> query;
  h.GetNodes(query, COUnit);
  int n = query.isize();
  svec<int> ref;
  GetNodes(ref, COUnit);
  int m = ref.isize();
  for (int i=0; i<query.isize(); i++)
      ref.push_back(query[i]);
  UniqueSort(ref);
  double ratio = (double)ref.isize()/(double)m;
  if (ratio < 1.1) { //TODO parameterise
    return false;
  }
  return true;
}

bool Search::IsNew(const SearchStack & test, const ConsensOverlapUnit & COUnit)
{
  int i, j;
  svec<int> query;

  test.GetNodes(query, COUnit);
  int n = query.isize();
  for (i=m_results.Size()-1; i>=0; i--) { 
    svec<int> ref;
    m_results.GetResult(i).GetNodes(ref, COUnit);
    int m = ref.isize();
    for (j=0; j<query.isize(); j++)
      ref.push_back(query[j]);

    UniqueSort(ref);
    double ratio = (double)ref.isize()/(double)m;
    if (ratio < 1.02) {
      return false;
    }
  }
  return true;
}

int Search::Evaluate(SearchStack & stack, int diffNodeCount, const ConsensOverlapUnit & COUnit)
{
    // TRY this to reduce processing tiny contigs or where enough difference doesn't exist from previously evaluated path
  //cout << "Call Evaluate with " << stack.Size() << " elements" << endl;

    if (m_results.Size() > 0 && (stack.Top().NodeCount() < 3 || diffNodeCount <0.02*stack.Top().NodeCount()) )
        return -1;

    SearchStack minimal;
    stack.Minimal(minimal);
    minimal.Length(COUnit);

    m_workHyp.clear();
    MakeHypothesis(m_workHyp, minimal, COUnit, (m_exhaust && !m_results.Override()));
    SetPairs(m_workHyp, COUnit);

    int to, from;
    int pairs = CountPairs(to, from, m_workHyp, COUnit);
    minimal.SetPairs(pairs);
    m_results.AddResult(minimal, m_exhaust);
    return to;
}

int Search::SelectLeftest(const ConsensOverlapUnit & COUnit, bool& seedReverse)
{
    if (m_results.Size() == 0) {
        return -1;
    }
    int hypLen = m_results.GetTopResult().Length();

    m_workHyp.clear();
    MakeHypothesis(m_workHyp, m_results.GetTopResult(), COUnit/*, true*/);
    m_workHyp.Reverse(hypLen);
    if(m_workHyp.Size()>0) {
        seedReverse = (m_workHyp[0].Ori() != 1);
        return m_workHyp[0].Read();
    }
    return -1;
}


void Search::SelectTopN(const ConsensOverlapUnit & COUnit)
{
    if (m_results.Size() == 0) { return; }

    if (m_exhaust) {
        svec<Hypothesis> raw;
        raw.resize(m_results.Size());
        for (int i=0; i<m_results.Size(); i++) {
            int hypLen = m_results.GetResult(i).Length();
            MakeHypothesis(raw[i], m_results.GetResult(i), COUnit, !m_results.Override());
        }
        for (int i=0; i<raw.isize(); i++) {
            if (raw[i].Size() == 0)
                continue;
            for (int j=i+1; j<raw.isize(); j++) {
                if (raw[j].Size() == 0)
                    continue;
                if(raw[i].IsNew(raw[j], COUnit)) { 
                    raw[j].clear();
                }
                if (raw[i].ContainsSubset(raw[j])) {
                    raw[j].clear();
                }
                if (raw[j].ContainsSubset(raw[i])) {
                    raw[i] = raw[j];
                    raw[j].clear();
                }
            }
        }

        bool bCont = false;
        for (int i=0; i<raw.isize(); i++) {
            if (raw[i].Size() == 0)
                continue;
            Commit(raw[i]);
            int to, from;
            SetPairs(raw[i], COUnit);
            int numOfPairs = CountPairs(to, from, raw[i], COUnit, true);
            raw[i].SetPairs(numOfPairs);
            raw[i].TrimRight(to);
            raw[i].TrimLeft(from);
            SetPairs(raw[i], COUnit);
            //raw[i].RemoveUnpaired();
      
            m_sink.Dump(raw[i], COUnit, bCont, m_minAltKeep);
            bCont = true;
        }
    } else { //None-exhaustive mode
        m_workHyp.clear();
        MakeHypothesis(m_workHyp, m_results.GetTopResult(), COUnit, !m_results.Override());
        int hypLen = m_results.GetTopResult().Length();

        Commit(m_workHyp);
        int to, from;
        SetPairs(m_workHyp, COUnit);
        int numOfPairs = CountPairs(to, from, m_workHyp, COUnit, true);
        m_workHyp.SetPairs(numOfPairs);
        m_workHyp.TrimRight(to);
        m_workHyp.TrimLeft(from);
        SetPairs(m_workHyp, COUnit);
        //m_workHyp.RemoveUnpaired();
        m_sink.Dump(m_workHyp, COUnit);
    }
}

void Search::SetPairs(Hypothesis & hyp, const ConsensOverlapUnit & COUnit)
{
    m_present.clear();
    m_present.resize(COUnit.GetNumReads(), -1);

    for (int i=0; i<hyp.Size(); i++) {
        HypothesisNode & h = hyp[i];
        int r = h.Read();
        int numPartner = COUnit.getNumOfPartners(r);
        if (numPartner > 0) {
            for (int x = 0; x<numPartner; x++) {
              int partner = COUnit.getPartner(r, x);
              if (x == 0 || m_present[partner] < 0)
                m_present.Set(partner, i);
            }
        } else {
            h.SetPaired(-2);      
        }
    }
    for (int i=0; i<hyp.Size(); i++) {
        HypothesisNode & h = hyp[i];
        int r = h.Read();
        if (h.Pair() != -2)
            h.SetPaired(m_present[r]);
        else
            h.SetPaired(-1);
    }
}

int Search::CountPairs(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint)
{
    if (m_pairDir == 0)
        return CountUnPairs(to, from, hyp, COUnit, bPrint);

    to          =  -1;
    from        =  0;

    // Stupid heuristics!!
    if (hyp.Size() < 2)
        return 0;

    to          =  hyp[hyp.Size()-1].Stop();
    from        =  hyp[0].Start();  

    int pairCnt = 0;
    for (int i=0; i<hyp.Size(); i++) {
        const HypothesisNode & currNode = hyp[i];
        if(currNode.Start()+m_discount>to-m_discount && pairCnt!=0) {
            return pairCnt; //Stop the extension as pair coverage is broken
        }
        if(currNode.Pair()<0) { continue; } //Node is not paired - do no count
        const HypothesisNode & currPair = hyp[currNode.Pair()];
        if (currPair.Start() < currNode.Start() ||
            currPair.Start() - currNode.Start() > m_libSize) { continue; } //Pair exceeds library size limit - do not count
        if (currPair.Ori() != m_pairDir*currNode.Ori())      { continue; }

        if(from > currNode.Start() || pairCnt==0) { from = currNode.Start(); } //Extend bracket start 
        if(to   < currPair.Stop()  || pairCnt==0) { to   = currPair.Stop();  } //Extend bracket stop
        pairCnt++; 
    }

  //TODO old code consistency!
  if(from>0) { from = from - 1; }

  return pairCnt;
}

int Search::CountPairs_fullStat(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint)
{
    if (m_pairDir == 0)
        return CountUnPairs(to, from, hyp, COUnit, bPrint);

    bPrint = false;

    if (bPrint) 
        cout << "Start CountPairs." << endl;

    to = -1;
    from = 0;

    // Stupid heuristics!!
    if (hyp.Size() < 2)
        return 0;

    //  bPrint = true;
    int i, j;
    // INEFFICIENT!!!!!!
    //svec<int> present;
    //svec<int> cov_seq;
    //svec<int> cov_pair;

    int max = 0;
    int havePartner = 0;
    for (i=0; i<hyp.Size(); i++) {
        const HypothesisNode & h = hyp[i];
        if (h.Stop() > max)
            max = h.Stop();
    
        if (COUnit.getNumOfPartners(h.Read()) > 0)
            havePartner++;
    }
    // If they are unpaired, there is no point in going further.
    //if (havePartner < 2)
    //return 0;

    m_cov_seq.clear();
    m_cov_pair.clear();
    m_cov_seq.resize(max, 0);
    m_cov_pair.resize(max, 0);
 
    //m_cov_seq_strict.clear();
    //m_cov_pair_strict.clear();
    //m_cov_seq_strict.resize(max, 0);
    // m_cov_pair_strict.resize(max, 0);

    // HARD CODED!!!
    //int estPairSize = 350;

 
    m_present.clear();
    m_present.resize(COUnit.GetNumReads(), -1);

    int pairs = 0;
    for (i=0; i<hyp.Size(); i++) {
        const HypothesisNode & h = hyp[i];
        int r = h.Read();
        int numPartner = COUnit.getNumOfPartners(r);
        for (int x=0; x<numPartner; x++) {
            int partner = COUnit.getPartner(r, x);
            if (x == 0 || m_present[partner] < 0)
                m_present.Set(partner, i);
        }
    }
    for (i=0; i<hyp.Size(); i++) {
        const HypothesisNode & h = hyp[i];

        // Sequence coverage
        for (j=h.Start(); j<h.Stop(); j++)
            m_cov_seq.Set(j, m_cov_seq[j]+1);

        //for (j=h.Start()+m_discount/2; j<h.Stop()-m_discount/2; j++)
        //m_cov_seq_strict.Set(j, m_cov_seq_strict[j]+1);
        //if (bPrint && h.Read() == 771)
        //cout << "READ " << h.Read() << endl;

        // NOTE: This can actually go...
        if (m_pairDir == -1 && h.Ori() != 1)
            continue;

        int r = h.Read();
        if (m_present[r] > -1) {
            //if (bPrint && h.Read() == 771) {
            //cout << "Has partner " << endl;
            //}
            const HypothesisNode & g = hyp[m_present[r]];

            if (g.Start() < h.Start())
                continue;
 
            int start = h.Start();
            int stop = g.Stop();
            //if (stop < start) {
            //int tmp = start;
            //start = stop;
            //stop = tmp;
            //}
            // WARNING: Hard-coded limit for library size!!!!
 
            // Enforce proper location of pairs, i.e. pointing towards each other

            if (g.Ori() != m_pairDir * h.Ori())
                continue;

            // HARD CODED!!!
            if (stop - start > 10000)
                continue;
            for (j=start + m_discount; j<stop-m_discount; j++)
                m_cov_pair.Set(j, m_cov_pair[j]+1);

       
        }
        //  pairs++;
    }

    to = m_cov_seq.isize();
    from = 0;
    double last = 0.;
    bool bYes = false;
    bool bYesGoodCov = false;
    if (bPrint) {
        cout << "Printing coverage per base:" << endl;
    }

    //int minCovSingle;

    for (i=m_discount; i<m_cov_seq.isize()-m_discount; i++) {
        double ratio = (double)m_cov_pair[i]/(double)(m_cov_seq[i]+1);
        if (ratio > 0.2)
            bYes = true;
        if (ratio > 0.0001) {
            if (from == 0)
                from = i - m_discount - 1;
        }

        if (bPrint) {
            cout << i << " " <<  m_cov_pair[i] << " " << m_cov_seq[i] << " " << ratio << " " << last << endl;
        }


        if (ratio <= 0.000001 && bYes) {
            to = i + m_discount;
            if (bPrint) {
                cout << "Found break point at pos " << to << endl;
            }
            break;
        }
        last = ratio;
    }

    if (from < 0)
        from = 0;
 
    // NO BREAKING!!!
    //if (bPrint) {
    //  to = m_cov_seq.isize();
    //  from = 0;
    //}

    //cout << "from: " << from << " to: " << to << endl;  

    for (i=0; i<hyp.Size(); i++) {
        const HypothesisNode & h = hyp[i];
        if (h.Start() < from)
            continue;
        if (to > 0 && h.Stop() > to)
            continue;

        // NOTE: This can go
        if (m_pairDir == -1 && h.Ori() != 1)
            continue;
        int r = h.Read();
        if (m_present[r] > -1) {
            const HypothesisNode & g = hyp[m_present[r]];
            bool bGood = true;
            if (g.Start() < h.Start())
                bGood = false;

            // WARNING: Hard-coded limit for library size!!!!
            if (g.Start() - h.Start() > 10000)
                bGood = false;

            // Note: this DOES NOT WORK for the small mouse example.
            // Enforce proper location of pairs, i.e. pointing towards each other
            if (g.Ori() != m_pairDir * h.Ori())
                bGood = false;
            if (bGood)
                pairs++;
        }
    }

    if (bPrint) {
        cout << "Returning coordinates: from=" << from << " to=" << to << endl;
    }

    // Not sure what to return here...
    //return pairs;
    //if (to < 0)
    //return 0;
    return pairs;
}



int Search::CountUnPairs(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint)
{
    bPrint = false;

    //return hyp.Size();
    
    if (bPrint) 
        cout << "Start CountUnPairs." << endl;

    to = -1;
    from = 0;

    // Stupid heuristics!!
    if (hyp.Size() < 2)
        return 0;

    int i, j;
 
    int max = 0;
    for (i=0; i<hyp.Size(); i++) {
        const HypothesisNode & h = hyp[i];
        if (h.Stop() > max)
            max = h.Stop();    
    }

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return max;
    
 
    m_cov_seq.clear();
    m_cov_seq.resize(max+1, 0);
 
    int total = 0;
    for (i=0; i<hyp.Size(); i++) {
        const HypothesisNode & h = hyp[i];

        // Sequence coverage
        for (j=h.Start()+m_discount/2; j<h.Stop()-m_discount/2; j++) {
            m_cov_seq.Set(j, m_cov_seq[j]+1);
            total++;
        }
    }
 
    int avg = total / (m_cov_seq.isize() + 1);
    int minSingle = 1 + avg / 20; // HARD CODED!!!!

    to = m_cov_seq.isize();
    from = 0;

    if (bPrint) {
        cout << "Printing coverage per base:" << endl;
    }
  
    //bool bAbove - false;
    if (bPrint) {
        cout << "Coverage: " << avg << " thresh: " << minSingle << endl;
    }

    if (minSingle < 1)
        minSingle = 1;

    for (i=m_discount; i<m_cov_seq.isize()-m_discount; i++) {
        if (bPrint) {
            cout << i << " " << m_cov_seq[i] << endl;
        }
        if (m_cov_seq[i-1] > minSingle && m_cov_seq[i] <= minSingle) {
            to = i;
            if (bPrint) {
                cout << "Found SINGLE break point at pos " << to << endl;
            }
            break;
      
        }
    }

    if (from < 0)
        from = 0;
 
    if (bPrint) {
        cout << "Returning coordinates: from=" << from << " to=" << to << endl;
    }

    return to;
}


// Needs work...
void Search::MakeHypothesis(Hypothesis & hyp, const SearchStack & ss, const ConsensOverlapUnit & COUnit, bool bAddIn)
{
    const svec<SearchNode> & st = ss.Raw();
    const ConsensReads & consReads = COUnit.getConsReads();
  
    hyp.Reserve(2*st.isize());

    int i, j;
    int off = 0;

    //******************** Inefficient!!!! *********************  
    m_localUsed.clear();

    if (bAddIn) {
        m_localUsed.resize(consReads.getNumOfReads(), 0);
        m_ids.resize(ss.Size());
        for (i=0; i<ss.Size(); i++) {
            m_ids[i] = st[i].Read();
            m_localUsed.Set(m_ids[i], 1);
        }
        sort(m_ids.begin(), m_ids.end());
    }
    int lowest = 0;
    for (i=0; i<ss.Size(); i++) {
        HypothesisNode node;
        node = st[i];
    
        int dSize = consReads.getSize(node.Read());
        off += node.Offset();
        node.SetCoords(off, off + dSize);
        hyp.Add(node);

	if (node.Start() < lowest)
	  lowest = node.Start();
	
        if (bAddIn && i+1 < ss.Size()) {
            int offLimit  = st[i+1].Offset();
            int numOfLaps = COUnit.GetNumLaps(node.Read());
            for (j=0; j<numOfLaps; j++) {
                const ReadOverlap & s = COUnit.GetLap(node.Read(), j);
                if (m_globalUsed[s.getOverlapIndex()] > 0)
                    continue;
                if (m_localUsed[s.getOverlapIndex()] > 0)
                    continue;
	
                // Only add in partnered COUnit & COUnit that have their partner placed already.
                int numPartner = COUnit.getNumOfPartners(s.getOverlapIndex());
                if (numPartner == 0)
                    continue;
                int partner = COUnit.getPartner(s.getOverlapIndex(), 0);
                if (!binary_search(m_ids.begin(), m_ids.end(), partner)) {
                    continue;
                }

                if (s.getDirection() != node.Ori()) {
                    continue;
                }
                if (s.getContactPos() > offLimit) {
                    continue;
                } else {
                }

                if (COUnit.HasLap(st[i+1].Read(), s.getOverlapIndex())) {
                    HypothesisNode t = SearchNode(s.getOverlapIndex(), -1, s.getOrient() * node.Ori(), s.getContactPos());
                    dSize = consReads.getSize(s.getOverlapIndex());
                    t.SetCoords(off + t.Offset(), off + t.Offset() + dSize);
                    hyp.Add(t);	
                    m_localUsed.Set(s.getOverlapIndex(), 1);
                } else {
                }   
            }
        }
    }
    hyp.Sort();
    for (i=0; i<hyp.Size(); i++)
      hyp[i].AddOffset(-lowest);
}

void Search::Commit(const Hypothesis & hyp)
{
    int i;
    for (i=0; i<hyp.Size(); i++) {
        m_globalUsed.Set(hyp[i].Read(), m_globalUsed[hyp[i].Read()]+1);
        m_report++;
        //cout << "Global used: " << hyp[i].Read() << endl;
    }
    if (m_pGlob != NULL /*&& m_report > 100*/) {
        //cout << "Sync'ing up with global usage handler." << endl;
        m_pGlob->Sync(m_globalUsed);
        m_report = 0;
    }
    m_globalUsed.Reset();
}

bool Search::DoSearchAll(const ConsensOverlapUnit & COUnit, int numAvailableReads, int startWithRead)
{
    srand(time(NULL));
    if (m_exhaust) {
        m_usage.Clear();
        m_usage.Resize(COUnit.GetNumReads(), numAvailableReads);
    }

    if (m_globalUsed.isize() == 0)
        m_globalUsed.resize(COUnit.GetNumReads(), 0);

    m_results.Reset();
    m_results.SetCapacity(m_exhaust);

    int i;
 
    //cout << "Processing all remaining COUnit (fw/rc)." << endl;
    for (i=startWithRead; i<COUnit.GetNumReads(); i++) {
        if (!IsUsedGlobal(i)
            && HasExtensions(COUnit, i)) {
	    //cout << "Left w/ " << i << endl;
            bool seedReversed = false;
            int start = SelectStartNode(COUnit, i, seedReversed);
	    //cout << "Yielded " << start << endl;
            if(start>=0) { DoSearch(COUnit, start, seedReversed); }
        }
    }
    return true;
}

bool Search::HasExtensions(const ConsensOverlapUnit & COUnit, int id) const
{
    int ori = 1;
    int numLaps = COUnit.GetNumDirLaps(id, ori);
    for (int i=0; i<numLaps; i++) {
        const ReadOverlap & l = COUnit.GetDirLap(id, i, ori);	
        if (m_globalUsed[l.getOverlapIndex()] == 0)
            return true;
    }
    return false; 
}

void Search::SearchCore(const ConsensOverlapUnit & COUnit, int index, bool seedReverse)
{
  //cout << "call DoSearch w/ " << index << " left " <<  COUnit.GetNumDirLaps(index, -1) << " ";
  //cout << " right " <<  COUnit.GetNumDirLaps(index, 1) << endl;
   
    m_lastNoPairs = -1;
    if (m_exhaust)
        m_usage.Clear();

    if (m_globalUsed.isize() == 0)
        m_globalUsed.resize(COUnit.GetNumReads(), 0);

    m_used.clear();
    m_used.resize(COUnit.GetNumReads(), 0);

    m_results.Reset();

    SearchStack stack;
    int diffNodeCount = 0; // Used to keep track of differences from latest stack

    SearchNode init(index);
    if (seedReverse) {
        init.SetOri(-1);
    }
    init.SetNodeCount(1);
    stack.Push(init);
    SetUsed(init);
    if (m_exhaust) {
        m_usage.SetUsed(init.Read(), init.Pos());
    }
    int i;

    SearchNode n;
    int pos = 0;
    int backoff = -1;

    while (!stack.Empty()) {
        SearchNode & n = stack.Top();
        int curr = stack.Curr();
        int pos = n.Pos();
        int index = n.Counter();
        int ori = n.Ori();
        int nodeCount = n.NodeCount();

        if (index == COUnit.GetNumDirLaps(n.Read(), ori)) {
            if (!n.Ext()) {
                int limit = Evaluate(stack, diffNodeCount, COUnit);
                diffNodeCount = 0;
                //break; //MGG: This should make the search greedy!!!

                if (limit < n.Pos()) {	  
                    // Do or do not not backoff all the way!!
                    backoff = limit;
                } else {
                    backoff = -1; 
                }
            }
            // Backtrack 
            SearchNode pop;
            do {
                stack.Pop(pop);
                // The current implementation is slightly too greedy
                if (m_exhaust) {
                    if (backoff > 0 && !(pop.Pos() <= backoff)) {
                    } else {
                        SetUsed(pop, false); // Let's allow for re-usage
                    }
                }
            } while(backoff >=0 && pop.Pos() > backoff);
            backoff = -1;
            continue;
        }
    
        int overlapCnt = COUnit.GetNumDirLaps(n.Read(), ori); 
        int limit      = COUnit.GetNumDirLaps(n.Read(), ori); //Limit the number of overlaps to consider
        if(m_exhaust) {
          limit = min(COUnit.getConsReadSize(n.Read()), overlapCnt);
        }
        for (i=index; i < overlapCnt; i++) {
            n.IncCounter();
            if(i>limit) {
                continue; 
            }
            const ReadOverlap & l = COUnit.GetDirLap(n.Read(), i, ori);
            SearchNode to_push(l.getOverlapIndex(), curr, ori*l.getOrient(), l.getContactPos(), pos + l.getContactPos());
            to_push.SetNodeCount(nodeCount+1); 

            if (!IsUsed(to_push) &&
                 (!m_exhaust || (m_exhaust && !m_usage.IsUsed(to_push.Read(), to_push.Pos())))) {
                n.SetExt();
                stack.Push(to_push);
                diffNodeCount++;
                SetUsed(to_push);
                if (m_exhaust)
                    m_usage.SetUsed(to_push.Read(), to_push.Pos());
                break;
            }
        }
    }
    m_results.Sort();
}

int Search::SelectStartNode(const ConsensOverlapUnit & COUnit, int index, bool& seedReverse) {
    SearchCore(COUnit, index, seedReverse);
    int z = SelectLeftest(COUnit, seedReverse);
    if (z < 0) {
        cout << "ERROR!!!" << endl;
    }
    return z;
}

void Search::DoSearch(const ConsensOverlapUnit & COUnit, int index, bool seedReverse) {
    SearchCore(COUnit, index, seedReverse);
    SelectTopN(COUnit);
}


