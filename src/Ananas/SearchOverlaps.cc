#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>
#include "base/FileParser.h"
#include "src/Ananas/SearchOverlaps.h"
#include "src/Ananas/GlobUsage.h"

bool Hypothesis::ContainsSubset(const Hypothesis & h)
{
    int i, j;
    int l1 = 0;
    int l2 = 0;

    if (h.isize() == 0)
        return true;

    // Stupid heuristics!!!!!
    int max = 3;
    int k = 0;
    for (int j=h.isize()-1; j>=0; j--) {
        const HypothesisNode & n = h[j];
        for (i=0; i<isize(); i++) {
            if (n.CloseEnough(m_main[i]))
                return true;
        }
        k++;
        if (k >= max)
            break;
    }
    return false;
}

bool Search::IsNew(const SearchStack & test, const ConsensOverlapUnit & COUnit)
{
  int i, j;
  svec<int> query;

  test.GetNodes(query, COUnit);
  int n = query.isize();
  for (i=m_results.isize()-1; i>=0; i--) { 
    svec<int> ref;
    m_results[i].GetNodes(ref, COUnit);
    int m = ref.isize();
    for (j=0; j<query.isize(); j++)
      ref.push_back(query[j]);

    UniqueSort(ref);
    double ratio = (double)ref.isize()/(double)m;
    //cout << "m=" << m << " merge=" << ref.isize() << " q=" << n << endl;
    if (ratio < 1.02) {
      //cout << "Toasted!" << endl;
      return false;
    }
  }
  return true;
}

int Search::Evaluate(SearchStack & stack, const ConsensOverlapUnit & COUnit)
{
    SearchStack minimal;
    stack.Minimal(minimal);
    //cout << "Solution." << endl;

 
    int len = minimal.Length(COUnit);
    //if (m_exhaust && len < m_minAltKeep)
    //return -1;

    // TRY this to reduce processing tiny contigs.
    if (m_results.isize() > 0 && minimal.Size() < 3)
        return -1;

    SearchStack to_test = minimal;
    to_test.Reverse();
 
    m_workHyp.clear();


    MakeHypothesis(m_workHyp, to_test, COUnit, (m_exhaust & !m_override) ); //TEST Last parameter.
    int to, from;
    int pairs = CountPairs(to, from, m_workHyp, COUnit);
    //cout << "Pairs: " << pairs << " to: " << to << endl;

    if (to >= 0) {
        //cout << "Trimming hypothesis: " << to << endl;
        //TESTING!!!!!!!
        //minimal.Trim(COUnit, to);
    }

 
    minimal.SetPairs(pairs);
  
    if (m_exhaust) {    
      if (m_results.isize() == 0) {
	m_results.push_back(minimal);
      } else {
	if (!IsNew(minimal, COUnit)) { // Not different enough??
	  return to;
	}
	if (m_results.isize() < m_maxResults) {
	  m_results.push_back(minimal);  
	} else {
	  m_results.Sort();
	  m_override = true;
	  if (m_results[0] < minimal) {     
	    m_results[0] = minimal;
	  }
	}
      }


      /*
        //cout << "Num hyps: " << m_results.isize()  << endl;
        if (m_override || m_results.isize() > m_maxResults) {
            if (!m_override)
                cout << "WARNING: Hypothesis buffer overflow!! Resetting to override." << endl;
            // NOTE: Do something more sensible here.
            m_override = true;
            if (m_results[0] < minimal) {
                m_results[0] = minimal;
            }
        } else {
            if (pairs > 0 || m_lastNoPairs == -1) {
                if (pairs == 0)
                    m_lastNoPairs = m_results.isize();
                m_results.push_back(minimal);
            } else {
                // No pairs, let's keep only one here.
                if (m_results[m_lastNoPairs] < minimal) {
                    m_results[m_lastNoPairs] = minimal;
                }
            }
	    }*/

    } else {
        // Keep only one top-N (for memory)
        if (m_results.isize() == 0) {
            m_results.push_back(minimal);
        } else {
            if (m_results[0] < minimal) {
                m_results[0] = minimal;
            }
        }
    }
    return to;
}

int Search::SelectLeftest(const ConsensOverlapUnit & COUnit, bool rc)
{
    m_results.Sort();
    if (m_results.isize() == 0) {
        //cout << "(none)" << endl;
        return -1;
    }

    int hypLen = m_results[m_results.isize()-1].Length();

    //cout << "Length: " << hypLen << endl;
    //m_results[m_results.isize()-1].Print(COUnit);

    SearchStack & result = m_results[m_results.isize()-1];
    result.Reverse();
 
    m_workHyp.clear();
    MakeHypothesis(m_workHyp, result, COUnit/*, true*/);
    if (rc) {
        m_workHyp.Reverse(hypLen);
    }

    for (int i=0; i<m_workHyp.isize(); i++) {
        if (m_workHyp[i].Ori() == 1)
            return m_workHyp[i].Read();
    }
    return -1;
}


void Search::SelectTopN(const ConsensOverlapUnit & COUnit, bool rc)
{
    //cout << "Printing longest hypothesis." << endl;
    if (m_results.isize() == 0) {
        //cout << "(none)" << endl;
        return;
    }
    m_results.Sort();
    int i, j;
    int hypLen = m_results[m_results.isize()-1].Length();
    int pairs = m_results[m_results.isize()-1].Pairs();

    //WARNING: Disable this for now (debugging!)
    if (m_exhaust) {
        cout << "Reporting ALL hypotheses." << endl;
 
        svec<Hypothesis> raw;
        raw.resize(m_results.isize());
        for (i=0; i<m_results.isize(); i++) {
            //cout << "# " << i << " pairs: " << pairs << " Length: " << hypLen << endl;
            SearchStack resultN = m_results[i];
            resultN.Reverse();
            //result.Print(COUnit);
            MakeHypothesis(raw[i], resultN, COUnit, !m_override);
            if (rc) {
                raw[i].Reverse(hypLen);
            }
        }

        cout << "Starting to toast..." << endl;
        // Toast hypotheses... (get rid of the ones that are too similar)
        for (i=raw.isize()-1; i>0; i--) {
            if (raw[i].isize() == 0)
                continue;
            //cout << "Checking hyp " << i << endl;
            for (j=i-1; j>=0; j--) {
                if (raw[j].isize() == 0)
                    continue;
                if (raw[i].ContainsSubset(raw[j])) {
                    cout << "Toasting hypothesis " << j << " because of " << i << endl;
                    raw[j].clear();
                }
                if (raw[j].ContainsSubset(raw[i])) {
                    cout << "Toasting hypothesis (reverse)" << i << " because of " << j << endl;
                    raw[i] = raw[j];
                    raw[j].clear();
                }
            }
        }

        bool bCont = false;
        cout << "Remaining hypotheses." << endl;
        for (i=raw.isize()-1; i>=0; i--) {
            if (raw[i].isize() == 0)
                continue;
            int to, from;
            CountPairs(to, from, raw[i], COUnit, true);

            raw[i].TrimRight(to);
            raw[i].TrimLeft(from);
            SetPairs(raw[i], COUnit);
            //raw[i].RemoveUnpaired();
      
            m_sink.Dump(raw[i], COUnit, bCont, m_minAltKeep);
            Commit(raw[i]);
            bCont = true;
        }
        cout << "done." << endl;
        return; // Done here.
    } else {
        //cout << "Pairs:" << pairs << " Length: " << hypLen << endl;
    }

    SearchStack & result = m_results[m_results.isize()-1];
    //result.Print(COUnit);
    result.Reverse();
 
    m_workHyp.clear();
    MakeHypothesis(m_workHyp, result, COUnit, !m_override);
    if (rc) {
        m_workHyp.Reverse(hypLen);
    }

    int to, from;
    //cout << "FINAL Check, " << m_sink.LastContig() << endl;
    CountPairs(to, from, m_workHyp, COUnit, true);
    m_workHyp.TrimRight(to);
    m_workHyp.TrimLeft(from);
    SetPairs(m_workHyp, COUnit);
    //hyp.RemoveUnpaired();
    m_sink.Dump(m_workHyp, COUnit);
    //cout << "Written to" << m_sink.LastContig() << endl;
    Commit(m_workHyp);
 
}

void Search::SetPairs(Hypothesis & hyp, const ConsensOverlapUnit & COUnit)
{
    int i, j;

    //cout << "Listing present (before)" << endl;
    //for (i=0; i<m_present.isize(); i++) {
    //  cout << i << " " << m_present[i] << endl;
    //}
    m_present.clear();
    m_present.resize(COUnit.GetNumReads(), -1);

    //cout << "Listing present (after)" << endl;
    //for (i=0; i<m_present.isize(); i++) {
    //  cout << i << " " << m_present[i] << endl;
    //}

    for (i=0; i<hyp.isize(); i++) {
        HypothesisNode & h = hyp[i];
        int r = h.Read();
        int numPartner = COUnit.getNumOfPartners(r);
        if (numPartner > 0) {
            for (int x = 0; x<numPartner; x++) {
                //cout << "Set partner " << r << " -> " << COUnit.getPartner(r, x) << endl;
                m_present.Set(COUnit.getPartner(r, x), i);
            }
        } else {
            h.SetPaired(-2);      
        }
    }
    for (i=0; i<hyp.isize(); i++) {
        HypothesisNode & h = hyp[i];
        int r = h.Read();
        //if (present[r] > -1) {
        if (h.Pair() != -2)
            h.SetPaired(m_present[r]);
        //}
    }
  
}

int Search::CountPairs(int & to, int & from, const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrint)
{
    if (m_pairDir == 0)
        return CountUnPairs(to, from, hyp, COUnit, bPrint);

    bPrint = false;

    if (bPrint) 
        cout << "Start CountPairs." << endl;

    to = -1;
    from = 0;

    // Stupid heuristics!!
    if (hyp.isize() < 2)
        return 0;

    //  bPrint = true;
    int i, j;
    // INEFFICIENT!!!!!!
    //svec<int> present;
    //svec<int> cov_seq;
    //svec<int> cov_pair;

    int max = 0;
    int havePartner = 0;
    for (i=0; i<hyp.isize(); i++) {
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
    for (i=0; i<hyp.isize(); i++) {
        const HypothesisNode & h = hyp[i];
        int r = h.Read();
        int numPartner = COUnit.getNumOfPartners(r);
        for (int x=0; x<numPartner; x++) {
            int partner = COUnit.getPartner(r, x);
            if (x == 0 || m_present[partner] < 0)
                m_present.Set(partner, i);
        }
    }
    for (i=0; i<hyp.isize(); i++) {
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
            //if (bPrint && h.Read() == 771) {
            //cout << "Located right " << endl;
            //}
 
            int start = h.Start();
            int stop = g.Stop();
            //if (stop < start) {
            //int tmp = start;
            //start = stop;
            //stop = tmp;
            //}
            // WARNING: Hard-coded limit for library size!!!!
 
            // Enforce proper location of pairs, i.e. pointing towards each other

            int plusminus = 1;
            //if (g.Ori() == h.Ori())
            //plusminus = -1;
            if (g.Ori() != m_pairDir * h.Ori())
                continue;

            //if (bPrint && h.Read() == 771) {
            //cout << "Adding pair info " << start + m_discount << " - " << stop - m_discount << " len: " << stop - start << endl;
            //}

            //if (stop - start > estPairSize) {
            //	for (j=start + m_discount; j<stop-m_discount; j++)
            //	  m_cov_pair_strict.Set(j, m_cov_pair_strict[j]+plusminus);
            //}
     
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
        /*if (havePartner == 0) { // Simple & stupid heuristics for unpaired data sets (should make this configurable!!)
          if (m_cov_seq[i-1] > 1 && m_cov_seq[i] <= 1) {
          to = i;
          if (bPrint) {
          cout << "Found SINGLE break point at pos " << to << endl;
          }
          break;

          }
          }*/
    


        //if (m_cov_seq_strict[i] > 5 || m_cov_pair_strict[i] > 0)
        // bYesGoodCov = true;

        /*if (bYesGoodCov && m_cov_seq_strict[i] == 0 && m_cov_pair_strict[i] == 0) {
          to = i + m_discount;
          if (bPrint) {
          cout << "Found break point (STRICT) at pos " << to << endl;
          }
          break;

          }*/
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

    for (i=0; i<hyp.isize(); i++) {
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

    if (bPrint) 
        cout << "Start CountUnPairs." << endl;

    to = -1;
    from = 0;

    // Stupid heuristics!!
    if (hyp.isize() < 2)
        return 0;

    int i, j;
 
    int max = 0;
    for (i=0; i<hyp.isize(); i++) {
        const HypothesisNode & h = hyp[i];
        if (h.Stop() > max)
            max = h.Stop();    
    }
 
    m_cov_seq.clear();
    m_cov_seq.resize(max, 0);
 
    int total = 0;
    for (i=0; i<hyp.isize(); i++) {
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
    //bAddIn = false;

    const svec<SearchNode> & st = ss.Raw();
    const ConsensReads & consReads = COUnit.getConsReads();
  
    int i, j;
    int off = 0;

    //cout << "START" << endl;

    //******************** Inefficient!!!! *********************  
    m_localUsed.clear();

    // HARD CODED!!!!!
    //if (ss.Size() > 300)
    //bAddIn = false;
  
    // WARNING: Dynamic allocation (get rid of it!!!)
    //svec<int> ids;

    if (bAddIn) {
        m_localUsed.resize(consReads.getNumOfReads(), 0);
        m_ids.resize(ss.Size());
        for (i=0; i<ss.Size(); i++) {
            m_ids[i] = st[i].Read();
            m_localUsed.Set(m_ids[i], 1);
            //if (m_ids[i] == 295)
            //cout << "Read 295 is alredy in HERE" << endl;
        }
        m_ids.Sort();
    }

    for (i=0; i<ss.Size(); i++) {
        HypothesisNode node;
        node = st[i];
    
        int dSize = consReads.getSize(node.Read());
        off += node.Offset();
        node.SetCoords(off, off + dSize);
        hyp.Add(node);

        if (bAddIn && i+1 < ss.Size()) {
            int offLimit = st[i+1].Offset();
            for (j=0; j<COUnit.GetNumLaps(node.Read()); j++) {
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
                if (m_ids.BinSearch(partner) < 0) {
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
                    t.SetCoords(off + t.Offset(), off + t.Offset() + dSize);
                    hyp.Add(t);	
                    m_localUsed.Set(s.getOverlapIndex(), 1);
                } else {
                }   
            }
        }
    }
    hyp.Sort();
}


void Search::Commit(const Hypothesis & hyp)
{
    int i;
    for (i=0; i<hyp.isize(); i++) {
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

bool Search::DoSearchAll(const ConsensOverlapUnit & COUnit, int startWithRead)
{
    m_override = false;
    //cout << "Enter main loop." << endl;
    m_usage.Resize(COUnit.GetNumReads());

    if (m_globalUsed.isize() == 0)
        m_globalUsed.resize(COUnit.GetNumReads(), 0);

    int i;

 
    //cout << "Processing all remaining COUnit (fw/rc)." << endl;
    for (i=startWithRead; i<COUnit.GetNumReads(); i++) {
        if (!IsUsedGlobal(i)
            && HasExtensions(COUnit, i)) {
            //DoSearch(COUnit, 17);
            //cout << "SEARCHING." << endl;
            //cout << "Finding start read: " << i << endl;
            int start = DoSearch(COUnit, i, true);
            //cout << "Searching for real: " << start << endl;
            DoSearch(COUnit, start);
        }
    }
    return true;
}

bool Search::HasExtensions(const ConsensOverlapUnit & COUnit, int id) const
{
    int i;
    int ori = 1;
    int count = 0;
    for (i=0; i<COUnit.GetNumDirLaps(id, ori); i++) {
        const ReadOverlap & l = COUnit.GetDirLap(id, i, ori);	
        if (m_globalUsed[l.getOverlapIndex()] == 0)
            count++;
    }
    return (count > 0);
}

int Search::DoSearch(const ConsensOverlapUnit & COUnit, int index, bool rc)
{
    //cout << "Start searching read " << index << endl;

    m_lastNoPairs = -1;
    m_usage.Clear();

    if (m_globalUsed.isize() == 0)
        m_globalUsed.resize(COUnit.GetNumReads(), 0);

    m_used.clear();
    m_used.resize(COUnit.GetNumReads(), 0);
    m_results.clear();
  
    if (m_exhaust)
        m_results.reserve(m_maxResults);
    else
        m_results.reserve(1);
  
 
    //m_usedFW.clear();
    //m_usedFW.resize(COUnit.GetNumReads(), 0);
    //m_usedRC.clear();
    //m_usedRC.resize(COUnit.GetNumReads(), 0);
 
    SearchStack stack;
    SearchNode init(index);
    if (rc) {
        init.SetOri(-1);
    }
    stack.Push(init);
    SetUsed(init);
    if (m_exhaust)
        m_usage.SetUsed(init.Read(), init.Pos());

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
  
        if (index == COUnit.GetNumDirLaps(n.Read(), ori)) {
            bool bReuse = false;
            if (!n.Ext()) {
                // Store hypothesis
                //cout << "End, evaluating." << endl;
	
                // NOTE: This would be mch faster, but it doesn't work quite yet.
                //cout << "Call Evaluate." << endl;
                bReuse = true;
                int limit = Evaluate(stack, COUnit);

                //break; //MGG: This should make the search greedy!!!

                if (limit < n.Pos()) {	  

                    // Do or do not not backoff all the way!!
                    backoff = limit;
                } else {
                    backoff = -1;
                }
            }
            // Go backwards
            SearchNode pop;
            int popCount = 0;
            do {
                stack.Pop(pop);
      
                // IMPORTANT NOTE: Revisit this condition!!!!
                // The current implementation is slightly too greedy
                if (m_exhaust) {
                    if (backoff > 0 && !(pop.Pos() <= backoff)) {
                        // cout << "Don't re-use
                    } else {
                        SetUsed(pop, false); // Let's allow for re-usage
                    }
                }
                popCount++;
                //cout << "Popped " << pop.Read() << " " << pop.Pos() << endl;
            } while(backoff >=0 && pop.Pos() > backoff);
            backoff = -1;
            continue;
        }
    
        //cout << "Trying to extend " << n.Read() << " " << n.Pos() << endl;
        // int ext = index;
        //for (i=ext; i>=0; i--) {
        for (i=index; i < COUnit.GetNumDirLaps(n.Read(), ori); i++) {
            const ReadOverlap & l = COUnit.GetDirLap(n.Read(), i, ori);	
            SearchNode to_push(l.getOverlapIndex(), curr, ori*l.getOrient(), l.getContactPos(), pos + l.getContactPos());
            n.IncCounter();

            if (!m_globalUsed[to_push.Read()] && !IsUsed(to_push) &&
                !m_usage.IsUsed(to_push.Read(), to_push.Pos())) {
                n.SetExt();
                stack.Push(to_push);
                //cout << " -> " << to_push.Read() << " " << to_push.Pos() << " " << COUnit.Sequences().Name(to_push.Read()) << endl; 
                SetUsed(to_push);
                if (m_exhaust)
                    m_usage.SetUsed(to_push.Read(), to_push.Pos());
                break;
            }
        }
    }

    if (rc) {
        int z = SelectLeftest(COUnit, rc);
        //cout << "Return: " <<  z << endl;
        if (z < 0) {
            cout << "ERROR!!!" << endl;
        }
        return z;
    }

    SelectTopN(COUnit, rc);

    return -1;
}


