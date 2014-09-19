#define FORCE_DEBUG
#include "src/Ananas/SequenceGraph.h"

// Stupid recursive function...
bool GraphEnumerate::Search(svec<SeqGraph> & linear, const SeqGraph & in, int index)
{
  int i, j;
  const SeqNode & s = in[index];
  m_work.push_back(s);

  for (i=0; i<s.NumOut(); i++) {
    int next = s.Out(i);
    Search(linear, in, next);
  }
  if (s.NumOut() == 0) {
    linear.push_back(m_work);
  }

  m_work.pop();
  return true;
}

void GraphEnumerate::Enumerate(svec<SeqGraph> & linear, const SeqGraph & in)
{
  int i, j, k;

  // Search for start nodes
  for (i=0; i<in.isize(); i++) {
    //SeqGraph graph;
    if (in[i].NumIn() > 0) // This does not work if it is circular...
      continue;

    // cout << 

    m_work.clear();
    Search(linear, in, i);
 
  }
}

void SequenceBuild::Build(DNAVector & out, const vecDNAVector & in, const SeqGraph & g)
{
  int i, j;
  out.resize(0);
  for (j=0; j<g.isize(); j++) {
    const SeqCoords & cc = g[j].Contrib(0); // One is sufficient
    DNAVector tmp;
    tmp.SetToSubOf(in[cc.Seq()], cc.First(), cc.Last() - cc.First() + 1);
    out += tmp;
  }
}

void SequenceBuild::Build(vecDNAVector & all, const vecDNAVector & in, const svec<SeqGraph> & graph)
{
  int i, j;
  all.resize(0);

  for (i=0; i<graph.isize(); i++) {
    const SeqGraph & g = graph[i];
    DNAVector out;
    for (j=0; j<g.isize(); j++) {
      const SeqCoords & cc = g[j].Contrib(0); // One is sufficient
      DNAVector tmp;
      tmp.SetToSubOf(in[cc.Seq()], cc.First(), cc.Last() - cc.First() + 1);
      out += tmp;
    }
    bool bTaken = false;
    for (j=0; j<all.isize(); j++) {
      double ident = all[j].FindIdent(out);
      if (ident >=0.995) {
	bTaken = true;
	break;
      }	
    }
    if (!bTaken)
      all.push_back(out);
  }
  
}


GraphConstructor::GraphConstructor(const vecDNAVector * pDna) {
  m_pDna = pDna;

 
}


// This method is a prime example for bad programming... on many levels :(
void GraphConstructor::Construct(SeqGraph & out)
{
  int i, j, k, l;


  // Do this pairwise.
  for (i=0; i<m_pDna->isize(); i++) {
    for (l=i+1; l<m_pDna->isize(); l++) {
      cout << "SIZE: " << (*m_pDna)[l].isize() << endl;
      vecDNAVector tmpBases;
      tmpBases.push_back((*m_pDna)[i]);
      vecDNAVector tmpBasesTarget;
      tmpBasesTarget.push_back((*m_pDna)[l]);
      svec<int> contigStarts;
      svec<int> contigDevs;
      contigStarts.push_back(0);
      contigDevs.push_back(0);

      // Super inefficient to have this allocated in each instance!
      KmerSuperAligner sa;
      sa.SetWordSize(12);
      sa.SetNumKmers(2);
      sa.SetLookAhead(0);
      sa.SetNewLookahead(12, 2);

      sa.SetRefBases(tmpBasesTarget);

      SuperAlign result;
      sa.Align(result, tmpBases, contigStarts, contigDevs);
    

      int lastStartQuery = 0;
      int lastStartTarget = 0;
      int lastQueryNodeIndex = -1;
      int lastTargetNodeIndex = -1;
      int lastJointNodeIndex = -1;
      
       //int lastEndQuery = m_pDna[i]->isize()-1;
      //int lastEndTarget = m_pDna[l]->isize()-1;

      for (j=0; j<result.GetMatchCount(); j++) {
	const SuperMatch & m = result.GetMatch(j);
	int contig = m.GetContig();
	int startQuery = m.GetFirstBase();
	int endQuery =  m.GetLastBase();
	
	int targetID = m.GetRefID();
	bool targetRC = m.GetRC();
	int startTarget = m.GetRefStart();
	int endTarget = m.GetRefEnd();
	int len = m.GetLastBase() - m.GetFirstBase() + 1;
	
	// Target (single contrib)
	
	if (startTarget > lastStartTarget) {
	  lastTargetNodeIndex = out.AddSeqNode();
	  out.AddCoords(lastTargetNodeIndex, l, lastStartTarget, startTarget-1);
	  out.Connect(lastJointNodeIndex, lastTargetNodeIndex);
	} else {
	  lastTargetNodeIndex = lastJointNodeIndex;
	}
	if (startQuery > lastStartQuery) {
	  lastQueryNodeIndex = out.AddSeqNode();
	  out.AddCoords(lastQueryNodeIndex, i, lastStartQuery, startQuery-1);
	  out.Connect(lastJointNodeIndex, lastQueryNodeIndex);
	} else {
	  lastQueryNodeIndex = lastJointNodeIndex;
	}
	// Add overlap
	lastJointNodeIndex = out.AddSeqNode();
	out.AddCoords(lastJointNodeIndex, l, startTarget, endTarget);
	out.AddCoords(lastJointNodeIndex, i, startQuery, endQuery);
 
	out.Connect(lastQueryNodeIndex, lastJointNodeIndex);
	out.Connect(lastTargetNodeIndex, lastJointNodeIndex);

	lastQueryNodeIndex = -1;
	lastTargetNodeIndex = -1;
	//int lastJointNodeIndex = -1;
 	
	lastStartTarget = endTarget + 1;
	lastStartQuery = endQuery + 1;


	cout << "Report match " << contig << " vs " << targetID << " ";
	cout << startQuery << "-" << endQuery << " and " << startTarget << "-" << endTarget << endl;

	cout << "SIZE Check (2): " << (*m_pDna)[l].isize() << endl;
  	
	const DNAVector & a = (*m_pDna)[i];
	const DNAVector & b = (*m_pDna)[l];
	for (k=startQuery; k<=endQuery; k++) {
	  cout << a[k];
	}
	cout << endl;
	for (k=startTarget; k<=endTarget; k++) {
	  cout << b[k];
	}
	cout << endl;
	cout << "SIZE Check (3): " << (*m_pDna)[l].isize() << endl;
      }
      cout << "SIZE Check (4): " << (*m_pDna)[l].isize() << endl;
       

      cout << "Try LAST Target " << l << " " << lastStartTarget << " " << (*m_pDna)[l].isize() << endl;
      if (lastStartTarget <  (*m_pDna)[l].isize()) {
	int index = out.AddSeqNode();
	cout << "LAST Target " << l << " " << lastStartTarget << " " << (*m_pDna)[l].isize() << endl;
	out.AddCoords(index, l, lastStartTarget, (*m_pDna)[l].isize()-1);
	out.Connect(lastJointNodeIndex, index);
      }
      if (lastStartQuery  < (*m_pDna)[i].isize()) {
	int index = out.AddSeqNode();
	cout << "LAST Query" << endl;
	out.AddCoords(index, i, lastStartQuery, (*m_pDna)[i].isize() - 1);
	out.Connect(lastJointNodeIndex, index);
      }
     
    }
  }
}

void GraphConstructor::Sequence(DNAVector & out, const SeqGraph & in)
{
}



