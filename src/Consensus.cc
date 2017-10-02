#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "Consensus.h"
#include "SearchOverlaps.h"
#include "ConsensOverlapUnit.h"
#include "ContScaff.h"

class Kmer
{
public:
    Kmer() {
        m_pos = -1;
        m_mul = -1;
    }
    Kmer(const DNAVector & d, int from, int size, int pos_base) {
        m_mul = -1;
        Set(d, from, size, pos_base);
    }

    void Set(const DNAVector & d, int from, int size, int pos_base) {
        d.Substring(m_kmer, from, size);
        m_pos = pos_base + from;
    }
 
    void IncMult() {
        m_mul++;
    }
 
    void SetMult(int m) {
        m_mul = m;
    }
    void SetPos(int i) {
        m_pos = i;
    }

    const string & GetKmer() const {return m_kmer;}
    int Pos() const {return m_pos;}
    int Mult() const {return m_mul;}

    bool operator < (const Kmer & k) const {
        //return m_kmer < k.m_kmer;
        if (m_kmer != k.m_kmer)
            return m_kmer < k.m_kmer;   
        return m_pos < k.m_pos;
    }
    bool operator != (const Kmer & k) const {
        return !operator == (k);
    }
    bool operator == (const Kmer & k) const {
        int diff = m_pos - k.m_pos;
        if (diff < 0)
            diff = -diff;
        // HARD CODED!!!!
        if (m_kmer == k.m_kmer && diff < 10)
            return true;
        return false;
    }

    void Extend(char l) {
        int i;
        for (i=1; i<m_kmer.size(); i++) {
            m_kmer[i-1] = m_kmer[i];
        }
        m_kmer[m_kmer.size()-1] = l;
    }

private:
    string m_kmer;
    int m_pos;
    int m_mul;
};


class KmerAssembly
{
public:
    KmerAssembly() {
        m_size = 128;
    }

    void Clear() {
        m_kmers.clear();
    }

    void Add(const DNAVector & d, int base_pos) {
        int i;
        for (i=0; i<=d.isize()-m_size; i++) {
            m_kmers.push_back(Kmer(d, i, m_size, base_pos));
        }
    }

    void Build(DNAVector & out) {
        Sort(m_kmers);

        int start = -1;
   
        int startMul = 0;
        int i, j;

        svec<Kmer> comp;
    
        Kmer tmp;
        for (i=0; i<m_kmers.isize(); i++) {
            tmp = m_kmers[i];
            tmp.SetMult(1);
            for (j=i+1; j<m_kmers.isize(); j++) {
                if (m_kmers[j] != m_kmers[i]) {
                    comp.push_back(tmp);
                    //cout << "Adding: " << tmp.GetKmer() << " " << tmp.Mult() << endl;
                    if (tmp.Pos() == 0 && tmp.Mult() > startMul) {
                        startMul = tmp.Mult();
                        start = comp.isize()-1;
                    } 
                    break;
                } else {
                    tmp.IncMult();
                    //if (m_kmers[j].Pos() < tmp.Pos()) {
                    //cout << "WERID!!" << endl;
                    //tmp.SetPos(m_kmers[j].Pos());
                    //}
                }
            }
            i = j-1;
        }

        // Last element (stupid!!!)
        comp.push_back(tmp);
        //cout << "Adding: " << tmp.GetKmer() << " " << tmp.Mult() << endl;
        if (tmp.Pos() == 0 && tmp.Mult() > startMul) {
            startMul = tmp.Mult();
            start = comp.isize()-1;
        } 


        // Now we have unique kmers w/ counts.
        comp.push_back(tmp);
        Sort(comp);
        svec<int> used;
        used.resize(comp.isize(), 0);
        if (start == -1) {
            cout << "ERROR: Start " << start << endl;
            ///for (i=0; i<comp.isize()
        }
        m_kmers = comp; // Stupid...


        string result = comp[start].GetKmer();
        int best = -1;
        int lastPos = 0;
        do {
            const Kmer & k = comp[start];
            //cout << "Kmer: " << k.GetKmer() << endl;
      
            svec<Kmer> ext;
            ext.resize(4);
            for (i=0; i<ext.isize(); i++) {
                ext[i] = k;
                ext[i].SetPos(0); // Get the first one
                ext[i].Extend(NucLetter(i));
            }

            best = -1;
            int bestBase = -1;
            int bestCount = 0;
            for (i=0; i<ext.isize(); i++) {
                //cout << "Lookup " << ext[i].GetKmer() << endl;
                int index = BinSearchFuzzy(comp, ext[i]);
                index = FindBest(comp, ext[i], index, lastPos);
                if (index < 0)
                    continue;
                //cout << "Index=" << index << " pos " <<  comp[index].Pos() << endl;
	
                if (used[index] == 0 && comp[index].Mult() > bestCount) {
                    bestCount = comp[index].Mult();
                    best = index;
                    bestBase = i;
                }
            }
            if (bestBase != -1) {
                //cout << "Best ext: " << NucLetter(bestBase) << endl;
                if (best == -1)
                    cout << "ERROR!" << endl;
                used[best] = 1;
                result += NucLetter(bestBase);
                lastPos = comp[best].Pos();
                //cout << result << endl;
            } else {
	      //cout << "Missing..." << endl;
	    }
            start = best;
        } while (best != -1);
        out.SetFromBases(result);
    }
    

private:
    int FindBest(const svec<Kmer> & comp, Kmer & k, int index, int pos) {
        if (index < 0)
            return index;
        int i;
        int min = m_size;
        int bestIndex = -1;
        for (i=index; i<comp.isize(); i++) {
            if (comp[i].GetKmer() != k.GetKmer())
                break;
            int p = comp[i].Pos();
            int diff = p - pos - 1;
            //cout << "Compare " << pos << " " << p << " " << diff << " i=" << i << endl;
            if (diff < 0)
                diff = -diff * 2;
            if (diff < min) {
                min = diff;
                bestIndex = i;
            }
        }
        return bestIndex;
    }


    svec<Kmer> m_kmers;
    int m_size;
};


//-----------------------------------------------------------
void ConsensusBuilder::BuildWithGaps(DNAVector & out, const Contig& cont, const ConsensOverlapUnit & COUnit)
{
  /*
  const ConsensReads & consReads = COUnit.getConsReads();

  
  KmerAssembly kmers;
  
  int n = 0;
  for (int i=0; i<cont.isize(); i++) {
    int r = cont[i].Read();
    DNAVector d = consReads[r];
    if (cont[i].Ori() == -1)
      d.ReverseComplement();
    if (cont[i].Start() + d.isize() > n)
      n = cont[i].Start() + d.isize();
    kmers.Add(d, cont[i].Start());
    
  }    
  kmers.Build(out);*/


  const ConsensReads & consReads = COUnit.getConsReads();

  int len = 0;
  for (int i=0; i<cont.isize(); i++) {
    if (cont[i].Stop() > len) {
      len = cont[i].Stop();
    }
  }

  Consensus cons;
  cons.resize(len);
  
  int n = 0;
  for (int i=0; i<cont.isize(); i++) {
    int r = cont[i].Read();
    DNAVector d = consReads[r];
    if (cont[i].Ori() == -1)
      d.ReverseComplement();
    if (cont[i].Start() + d.isize() > n)
      n = cont[i].Start() + d.isize();
    for (int x=0; x<d.isize(); x++) {
      cons.Add(cont[i].Start() + x, d[x]);
    }
  }    
  out.resize(len);
  for (int i=0; i<len; i++) 
    out[i] = cons.GetFirst(i);
  
}

void ConsensusBuilder::Build(DNAVector & out, const Contig& cont, const ConsensOverlapUnit & COUnit)
{
  const ConsensReads & consReads = COUnit.getConsReads();

  int len = 0;
  for (int i=0; i<cont.isize(); i++) {
    if (cont[i].Stop() > len) {
      len = cont[i].Stop();
    }
  }

  Consensus cons;
  cons.resize(len);
  
  int n = 0;
  for (int i=0; i<cont.isize(); i++) {
    int r = cont[i].Read();
    DNAVector d = consReads[r];
    if (cont[i].Ori() == -1)
      d.ReverseComplement();
    if (cont[i].Start() + d.isize() > n)
      n = cont[i].Start() + d.isize();
    for (int x=0; x<d.isize(); x++) {
      cons.Add(cont[i].Start() + x, d[x]);
    }
  }    
  out.resize(len);
  for (int i=0; i<len; i++) 
    out[i] = cons[i];
}


LayoutSink::LayoutSink()
{
  m_counter = -1;
  m_minor = 0;
  m_pLayout = NULL;
  m_minIdent = 0.99;
  m_index  = 0;
  m_prefix = "Sample1";
}

LayoutSink::~LayoutSink()
{
  if (m_pLayout != NULL) {
    fprintf(m_pLayout, "<SCAFFOLD_READCOUNT> %s %d </SCAFFOLD_READCOUNT>\n", m_lastScaffName.c_str(), (int)m_currReads.size());
    fprintf(m_pLayout, "<SCAFFOLD_PAIRCOUNT> %s %d </SCAFFOLD_PAIRCOUNT>\n", m_lastScaffName.c_str(), (int)m_currPairedReads.size()/2);
    fprintf(m_pLayout, "</SCAFFOLD> %s\n\n", m_lastScaffName.c_str());
    fclose(m_pLayout);
  }
}
 
void LayoutSink::SetLayoutFile(const string & layout)
{
  if (m_pLayout != NULL) { fclose(m_pLayout); }
  m_pLayout = fopen(layout.c_str(), "w");
  if (m_pLayout == NULL) 
    cout << "ERROR: Could not open file " << layout << " for writing!!" << endl;
}

bool LayoutSink::IsSame(const DNAVector & d)
{
  int n = m_lastCons.isize();
  int m = d.isize();
  
  if (m < n)
    n = m;

  if (n == 0)
    return false;

  int i;
  int match = 0;
  
  for (i=0; i<n; i++) {
    if (m_lastCons[i] == d[i])
      match++;
  }
  double ident = (double)match/(double)n;
  if (ident >= m_minIdent)
    return true;
  else
    return false;

}

void LayoutSink::Dump(const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrev, int minLen)
{
  if (!bPrev || m_counter < 0) {
    m_counter++;
    m_minor = 0;
    m_lastCons.resize(0);
  } 

  if (bPrev)
    m_minor++;

  int currNumReads_contig       = 0;
  int currNumPairedReads_contig = 0;

  int i;
  char name[512];
  sprintf(name, ">Contig_%s_%3d_%7d_%3d", m_prefix.c_str(), m_index, m_counter, m_minor);
  for (i=0; i<(int)strlen(name); i++) {
    if (name[i] == ' ')
      name[i] = '0';
  }
  m_lastContigName = name;
  
  const ConsensReads & rr = COUnit.getConsReads();
  if(!bPrev) {
    if(m_counter>0) { 
      fprintf(m_pLayout, "<SCAFFOLD_READCOUNT> %s %d </SCAFFOLD_READCOUNT>\n", m_lastScaffName.c_str(), (int)m_currReads.size());
      fprintf(m_pLayout, "<SCAFFOLD_PAIRCOUNT> %s %d </SCAFFOLD_PAIRCOUNT>\n", m_lastScaffName.c_str(), (int)m_currPairedReads.size()/2);
      m_currPairedReads.clear();
      m_currReads.clear();
      fprintf(m_pLayout, "</SCAFFOLD> %s\n\n", m_lastScaffName.c_str());
    }
    char scaffName[1024];
    sprintf(scaffName, ">Scaffold_%5d", m_counter);
    for (unsigned int x = 0; x<strlen(scaffName); x++) {
      if (scaffName[x] == ' ')
      scaffName[x] = '0';
    }
    fprintf(m_pLayout, "<SCAFFOLD> %s\n", scaffName);
    m_lastScaffName = scaffName;
  }
  fprintf(m_pLayout, "<CONTIG> %s %d \n", name, hyp.Pairs());
  for (i=0; i<hyp.Size(); i++) {
    int rIdx = hyp[i].Read();
    int numPartner = COUnit.getNumOfPartners(rIdx);
    fprintf(m_pLayout, "%d\t%d\t%d - %d\t%d\t%d\n", rIdx, hyp[i].Ori(), hyp[i].Start(), hyp[i].Stop(), 
            hyp[i].Pair()>=0?hyp[hyp[i].Pair()].Read():hyp[i].Pair(),  
            hyp[i].Pair()>=0?hyp[hyp[i].Pair()].Ori():hyp[i].Pair());

    const svec<int>& consMemIds = COUnit.getConsMembers(rIdx);
    for(int cmId:consMemIds) { //Add all the raw reads
      m_currReads[cmId] = true;
    }
    currNumReads_contig += COUnit.getConsensCount(rIdx);
    if (hyp[i].Pair() >= 0) {
      for(int cmId:consMemIds) { //Add all the raw read pairs
        m_currPairedReads[cmId] = true;
      }
      currNumPairedReads_contig += COUnit.getConsensCount(rIdx);
    }
  }
  fprintf(m_pLayout, "<CONTIG_READCOUNT> %s %d </CONTIG_READCOUNT>\n", name, currNumReads_contig);
  fprintf(m_pLayout, "<CONTIG_PAIRCOUNT> %s %d </CONTIG_PAIRCOUNT>\n", name, currNumPairedReads_contig/2);
  fprintf(m_pLayout, "</CONTIG> %s\n", name);
  fflush(m_pLayout);
}

void LayoutSink::fastaFromAssembly(const string& fastaFile, Assembled& asmb, const ConsensOverlapUnit & COUnit, int minLen, bool bUseGaps)
{
  FILE* pFasta = fopen(fastaFile.c_str(), "w");
  if (pFasta == NULL) {
    cout << "ERROR: Could not open file " << fastaFile << " for writing!!" << endl;
  }
  for(int scafCnt=0; scafCnt<asmb.isize(); scafCnt++) {
    Scaffold& currScaff = asmb[scafCnt];
    for(int contCnt=0; contCnt<currScaff.isize(); contCnt++) {
      Contig& currCont = currScaff[contCnt];
      ConsensusBuilder cons;
      DNAVector dna;
      if (bUseGaps)
	cons.BuildWithGaps(dna, currCont, COUnit);
      else
	cons.Build(dna, currCont, COUnit);
      if (IsSame(dna)) {
        cout << "Deep toasting sequence " << currCont.Name() << endl;
        currCont.SetDiscard(true);
        continue;
      }
      if(dna.isize()<minLen) { 
        currCont.SetDiscard(true);
        continue; 
      }
      fprintf(pFasta, "%s\n", currCont.Name().c_str());
      for (int i=0; i<dna.isize(); i++) {
        fprintf(pFasta, "%c", dna[i]);
      }
      fprintf(pFasta, "\n");
      fflush(pFasta); 
    }
  }
}
