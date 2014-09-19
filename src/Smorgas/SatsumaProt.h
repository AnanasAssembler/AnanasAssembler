#ifndef _SATSUMA_PROT_H_
#define _SATSUMA_PROT_H_

#include "aligns/KmerAlignCore.h"
#include <string>

#include "src/CrossCorr.h"
#include "src/Satsuma/XCorrDynProg.h"
#include "util/SysTime.h"
#include "src/MakeSixFrames.h"
#include "src/MultiProtein.h"
#include "src/MultXCorr.h"
#include "src/SequenceStream.h"

class ProteinXCorr : public MultiSizeXCorr
{
};

class Hit
{
public:
  Hit(int contig, int matches) {
    m_contig = contig;
    m_matches = -matches;
  }
  Hit() {
    m_contig = -1;
    m_matches = 0;
  }

  int Contig() const {return m_contig;}
  int Seeds() const {return m_matches;}

  bool operator < (const Hit & h) const {
    return (m_matches < h.m_matches);
  }

  void operator = (int i) {
    m_contig = i;
  }

private:
  int m_contig;
  int m_matches;
};

//TODO Needs commenting and access level corrections
class SatsumaProt; //Forward declaration
class SatsumaProtParams 
{
friend class SatsumaProt;
public:
  SatsumaProtParams(bool bExhaust,  bool bSelf, bool bSame, int protein_K, int kmerStep,
                    int filter, int failAllowed, int block, int n_blocks, int limit, 
                    bool bRNA, bool bQuiet, double cutoff, int slide, float eThresh
                    ):bExhaust(bExhaust), bSelf(bSelf), bSame(bSame), protein_K(protein_K), 
                      kmerStepSize(kmerStep), prefilterType(filter), allowFailHits(failAllowed), 
                      block(block),  n_blocks(n_blocks), limit(limit), bRNA(bRNA), bQuiet(bQuiet),
                      ungapped_cutoff(cutoff), filterWSlide(slide), eValThresh(eThresh) {}
protected:
  bool bExhaust;          ///Do not prefilter and do exhaustive search
  bool bSelf;             ///Include self alignments
  bool bSame;
  int  protein_K;         ///K-mer size
  int  kmerStepSize;      ///Step size for generating kmers in the KmerAlignCore, should be set to 1 for exhaustive set
  int  prefilterType;     ///1:max fixed distance k-mer 2:maximum number of k-mers 3.pre-filter 2 and filter with 1
  int  allowFailHits;     ///Number of hit failures to allow before stopping search
  int  block;             ///Search only this subset 
  int  n_blocks;          ///Number of blocks
  int  limit;             ///Number of results to display
  bool bRNA;              ///Do RNA alignment
  bool bQuiet;            ///Verbosity
  double ungapped_cutoff;
  int filterWSlide;       ///Filter 2 Window slide, if set to 1 the window sliding will cover all kmers.
  double eValThresh;      ///Threshold of Evalues to accept
};

class SatsumaProt 
{
public:
  SatsumaProt(string& refFile , string& dbName, string& queryFile, SatsumaProtParams& p):params(p) {
    init(refFile, dbName, queryFile);
  }

  void init(string& refFile, string& dbName, string& queryFile);
 
  void readTarget(string& refFile, string& dbName); 

  void readQuery(string& queryFile);

  void setQuery(const vecDNAVector & d);

  void alignAll(ostream & o = cout);  

 void filterTopHits_1(DNAVector d); 

 void filterTopHits_2(DNAVector d); 

 void filterTopHits_3(DNAVector d); 

private:
  void ResetQuals(vecDNAVector & d) 
  {
    int i, j;
    for (i=0; i<d.isize(); i++) {
      for (j=0; j<d[i].isize(); j++)
        d[i].SetQual(j, 0);
    }
  }

  void ResetOneQual(DNAVector & d, int diff) {
    if (diff < 0)
      diff = -diff;
    diff = (diff % d.isize());

    d.SetQual(diff, 0);  
  }

  inline int SetQual(DNAVector & d, int diff) {
    if (diff < 0)
      diff = -diff;
    diff = (diff % d.isize());
    int v = d.Qual(diff);
    if (v < 255) {
      v++;
      if (v < 0)
        cout << "ERROR!" << endl;
      d.SetQual(diff, v);
    }
    return v;
  }

  vecDNAVector   ref;                  ///target protein k-mer map
  SequenceStream queryStream;          ///query proteins as stream
  TranslateBasesToNumberProtein     trans; 
  KmerAlignCore<KmerAlignCoreRecord> core;
  ProteinXCorr xc;
  svec<Hit> allhits;
  svec<int> hit_ref;
  svec<int> hit_diff;
  svec<int> contrib;
  svec<int> filter1Passed;
  SatsumaProtParams params;
};
#endif //_SATSUMA_PROT_H_
