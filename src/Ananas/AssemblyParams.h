#ifndef _ASSEMBLYPARAMS_H_
#define _ASSEMBLYPARAMS_H_
//======================================================
class AssemblyParams 
{
public:
  AssemblyParams(bool singleStrand=true, int stepSize=10, int seedSize=15, 
                  float minIdent=0.75, float minEndCoverage=0.98,  int minOverlap=40, 
                  int alignmentBound=3, int minBasePerScaf=200)
                 :m_singleStrand(singleStrand), m_subReadStep(stepSize), m_seedSize(seedSize),
                  m_minIdent(minIdent), m_minEndCoverage(minEndCoverage), m_minOverlap(minOverlap),
                  m_alignmentBound(alignmentBound), m_minBasePerScaf(minBasePerScaf) { }

    bool  isSingleStrand() const    { return m_singleStrand;   }
    int   getSubreadStep() const    { return m_subReadStep;    }  
    int   getSeedSize()  const      { return m_seedSize;       } 
    float getMinIdentity() const    { return m_minIdent;       }
    float getMinEndCover() const    { return m_minEndCoverage; }  
    int   getMinOverlap() const     { return m_minOverlap;     }
    int   getAlignBand() const      { return m_alignmentBound; }
    int   getMinBasePerScaf() const { return m_minBasePerScaf; } 

    void  setSubreadStep(int sst)   { m_subReadStep    = sst;  }  
    void  setSeedSize(int ss)       { m_seedSize       = ss;   } 
    void  setMinIdentity(float idt) { m_minIdent       = idt;  }
    void  setMinEndCover(int mec)   { m_minEndCoverage = mec;  }  
    void  setMinOverlap(int mo)     { m_minOverlap     = mo;   }
    void  setAlignBand(int ab)      { m_alignmentBound = ab;   }
    void  setMinBasePerScaf(int mb) { m_minBasePerScaf = mb;   } 


private: 
  bool    m_singleStrand;   /// Flag specifying whether the reads are single or double strand
  int     m_subReadStep;    /// Block step size used for constructing subreads
  int     m_seedSize;       /// Seed size for selecting candidate reads for assembly extension
  float   m_minIdent;       /// Minimum identity for accepting a candidate read as an assembly extension
  float   m_minEndCoverage; /// Minimum coverage of the ends of the read for accepting an assembly extension
  int     m_minOverlap;     /// Minimum overlap of a read for accepting as an assembly extension
  int     m_alignmentBound; /// Alignment bandwidth used for local alignment to decide on choosing candidate reads
  int     m_minBasePerScaf; /// Minimum number of reads per consstructed scaffold that will be acceptable 
};
//======================================================

#endif // _ASSEMBLYPARAMS_H_
