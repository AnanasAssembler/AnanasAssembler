#ifndef _READOVERLAP_H_
#define _READOVERLAP_H_

#include <stdint.h>
#include "base/SVector.h"
#include "src/Ananas/Reads.h"

//======================================================

class ReadOverlap 
{
//TODO changing bool value of direction and orientation to int to accord with Search code but this is NOT correct
public:
  ReadOverlap():m_overlapIndex(-1), m_contactPos(-1), m_score(-1.), m_direction(-1),  m_orient(-1) {};
  ReadOverlap(int oI, int cP, float s, int d, int o)
                 :m_overlapIndex(oI), m_contactPos(cP), m_score(s), m_direction(d), m_orient(o) { }

  void set(int oI, int cP, float s, int d, int o); 

  int     getOverlapIndex() const       { return m_overlapIndex;  }  
  int     getContactPos() const         { return m_contactPos;    }  
  float   getScore() const              { return m_score;         }
  int     getDirection() const          { return m_direction;     }
  int     getOrient() const             { return m_orient;        }

  bool operator < (const ReadOverlap & rO) const; 

private: 
  int     m_overlapIndex;  /// The index of the read to which this item refers to
  int     m_contactPos;    /// The index in the read where this overlap occurs from 
  float   m_score;         /// Identity score of the overlap
  int     m_direction;     /// The direction of the overlap (true: right  false: left)
  int     m_orient;        /// The orientation of the overlap (Whether the reads are of the same strand)
};
//======================================================


//======================================================
class ReadInfo
{
public:
  ReadInfo(): m_rightOverlaps(), m_leftOverlaps(), m_overlapIds() {}

  const svec<ReadOverlap>& getOverlaps(int isRight) const { 
    if(isRight==1) {  return m_rightOverlaps; } 
    else           {  return m_leftOverlaps;  } 
  }

  void addOverlap(int overlapIndex, int contactPos, float score, int direction, int orient) {
    addOverlap(ReadOverlap(overlapIndex, contactPos, score, direction, orient)); 
  }

  void addOverlap(const ReadOverlap & oL) {
    if(oL.getDirection()==1) { m_rightOverlaps.push_back(oL); }
    else                     { m_leftOverlaps.push_back(oL);  }
    m_overlapIds.push_back(oL.getOverlapIndex());
  }
    
  void sortLaps(); 

  bool organizeLaps(const ConsensReads & allReads, int id);  

  int getNumLaps() const { return m_rightOverlaps.isize() + m_leftOverlaps.isize(); }
  const ReadOverlap & getLap(int i) const {
    if(i<getNumRightLaps()) { return m_rightOverlaps[i]; }
    else { return m_leftOverlaps[i-getNumRightLaps()]; }
  }
  int getNumRightLaps() const {return m_rightOverlaps.isize();}
  const ReadOverlap & getRightLap(int i) const {return m_rightOverlaps[i];}
  int getNumLeftLaps() const {return m_leftOverlaps.isize();}
  const ReadOverlap & getLeftLap(int i) const {return m_leftOverlaps[i];}
  
  int getNumDirLaps(int dir) const { return (dir==1? getNumRightLaps(): getNumLeftLaps()); }

  const ReadOverlap & getDirLap(int i, int dir) const { return (dir==1? m_rightOverlaps[i]: m_leftOverlaps[i]); }

  bool hasLap(int id) const { return (binary_search(m_overlapIds.begin(), m_overlapIds.end(), id)); }
  
protected:
  bool isChimera(int len);

private:

  //int m_partner;  //TODO decide where this belongs
  svec<ReadOverlap> m_rightOverlaps;       /// Overlaps that extend from the right side
  svec<ReadOverlap> m_leftOverlaps;        /// Overlaps that extend from the left side
  svec<int>         m_overlapIds;          /// For fast searching
};
//======================================================


//======================================================
class AllReadOverlaps
{
public:

  AllReadOverlaps():m_overlaps(), m_chimera() {}

  //  TODO Can only be constructed with size as parameter so that user is aware that size should be set 
  AllReadOverlaps(int size):m_overlaps(), m_chimera() {resize(size);}
   
  const ReadInfo& operator[](int i) const                  { return m_overlaps[i];                      }
  const svec<ReadOverlap>& getRightOverlaps(int idx) const { return m_overlaps[idx].getOverlaps(1);     }
  const svec<ReadOverlap>& getLeftOverlaps(int idx)  const { return m_overlaps[idx].getOverlaps(-1);    }
 
  void addOverlap(int readIndex, int overlapIndex,
                  int contactPos, float score, int direction, int orient) {
    m_overlaps[readIndex].addOverlap(overlapIndex, contactPos, score, direction, orient);
  }

  void resize(int sz) { 
    m_overlaps.resize(sz);
  }

  void write(const string& overlapFile, int mode) const; 
  void writeAsc(const string& readOverlapFile) const; 
  void writeBin(const string& readOverlapFile) const; 
  void read(const string& readOverlapFile, int mode); 
  void readAsc(const string& readOverlapFile, const svec<int> & good, const ConsensReads& consReads); 
  void readAsc(const string& readOverlapFile); 
  void readBin(const string& readOverlapFile); 

  bool isChimeric(int i) const { return m_chimera[i]; }

private:
  const svec<ReadOverlap>& getReadOverlaps(int readIndex, int isRight) const { 
    return m_overlaps[readIndex].getOverlaps(isRight);
  }

  void postReadActions(const ConsensReads& consReads); 

  /** Internal helper class for writing/serialization */
  string getOverlapString(int index, int dir) const; 
  void addOverlapFromString(const string& strIn);

  svec<ReadInfo> m_overlaps;     /// All Overlap info , one item per read
  svec<int> m_chimera;           /// TODO comment
};
//======================================================

#endif // _READOVERLAP_H_
