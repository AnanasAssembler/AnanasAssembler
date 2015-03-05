#ifndef _READOVERLAP_H_
#define _READOVERLAP_H_

#include <stdint.h>
#include "base/ThreadHandler.h"
#include "base/SVector.h"
#include "src/Ananas/Reads.h"

//======================================================

class ReadOverlap 
{
public:
    ReadOverlap():m_overlapIndex(-1), m_contactPos(-1), m_direction(false),  m_orient(false) {};
    ReadOverlap(int oI, int cP, int d, int o) {
      set(oI, cP, d, o);
    }

    void set(int oI, int cP, int d, int o); 

    string toString() const;

    int     getOverlapIndex() const       { return m_overlapIndex;  }  
    int     getContactPos() const         { return m_contactPos;    }  

    int getDirection() const { 
      return (m_direction)?1:-1;
    }
    int getOrient() const { 
      return (m_orient)?1:-1;
    }

    bool operator < (const ReadOverlap & rO) const; 

private: 
    int     m_overlapIndex;  /// The index of the read to which this item refers to
    int     m_contactPos;    /// The index in the read where this overlap occurs from 
    bool    m_direction;     /// The direction of the overlap (true (+1): right  false (-1): left) 
    bool    m_orient;        /// The overlap orientat (Whether the reads are of the same strand) - same logic as direction
};
//======================================================


//======================================================
class ReadInfo
{
public:
    ReadInfo(): m_rightOverlaps(), m_leftOverlaps(), m_overlapIds(), m_numRightOL(0), m_numLeftOL(0) {
        m_rightOverlaps.reserve(1000);
        m_leftOverlaps.reserve(1000);
        m_overlapIds.reserve(2000);
    }

    const svec<ReadOverlap>& getOverlaps(int isRight) const { 
        if(isRight==1) {  return m_rightOverlaps; } 
        else           {  return m_leftOverlaps;  } 
    }

    void addOverlap(int overlapIndex, int contactPos, int direction, int orient) {
        addOverlap(ReadOverlap(overlapIndex, contactPos, direction, orient)); 
    }

    void addOverlap(const ReadOverlap & oL) {
        if(oL.getDirection()==-1) { m_leftOverlaps.push_back(oL);  m_numLeftOL++;     }
        else                      { m_rightOverlaps.push_back(oL); m_numRightOL++;    }
        m_overlapIds.push_back(oL.getOverlapIndex());
    }
    
    void sortLaps(); 
    void sortOverlapIndexes(); 

    bool organizeLaps(const ConsensReads & allReads, int id);  

    const ReadOverlap & getLap(int i) const {
        if(i<getNumRightLaps()) { return m_rightOverlaps[i]; }
        else { return m_leftOverlaps[i-getNumRightLaps()]; }
    }
    int getNumLaps() const                       { return getNumRightLaps()+getNumLeftLaps(); }
    int getNumRightLaps() const                  { return m_numRightOL;                       }
    int getNumLeftLaps() const                   { return m_numLeftOL;                        }
    const ReadOverlap & getRightLap(int i) const { return m_rightOverlaps[i];                 }
    const ReadOverlap & getLeftLap(int i) const  { return m_leftOverlaps[i];                  }
  
    int getNumDirLaps(int dir) const                    { return (dir==1? getNumRightLaps(): getNumLeftLaps());                 }
    const ReadOverlap & getDirLap(int i, int dir) const { return (dir==1? m_rightOverlaps[i]: m_leftOverlaps[i]);               }
    bool hasLap(int id) const                           { return (binary_search(m_overlapIds.begin(), m_overlapIds.end(), id)); }
  
protected:
    bool isChimera(int len);

private:
    svec<ReadOverlap> m_rightOverlaps;       /// Overlaps that extend from the right side
    svec<ReadOverlap> m_leftOverlaps;        /// Overlaps that extend from the left side
    svec<int>         m_overlapIds;          /// For fast searching
    int               m_numRightOL;          /// Kept for saving on computing size of vector (freq usage)
    int               m_numLeftOL;           /// Number of left overlaps (as above for right overlaps)
};
//======================================================


//======================================================
class AllReadOverlaps
{
public:
    AllReadOverlaps():m_overlaps(), m_chimera(), m_mutex() {}

    //  TODO Can only be constructed with size as parameter so that user is aware that size should be set 
    AllReadOverlaps(int size):m_overlaps(), m_chimera(), m_mutex() {resize(size);}
   
    const ReadInfo& operator[](int i) const                  { return m_overlaps[i];                      }
    const svec<ReadOverlap>& getRightOverlaps(int idx) const { return m_overlaps[idx].getOverlaps(1);     }
    const svec<ReadOverlap>& getLeftOverlaps(int idx)  const { return m_overlaps[idx].getOverlaps(-1);    }
 
    void addOverlap(int readIndex, int overlapIndex,
                    int contactPos, int direction, int orient) {
        m_overlaps[readIndex].addOverlap(overlapIndex, contactPos, direction, orient);
    }

    void addOverlapSync(int readIndex, int overlapIndex,
                    int contactPos, int direction, int orient) {
        m_mutex.Lock();
        m_overlaps[readIndex].addOverlap(overlapIndex, contactPos, direction, orient);
        m_mutex.Unlock();
    }

    void resize(int sz)  { m_overlaps.resize(sz);     }
    int  getSize() const { return m_overlaps.isize(); }

    bool hasOverlap(int readIndex, int overlapIndex) const { 
        return m_overlaps[readIndex].hasLap(overlapIndex); 
    }
   
    void write(const string& overlapFile, int mode) const; 
    void writeAsc(const string& readOverlapFile) const; 
    void writeBin(const string& readOverlapFile) const; 
    void load(const string& readOverlapFile, const svec<int> & good, const ConsensReads& consReads, int mode); 
    void loadAsc(const string& readOverlapFile, const svec<int> & good, const ConsensReads& consReads); 
    void loadBin(const string& readOverlapFile, const svec<int> & good, const ConsensReads& consReads); 

    bool isChimeric(int i) const { return m_chimera[i]; }

    const svec<ReadOverlap>& getReadOverlaps(int readIndex, int isRight) const { 
      return m_overlaps[readIndex].getOverlaps(isRight);
    }
  
    void actionsAfterOverlapSet();

private:
    void postReadActions(const ConsensReads& consReads); 

    /** Internal helper class for writing/serialization */
    string getOverlapString(int index, int dir) const; 
    void addOverlapFromString(const string& strIn);

    svec<ReadInfo> m_overlaps;     /// All Overlap info , one item per read
    svec<int> m_chimera;           /// TODO comment
    ThreadMutex  m_mutex;          /// To use for locking while assigning overlaps
};
//======================================================

#endif // _READOVERLAP_H_
