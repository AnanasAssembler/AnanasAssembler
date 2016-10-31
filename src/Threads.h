#ifndef _THREADS_H_
#define _THREADS_H_

#include "ryggrad/src/base/ThreadHandler.h"
#include "ThreadQueueVec.h"
#include "Reads.h"
#include "ReadOverlap.h"

template <class SuffixType>
class FindOverlapsThread : public IOneThread
{
public:
  FindOverlapsThread(const SuffixType& sr,
                           AllReadOverlaps& ol, int mode,
                           int from, int to, 
                           int rThresh, int tn): m_subreads(sr), m_overlaps(ol), m_mode(mode),
                                                 m_fromIdx(from), m_toIdx(to), m_runThresh(rThresh), 
                                                 m_threadIdx(tn) {}

protected:
  virtual bool OnDie() { return true; }
  virtual bool OnDo(const string & msg); 
  virtual bool OnInitialize(const string & msg) {    return true; }

  const SuffixType&  m_subreads;              /// Subreads to find overlaps
  AllReadOverlaps& m_overlaps;                /// All overlaps among consensus reads
  int m_mode;                                 /// Specify if none-extending overlaps should also be returned (1 and if not 0) 
  int m_fromIdx;
  int m_toIdx;
  int m_runThresh;                            /// Threshold where if number of existing overlaps exceed this, will skip finding more overlaps
  int m_threadIdx;
};

template <class SuffixType>
class FindOverlapsSingleThread : public FindOverlapsThread<SuffixType>
{
public:
  FindOverlapsSingleThread(ThreadQueueVec& tQueue, const SuffixType& sr,
                     AllReadOverlaps& ol, int mode,
                     int rThresh, int tn): FindOverlapsThread<SuffixType>(sr, ol, mode, 0, 0, rThresh, tn),
                                           m_threadQueue(tQueue) {}
protected:
  virtual bool OnDo(const string & msg); 

private:
  ThreadQueueVec& m_threadQueue;              /// Queue for picking up jobs from the queue
};

//======================================================
template <class SuffixType>
bool FindOverlapsThread<SuffixType>::OnDo(const string & msg) {
  int progCount = 0;
  int totSize   = (m_toIdx - m_fromIdx);
  int inc       = totSize/1000;
  if (inc < 1)
      inc = 1;

  for(int i=m_fromIdx; i<m_toIdx; i++) { 
    if(m_overlaps[i].getNumLaps()<m_runThresh) {
       m_subreads.findOverlaps(i, m_overlaps, m_mode, m_runThresh); 
     }
     progCount++;
     if (progCount % inc == 0) 
       cout << "\r===================== " << 100.0*progCount/totSize 
            << "%  " << flush; 
  }
  return true;
}


//======================================================

//======================================================
template <class SuffixType>
bool FindOverlapsSingleThread<SuffixType>::OnDo(const string & msg) {
  int totSize = this->m_threadQueue.getSize();
  int inc     = totSize/1000;
  if (inc < 1)
    inc = 1;
  int currIdx = this->m_threadQueue.getNext();
   while(currIdx>=0) {
     if(this->m_overlaps[currIdx].getNumLaps()<this->m_runThresh) {
       FILE_LOG(logDEBUG3) << "Finding overlaps for read idx: " << currIdx; 
       this->m_subreads.findOverlaps(currIdx, this->m_overlaps, this->m_mode, this->m_runThresh); 
     }
     if (currIdx  % inc == 0) 
            cout << "\r===================== " << 100.0*currIdx/totSize 
                 << "%  " << flush; 
     currIdx = this->m_threadQueue.getNext();
   }
   return true;
}

//======================================================


#endif // _THREADS_H_
