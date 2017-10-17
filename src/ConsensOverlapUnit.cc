#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/StringUtil.h"
#include "ryggrad/src/base/Logger.h"
#include "ryggrad/src/base/RandomStuff.h"
#include "ConsensOverlapUnit.h"
#include "ThreadQueueVec.h"

//======================================================

void ConsensOverlapUnit::findOverlaps(int numOfThreads, int mode, double identThresh, int maxOverlapPerRead, string groupedReadInfo) {
    if(groupedReadInfo=="") {
        createConsensReads(identThresh); 
    } else {
        m_consReads.load(groupedReadInfo, 1); //TODO allow for binary file
    }

    SubReads<ConsensReads>  subreads(m_consReads, m_params, false); // Object handling the subreads 
    int totSize   = m_consReads.getNumOfReads();
    if(numOfThreads>totSize) { numOfThreads = totSize; }

    FILE_LOG(logDEBUG1) << "Finding Overlaps";
    cout << "Finding All Overlaps..." << endl;

    m_overlaps.resize(totSize);          // Make sure enough memory is declared 
    ThreadQueueVec threadQueue(totSize); // Use for queueing instances threads should handle
    ThreadHandler th;

        for (int i=0; i<numOfThreads; i++) {
            char tmp[256];
            sprintf(tmp, "%d", i);
            string init = "init_";
            init += tmp;
            // set the limiting threshold for number of overlaps to 0 so that it is set to readsize*2 in the overlap finder 
            th.AddThread(new FindOverlapsSingleThread< SubReads<ConsensReads> >(threadQueue, subreads, m_overlaps, mode, maxOverlapPerRead, i));
            th.Feed(i, init);
        }
        while (!th.AllDone()) {
            usleep(10000);
        }

    m_overlaps.actionsAfterOverlapSet(); //Sorts Overlaps

    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed finding Overlaps." << endl;
}

void ConsensOverlapUnit::createConsensReads(float minMatchScore_p) {
    SubReads<RawReads>  subreads(m_rawReads, m_params, true); // Object handling the subreads 
    ReadGroups nIdentGroupInfo(m_rawReads.getNumOfReads());   // Make sure enough memory is declared 
    FILE_LOG(logDEBUG1) << "Finding Consensus Reads based on Near Ident Groupings";
    cout << "Finding Near Ident Groupings..." << endl;
    int progCount = 0;
    int totSize   = subreads.getSize();
    int inc = totSize/10000;
    if (inc == 0)
        inc = 1; //Let's not crash if we have too few reads 
    if(minMatchScore_p  > 0 && minMatchScore_p <= 1) {
      for(int i=0; i<totSize-1; i++) { 
          SubRead a = subreads.getByIndex(i); 
          SubRead b = subreads.getByIndex(i+1);
          if(!nIdentGroupInfo.isGrouped(a.getIndex(), b.getIndex()) 
             && a.getOffset() == b.getOffset()) {
              DNAVector aD, bD;
              aD.SetFromBases(subreads.getSeq(a, 0, m_rawReads.getSize(a.getIndex())-1));
              bD.SetFromBases(subreads.getSeq(b, 0, m_rawReads.getSize(b.getIndex())-1));
              float matchScore = aD.FindIdent(bD);  // Direct match identity without alignment 
              if (matchScore>=minMatchScore_p) { 
                  nIdentGroupInfo.group(a.getIndex(), b.getIndex()); 
              }
          } 
          progCount++;
          if (progCount % inc == 0) 
              cout << "\r===================== " << 100.0*progCount/totSize << "% " << flush; 
      }
    }
    nIdentGroupInfo.assignSingleGroups();

    FILE_LOG(logDEBUG1) << "Finished finding consensus reads";
    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed." << endl;
    m_consReads.reserveMem(nIdentGroupInfo.getNumOfGroups());
    for (int i=0; i<nIdentGroupInfo.getNumOfGroups(); i++) {
        const svec<int >& idxs = nIdentGroupInfo.getGroup(i);
        if(idxs.isize()==0) { continue; }
        m_consReads.addConsRead(idxs);
    }
    findPartners(); //Set the partners for the consensus reads
} 

//TODO determine what should be done
void ConsensOverlapUnit::Prune(const svec<int> & good) {
    for (int i=0; i<good.isize(); i++) {
        if (good[i] == 0)
            m_consReads.clearRead(i);
    }  
}

void ConsensOverlapUnit::findPartners() {
    m_partners.resize(getNumOfConsReads());
    for(int i=0; i<getNumOfConsReads(); i++) {
        const svec<int>& consMemIds = m_consReads.getConsMembers(i);
        for(int j=0; j<consMemIds.isize(); j++) {
            int consIdx = m_consReads.whichGroup(m_rawReads.getPairId(consMemIds[j]));
            if(consIdx!=-1) { m_partners.addPartner(i, consIdx); }
        }
    }
}


//======================================================


