#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/StringUtil.h"
#include "extern/logger/log.h"
#include "base/RandomStuff.h"
#include "src/Ananas/ContigClusterUnit.h"

//======================================================
void ContigClusterUnit::clusterContigs(int numOfThreads) {
    SubReads<RawReads>  subreads(m_rawReads, m_params, false); // Object handling the subreads 
    int totSize   = m_rawReads.getNumOfReads();

    cout << "Clustering Contigs..." << endl;

    AllReadOverlaps overlaps; 
    overlaps.resize(totSize); // Make sure enough memory is declared 

    ThreadQueueVec threadQueue(totSize);
    ThreadHandler th;
    if(numOfThreads>totSize) { numOfThreads = totSize; }
    for (int i=0; i<numOfThreads; i++) {
        char tmp[256];
        sprintf(tmp, "%d", i);
        string init = "init_";
        init += tmp;
        th.AddThread(new FindOverlapsSingleThread< SubReads<RawReads> >(threadQueue, subreads, overlaps, 1, i));    
        th.Feed(i, init);
    }

    while (!th.AllDone()) {
        usleep(10000);
    }
    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed finding Overlaps." << endl;

    ReadGroups nIdentGroupInfo(m_rawReads.getNumOfReads());   // Make sure enough memory is declared 
    for(int i=0; i<totSize; i++) {
      for(int dir=-1; dir<2; dir+=2) {
        const svec<ReadOverlap>& currOverlaps  = overlaps.getReadOverlaps(i, dir);
        for(int j=0; j<currOverlaps.isize(); j++) {
          if(!nIdentGroupInfo.isGrouped(i, currOverlaps[j].getOverlapIndex())) { 
            nIdentGroupInfo.group(i, currOverlaps[j].getOverlapIndex()); 
          }
        }
      }
    }         
    nIdentGroupInfo.assignSingleGroups();
    nIdentGroupInfo.setTags<RawReads>(m_rawReads);
    nIdentGroupInfo.write("temp1.tmp");
    writeContigClusters("temp.tmp", overlaps);
}

void ContigClusterUnit::writeContigClusters(const string& clusterFile, const AllReadOverlaps& overlaps) const {
    ofstream sout;
    sout.open(clusterFile.c_str());
    int totNumOfReads = overlaps.getSize();
    sout << totNumOfReads << endl;
    for(int index=0; index<totNumOfReads; index++) {
        for(int dir=-1; dir<2; dir+=2) {
           const svec<ReadOverlap>& currOverlaps  = overlaps.getReadOverlaps(index, dir);
           stringstream ss;
            for(int j=0; j<currOverlaps.isize(); j++) {
               ss << m_rawReads[index].Name()  << "\t" << m_rawReads[currOverlaps[j].getOverlapIndex()].Name() << "\t" 
                  << m_rawReads[index].size()  << "\t" << m_rawReads[currOverlaps[j].getOverlapIndex()].size() << "\t"
                  << currOverlaps[j].getContactPos() << "\t" << currOverlaps[j].getScore() << "\t" << "\t" << (dir==1?">":"<") << "\t"
                  << (currOverlaps[j].getOrient()==1?"+":"-") << endl;
            } 
            sout << ss.str();
        }
    }
    sout.close();
}
//======================================================
