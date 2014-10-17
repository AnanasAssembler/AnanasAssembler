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

    m_overlaps.resize(totSize); // Make sure enough memory is declared 

    ThreadQueueVec threadQueue(totSize);
    ThreadHandler th;
    if(numOfThreads>totSize) { numOfThreads = totSize; }
    for (int i=0; i<numOfThreads; i++) {
        char tmp[256];
        sprintf(tmp, "%d", i);
        string init = "init_";
        init += tmp;
        th.AddThread(new FindOverlapsSingleThread< SubReads<RawReads> >(threadQueue, subreads, m_overlaps, 1, i));    
        th.Feed(i, init);
    }

    while (!th.AllDone()) {
        usleep(10000);
    }
    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed finding Overlaps." << endl;

    for(int i=0; i<totSize; i++) {
      for(int dir=-1; dir<2; dir+=2) {
        const svec<ReadOverlap>& currOverlaps  = m_overlaps.getReadOverlaps(i, dir);
        for(int j=0; j<currOverlaps.isize(); j++) {
          if(!m_clusters.isGrouped(i, currOverlaps[j].getOverlapIndex())) { 
            m_clusters.group(i, currOverlaps[j].getOverlapIndex()); 
          }
        }
      }
    }         
    m_clusters.assignSingleGroups();
    m_clusters.setTags<RawReads>(m_rawReads);
}

void ContigClusterUnit::writeContigClusters(const string& clusterFile) const {
    m_clusters.write(clusterFile);
}

void ContigClusterUnit::writeContigPairs(const string& clusterFile) const {
    ofstream sout;
    sout.open(clusterFile.c_str());
    int totNumOfReads = m_overlaps.getSize();
    sout << totNumOfReads << endl;
    for(int index=0; index<totNumOfReads; index++) {
        for(int dir=-1; dir<2; dir+=2) {
           const svec<ReadOverlap>& currOverlaps  = m_overlaps.getReadOverlaps(index, dir);
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
