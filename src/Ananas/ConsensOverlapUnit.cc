#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/StringUtil.h"
#include "extern/logger/log.h"
#include "base/RandomStuff.h"
#include "src/Ananas/ConsensOverlapUnit.h"

//======================================================

void ConsensOverlapUnit::findOverlaps(int numOfThreads, string groupedReadInfo, double identThresh) {
    if(groupedReadInfo=="") {
        createConsensReads(identThresh); 
    } else {
        m_consReads.load(groupedReadInfo, 1); //TODO allow for binary file
    }

    SubReads<ConsensReads>  subreads(m_consReads, m_params, false); // Object handling the subreads 
    int totSize   = m_consReads.getNumOfReads();

    FILE_LOG(logDEBUG1) << "Finding Overlaps";
    cout << "Finding All Overlaps..." << endl;
    int progCount = 0;

    m_overlaps.resize(totSize); // Make sure enough memory is declared 

    ThreadHandler th;
    for (int i=0; i<numOfThreads; i++) {
        char tmp[256];
        sprintf(tmp, "%d", i);
        string init = "init_";
        init += tmp;
        int from = i*(totSize/numOfThreads);
        int to   = (i+1)*(totSize/numOfThreads);
        th.AddThread(new FindOverlapsThread(subreads, m_overlaps,from, to, i));    
        th.Feed(i, init);
    }

    while (!th.AllDone()) {
        usleep(10000);
    }
    cout << "\r===================== " << "100.0% " << flush; 
    cout << "Completed finding Overlaps." << endl;
}

void ConsensOverlapUnit::createConsensReads(float minMatchScore_p) {
    //int stepSize = 1/(1-(minMatchScore_p-0.01));
    //AssemblyParams params(m_params);
    //params.setSubreadStep(stepSize);
    //SubReads<RawReads>  subreads(m_rawReads, params, true); // Object handling the subreads 
    SubReads<RawReads>  subreads(m_rawReads, m_params, true); // Object handling the subreads 
    ReadGroups nIdentGroupInfo(m_rawReads.getNumOfReads());   // Make sure enough memory is declared 
    FILE_LOG(logDEBUG1) << "Finding Consensus Reads based on Near Ident Groupings";
    cout << "Finding Near Ident Groupings..." << endl;
    int progCount = 0;
    int totSize   = subreads.getSize();
    int inc = totSize/10000;
    if (inc == 0)
        inc = 1; //Let's not crash if we have too few reads 
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

//======================================================
bool FindOverlapsThread::OnDo(const string & msg) {
    int progCount = 0;
    int totSize   = (m_toIdx - m_fromIdx);
    int inc       = totSize/100;
    if (inc < 1)
        inc = 1;

    for(int i=m_fromIdx; i<m_toIdx; i++) { 
        m_subreads.findOverlaps(i, m_overlaps); 
        progCount++;
        if (progCount % inc == 0) 
            cout << "\r===================== " << 100.0*progCount/totSize 
                 << "%  " << flush; 
    }
    return true;
}


//======================================================
