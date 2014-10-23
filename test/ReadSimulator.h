#ifndef _READSIMULATOR_H_
#define _READSIMULATOR_H_

#include <string>
#include "extern/logger/log.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"
#include "base/StringUtil.h"
#include "base/RandomStuff.h"
#include "src/Ananas/ReadOverlap.h"

class ReadSimulator {
public:
    ReadSimulator():m_reads() {}

    /** Generated Paired end reads with given error rates and length */
    void generatePEReads(const vecDNAVector& inputSeqs, int interval, 
                       double sub, double indel); 

    int getReadCnt() { return m_reads.size(); }

    void addRead(const DNAVector& read, int offset) { m_reads.push_back(read); }

    void findAllOverlaps(AllReadOverlaps& allOverlaps, int minOverlap);

    void writeReads(ostream& sout); 

private:
    void mutate(const DNAVector & in, double sub, double indel, DNAVector & out); 
    static void getInfo(const DNAVector& d, string& origName, int& offset, int& length, int& index, bool& isLeft); 
    //Compare reads for sorting based on offset 
    static bool cmp(const DNAVector& d1, const DNAVector& d2);  
    void sortReads() { sort(m_reads.begin(), m_reads.end(), cmp); } 
 
    vector<DNAVector> m_reads;     /// The sequences for the simulated reads
};

#endif // _READSIMULATOR_H_


