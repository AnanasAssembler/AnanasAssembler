#ifndef _CONTIGCLUSTERUNIT_H_
#define _CONTIGCLUSTERUNIT_H_

#include "src/Ananas/SubReads.h"
#include "src/Ananas/AssemblyParams.h"
#include "src/Ananas/Reads.h"
#include "src/Ananas/ReadOverlap.h"
#include "src/Ananas/Threads.h"


//======================================================
class ContigClusterUnit
{
public:

    // Basic Constructor used for finding clusters 
    ContigClusterUnit(const AssemblyParams& params, const string& inputFile)
              :m_params(params), m_rawReads(inputFile, 0) {}

    int getNumOfRawReads() const              { return m_rawReads.getNumOfReads();  } 

    void clusterContigs(int numOfThreads);
    void writeContigClusters(const string& clusterFile, const AllReadOverlaps& overlaps) const;

private:

    AssemblyParams         m_params;          /// Object containing the various parameters required for assembly
    RawReads               m_rawReads;        /// A list of reads from which subreads where constructed
};

template <class SuffixType>
class ContigClusterThread: public FindOverlapsSingleThread<SuffixType>
{
};



#endif // _CONTIGCLUSTERUNIT_H_
