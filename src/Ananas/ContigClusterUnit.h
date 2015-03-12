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
              :m_params(params), m_rawReads(inputFile, 0, m_params.isSingleStrand()),
               m_overlaps(), m_clusters(m_rawReads.getNumOfReads()) {}

    int getNumOfRawReads() const              { return m_rawReads.getNumOfReads();  } 

    void clusterContigs(int numOfThreads);
    void writeContigClusters(const string& clusterFile) const;
    void writeContigPairs(const string& clusterFile) const;

private:

    AssemblyParams         m_params;          /// Object containing the various parameters required for assembly
    RawReads               m_rawReads;        /// A list of reads from which subreads where constructed
    AllReadOverlaps        m_overlaps;        /// Pairwise relations between contigs
    ReadGroups             m_clusters;        /// The contigs grouped into clusters

};

template <class SuffixType>
class ContigClusterThread: public FindOverlapsSingleThread<SuffixType>
{
};



#endif // _CONTIGCLUSTERUNIT_H_
