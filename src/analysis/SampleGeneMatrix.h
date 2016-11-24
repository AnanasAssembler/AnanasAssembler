#ifndef SAMPLEGENEMATRIX_H
#define SAMPLEGENEMATRIX_H

#include <map>
#include "ryggrad/src/base/SVector.h"

class SampleGeneMatrix
{
public:
    SampleGeneMatrix(): m_matrix(), m_sampleIds(), m_geneIds() {}
    
    int readCount(const string& geneId, const string& sampleId); 

    void load(const string& fileName);

    void write(const string& fileName);


 private:
    map < string, map <string, int> > m_matrix;      /// Matrix mapping geneId onto total count in every sample
    map<string, bool> m_sampleIds;                   /// A uniqe list of all sample Ids
    map<string, bool> m_geneIds;                     /// A uniqe list of all gene Ids
};

#endif //SAMPLEGENEMATRIX_H

