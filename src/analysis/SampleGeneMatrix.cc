#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/base/StringUtil.h"
#include "src/analysis/SampleGeneMatrix.h"

void SampleGeneMatrix::load(const string& fileName) {
    FlatFileParser parser;
    parser.Open(fileName);

    while (parser.ParseLine()) {
        if (parser.GetItemCount() != 3)
            continue;
        string sampleId = parser.AsString(0);
        int underscorePos1 = sampleId.find("_", 0);
        int underscorePos2 = sampleId.find("_", underscorePos1+1);
        if(underscorePos1>0 && underscorePos2>0) {
            sampleId           = sampleId.substr(underscorePos1+1, underscorePos2-underscorePos1-1);
        }
   
        int readCnt    = parser.AsInt(1);

        string geneId  = parser.AsString(2);

       // Add the matrix value
       if(m_matrix.count(geneId)==1 && m_matrix[geneId].count(sampleId)==1) {
           m_matrix[geneId][sampleId] += readCnt; 
       } else {
           m_matrix[geneId][sampleId] = readCnt;
           m_sampleIds[sampleId]      = true;   /// Register this sample Id in case it has not been registered before
           m_geneIds[geneId]          = true;   /// Regiester this gene Id in case it has not been registered before
       }
    }
}

void SampleGeneMatrix::write(const string& fileName) {
    FILE * pOut = fopen(fileName.c_str(), "w");
    if (pOut == NULL) {
        cout << "ERROR: Could not open file " << fileName << " for writing." << endl;
    }

    int sampleCnt = 0;
    fprintf(pOut, "GeneId\t");
    for(auto const& ent1 : m_sampleIds) {
        fprintf(pOut, "%s\t", ent1.first.c_str());
        sampleCnt++;
    }
    fprintf(pOut, "\n");

    int geneCnt   = 0;
    for(auto const& ent1 : m_geneIds) {
        fprintf(pOut, "%s\t", ent1.first.c_str());
        geneCnt++;
        for(auto const& ent2 : m_sampleIds) {
            fprintf(pOut, "%d\t", readCount(ent1.first, ent2.first));
        }
        fprintf(pOut, "\n");
    }
    cout << "Matrix contains: " << sampleCnt << " samples and " << geneCnt << " genes." << endl;

}

int SampleGeneMatrix::readCount(const string& geneId, const string& sampleId) {
    if(m_matrix.count(geneId)==1 && m_matrix[geneId].count(sampleId)==1) {
        return m_matrix[geneId][sampleId];
    } else {
        return 0;
    }
} 
