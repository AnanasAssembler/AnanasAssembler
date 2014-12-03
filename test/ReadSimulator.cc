#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "test/ReadSimulator.h"
 

//======================================================

void ReadSimulator::generatePEReads(const vecDNAVector& inputSeqs, int interval, 
                     double sub, double indel) {
   int i, j;
    int len = 100;
    int sep = 100;
    
    int plusminus = 60;
    
    int totCnt = 0;
    for (j=0; j<(int)inputSeqs.size(); j++) {
        const DNAVector & d = inputSeqs[j];
        for (i=0; i<d.isize()-2*len-sep; i+= interval) {
            int pos1 = i + RandomInt(plusminus) - plusminus/2;
            if (pos1 < 0)
	        pos1 = 0;
	    int pos2 = pos1 + len + sep;
	    if (pos2 + len >= d.isize())
	        pos2 = d.isize() - len;
	    DNAVector left, right;
	    left.SetToSubOf(d, pos1, len);   
	    right.SetToSubOf(d, pos2, len);
	    DNAVector left_err, right_err;
	    mutate(left, sub, indel, left_err);
	    mutate(right, sub, indel, right_err);
	    right_err.ReverseComplement();
            char strLeft[1024], strRight[1024];
            // Index added after underscore is for identifying same base name for a pair
	    sprintf(strLeft, "%s#%d/1\\%d\\%d\\%d", inputSeqs.Name(j).c_str(), totCnt, pos1, len, totCnt);
	    sprintf(strRight, "%s#%d/2\\%d\\%d\\%d", inputSeqs.Name(j).c_str(), totCnt, pos2, len, totCnt+1);
	    totCnt+=2;
            left_err.setName(strLeft);
            right_err.setName(strRight);
            addRead(left_err, pos1);
            addRead(right_err, pos2);
        }
    }
    sortReads();
}

void ReadSimulator::findAllOverlaps(AllReadOverlaps& allOverlaps, int minOverlap) {
     FILE_LOG(logINFO) << "Finding overlaps";
    for (int i=0; i<(int)m_reads.size(); i++) {
        FILE_LOG(logDEBUG1) << "Index i: " << i;
        string origName_i;
        int offset_i, length_i, index_i;
        bool isLeft_i;
        getInfo(m_reads[i], origName_i, offset_i, length_i, index_i, isLeft_i);
        for(int j=i+1; j<(int)m_reads.size();j++) {
            FILE_LOG(logDEBUG2) << "Index j: " << j;
            string origName_j;
            int offset_j, length_j, index_j;
            bool isLeft_j;
            getInfo(m_reads[j], origName_j, offset_j, length_j, index_j, isLeft_j);
            int contactPos  = offset_j - offset_i;
            FILE_LOG(logDEBUG2) << "Offset_j: " << offset_j << " Offset_i: " << offset_i;
            if(origName_j!=origName_i || contactPos>(length_i-minOverlap)) { break; }
            int overlapDir = 1;
            allOverlaps.addOverlap(i, j, contactPos, overlapDir, true);
        }
        for(int k=i-1; k>=0; k--) {
            FILE_LOG(logDEBUG2) << "Index k: " << k;
            string origName_k;
            int offset_k, length_k, index_k;
            bool isLeft_k;
            getInfo(m_reads[k], origName_k, offset_k, length_k, index_k, isLeft_k);
            int contactPos  = offset_i - offset_k;
            FILE_LOG(logDEBUG2) << "Offset_i: " << offset_i << " Offset_k: " << offset_k;
            if(origName_k!=origName_i || contactPos>(length_i-minOverlap)) { break; }
            int overlapDir = -1;
            allOverlaps.addOverlap(i, k, contactPos, overlapDir, true);
        }
    }
}

void ReadSimulator::writeReads(ostream& sout) {
    vector<char> separators;
    separators.push_back('\\');
    for(int i=0; i<(int)m_reads.size(); i++) {
        string name = m_reads[i].getName();
        vector<string> tokens;
        Tokenize(name, separators, tokens);
        string printName = tokens[0]; //Removing extra tags to standardize for pair association
        sout << printName << endl
             << m_reads[i].AsString()   << endl;
    }  
}

void ReadSimulator::mutate(const DNAVector & in, double sub, double indel, DNAVector & out) {
    int i, j;
    string n;
    for (i=0; i<in.isize(); i++) {
        if (RandomFloat(1.) < indel) {
            if (RandomFloat(1.) < 0.5) {
	        // Deletion
	        i += RandomInt(3) + 1;
            } else {
	        // Insertion
	        int plus = RandomInt(3) + 1;
	        for (j=0; j<plus; j++) {
	            n += NucLetter(RandomInt(4));
	        }
            }
        }
    
        if (RandomFloat(1.) < sub) {
            int s = NucIndex(in[i]);
            s += RandomInt(3) + 1;
            s = s % 4;
            n += NucLetter(s);
        } else {
            n += in[i];
        }
    }
    out.SetFromBases(n);
}

void ReadSimulator::getInfo(const DNAVector& d, string& origName, int& offset, int& length, int& index, bool& isLeft) {
    vector<char> separators(3);
    separators[0] =  '\\';
    separators[1] = '/';
    separators[2] = '#';
    string name = d.getName();
    vector<string> tokens;
    Tokenize(name, separators, tokens);
    origName = tokens[0];
    isLeft   = ((atoi(tokens[2].c_str())==1)?true:false);
    offset   = atoi(tokens[3].c_str());
    length   = atoi(tokens[4].c_str());
    index    = atoi(tokens[5].c_str());
}

//Compare reads for sorting based on offset 
bool ReadSimulator::cmp(const DNAVector& d1, const DNAVector& d2) { 
    string origName_1, origName_2;
    int offset_1, length_1, index_1, offset_2, length_2, index_2;
    bool isLeft_1, isLeft_2;
    getInfo(d1, origName_1, offset_1, length_1, index_1, isLeft_1);
    getInfo(d2, origName_2, offset_2, length_2, index_2, isLeft_2);
    if(origName_1==origName_2) {
        return (offset_1<offset_2);
    } else {
        return (origName_1<origName_2);
    }
}

//======================================================
