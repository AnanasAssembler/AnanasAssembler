#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <stdlib.h>
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/util/mutil.h"
#include "ryggrad/src/base/Logger.h"
#include "ryggrad/src/base/StringUtil.h"
#include "Reads.h"
#include "Consensus.h"

//======================================================
string ConsensRead::toString() const{
    string temp;
    toString(temp);
    return temp;
}

void ConsensRead::toString(string& outString) const{
    stringstream sout;
    for(int i=0; i<m_readIdxs.isize(); i++) {
        sout << m_readIdxs[i] << "\t";
    }
    outString = sout.str();
}

string ConsensRead::genSeqName() {
    stringstream  ss;
    ss << ">CSeq:#" << getNumOfReads() << "_"<< m_readIdxs[0]; 
    return ss.str();
}
//======================================================


//======================================================
void ConsensReads::setConsensus(ConsensRead& cRead) {
    Consensus cons;
    bool setConsFlag = false;
    svec<int> readIdxs = cRead.getReads();
    int totLen = m_rawReads.getSize(readIdxs[0]);
    for(int i=1; i<readIdxs.isize(); i++) {  // Start from next index (1) as 0 has alread been looked at
      int tempSz = m_rawReads.getSize(readIdxs[i]);
      if(tempSz>totLen) { totLen = tempSz; }
    }
    if(cRead.isSingle()) { 
        cRead.setSize(totLen);
    } else {
        cons.resize(totLen); //Resize to length of the reads
        for(int i=0; i<readIdxs.isize(); i++) {
            int idx             = readIdxs[i];
            if(m_rawReads.hasSeq(idx)) { // Might be in memory eff mode where sequences don't exist
                DNAVector read    = m_rawReads[idx];
                for(int j=0; j<read.isize(); j++) {
                    cons.Add(j, read[j]); 
                }
                setConsFlag = true;
            } 
        }

        if(setConsFlag) {
            DNAVector out;
            out.resize(totLen);
            for (int i=0; i<totLen; i++) {
                out[i] = cons[i];
            }
            cRead.setSeq(out, m_isSingleStrand_p);
        } else {
            cRead.setSize(totLen);
        }
    }
}

void ConsensReads::addConsRead(const svec<int>& rawIdxs) { 
    //Log warning as empty set should not be passed here 
    if(rawIdxs.size()<1) { FILE_LOG(logWARNING) << "Empty set of ids cannot be used for consensus"; } 
    if(hasGroupInfo(rawIdxs)) { return; } //Group has already been created
    ConsensRead cRead;
    //Set group info for accessing for each raw read
    for(int i=0; i<rawIdxs.isize(); i++) {
        //Each index is only allowed once and only in one group
        if(!hasGroupInfo(rawIdxs[i])) {
            cRead.addRead(rawIdxs[i]);
            m_groupInfo[rawIdxs[i]]=m_consReads.size(); //Set the group index to the most current group addition
        }
    }
    setConsensus(cRead);
    cRead.setName();
    m_consReads.push_back(cRead); 
}

bool ConsensReads::hasGroupInfo(const svec<int>& idxs) const  {
    for(int i=0; i<idxs.isize(); i++) {
        if(hasGroupInfo(idxs[i])) { return true; }
    }
    return false;
}

void ConsensReads::write(const string& consReadFile, int mode) const {
    if(mode==0) { writeBin(consReadFile); } 
    else        { writeAsc(consReadFile); }
}

void ConsensReads::writeAsc(const string& consReadFile) const {
    ofstream sout;
    sout.open(consReadFile.c_str());
    for(int i=0; i<m_consReads.isize(); i++) {
        sout << m_consReads[i].toString() << endl;
    } 
    sout.close();
}

void ConsensReads::writeSeqsAsc(const string& fastaReadsFile) const {
    ofstream sout;
    sout.open(fastaReadsFile.c_str());
    for(int i=0; i<m_consReads.isize(); i++) {
        sout << (*this)[i].getName() << endl << (*this)[i].AsString() << endl;
    } 
    sout.close();
}

void ConsensReads::writeNamesAsc(const string& fastaReadsFile) const {
    ofstream sout;
    sout.open(fastaReadsFile.c_str());
    for(int i=0; i<m_consReads.isize(); i++) {
        sout << (*this)[i].getName() << " (" << i << ") - \t";
        svec<int> readIdxs = m_consReads[i].getReads();
        for(int j=0; j<readIdxs.isize(); j++) {
          sout << m_rawReads[readIdxs[j]].Name() << " (" << readIdxs[j] << ") \t";
        }
        sout << endl;
    } 
    sout.close();
}

void ConsensReads::writeBin(const string& consReadFile) const {
    CMWriteFileStream fs;
    fs.Open(consReadFile.c_str());
    for(int i=0; i<m_consReads.isize(); i++) {
        fs.Write((m_consReads[i].toString()+"\n").c_str());
    } 
    fs.Close();
}

void ConsensReads::load(const string& consReadFile, int mode) {
    if(mode==0) { loadBin(consReadFile); } 
    else        { loadAsc(consReadFile); }
}

void ConsensReads::loadAsc(const string& consReadFile) {
    ifstream sin;
    sin.open(consReadFile.c_str());
    string line;
    while(getline(sin, line)) {
        addConReadFromString(line);
    } 
    sin.close();
}
 
void ConsensReads::loadBin(const string& consReadFile) {
    CMReadFileStream fs;
    fs.Open(consReadFile.c_str());
    CMString line;
    while(!fs.IsEnd()) {
        fs.Read(line);
        addConReadFromString((const char*) line);
    }
    fs.Close();
}

void ConsensReads::addConReadFromString(const string& strIn){
    vector<string> tokens;
    Tokenize(strIn, tokens);
    svec<int> idxs;
    int tot = tokens.size();
    idxs.resize(tot);
    for(int i=0; i<tot; i++) {
        idxs[i] = (atoi(tokens[i].c_str()));
    }
    addConsRead(idxs);
}

//======================================================


//======================================================
void RawReads::write(const string& rawReadsFile) const {
    m_reads.Write(rawReadsFile);
}

void RawReads::load(const string& rawReadsFile) {
    m_reads.Read(rawReadsFile);
    if(!m_isSingleStrand_p) {
      m_readsRC = m_reads;
      m_readsRC.ReverseComplement();
    }
    // Set pair info
    m_pairInfo.resize(m_reads.isize());
    int cc = 0;
    for (int i=0; i<m_reads.isize(); i++) {
        char name[1024];
        strcpy(name, m_reads.Name(i).c_str());
        int n = strlen(name);
        bool b = false;
        if (name[n-2] != '/') continue;
        if (name[n-1] == '1') {
            name[n-1] = '2';
            b = true;
        } else {
            if (name[n-1] == '2') {
                name[n-1] = '1';
                b = true;
            }
        }
        if (!b)
            continue;
        int index = m_reads.NameIndex(name);
        if (index < 0) continue;
        m_pairInfo[i] = index;
        cc++;
    }
}

void RawReads::writePairSzInfo(const string& pairSzInfoFile) const {
    ofstream sout;
    sout.open(pairSzInfoFile.c_str());
    sout << m_pairInfo.isize() << endl;
    for(int i=0; i<m_pairInfo.isize(); i++) {
        sout << m_pairInfo[i] << "\t" << getSize(i) <<  endl; 
    }
    sout.close();
}

void RawReads::loadPairSzInfo(const string& pairSzInfoFile) {
    ifstream sin;
    sin.open(pairSzInfoFile.c_str());
    string line;
    getline(sin, line);
    m_pairInfo.resize(atoi(line.c_str()));
    m_sizeInfo.resize(atoi(line.c_str()));
    int idxCnt = 0;
    while(getline(sin, line)) {
        vector<string> tokens;
        Tokenize(line, tokens);
        m_pairInfo[idxCnt] = atoi(tokens[0].c_str());
        m_sizeInfo[idxCnt] = atoi(tokens[1].c_str());
        idxCnt++;
    }
    sin.close();
} 

//======================================================


//======================================================
void ReadGroups::group(int rIdx1, int rIdx2) {
    //If already grouped with each other do nothing
    if(isGrouped(rIdx1, rIdx2)) { return; }

    //If neither have groups, create new group
    if(!hasGroup(rIdx1) && !hasGroup(rIdx2)) {
        svec<int> newGroup;
        newGroup.resize(2);
        newGroup[0] = rIdx1;
        newGroup[1] = rIdx2;
        m_groups.push_back(newGroup);
        m_groupInfo[rIdx1] = m_groups.size()-1; 
        m_groupInfo[rIdx2] = m_groups.size()-1;
        return;
    }

    //If one has a group and the other doesn't, assign the group
    int g1 = whichGroup(rIdx1);
    if(!hasGroup(rIdx2)) {
        m_groupInfo[rIdx2] = g1; 
        m_groups[g1].push_back(rIdx2);
        return;
    }
    int g2 = whichGroup(rIdx2);
    if(!hasGroup(rIdx1)) {
        m_groupInfo[rIdx1] = g2; 
        m_groups[g2].push_back(rIdx1);
        return;
    }

    //If they both have groups but different to one another, merge
    if(g1<g2) {
        const svec<int>& group2 = getGroup(g2);
        for(int i=0; i<group2.isize(); i++) {
            m_groupInfo[group2[i]] = g1;
            m_groups[g1].push_back(group2[i]);
        }
        m_groups[g2].clear();
        return;
    } else {
        const svec<int>& group1 = getGroup(g1);
        for(int i=0; i<group1.isize(); i++) {
            m_groupInfo[group1[i]] = g2;
            m_groups[g2].push_back(group1[i]);
        }
        m_groups[g1].clear();
        return;
    }
}

void ReadGroups::assignSingleGroups() {
    for(int i=0; i<m_groupInfo.isize(); i++) {
        if(!hasGroup(i)) {
            svec<int> newGroup;
            newGroup.resize(1);
            newGroup[0] = i;
            m_groups.push_back(newGroup);
            m_groupInfo[i] = m_groups.size()-1;
        }
    }
}

void ReadGroups::write(const string& readGroupFile) const {
    ofstream sout;
    sout.open(readGroupFile.c_str());
    write(sout);
    sout.close();
}

void ReadGroups::write(ostream& sout) const {
    for(int i=0; i<m_groups.isize(); i++) {
        for(int j=0; j<m_groups[i].isize(); j++) {
            int idx = m_groups[i][j];
            if(m_tags.isize()>idx) { 
                sout << m_tags[idx] << "\t";
            } else {
                sout << m_groups[i][j] << "\t" ; 
            }
        }
        if(m_groups[i].isize()>0) { sout << endl; }
    }
}
//======================================================


//======================================================

void SubvecDNAVector::Read(const string & fileName, const svec<int>& idsToKeep) {
    FlatFileParser parser;
    parser.Open(fileName);
    int i;
    m_data.clear();

    int chunk = 20000; // reserve some space
    int bigChunk = 200000000; // 200 Megs?
    m_data.resize(chunk);

    DNAVector * pVec = NULL;
    DNAVector tmpVec;

    int j = 0;
    int regReadCount = 0;
    int totReadCount = 0;
    svec<int> idsToKeepSrt = idsToKeep;
    sort(idsToKeepSrt.begin(), idsToKeepSrt.end());
    svec<int>::iterator iter = idsToKeepSrt.begin();
    while (parser.ParseLine() && (iter!=idsToKeepSrt.end() || idsToKeepSrt.isize()==0)) {
        if (parser.GetItemCount() == 0) { continue; }
        const char * p = parser.AsString(0).c_str();
        if (regReadCount == 0 && p[0] == '@') { // It's a fastq file!!!
            m_data.resize(0);
            ReadQ(fileName);
            return;
        }
        if (p[0] == '>') {
            if(idsToKeepSrt.isize()!=0 && totReadCount!=(*iter)) {
                totReadCount++;
                parser.ParseLine(); //Skip to next FASTA entry
                continue;
            }

            string tmpName = parser.AsString(0);
            for (int x=1; x<parser.GetItemCount(); x++) {
                tmpName += "_";
                tmpName += parser.AsString(x);
            }
            if (regReadCount >= m_data.isize())
                m_data.resize(regReadCount + chunk);

            pVec = &(*this)[regReadCount];
            pVec->SetName(tmpName);

            if(parser.ParseLine() && parser.GetItemCount()!=0) {
                const char * p = parser.AsString(0).c_str();
                int n = strlen(p);
                for (i=0; i<n; i++) {
                    if (j >= tmpVec.isize())
                        tmpVec.resize(j + bigChunk);
                    tmpVec[j] = p[i];
                    j++;
                }
                if (pVec != NULL) {
                    pVec->SetToSubOf(tmpVec, 0, j);
                    pVec->ToUpper();
                    j = 0;
                }
            }

            if(idsToKeepSrt.isize()!=0) { // idsToKeep is none empty
                linkIndex(totReadCount, regReadCount);
                iter++;
            }

            regReadCount++;
            totReadCount++;
        }
    }

    m_data.resize(regReadCount);
    setupMap();
}
