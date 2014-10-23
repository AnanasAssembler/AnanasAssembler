#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <map>
#include <sstream>
#include <stdint.h>
#include "base/SVector.h"
#include "util/mutil.h"
#include "extern/logger/log.h"
#include "src/DNAVector.h"
#include "src/Ananas/ReadOverlap.h"

//======================================================
void AllReadOverlaps::write(const string& overlapFile, int mode) const {
    if(mode==0) { writeBin(overlapFile); } 
    else        { writeAsc(overlapFile); }
}

void AllReadOverlaps::writeAsc(const string& readOverlapFile) const {
    ofstream sout;
    sout.open(readOverlapFile.c_str());
    int totNumOfReads = m_overlaps.isize();
    sout << totNumOfReads << endl;
    for(int i=0; i<totNumOfReads; i++) {
        sout << getOverlapString(i, 1);
        sout << getOverlapString(i, -1);
    }
    sout.close();
}

void AllReadOverlaps::writeBin(const string& readOverlapFile) const {
    CMWriteFileStream fs;
    fs.Open(readOverlapFile.c_str());
    int totNumOfReads = m_overlaps.isize();
    fs.Write(totNumOfReads); // To use for allocate memory when reading
    for(int i=0; i<totNumOfReads; i++) {
        // Loop for false (0) and true (1)
        for(int dir=0; dir<2; dir++) {
            const svec<ReadOverlap>& overlaps  = getReadOverlaps(i, dir);
            for(int j=0; j<overlaps.isize(); j++) {
                fs.Write(i);
                fs.Write(overlaps[j].getOverlapIndex());
                fs.Write(overlaps[j].getContactPos());
                fs.Write(overlaps[j].getScore());
                fs.Write(overlaps[j].getDirection());
                fs.Write(overlaps[j].getOrient());
            }
        }
    }
    fs.Close();
}

string AllReadOverlaps::getOverlapString(int index, int dir) const {
    stringstream ss;
    const svec<ReadOverlap>& overlaps  = getReadOverlaps(index, dir);
    for(int j=0; j<overlaps.isize(); j++) {
        ss << index  << "\t" << overlaps[j].getOverlapIndex() << "\t" 
           << overlaps[j].getContactPos() << "\t" << overlaps[j].getScore() << "\t" << "\t" << (dir==1?">":"<") << "\t"
           << (overlaps[j].getOrient()==1?"+":"-") << endl;
    }
    return ss.str();
}

void AllReadOverlaps::load(const string& readOverlapFile, int mode) {
    if(mode==0) { loadBin(readOverlapFile); } 
    else        { loadAsc(readOverlapFile); }
}


//TODO this is temporary - remove
void AllReadOverlaps::loadAsc(const string& readOverlapFile, const svec<int> & good, const ConsensReads& consReads) {
    ifstream sIn;
    sIn.open(readOverlapFile.c_str());
    string line;
    getline(sIn, line);
    int totNumOfReads = atoi(line.c_str());
    resize(totNumOfReads);
    while(getline(sIn, line)) {
        CMTokenizer tokenizer;
        tokenizer.AddDelimiter("\t");
        CMPtrStringList tokens;
        tokenizer.Tokenize(tokens, line.c_str());
        if(tokens.length()<6) { 
            FILE_LOG(logERROR) << "Wrong overlap file format - six columns required"; 
            return;
        }
        string dir    = (const char*) *tokens[4];
        string strand = (const char*) *tokens[5];
        if(good.isize() > 0 && (!good[atoi((const char*)*tokens[0])] && !good[atoi((const char*)*tokens[1])])) {
            continue;
        }
        addOverlap(atoi((const char*)*tokens[0]), 
                   atoi((const char*)*tokens[1]), atoi((const char*)*tokens[2]), 
                   atof((const char*)*tokens[3]), (dir==">"?1:-1), (strand=="+"?1:-1));

    }
    sIn.close();
    postReadActions(consReads); 
}

void AllReadOverlaps::loadAsc(const string& readOverlapFile) {
    ifstream sIn;
    sIn.open(readOverlapFile.c_str());
    string line;
    getline(sIn, line);
    int totNumOfReads = atoi(line.c_str());
    resize(totNumOfReads);
    while(getline(sIn, line)) {
        addOverlapFromString(line);
    } 
    sIn.close();
    //  postReadActions();  TODO
}

void AllReadOverlaps::loadBin(const string& readOverlapFile) {
    CMReadFileStream fs;
    fs.Open(readOverlapFile.c_str());
    int totNumOfReads;
    fs.Read(totNumOfReads);  //To use for allocate memory when reading
    resize(totNumOfReads);
    while(!fs.IsEnd()) {
        int dir, orient;
        int readIdx, overlapIdx;
        int contactPos;
        float score; 
        fs.Read(dir);
        fs.Read(readIdx);
        fs.Read(overlapIdx);
        fs.Read(contactPos);
        fs.Read(score);
        fs.Read(orient);
     
        addOverlap(readIdx, overlapIdx, contactPos, score, dir, orient);
    }
    fs.Close();
}

void AllReadOverlaps::addOverlapFromString(const string& strIn){
    CMTokenizer tokenizer;
    tokenizer.AddDelimiter("\t");
    CMPtrStringList tokens;
    tokenizer.Tokenize(tokens, strIn.c_str());
    if(tokens.length()<6) { 
        FILE_LOG(logERROR) << "Wrong overlap file format - six columns required"; 
        return;
    }
    string dir    = (const char*) *tokens[4];
    string strand = (const char*) *tokens[5];
    addOverlap(atoi((const char*)*tokens[0]), 
               atoi((const char*)*tokens[1]), atoi((const char*)*tokens[2]), 
               atof((const char*)*tokens[3]), (dir==">"?1:-1), (strand=="+"?1:-1));
}
 
void AllReadOverlaps::actionsAfterOverlapSet() {
    for (int i=0; i<m_overlaps.isize(); i++) {
        m_overlaps[i].sortLaps();
    }
}

void AllReadOverlaps::postReadActions(const ConsensReads& consReads) {
    actionsAfterOverlapSet();
    m_chimera.resize(m_overlaps.isize(), 0);
    int chimeras = 0;
    for (int i=0; i<m_overlaps.isize(); i++) {
        bool b = m_overlaps[i].organizeLaps(consReads, i);
        if (b) {
            m_chimera[i] = 1;
            chimeras++;
        }
    }
    cout << "Reported " << chimeras << " reads as possibly chimeric." << endl;
}
//======================================================



//======================================================
void ReadOverlap::set(int oI, int cP, float s, int d, int o) {
    m_overlapIndex = oI;
    m_contactPos   = cP;
    m_score        = s;
    m_direction    = d;
    m_orient       = o;
}


bool ReadOverlap::operator < (const ReadOverlap & rO) const {
    if (getDirection() != rO.getDirection()) { 
        return (getDirection() < rO.getDirection());
    }
    if (getContactPos() != rO.getContactPos()) {
        return (getContactPos() < rO.getContactPos());
    }
    return (getOverlapIndex() < rO.getOverlapIndex()); 
}

string ReadOverlap::toString() const {
    stringstream ss;
    ss << getOverlapIndex() << "\t" 
       << getContactPos() << "\t" << getScore() << "\t" 
       << "\t" << (getDirection()==1?">":"<") << "\t"
       << (getOrient()==1?"+":"-");
    return ss.str();
}


//======================================================
void ReadInfo::sortLaps() {
    Sort(m_rightOverlaps);
    Sort(m_leftOverlaps);
    Sort(m_overlapIds);
}

bool ReadInfo::organizeLaps(const ConsensReads& allReads, int id) {
    if (isChimera(allReads.getSize(id))) {
        m_rightOverlaps.clear();
        m_leftOverlaps.clear();
        return true;
    }
    return false;
}

bool ReadInfo::isChimera(int len)
{
    int iRight = -1;
    int iLeft = -1;

    if (getNumDirLaps(1) < 5 && getNumDirLaps(-1) < 5) // HARD CODED
        return false;

    int i;
    for (i=0; i<getNumDirLaps(1); i++) {
        const ReadOverlap & l = getDirLap(i, 1);
        if (l.getContactPos() > 3) { // HARD CODED
            iRight = i;
            break;
        }
    }
    for (i=0; i<getNumDirLaps(-1); i++) {
        const ReadOverlap & l = getDirLap(i, -1);
        if (l.getContactPos() > 3) { //HARD CODED
            iLeft = i;
            break;
        }
    }
    if (iRight == -1 || iLeft == -1)
        return false;
 
    const ReadOverlap & left = getDirLap(iLeft, -1);
    const ReadOverlap & right = getDirLap(iRight, 1);
  
    int minRight = right.getContactPos();
    int maxLeft = len - left.getContactPos();
    int slack = 20; // HARD CODED
    int diff = maxLeft - minRight; 
    if (diff < slack)
        return true;

    return false;

}
//======================================================

