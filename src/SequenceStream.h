#ifndef SEQUENCESTREAM_H
#define SEQUENCESTREAM_H

#include "src/DNAVector.h"
#include "base/FileParser.h"


class SequenceStream
{
public:
  SequenceStream() {
    m_frame = 0;
    m_trans = false;
    m_k = 0;
    m_wantTrans = true;
    m_fq = false;
    m_extCount = 0;
  }

  void Open(const string & fileName);

  //True if there is one, false if EOF
  bool GetNext(DNAVector & d, string & name);

  void SetProtein(bool b) {
    m_wantTrans = b;
  }
  
  void SetSource(const vecDNAVector & s) {
    m_extCount = 0;
    m_extSource = s;
  }

private:
  bool NeedsTrans(const DNAVector & d);

  FlatFileParser m_parser;
  bool m_trans;
  bool m_wantTrans;
  int m_frame;
  vecDNAVector m_cache;
  int m_k;
  string m_name;
  bool m_fq;
  vecDNAVector m_extSource;
  int m_extCount;
};




#endif //SEQUENCESTREAM

