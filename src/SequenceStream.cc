#include "src/SequenceStream.h"
#include "src/MakeSixFrames.h"

const char * NameClean(const string & n)
{
  //return n.c_str();
  const char * p = n.c_str();
  return &p[1];
}

void SequenceStream::Open(const string & fileName)
{
  m_parser.Open(fileName);
  m_parser.ParseLine(); 
  m_trans = false;
  m_name = NameClean(m_parser.AsString(0));
  if ((m_parser.AsString(0))[0] != '>')
    m_fq = true;
}

  //True if there is one, false if EOF
bool SequenceStream::GetNext(DNAVector & d, string & name)
{
  if (m_trans) {
    if (m_k < 6) {
      name = m_cache.Name(m_k);
      d = m_cache[m_k];
      m_k++;
      return true;
    }
  }

  MakeSixFrames frames;

  if (m_extSource.isize() == 0) {
    bool b = false;
    do {
      b = m_parser.ParseLine();    
    } while(b && m_parser.GetItemCount() == 0);
    
    if (!b)
      return false;
    bool bCont = false;
    
    name = m_name;
    if (m_fq) {
      d.SetFromBases(m_parser.Line());    
      m_parser.ParseLine();
      m_parser.ParseLine();
      m_parser.ParseLine();
      if (m_parser.GetItemCount() > 0)   
	m_name = m_parser.AsString(0);
      bCont = true;
    } else {
      string seq = m_parser.Line();
      while (m_parser.ParseLine()) {
	if (m_parser.GetItemCount() == 0)
	  continue;
	if ((m_parser.Line())[0] == '>') {
	  bCont = true;
	  m_name = NameClean(m_parser.Line());
	  break;
	}
	seq += m_parser.Line();
      }
      d.SetFromBases(seq);
    }
  } else { // Have external source
    if (m_extCount == m_extSource.isize())
      return false;
    d = m_extSource[m_extCount];
    name = m_extSource.Name(m_extCount);
    m_extCount++;
    
  }

  

  if (m_wantTrans) {
    m_trans = NeedsTrans(d);
    // if (m_trans) {
    //  cout << "Yes!" << endl;
    //} else {
    //  cout << "No!" << endl;
    // }
  }


  if (m_trans) {
    m_k = 1;
    vecDNAVector tmp;
    tmp.push_back(d, name);
    frames.AllSixFrames(m_cache, tmp);
    m_cache.DoProteins(true);
    d = m_cache[0];
    name = m_cache.Name(0);
  }

  return true;
  
}


bool SequenceStream::NeedsTrans(const DNAVector & d) 
{
  int i;
  for (i=0; i<d.isize(); i++) {
    if (d[i] != 'A' && d[i] != 'C' && d[i] != 'G' && d[i] != 'T' && d[i] != 'N')
      return false;
  }
  return true;

}


