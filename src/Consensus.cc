#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "Consensus.h"
#include "SearchOverlaps.h"
#include "ConsensOverlapUnit.h"
#include "ContScaff.h"

void ConsensusBuilder::Build(DNAVector & out, const Contig& cont, const ConsensOverlapUnit & COUnit)
{
  const ConsensReads & consReads = COUnit.getConsReads();

  int len = 0;
  for (int i=0; i<cont.isize(); i++) {
    if (cont[i].Stop() > len) {
      len = cont[i].Stop();
    }
  }

  Consensus cons;
  cons.resize(len);
  
  int n = 0;
  for (int i=0; i<cont.isize(); i++) {
    int r = cont[i].Read();
    DNAVector d = consReads[r];
    if (cont[i].Ori() == -1)
      d.ReverseComplement();
    if (cont[i].Start() + d.isize() > n)
      n = cont[i].Start() + d.isize();
    for (int x=0; x<d.isize(); x++) {
      cons.Add(cont[i].Start() + x, d[x]);
    }
  }    
  out.resize(len);
  for (int i=0; i<len; i++) 
    out[i] = cons[i];
}

LayoutSink::LayoutSink()
{
  m_counter = -1;
  m_minor = 0;
  m_pLayout = NULL;
  m_minIdent = 0.99;
  m_index  = 0;
  m_prefix = "Sample1";
}

LayoutSink::~LayoutSink()
{
  if (m_pLayout != NULL) {
    fprintf(m_pLayout, "<SCAFFOLD_READCOUNT> %s %d </SCAFFOLD_READCOUNT>\n", m_lastScaffName.c_str(), (int)m_currReads.size());
    fprintf(m_pLayout, "<SCAFFOLD_PAIRCOUNT> %s %d </SCAFFOLD_PAIRCOUNT>\n", m_lastScaffName.c_str(), (int)m_currPairedReads.size()/2);
    fprintf(m_pLayout, "</SCAFFOLD> %s\n\n", m_lastScaffName.c_str());
    fclose(m_pLayout);
  }
}
 
void LayoutSink::SetLayoutFile(const string & layout)
{
  if (m_pLayout != NULL) { fclose(m_pLayout); }
  m_pLayout = fopen(layout.c_str(), "w");
  if (m_pLayout == NULL) 
    cout << "ERROR: Could not open file " << layout << " for writing!!" << endl;
}

bool LayoutSink::IsSame(const DNAVector & d)
{
  int n = m_lastCons.isize();
  int m = d.isize();
  
  if (m < n)
    n = m;

  if (n == 0)
    return false;

  int i;
  int match = 0;
  
  for (i=0; i<n; i++) {
    if (m_lastCons[i] == d[i])
      match++;
  }
  double ident = (double)match/(double)n;
  if (ident >= m_minIdent)
    return true;
  else
    return false;

}

void LayoutSink::Dump(const Hypothesis & hyp, const ConsensOverlapUnit & COUnit, bool bPrev, int minLen)
{
  if (!bPrev || m_counter < 0) {
    m_counter++;
    m_minor = 0;
    m_lastCons.resize(0);
  } 

  if (bPrev)
    m_minor++;

  int currNumPairedReads = 0;

  int i;
  char name[512];
  sprintf(name, ">Contig_%s_%3d_%7d_%3d", m_prefix.c_str(), m_index, m_counter, m_minor);
  for (i=0; i<(int)strlen(name); i++) {
    if (name[i] == ' ')
      name[i] = '0';
  }
  m_lastContigName = name;
  
  const ConsensReads & rr = COUnit.getConsReads();
  char scaffName[1024];
  if(!bPrev) {
    if(m_counter>0) { 
      fprintf(m_pLayout, "<SCAFFOLD_READCOUNT> %s %d </SCAFFOLD_READCOUNT>\n", m_lastScaffName.c_str(), (int)m_currReads.size());
      fprintf(m_pLayout, "<SCAFFOLD_PAIRCOUNT> %s %d </SCAFFOLD_PAIRCOUNT>\n", m_lastScaffName.c_str(), (int)m_currPairedReads.size()/2);
      m_currPairedReads.clear();
      m_currReads.clear();
      fprintf(m_pLayout, "</SCAFFOLD> %s\n\n", m_lastScaffName.c_str());
    }
    char tmp[1024];
    sprintf(scaffName, ">Scaffold_%5d", m_counter);
    for (unsigned int x = 0; x<strlen(scaffName); x++) {
      if (scaffName[x] == ' ')
      scaffName[x] = '0';
    }
    fprintf(m_pLayout, "<SCAFFOLD> %s\n", scaffName);
    m_lastScaffName = scaffName;
  }
  fprintf(m_pLayout, "<CONTIG> %s %d \n", name, hyp.Pairs());
  for (i=0; i<hyp.Size(); i++) {
    int r = hyp[i].Read();
    int numPartner = COUnit.getNumOfPartners(r);
    fprintf(m_pLayout, "%d\t%d\t%d - %d\t%d\t%d\n", r, hyp[i].Ori(), hyp[i].Start(), hyp[i].Stop(), 
            hyp[i].Pair()>=0?hyp[hyp[i].Pair()].Read():hyp[i].Pair(),  
            hyp[i].Pair()>=0?hyp[hyp[i].Pair()].Ori():hyp[i].Pair());
    m_currReads[r] = true;
    if (hyp[i].Pair() >= 0) {
      m_currPairedReads[r] = true;
      currNumPairedReads++;
    }
  }
  fprintf(m_pLayout, "<CONTIG_READCOUNT> %s %d </CONTIG_READCOUNT>\n", name, hyp.Size());
  fprintf(m_pLayout, "<CONTIG_PAIRCOUNT> %s %d </CONTIG_PAIRCOUNT>\n", name, currNumPairedReads/2);
  fprintf(m_pLayout, "</CONTIG> %s\n", name);
  fflush(m_pLayout);
}

void LayoutSink::fastaFromAssembly(const string& fastaFile, const Assembled& asmb, const ConsensOverlapUnit & COUnit, int minLen)
{
  FILE* pFasta = fopen(fastaFile.c_str(), "w");
  if (pFasta == NULL) {
    cout << "ERROR: Could not open file " << fastaFile << " for writing!!" << endl;
  }
  for(int scafCnt=0; scafCnt<asmb.isize(); scafCnt++) {
    Scaffold currScaff = asmb[scafCnt];
    for(int contCnt=0; contCnt<currScaff.isize(); contCnt++) {
      Contig currCont = currScaff[contCnt];
      ConsensusBuilder cons;
      DNAVector dna;
      cons.Build(dna, currCont, COUnit);
      if (IsSame(dna)) {
        cout << "Deep toasting sequence " << currCont.Name() << endl;
        return;
      }
      if(dna.isize()<minLen) { continue; }
      fprintf(pFasta, "%s\n", currCont.Name().c_str());
      for (int i=0; i<dna.isize(); i++) {
        fprintf(pFasta, "%c", dna[i]);
      }
      fprintf(pFasta, "\n");
      fflush(pFasta); 
    }
  }
}
