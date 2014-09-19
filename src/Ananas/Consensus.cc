#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Ananas/Consensus.h"
#include "src/Ananas/SearchOverlaps.h"
#include "src/Ananas/ConsensOverlapUnit.h"
#include "src/Ananas/ContScaff.h"

void ConsensusBuilder::Build(DNAVector & out, const Contig& cont, const ConsensOverlapUnit & COUnit)
{
  int len = 0;
  for (int i=0; i<cont.isize(); i++) {
    if (cont[i].Stop() > len) {
      len = cont[i].Stop()+1;
    }
  }

  Consensus cons;
  cons.resize(len);
  const ConsensReads & consReads = COUnit.getConsReads();
  
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
  //NOTE: There is a potential bug here where garbage bases are appended at the end.
  out.resize(out.isize()-1);
}

LayoutSink::LayoutSink()
{
  m_counter = -1;
  m_minor = 0;
  m_pLayout = NULL;
  m_minIdent = 0.99;
  m_prefix = 0;
}

LayoutSink::~LayoutSink()
{
  if (m_pLayout != NULL) {
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

  int i;
  char name[512];
  sprintf(name, ">Contig_%3d_%7d_%3d", m_prefix, m_counter, m_minor);
  for (i=0; i<(int)strlen(name); i++) {
    if (name[i] == ' ')
      name[i] = '0';
  }
  m_lastContigName = name;
  
  const ConsensReads & rr = COUnit.getConsReads();
  char scaffName[1024];
  if(!bPrev) {
    if(m_counter>0) { 
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
  fprintf(m_pLayout, "<CONTIG> %s\n", name);
  for (i=0; i<hyp.isize(); i++) {
    int r = hyp[i].Read();
    
    int numPartner = COUnit.getNumOfPartners(r);
    fprintf(m_pLayout, "%d\t%d\t%d - %d\t", r, hyp[i].Ori(), hyp[i].Start(), hyp[i].Stop());

    if (hyp[i].Pair() >= 0) {
      fprintf(m_pLayout, "\tpaired\t%d\t%d", hyp[hyp[i].Pair()].Read(),  hyp[hyp[i].Pair()].Ori());
    } else {
      if (hyp[i].Pair() == -1) {
	fprintf(m_pLayout, "\tsingle\tn/a");
      } else {
	fprintf(m_pLayout, "\tfragment\tn/a");
      }
    }
    fprintf(m_pLayout, "\n");
  }
  fprintf(m_pLayout, "</CONTIG> %s\n", name);
  fflush(m_pLayout);
}

void LayoutSink::fastaFromAssembly(const string& fastaFile, const Assembled& asmb, const ConsensOverlapUnit & COUnit)
{
  FILE* pFasta = fopen(fastaFile.c_str(), "w");
  if (pFasta == NULL) {
    cout << "ERROR: Could not open file " << fastaFile << " for writing!!" << endl;
  }
  char name[512];
  for(int scafCnt=0; scafCnt<asmb.isize(); scafCnt++) {
    Scaffold currScaff = asmb[scafCnt];
    for(int contCnt=0; contCnt<currScaff.isize(); contCnt++) {
      Contig currCont = currScaff[contCnt];
      sprintf(name, ">Contig_%3d_%7d_%3d", m_prefix, scafCnt, contCnt);
      for (int i=0; i<(int)strlen(name); i++) {
        if (name[i] == ' ')
          name[i] = '0';
      }
      ConsensusBuilder cons;
      DNAVector dna;
      cons.Build(dna, currCont, COUnit);
      if (IsSame(dna)) {
        cout << "Deep toasting sequence " << name << endl;
        return;
      }
      fprintf(pFasta, "%s\n", name);
      for (int i=0; i<dna.isize(); i++) {
        fprintf(pFasta, "%c", dna[i]);
      }
      fprintf(pFasta, "\n");
      fflush(pFasta); 
    }
  }
}
