#include "src/Papaya/NonRedKmerTable.h"

bool Regular(char c) {
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
    return true;
  else
    return false;
}


// Restricts all k-mers to what's in here.
void NonRedKmerTable::SetUp(const vecDNAVector & templ, bool noNs)
{
  int i, j;

  int l = 0;
  for (i=0; i<templ.isize(); i++)
    l += templ[i].isize()-m_k+1;

  cout << "Allocating: " << l << endl;
  m_data.resize(l);

  l = 0;
  for (i=0; i<templ.isize(); i++) {
    const DNAVector & d = templ[i];
    for (j=0; j<=d.isize()-m_k; j++) {      
      d.Substring(m_data[l], j, m_k);
      
      if (noNs) {
	bool bGood = true;
	for (int x=0; x<m_k; x++) {
	  if (!Regular((m_data[l])[x]))
	    bGood = false;
	}
	if (!bGood)
	  continue;
      }
      


      //cout << m_data[l] << endl;
      l++;
    }
  }
  cout << "Resizing: " << l << endl;
  m_data.resize(l);

  UniqueSort(m_data);
  m_counts.resize(m_data.isize(), 0);
  cout << "Done." << endl;
  
 
}

void NonRedKmerTable::AddData(const vecDNAVector & data)
{
  int i, j;

  cout << "Reads: " << data.isize() << endl;
  for (i=0; i<data.isize(); i++) {
    
    if (i % 100000 == 0) {
      cout << 100.*(double)i/(double)data.isize() << " %" << endl;
    }
    const DNAVector & d = data[i];
    for (j=0; j<=d.isize()-m_k; j++) {      
      string s; 
      d.Substring(s, j, m_k);
	

      int k = BinSearch(m_data, s);
      if (k < 0)
	continue;
      m_counts[k]++;
    }
  }

}

void NonRedKmerTable::AddData(vecDNAVectorStream & data)
{
  int i, j;

  cout << "Reads: (unknown)" << endl;
  while (true) {
    const DNAVector & d = data.Next();

    if (d.isize() == 0)
      break;

    for (j=0; j<=d.isize()-m_k; j++) {      
      string s; 
      d.Substring(s, j, m_k);
	

      int k = BinSearch(m_data, s);
      if (k < 0)
	continue;
      m_counts[k]++;
    }
  }

}

