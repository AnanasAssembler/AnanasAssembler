

#include "src/Papaya/KmerHash.h"


KmerHash::KmerHash(int size)
{
  int i;
  int m = 1;

  m_origSize = size;
  m_halfK = (size+1)/2;

  
  for (i=0; i<m_halfK; i++) {
    m *= 4;
  }
  cout << "k/2=" << m_halfK << "  table: " << m << endl;
  m_table.resize(m);
}


void KmerHash::Build(const vecDNAVector & d)
{
  int i, j;
  long long bases = 0;

  for (i=0; i<d.isize(); i++) {
    const DNAVector & s = d[i];
    for (j=0; j<=s.isize()-2*m_halfK; j++) {
      int a = KmerToNumberOne(s, j);
      int b = KmerToNumberTwo(s, j);

      if (a < 0 || b < 0)
	continue;

      //cout << "a=" << a << " b=" << b << endl;
      int index = Hash(a, b);
      //cout << "Hash=" << index << endl;
      //cout << s[j] << s[j+1] << s[j+2] << s[j+3] << "  ";
      //cout << "Index: " << index << endl;
      m_table[index].Inc(b);

      bases ++;

      if (bases % 1000000 == 0) {
	cout << "Processed: " << (double)bases/1000000 << " Mb" << endl;
      }
    }
  }
}

int KmerHash::Count(const DNAVector & d, int pos) const
{
 
  int a = KmerToNumberOne(d, pos);
  int b = KmerToNumberTwo(d, pos);
  if (a < 0 || b < 0)
    return 0;
  
  int index = Hash(a, b);
  return m_table[index].Count(b);
}

