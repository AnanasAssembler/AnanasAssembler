#ifndef KMERHASH_H
#define KMERHASH_H

#include "src/DNAVector.h"


class KmerHashElement
{
 public:
  KmerHashElement(int halfKey = -1, int count = 1) {
    m_halfKey = halfKey;
    m_count = count;
  }
  void Inc() {
    m_count++;
  }
  
  int HalfKey() const {return m_halfKey;}
  int Count() const {return m_count;}

 private:
  int m_halfKey;
  int m_count;
};


class KmerHashEntry
{
 public:
  KmerHashEntry() {}

  int Count(int halfKey) const {
    int i;
    for (i=0; i<m_elements.isize(); i++) {
      if (halfKey == m_elements[i].HalfKey()) {
	return m_elements[i].Count();
      }
    }
    return 0;
  }

  void Add(int halfKey) {
    m_elements.push_back(halfKey);
  }

  void Inc(int halfKey) {
    int i;
    for (i=0; i<m_elements.isize(); i++) {
      if (halfKey == m_elements[i].HalfKey()) {
	m_elements[i].Inc();
	//cout << "Found." << endl;
	return;
      }
    }
    //cout << "Adding." << endl;
    Add(halfKey);
  }

 private:
  svec<KmerHashElement> m_elements;
};


class KmerHash
{
 public:
  KmerHash(int size = 24);


  void Build(const vecDNAVector & d);

  int Count(const DNAVector & d, int pos) const;


 private:
  int Index(const DNAVector & d, int pos) const {
    return Hash(KmerToNumberOne(d, pos), KmerToNumberTwo(d, pos));
  }

  int Hash(int n1, int n2) const {
    int h = n1 ^ n2;
    return h;
  }

  int KmerToNumberOne(const DNAVector & d, int pos) const {
    return HalfKmerToNumber(d, pos);
  }

  int KmerToNumberTwo(const DNAVector & d, int pos) const {
    return HalfKmerToNumber(d, pos+m_halfK);
  }


  int HalfKmerToNumber(const DNAVector & d, int pos) const {
    int i;
    int num = 0;
    for (i=pos; i<pos+m_halfK; i++) {
      num = (num << 2);
      if (i >= d.isize())
	break;
      if (d[i] == 'N')
	return -1;
      int n = NucIndex(d[i]);
      //cout << d[i] << " " << n << endl;
      num += n;
      //cout << "n=" << n << " m=" << num << endl;
   }
    return num;
  }

  int m_halfK;
  int m_origSize;

  svec<KmerHashEntry> m_table;
  
};







#endif

