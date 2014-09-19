#ifndef FORCE_DEBUG
#define NDEBUG
#endif



#include "src/Saguaro/HMMDistance.h"
#include "src/DNAVector.h"
#include "util/mutil.h"
#include "base/FileParser.h"


double SimpleDistLocal(char v1, char v2) 
{
  if(v1=='N' || v2=='N' || v1=='X' || v2=='X') return -1.;
  else return DNA_Diff(v1,v2);
}

double SimpleDistLocalFloat(double a, double b) 
{
  if (a < 0. || b < 0.0)
    return -1.;

  double v = 1. - (a * b + (1.- a) * (1. - b));
  return v;
}



void HMMFeature::Read(CMReadFileStream & f, int ver)
{
  int n;
  f.Read(n);
  m_data.resize(n);
  int i;
  CMString name;
  f.Read(name);
  m_chr = (const char*)name;
  f.Read(m_pos);
  for (i=0; i<m_data.isize(); i++) {
    double d;
    f.Read(d);
    m_data[i] = d;
  }
  if (ver > 1) {
    f.Read(m_weight);
  }
}

void HMMFeature::Write(CMWriteFileStream & f)
{
  int n = m_data.isize();
  f.Write(n);

  int i;
  CMString name = (const char*)m_chr.c_str();
  f.Write(name);
 
  f.Write(m_pos);
  for (i=0; i<m_data.isize(); i++) {
    double d = m_data[i];
    f.Write(d);
  }
  f.Write(m_weight);
}

//============================================================


void CompressedFeature::SetSize(int n)
{
  m_data.resize(n, 0);
}

void CompressedFeature::Set(const char * p) 
{
  for (int i=0; i<m_data.isize(); i++)
    m_data[i] = p[i];
}

void CompressedFeature::SetFloat(const svec<double> & d) 
{
  m_floats = d;
  m_data.resize(0);
}


void CompressedFeature::SimpleFeature(svec<double> & f)
{
  int i, j;
  int n = m_data.isize();

  if (m_floats.isize() > 0) {
    f = m_floats;
  } else {
    f.resize(n);
    char cmp = m_data[n-1];
    i = n-1;
    while (cmp == -1) {
      cmp = m_data[i];
      i--;
    }
    for (i=0; i<n; i++) {
      f[i] = SimpleDistLocal(m_data[n-1], m_data[i]);
    }
  } 
  /*
  if (f[0] == 1) {
    for (i=0; i<n; i++) {
      if (f[i] < 0.)
	continue;
      f[i] = 1. - f[i];
    }
    }*/
}


void CompressedFeature::MakeFeature(HMMFeature & f) const
{
  int i, j;
  int n = m_data.isize();

  if (m_floats.isize() > 0) {
    n = m_floats.isize();
    f.resize(n*n);
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	f[i*n + j] = SimpleDistLocalFloat(m_floats[i], m_floats[j]);
      }
    }  
  } else {
    f.resize(n*n);
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	f[i*n + j] = SimpleDistLocal(m_data[i], m_data[j]);
      }
    }  
  }
}



void CompressedFeature::Read(CMReadFileStream & f, int ver)
{
  int n = 0;
  f.Read(n);
  SetSize(n);
  for (int i=0; i<m_data.isize(); i++)
    f.Read(m_data[i]);
  
  if (ver >= 4) {
    f.Read(n);
    m_floats.resize(n, 0);
    for (int i=0; i<n; i++)
      f.Read(m_floats[i]);
  }
  
}

void CompressedFeature::Write(CMWriteFileStream & f)
{
  f.Write(m_data.isize());
  for (int i=0; i<m_data.isize(); i++)
    f.Write(m_data[i]);

  f.Write(m_floats.isize());

  for (int i=0; i<m_floats.isize(); i++)
    f.Write(m_floats[i]);
}


//=========================================================================

void HMMFeatureVector::MergeRead(const string & name)
{
  int n;
  int ver;
  CMReadFileStream f;
  f.Open(name.c_str());
  
  f.Read(ver);

  int m = m_data.isize();

  f.Read(n);
  m_data.resize(n+m);
  int i;
  
  for (i=0; i<n; i++)
    m_data[i+m].Read(f, ver);
  
  int nn;
  f.Read(nn);
  m_names.resize(nn);
  for (i=0; i<nn; i++) {
    CMString s;
    f.Read(s);
    m_names[i] = (const char*)s;
  }

  if (ver > 2) {
    m_comp.resize(n+m);
    for (i=0; i<n; i++) {
      m_comp[i+m].Read(f, ver);
    }
  }

  f.Close();
}


void HMMFeatureVector::Read(const string & name)
{
  m_data.clear();
  MergeRead(name);
}
	   
void HMMFeatureVector::Write(const string & name)
{
  int n = m_data.isize();
  int ver = 4;
  CMWriteFileStream f;
  f.Open(name.c_str());
  
  f.Write(ver);
  f.Write(n);
  
  int i;
  
  for (i=0; i<n; i++)
    m_data[i].Write(f);

  int nn = m_names.isize();
  f.Write(nn);
  
  for (i=0; i<nn; i++) {
    CMString s = (const char*)m_names[i].c_str();
    f.Write(s);
  }
  

  for (i=0; i<n; i++) {
    m_comp[i].Write(f);
  }
  
  f.Close();
}


	   
//===========================================================

double HMMEuclidianDistance::Distance(const HMMFeature & v1, const HMMFeature & v2)
{
  int i;
  double d = 0.;
  for (i=0; i<v1.isize(); i++) {
    if (v1[i] < -0.5 || v2[i] < -0.5)
      continue;
    d += (v1[i] - v2[i])*(v1[i] - v2[i]);
  }
  d *= v1.Weight();
  d *= v2.Weight();

  return d;
}


HMMTreeDistance::HMMTreeDistance()
{
  m_shift = 8;

  m_n = (1 << m_shift);
  
  m_cache.resize(m_n*m_n, 0.);
 
  int i, j;
  for (i=1; i<m_n; i++) {
    double d1 = (double)i/(double)(m_n-1);
    for (j=0; j<m_n; j++) {
      double d2 = (double)j/(double)(m_n-1);
      int index = Index(d1, d2);
      m_cache[index] = log(exp(-d1)+d2-2*d2*exp(-d1));
    }
  }

}

int HMMTreeDistance::Index(double d1, double d2)
{
  int i1 = (int)(d1 * (m_n-1) + 0.5);
  int i2 = (int)(d2 * (m_n-1) + 0.5);
  
  if (i1 < 0 || i2 < 0 || i1 >= m_n || i2 >= m_n)
    return -1;


  int i = (i1 << m_shift) + i2;

  //if (i >=  65536) {
  //cout << d1 << "  " << d2 << endl;
  //}
  return i;
}

// Distance matrix score based on tree scoring calculation
// v1 is hypothesis vector, v2 is observation vector
double HMMTreeDistance::Distance(const HMMFeature & v1, const HMMFeature & v2)
{
  int i;
  double d = 0.;
  for (i=0; i<v1.isize(); i++) {
    if (v1[i] == 0 || v1[i] < -0.5 || v2[i] < -0.5)
      continue;
    
    int index = Index(v1[i], v2[i]);
    double val2 = 0.;
    if (index >= 0)
      val2 = m_cache[Index(v1[i], v2[i])];
    else
      val2 = log(exp(-v1[i])+v2[i]-2*v2[i]*exp(-v1[i]));

    //double diff = val2 - val;
    //if (diff > 0.05 || diff < -0.05)
    //cout << "val=" << val << "  cache=" << val2 << " v1=" << v1[i] << endl; 
    d += val2;
  }
  d *= v1.Weight();
  d *= v2.Weight();
  if (d < 0.) {
    //cout << d << endl;
    d = -d;
  }
  return d;
}

//HMMTrees::HMMTrees()
//{
//  m_pDist = new HMMEuclidianDistance;
//}

HMMTrees::HMMTrees()
{
  m_pDist = new HMMTreeDistance;
  m_bSorted = false;
}

HMMTrees::~HMMTrees()
{
  delete m_pDist;
}


void Swap(HMMFeatureVector &v, int i, int j) 
{
  HMMFeature tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
}

void HMMTrees::Sort()
{
  m_bSorted = true;
  cout << "Sorting models..." << endl;
  
  /*m_searchNames.resize(m_names.isize());

  for (int i=0; i<m_names.isize(); i++) {
    m_searchNames[i].SetIndex(m_names[i].Index());
    m_searchNames[i].SetName(m_names[i].String());
    }*/
  m_searchNames = m_names;

  ::Sort(m_searchNames);

  /*
  int i;
  for (i=0; i<m_names.isize(); i++) {
    int j = m_names[i].Index();

    Swap(m_model, i, j);
    Swap(m_weights, i, j);
    Swap(m_same, i, j);
    Swap(m_diff, i, j);

    m_names[i].SetIndex(i);
    }*/
}

// Read and add from ASCII file
void HMMTrees::AddRead(const string & file)
{
 
  FlatFileParser parser;
  parser.Open(file);
  int i, j;

 
  while (parser.ParseLine()) {
    
    if (parser.GetItemCount() == 0)
      continue;

    string name = parser.AsString(0);
    HMMFeature f;

    parser.ParseLine();
    int n = parser.GetItemCount();
    f.resize(n*n);

    if(parser.IsEndOfFile()) continue;

    for (i=0; i<n; i++) 
      m_model.SetName(i, parser.AsString(i));

    for (i=0; i<n; i++) {

      parser.ParseLine();

      for (j=0; j<n; j++) {
	double d = parser.AsFloat(j+1);
	int index = i*n +j;
	f[index] = d;
      }

    }

    Add(f, name);

  }

}

void HMMTrees::PrettyPrint()
{
  int i, j, k;
  int n = m_model.GetNameCount();
  //cout << "Names in models: " << n << endl << endl;
  for (k=0; k<m_model.isize(); k++) {
    //cout << "------------------------------------------------------" << endl;
    printf("%s\n", m_names[k].String().c_str());

    for (i=0; i<n; i++) {
      if (i > 0)
	printf("\t");
      printf("%s", m_model.GetName(i).c_str());     
    }
    printf("\n");

    const HMMFeature & feat = m_model[k];

    for (i=0; i<n; i++) {
      printf("%s", m_model.GetName(i).c_str());
      for (j=0; j<n; j++) {
	int index = i * n + j;
	printf("\t%2.2f", feat[index]);	
      }
      printf("\n");
    }
  }
  //cout << "Done printing." << endl;
}

void HMMTrees::PrettyPrintOne(const string & name) 
{
  int i, j, k;
  int n = m_model.GetNameCount();
  svec<int> bad;
  bad.resize(n, 0);
  //cout << "Names in models: " << n << endl << endl;
  for (k=0; k<m_model.isize(); k++) {
    if (m_names[k].String() != name)
      continue;
    //cout << "------------------------------------------------------" << endl;
    //printf("%s\n", m_names[k].String().c_str());
    const HMMFeature & feat = m_model[k];

    for (i=0; i<n; i++) {
      double val = feat[i*n+i];
      if (val > 0.5) {
	bad[i] = 1;
	continue;
      }
      if (i > 0)
	printf("\t");
      printf("%s", m_model.GetName(i).c_str());     
    }
    printf("\n");


    for (i=0; i<n; i++) {
    
      if (bad[i])
	continue;
      printf("%s", m_model.GetName(i).c_str());
      for (j=0; j<n; j++) {
	if (bad[j])
	  continue;
	int index = i * n + j;
	printf("\t%2.2f", feat[index]);	
      }
      printf("\n");
    }
  }
}

void HMMTrees::PrettyPrint(const string & file)
{
  FILE * p = fopen(file.c_str(), "w");

  int i, j, k;
  int n = m_model.GetNameCount();
  //cout << "Names in models: " << n << endl << endl;
  for (k=0; k<m_model.isize(); k++) {
    //cout << "------------------------------------------------------" << endl;
    fprintf(p, "%s\n", m_names[k].String().c_str());

    for (i=0; i<n; i++) {
      if (i > 0)
	fprintf(p, "\t");
      fprintf(p, "%s", m_model.GetName(i).c_str());     
    }
    fprintf(p, "\n");

    const HMMFeature & feat = m_model[k];

    for (i=0; i<n; i++) {
      fprintf(p, "%s", m_model.GetName(i).c_str());
      for (j=0; j<n; j++) {
	int index = i * n + j;
	fprintf(p, "\t%f", feat[index]);	
      }
      fprintf(p, "\n");
    }
  }
  fclose(p);
}

