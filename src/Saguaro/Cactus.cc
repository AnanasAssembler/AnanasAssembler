#include "src/Saguaro/Cactus.h"
#include "base/CommandLineParser.h"
#define FORCE_DEBUG



void SaveCacti(svec<Cactus> & all, const string & name)
{
  FILE * p = fopen(name.c_str(), "w");
  for (int i=0; i<all.isize(); i++)
    all[i].Write(p);
  fclose(p);
}


void LoadCacti(svec<Cactus> & all, const string & name)
{
  FlatFileParser f;
  f.Open(name);
  bool b = false;
  do {
    Cactus cc;
    if (f.IsEndOfFile())
      break;
    //cout << "reading." << endl;
    b = cc.Read(f);
    if (b) {
      all.push_back(cc);
    }
  } while(b);
}

void Cactus::Merge(const Cactus & cc)
{
  //out = all[indices[0]]; // Copy over labels etc.
  int i, j, k;
  int n = cc.Size();
  double div = Weight() + cc.Weight();
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      double val = Get(i, j)*Weight() + cc.Get(i, j)*cc.Weight();
      //for (k=0; k<indices.isize(); k++) {
      //	val += all[indices[k]].Get(i, j);
      //div += 1.;
      //}
      Set(i, j, val/div);
    }
  }
  m_weight += cc.Weight();
}

bool Cactus::Read(FlatFileParser & f)
{
  int i;
  if (!f.ParseLine())
    return false;
  if (f.GetItemCount() == 0)
    return false;
  m_name = f.AsString(0);
  if (!f.ParseLine())
    return false;
    
  //cout << "Read successful!!" << endl;
  for (i=0; i<f.GetItemCount(); i++) {
    m_label.push_back(f.AsString(i));    
  }
  m_data.resize(m_label.isize());
  i = 0;
  for (i=0; i<m_label.isize(); i++) {
    f.ParseLine();
    svec<double> & d = m_data[i];
    for (int j=1; j<f.GetItemCount(); j++) {
      //cout << "Reading value " << f.AsFloat(j) << endl;
      d.push_back(f.AsFloat(j));
    }    
  }

  if (m_data.isize() > 0)
    return true;
  return false;
}

void Cactus::Print()
{
  int i, j;
  cout << m_name << endl;
  for (i=0; i<m_label.isize(); i++) {
    cout << m_label[i] << "\t";
  }
  cout << endl;
  for (i=0; i<m_label.isize(); i++) {
    cout << m_label[i];
    svec<double> & d = m_data[i];
    for (j=0; j<m_label.isize(); j++) {
      cout << "\t" << d[j];
    }
    cout << endl;
  }
  
}

double Cactus::Distance(const Cactus & c)
{
  double sum = 0.;
  double div = 0.;
  int i, j;
  for (i=0; i<m_label.isize(); i++) {
    svec<double> & d = m_data[i];
    for (j=0; j<m_label.isize(); j++) {
      div += 1.;
      double dd = d[j] - c.Get(i, j);
      sum += dd * dd;
    }
  }
  return sum / div;
}

void Cactus::Write(FILE * p)
{
  int i, j;
  fprintf(p, "%s\n", m_name.c_str());
  for (i=0; i<m_label.isize(); i++) {
    if (i > 0)
      fprintf(p, "\t");
    fprintf(p, "%s", m_label[i].c_str());
  }
  fprintf(p, "\n");
 
  for (i=0; i<m_label.isize(); i++) {
    //cout << m_labels[i];
    fprintf(p, "%s", m_label[i].c_str());
    svec<double> & d = m_data[i];
    for (j=0; j<m_label.isize(); j++) {
      fprintf(p, "\t%1.6f", d[j]); 
	      //cout << "\t" << d[j];
    }
    fprintf(p, "\n");
  //cout << endl;
  }
  
}
