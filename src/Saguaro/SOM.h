#ifndef SOM_H_
#define SOM_H_

#include "base/SVector.h"
#include "base/RandomStuff.h"


class SOMNode;

class CoordsList
{
 public:
  CoordsList() {}

  void Read(const string & file);
  void Write(const string & file);

  int isize() const {return m_data.isize();}

  void push_back(svec<double> & d) {
    m_data.push_back(d);
  }

  svec<double> & operator [] (int i) {return m_data[i];}
  const svec<double> & operator [] (int i) const {return m_data[i];}

 private:
  svec< svec<double> > m_data;
};



class SOMConnection
{
 public:
  SOMConnection() {}

  void Add(SOMNode * p, double dist) {
    m_nodes.push_back(p);
    m_distance.push_back(dist);
  }


  int GetNodeCount() const {return m_nodes.isize();}
  SOMNode * GetNode(int i) {return m_nodes[i];}
  double GetDist(int i) {return m_distance[i];}

 private:
  svec<SOMNode *> m_nodes;
  svec<double> m_distance;
};


class SOMNode
{
 public:
  SOMNode() {
    m_hits = 0.;
  }

  void SetDimRange(int nExt, double fromExt, double toExt) {
    m_extCoords.resize(nExt);
    //m_mapCoords.resize(nMap, 0.);
    int i;
    for (i=0; i<nExt; i++) {
      double r = RandomFloat(toExt-fromExt);
      m_extCoords[i] = fromExt + r;
    }
  }
  
  SOMConnection & Connect() {return m_connect;}

  double Similarity(svec<double> & feat) const {
    double sum = 0.;
    double n = 0.;
    for (int i=0; i<m_extCoords.isize(); i++) {
      double d = feat[i] - m_extCoords[i];
      if (d < 0.)
	d = -d;
      n += 1.;
      sum += exp(-d);
      
    }
    return (sum/n);
  }

  double PhysDistance(const SOMNode * p) {
    for (int i=0; i<m_connect.GetNodeCount(); i++) {
      if (m_connect.GetNode(i) == p) {
	return m_connect.GetDist(i);
      }
    }
    return -1;
  }

  double Distance(const svec<double> & feat) const {
    svec<int> valid;
    valid.resize(feat.isize(), 1);
    return Distance(feat, valid);
  }


  double Distance(const svec<double> & feat, const svec<int> & valid) const {
    double n = 0.;
    double d = 0.;
    for (int i=0; i<m_extCoords.isize(); i++) {
      if (valid[i] > 0) {
	n += 1.;
	d += (feat[i] - m_extCoords[i]) * (feat[i] - m_extCoords[i]);
      }
    }
    return (d/n);
  }

  void IncHits() {m_hits += 1.;}

  void Update(svec<double> & feat, svec<int> & valid, double distance, double weight = 1.) {
    //double w = exp(-distance*distance) * weight;
    double w = 1/(distance*distance + 1.) * weight;
    for (int i=0; i<m_extCoords.isize(); i++) {
      if (valid[i] > 0) {
	//cout << distance << " " << m_extCoords[i] << " -> ";
	m_extCoords[i] = m_extCoords[i]*(1.-w) + feat[i]*w;
	//cout << m_extCoords[i] << endl;
      }
    }    
  }

  void UpdateDepend(svec<double> & feat, svec<int> & valid, double weight) {
    for (int i=0; i<m_connect.GetNodeCount(); i++) {
      double distance = m_connect.GetDist(i);
      SOMNode * p = m_connect.GetNode(i);
      p->Update(feat, valid, distance, weight);           
    }    
  }

  double DensityScore() {
    //===================================================
    return m_hits;
    //===================================================

    double sum = m_hits;
    svec<int> valid;
    valid.resize(m_extCoords.isize(), 1);
    for (int i=0; i<m_connect.GetNodeCount(); i++) {
      double distance = m_connect.GetDist(i);
      SOMNode * p = m_connect.GetNode(i);
      double d = p->Distance(m_extCoords, valid);    
      //d = log(d);
     
      sum += p->Hits() * exp(-distance)/d;
      //cout << sum << "\t" << d << "\t" << exp(-distance) << "\t" <<  p->Hits() << endl;
      //sum += d * 1./(distance*distance + 1.);
    }
    return sum;
  }

  void Print() {
    int n1 = 0;
    int n2 = 0;
    cout << "(";
    for (int i=0; i<m_extCoords.isize(); i++) {
      if (i > 0)
	cout << ", ";
      cout << m_extCoords[i];
      if (i < 10 && m_extCoords[i] > 0.5)
	n1++;
      if (i >= 10 && m_extCoords[i] < 0.5)
	n2++;
    }
    cout << ")";

    cout << " hits: " << m_hits;
    if (n1 >8 && n2 > 8)
      cout << " SEPARATION";
    cout << endl;
  }

  double Hits() const {return m_hits;}
  svec<double> & Coords() {return m_extCoords;}

 private:
  svec<double> m_extCoords;
  //svec<double> m_mapCoords;

  SOMConnection m_connect;
  double m_hits;
  
};






class SOMOrganizer
{
 public:
  SOMOrganizer() {}
  virtual ~SOMOrganizer() {}

  virtual void Organize(svec<SOMNode> & out, int n) = 0;

  void SetDimRange(int dim, double from, double to) {
    m_dim = dim;
    m_from = from;
    m_to = to;
  }

 protected:
  int m_dim;
  double m_from;
  double m_to;
};


class OneDimSOMOrganizer : public SOMOrganizer
{
  
  virtual void Organize(svec<SOMNode> & out, int n);

};

class TwoDimSOMOrganizer : public SOMOrganizer
{
  
  virtual void Organize(svec<SOMNode> & out, int n);

};





class SelfOrgFeatMap
{
 public:
  SelfOrgFeatMap() {
    m_weight = 1.;
    m_count = 0.;
    m_topIndex = -1;
    m_secondIndex = -1;
    m_maxDist = 0.;
  }

  void Organize(SOMOrganizer & o, int n) {
    o.Organize(m_nodes, n);
  } 

  SOMNode & GetNode(int i) {return m_nodes[i];}

  void Train(svec<double> & feat, svec<int> & valid);

  void IncCounter() {m_count++;}

  int DensestNode();
  int DensestNode(const CoordsList & l);

  int BestMatch(svec<double> & feat, svec<int> & valid);

  int GetTopContenders(svec<int> & idx, int n);
  
  void Print();
  double Weight() const {return m_weight;}

  int GetTop() const {return m_topIndex;}
  int GetSecond() const {return m_secondIndex;}
  double GetDist() const {return m_maxDist;}

 private:
  svec<SOMNode> m_nodes;
  double m_weight;
  double m_count;

  int m_topIndex;
  int m_secondIndex;
  double m_maxDist;

};



#endif //SOM_H_

