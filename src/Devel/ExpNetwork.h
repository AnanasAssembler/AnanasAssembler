#ifndef EXP_NETWORK_H
#define EXP_NETWORK_H

#include "base/SVector.h"
#include <string>

#define INFINITE 999999999.

// Some expression vector that knwos distance
class ExpVector
{
 public:
  ExpVector() {}

  int isize() const {return m_val.isize();}
  double & operator [] (int i) {return m_val[i];}
  const double & operator [] (int i) const {return m_val[i];}
  void resize(int s, int d = 0) {
    m_val.resize(s, d);
  }
  void push_back(double v) {
    m_val.push_back(v);
  }


  double Dist(const ExpVector & s) const {
    int i;
    double d = 0.;
    for (i=0; i<m_val.isize(); i++) {
      d += (m_val[i] - s.m_val[i])*(m_val[i] - s.m_val[i]); 
    }
    d /= (double)m_val.isize();
    return d;
  }
 
 private:
  svec<double> m_val;
};


class Connector
{
 public:
  Connector() {
    m_id = -1;
    m_from = -1;
    m_to = -1;
    m_penSame = 0.;
    m_penDiff = 0.;
  }

  void SetID(int i) {m_id = i;}
  void Set(int from, int to) {
    m_from = from;
    m_to = to;
  }
  
  int From() const {return m_from;}
  int To() const {return m_to;}

 private:
  int m_id;
  int m_from;
  int m_to;
  double m_penSame;
  double m_penDiff;
};


//--------------------------------------
class Node
{
 public:
  Node() {
    m_id = -1;
    m_in = 0;
    m_out = 0;
  }

  void SetIn(int c = 1) {
    m_in += c;
  }
  void SetOut(int c = 1) {
    m_out += c;
  }

  void SetID(int i) {
    m_id = i;
  }

  void SetData(const ExpVector & d) {
    m_data = d;
  }

  void SetName(const string & name) {
    m_name = name;
  }

  void SetDataSize(int n) {
    m_data.resize(n);
  }
  void SetData(int i, double v) {
    m_data[i] = v;
  }

  void SetDistScores(const svec<ExpVector> & model);

  double GetScore(int i) const {return m_scores[i];}
  void Propagate(Node & n, double penDiff = 12., double penSame = 0.) const;
  
  bool SetMin(int i, double s, int from, int model) {
    if (s < m_scores[i]) {
      m_scores[i] = s;
      m_from[i] = from;
      m_fromModel[i] = model;
      return true;
    }
    return false;
  }

  bool HasIn() const {return m_in > 0;}
  bool HasOut() const {return m_out > 0;}

  int ID() const {return m_id;}
  const string & Name() const {return m_name;}

  void Print() const {
    cout << "ID: " << m_id << " " << m_name << endl;
    for (int i=0; i<m_scores.isize(); i++) {
      cout << "  score: " << i << " = " << m_scores[i] + m_dist[i] << " -> " << m_from[i] << ", " << m_fromModel[i] << " dist: " << m_dist[i] << endl;
    }
  }

  bool Best(int & index, int & from, int & model) const {
    int i;
    if (model != -1) {
      index = model;
      model = m_fromModel[index];
      from = m_from[index];
      if (model != -1)
	return true;
      else
	return false;
    }

    int min = INFINITE;
    index = -1;
    from = -1;
    model = -1;
    for (i=0; i<m_scores.isize(); i++) {
      if (m_scores[i] + m_dist[i] < min) {
	min = m_scores[i] + m_dist[i];
	index = i;
	model = m_fromModel[i];
	from = m_from[i];
      }
    }
    if (index != -1)
      return true;
    else
      return false;
  }
  const ExpVector & Data() const {return  m_data;}
  double Dist(const Node & n) const {
    return m_data.Dist(n.Data());
  }
  double Dist(const ExpVector & d) const {
    return m_data.Dist(d);
  }

 private:
  int m_id;  
  string m_name;
  ExpVector m_data;
  svec<double> m_dist;
 
  svec<double> m_scores;
  svec<int> m_from;
  svec<int> m_fromModel;
  bool m_in;
  bool m_out;
};


class Network
{
 public:
  Network() {}

  void AddNode(const string & s, const ExpVector & v) {
    Node n;
    n.SetName(s);
    n.SetID(m_nodes.isize());
    n.SetData(v);
    m_nodes.push_back(n);
  }

  void Connect(const string & a, const string & b);

  void Process(svec<ExpVector> & model);
  void Print() const;
  
  int isize() const {return m_nodes.isize();}
  const Node & operator[] (int i) {return m_nodes[i];}
  int Model(int i) const {return m_bestModel[i];}
 private:
  bool TraceBack();

  int Find(const string & n);

  svec<Node> m_nodes;
  svec<Connector> m_conn;
  svec<int> m_bestModel;
};



#endif //EXP_NETWORK_H
