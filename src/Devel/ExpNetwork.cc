#define FORCE_DEBUG

#include "src/Devel/ExpNetwork.h"


void Node::SetDistScores(const svec<ExpVector> & model)
{
  m_dist.resize(model.isize());

  m_scores.clear();
  m_from.clear();
  m_fromModel.clear();

  double def = INFINITE;
  if (!m_in)
    def = 0.;

  m_scores.resize(model.isize(), def);
  m_from.resize(model.isize(), -1);
  m_fromModel.resize(model.isize(), -1);

  int i;
  double sum = 0.;
  double min = INFINITE;
  for (i=0; i<model.isize(); i++) {
    double d = m_data.Dist(model[i]);
    m_dist[i] = d;
    if (d < min)
      min = d;
    sum += d;
  }
  
  // HARD CODED
  double discount = 0.;

  sum /= (double)model.isize();

  for (i=0; i<m_dist.isize(); i++) {
    m_dist[i] -= min - discount;   
  }
 
}

void Node::Propagate(Node & n, double penDiff, double penSame) const
{
  int i, j;
  for (i=0; i<m_dist.isize(); i++) {
    for (j=0; j<m_scores.isize(); j++) {
      double pen = penDiff;
      if (i == j)
	pen = penSame;
      double s = m_scores[i] + m_dist[i] + pen;
      n.SetMin(j, s, m_id, i);
    }
  }
}



void Network::Connect(const string & a, const string & b)
{
  Connector cc;
  cc.SetID(m_conn.isize());
  int x = Find(a);
  int y = Find(b);
  cc.Set(x, y);
  m_conn.push_back(cc);
  m_nodes[y].SetIn();
  m_nodes[x].SetOut();
}


int Network::Find(const string & n)
{
  for (int i=0; i<m_nodes.isize(); i++) {
    if (m_nodes[i].Name() == n)
      return i;
  }
  return -1;
}

void Network::Process(svec<ExpVector> & model)
{
  int i, j;
  for (i=0; i<m_nodes.isize(); i++) {
    m_nodes[i].SetDistScores(model);
  }

  for (i=0; i<m_nodes.isize(); i++) {
    cout << "Iteration " << i << " of " << m_nodes.isize() << endl;
    for (j=0; j<m_conn.isize(); j++) {
      int a = m_conn[j].From();
      int b = m_conn[j].To();
      m_nodes[a].Propagate(m_nodes[b]);
    }
    //Print();
  }

  TraceBack();

  //while (TraceBack()) {
    //RemoveUsed();
  //}
}

void Network::Print() const
{
  int i;
  for (i=0; i<m_nodes.isize(); i++) 
    m_nodes[i].Print();
}

bool Network::TraceBack()
{
  cout << "========= TRACEBACK ===========" << endl;
  m_bestModel.clear();
  m_bestModel.resize(m_nodes.isize(), -1);
  svec<int> callDepth;
  callDepth.resize(m_nodes.isize(), -1);

  int i, j, k;
  int mod = -1;
  int index = -1;
  bool bYes = false;
  for (i=0; i<m_nodes.isize(); i++) {
    //if (m_nodes[i].HasOut())
    //continue;
    
    cout << "------------------------" << endl;
    bYes = true;
    k = i;
    bool b = true;
    int depth = 0;
    while (b) {
      int index;
      int from;
      b = m_nodes[k].Best(index, from, mod);    
      
      cout << "  " << k << " " << m_nodes[k].Name() << " Model " << index << endl;
      if (depth > callDepth[k]) {
	m_bestModel[k] = index;
	callDepth[k] = depth;
      }
      k = from;
      depth++;
      if (k < 0)
	break;
      //m_nodes[i].Print();
    }
  }

  cout << "======== FINAL =============" << endl;
  for (i=0; i<m_nodes.isize(); i++) {
    cout << "NODE: " << i << " " << m_nodes[i].Name() << " Model: " << m_bestModel[i] << endl;
  }

  return bYes;
}
