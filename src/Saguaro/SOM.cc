#include "src/Saguaro/SOM.h"
#include "base/FileParser.h"


void CoordsList::Read(const string & file)
{
  FILE * pTest = fopen(file.c_str(), "r");
  if (pTest == NULL) { 
    cout << "File not available" << endl;
    return;
  }
  fclose(pTest);

  FlatFileParser parser;
  
  parser.Open(file);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    svec<double> d;
    for (int i=0; i<parser.GetItemCount(); i++) {
      d.push_back(parser.AsFloat(i));
    }
    push_back(d);
  }

}

void CoordsList::Write(const string & file)
{
  FILE * p = fopen(file.c_str(), "w");
  int i, j;
  for (i=0; i<isize(); i++) {
    const svec<double> & d = m_data[i];
    for (j=0; j<d.isize(); j++)
      fprintf(p, "%f  ", d[j]);
    fprintf(p, "\n");
  }
  
  fclose(p);
}




class NodeMap
{
public:
  NodeMap() {
    m_index = -1;
    m_score = 0.;
  }
  void Set(int index, double score) {
    m_index = index;
    m_score = score;
  }

  bool operator<(const NodeMap & m) const {
    return (m_score < m.m_score);
  }
  int Index() const {return m_index;}

private:
  int m_index;
  double m_score;

};




void OneDimSOMOrganizer::Organize(svec<SOMNode> & out, int n)
{
  out.resize(n);
  int i, j;
  for (i=0; i<out.isize(); i++) {
    out[i].SetDimRange(m_dim, m_from, m_to);
  }
  
  double half = n/2;

  for (i=0; i<n; i++) {
    SOMNode & node = out[i];
    for (j=0; j<n; j++) {
      if (i == j)
	continue;
      double d = i-j;
      if (d < 0)
	d = -d;
      if (d > half)
	d = n-d;

      out[i].Connect().Add(&out[j], d);
    }
  }
}


void TwoDimSOMOrganizer::Organize(svec<SOMNode> & out, int n)
{
  out.resize(n*n);
  int i, j;

  for (i=0; i<out.isize(); i++) {
    out[i].SetDimRange(m_dim, m_from, m_to);
  }
  for (i=0; i<out.isize(); i++) {
    double x1 = i / n;
    double y1 = i - i*(i/n);
    SOMNode & node = out[i];
    for (j=0; j<out.isize(); j++) {
      if (i == j)
	continue;
      double x2 = j / n;
      double y2 = j - j*(j/n);
      
      double d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
      out[i].Connect().Add(&out[j], d);
    }
  }
}



void SelfOrgFeatMap::Train(svec<double> & feat, svec<int> & valid) 
{
  m_weight = 0.2/(m_count + 1.);
  m_count += 0.0001;
  int best = BestMatch(feat, valid);
  //cout << "best=" << best << " w=" << m_weight << endl;
 
  //cout << "Training the feature: ";
  //for (int i=0; i<feat.isize(); i++) {
  //  cout << feat[i] << "\t" << valid[i] << "\t";
  //}
  //cout << endl;
  
  if (best == -1)
    return;

  SOMNode & node = m_nodes[best];
  //if (m_weight < 0.25)
  node.IncHits();

  node.Update(feat, valid, 0., m_weight);
  node.UpdateDepend(feat, valid, m_weight);
  
}

int SelfOrgFeatMap::GetTopContenders(svec<int> & idx, int n)
{
  throw;
  return 0;
}


int SelfOrgFeatMap::DensestNode(const CoordsList & used)
{

  int i, j;
  double score = 99999999999.;
  int index = -1;
  svec<NodeMap> nm;
  nm.resize(m_nodes.isize());
  m_topIndex = -1;
  m_secondIndex = -1;

  for (i=0; i<m_nodes.isize(); i++) {
    double d = m_nodes[i].DensityScore();
    
    nm[i].Set(i, -d);
    cout << "density " << i << "\t" << d << endl;
    if (d < score) {
      score = d;
      index = i;
      
    }
      
  }
  cout << "Function(1)" << endl;

  int top = 100;
  if (m_nodes.isize() < top)
    top = m_nodes.isize();
  
  
  Sort(nm);

  cout << "nm.isize: " << nm.isize() << endl;
  m_maxDist = 0.;
  for (i=0; i<top; i++) {
    int index = nm[i].Index();

    cout << "index=" << index << endl;
    double usedsumOne = 0.;
    for (int x=0; x<used.isize(); x++) {
      cout << "x=" << x << endl;
      double useddist = m_nodes[index].Distance(used[x]);
      svec<double> dummy = used[x];
      cout << "useddist=" << useddist << endl;
      for (int y=0; y<dummy.isize(); y++) {
	cout << "y=" << y << endl;
	dummy[y] = 1. - dummy[y];
      }
      cout << "Uff..."  << endl;
      double useddist2 = m_nodes[index].Distance(dummy);
      if (useddist2 < useddist)
	useddist = useddist2;

      useddist = sqrt(useddist);

      usedsumOne += useddist;
    }
    cout << "Check(1)" << endl;
    for (j=i+1; j<top; j++) {
      int index2 = nm[j].Index();
      double dist = m_nodes[index].Distance(m_nodes[index2].Coords());
      SOMNode dummy = m_nodes[index2];
      for (int x=0; x<dummy.Coords().isize(); x++) {
	(dummy.Coords())[x] = 1. - (dummy.Coords())[x];
      }
      double dist2 = m_nodes[index].Distance(dummy.Coords());
      if (dist2 < dist)
	dist = dist2;

      dist = sqrt(dist);
      
      double usedsum = 0.;
      for (int x=0; x<used.isize(); x++) {
	double useddist = m_nodes[index2].Distance(used[x]);
	svec<double> dummy = used[x];
	for (int y=0; y<dummy.isize(); y++) {
	  dummy[y] = 1. - dummy[y];
	}
	double useddist2 = m_nodes[index2].Distance(dummy);
	if (useddist2 < useddist)
	  useddist = useddist2;
	
	useddist = sqrt(useddist);
	usedsum += useddist;
      }
      cout << "Check(2)" << endl;
 
      cout << i << "\t" << j << "\t" << dist << "\t" << usedsumOne << "\t" << usedsum << endl;
      dist += usedsum;
      dist += usedsumOne;

      if (dist > m_maxDist) {
	cout << "New best: " << dist << endl;
	m_maxDist = dist;
	m_topIndex = index;
 	m_secondIndex = index2;
      }
    }

  }

  cout << "Printing all" << endl;
  for (i=0; i<m_nodes.isize(); i++) {
    int index = nm[i].Index();

    double d = m_nodes[index].DensityScore();

    cout << index << "\t" << d << "\t" << "\t";
    m_nodes[index].Print();
  }



  return m_secondIndex;
}


int SelfOrgFeatMap::DensestNode()
{
  int i, j;
  double score = 99999999999.;
  int index = -1;
  svec<NodeMap> nm;
  nm.resize(m_nodes.isize());
  m_topIndex = -1;
  m_secondIndex = -1;

  for (i=0; i<m_nodes.isize(); i++) {
    double d = m_nodes[i].DensityScore();
    
    nm[i].Set(i, -d);
    cout << "density " << i << "\t" << d << endl;
    if (d < score) {
      score = d;
      index = i;
      
    }
      
  }
  cout << "Function(2)" << endl;

  svec<int> dont;
  cout << "Removing close neighbors..." << endl;
  dont.resize(m_nodes.isize(), 0);
  
  int top = 100;
  if (m_nodes.isize() < top)
    top = m_nodes.isize();
  
  
  Sort(nm);

  for (i=0; i<top; i++) {
    cout << "Node " << nm[i].Index() << endl;
    for (j=i+1; j<top; j++) {
      
      //if (i == j)
      //continue;
      int x = nm[i].Index();
      int y = nm[j].Index();
      double d = m_nodes[x].PhysDistance(&m_nodes[y]);
      if (d < 1.5) {
	dont[y]++;
      }
    }
  }

  cout << "Printing sorted: " << endl;
  int first = nm[0].Index();
  m_topIndex = first;

  SOMNode dummy = m_nodes[first];
  for (i=0; i<dummy.Coords().isize(); i++) {
    (dummy.Coords())[i] = 1. - (dummy.Coords())[i];
  }

  m_maxDist = 0.;
  for (i=0; i<top; i++) {
    int index = nm[i].Index();
    for (j=i+1; j<top; j++) {
      int index2 = nm[j].Index();
      double dist = m_nodes[index].Distance(m_nodes[index2].Coords());
      SOMNode dummy = m_nodes[index2];
      for (int x=0; x<dummy.Coords().isize(); x++) {
	(dummy.Coords())[x] = 1. - (dummy.Coords())[x];
      }
      double dist2 = m_nodes[index].Distance(dummy.Coords());
      if (dist2 < dist)
	dist = dist2;
      
      if (dist > m_maxDist && i < top && dont[index2] == 0) {
	m_maxDist = dist;
	//m_topIndex = index;
 	m_secondIndex = index2;
     }
    }
  }

  for (i=0; i<m_nodes.isize(); i++) {
    int index = nm[i].Index();

    double d = m_nodes[index].DensityScore();

    cout << index << "\t" << d << "\t" << "\t";
    if (dont[index] == 0)
      cout << "OK\t";
    else
      cout << "nope\t";
    m_nodes[index].Print();
  }


 return index;
}

int SelfOrgFeatMap::BestMatch(svec<double> & feat, svec<int> & valid)
{
  int i;
  double score = 99999999999.;
  int index = -1;


  for (i=0; i<m_nodes.isize(); i++) {
    double d = m_nodes[i].Distance(feat, valid);
    //cout << d << endl;
    if (d < score) {
      score = d;
      index = i;
    }
      
  }
  return index;
}

void SelfOrgFeatMap::Print() 
{
  int i, j;
  for (i=0; i<m_nodes.isize(); i++) {
    cout << i << "\t";
    m_nodes[i].Print();
  }

}
