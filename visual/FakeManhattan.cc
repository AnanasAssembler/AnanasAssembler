
#define FORCE_DEBUG
#include "visual/Whiteboard.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "visual/Axes.h"
#include "visual/Geometry.h"
#include "visual/Compounds.h"

#include <iostream>
#include "base/RandomStuff.h"

class Data
{
public:
  Data() {}

  void Generate() {
    m_data.resize(20);
    int i, j;
    
    int size = 300;
    m_colId.resize(size);
    int c = 6;
    
    svec<int> spikes;
    spikes.push_back(55);
    spikes.push_back(99);
    spikes.push_back(240);
    //spikes.push_back(133);
    

    for (i=0; i<m_data.isize(); i++) {
      svec<double> & d = m_data[i];
      d.resize(size, 0);
      for (j=0; j<size; j++) {
	if (j == 80)
	  c = 5;
  	if (j == 112)
	  c = 7;
	if (j == 140)
	  c = 9;
 	if (j == 166)
	  c = 3;
	if (j == 190)
	  c = 11;
	if (j == 221)
	  c = 14;
	if (j == 249)
	  c = 8;
	if (j == 275)
	  c = 1;

	m_colId[j] = c;
	d[j] = RandomFloat(30.);
	for (int x=0; x<spikes.isize(); x++) {
	  double plus = spikes[x] - j;
	  if (plus < 0)
	    plus = -plus;
	  plus = 100./(plus + 1);
	  d[j] += plus;
	}
	cout << j << "\t" << d[j] << endl;
      }
    }
  }

  int Num() const {return m_data[0].isize();}
  int Total() const {return m_data.isize();}
  double Get(int x, int y) const {
    return (m_data[x])[y];
  }
  
  color GetColor(int i) const {
    return MakeUpColor(m_colId[i]);
  }

private:
  svec< svec < double > > m_data;
  svec<int> m_colId;
};

void OneFrame(const string & o, const Data & c, double phi, double theta)
{
  double x_offset = 20;
  double y_offset = 20;


  int i, j, k;
  ns_whiteboard::whiteboard board;

 
  double x_max = 300.;
  double y_max = 300.;


  //double x = 0.;
  //double y = 0.;
  board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(0, 0), 
  				      ns_whiteboard::xy_coords(x_max + 2*x_offset, y_max + 2*y_offset),
  				      color(0.99,0.99,0.99)));
  
  Geometry3D g;
  g.SetOffset(120.);
  g.SetRotation(phi, theta);
  //g.SetDist(160.);

  double d = 300.;

  
  for (k=0; k<c.Total(); k++) {
    for (j=0; j<c.Num(); j++) {
      double v = c.Get(k, j);
      color col = c.GetColor(j);
	
      double x = 4*k;
      double y = j;
      double z = v;
	
      Box box(1., col);
      double w = 2.;
      box.Draw(board, 
	       g.Coords(x, y, z),
	       g.Coords(x, y, z+w),
	       g.Coords(x, y+w, z+w),
	       g.Coords(x, y+w, z));
      
      //board.Add( new ns_whiteboard::rect( g.Coords(x, y, z), 
      //				    g.Coords(x+3, y+3, z+3),
      //				    col));
      
    }
  }
  

  board.Add( new ns_whiteboard::line( g.Coords(0, 0, 0), 
                                      g.Coords(d, 0, 0),
                                      1.,
                                      color(0.5,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(0, 0, 0), 
                                      g.Coords(0, d, 0),
                                      1.,
                                      color(0.,0.5,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(0, 0, 0), 
                                      g.Coords(0, 0, d),
                                      1.,
                                      color(0.,0.,0.5)));

  ofstream out(o.c_str());
  
  ns_whiteboard::ps_display display(out, x_max + 2 * x_offset, y_max + 2 * y_offset);
  board.DisplayOn(&display);
 
}


int main( int argc, char** argv )
{
 
  //commandArg<string> iStringO("-i","input");
  commandArg<string> aStringO("-o","outfile (post-script)");
 
  
  commandLineParser P(argc,argv);
  P.SetDescription("Animation example");

  //P.registerArg(iStringO);
  P.registerArg(aStringO);

  P.parse();

  string o = P.GetStringValueFor(aStringO);
  //string in = P.GetStringValueFor(iStringO);
   
  Data c;
  c.Generate();
 
  
  //OneFrame(o, c, -0.3, 0.2);
  //return 0;

  int i;
  int k = 1000;
 
  double phi = 0.;
  for (i=0; i<180; i++) {
    char name[256];
    sprintf(name, "%s%d.ps", o.c_str(), k);
    
    phi = 2*3.1415/180.*(double)i;
    
    OneFrame(name, c, phi, phi);
    k++;
  }
 
  return 0;
}
