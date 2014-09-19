#include "visual/Whiteboard.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "visual/Axes.h"
#include "visual/Geometry.h"

#include <iostream>




void OneFrame(const string & o, double x, double y)
{
  double x_offset = 20;
  double y_offset = 20;


  int i, j;
  ns_whiteboard::whiteboard board;

 
  double x_max = 300.;
  double y_max = 300.;

  //double x = 0.;
  //double y = 0.;
  board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(0, 0), 
  				      ns_whiteboard::xy_coords(x_max + 2*x_offset, y_max + 2*y_offset),
  				      color(0.99,0.99,0.99)));
  
  Geometry3D g;
  g.SetOffset(150.);
  g.SetRotation(x, y);

  double d = 100.;
  y = 0.;


  board.Add( new ns_whiteboard::line( g.Coords(0, 0, 0), 
				      g.Coords(d, 0, 0),
				      1.,
				      color(0.99,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(0, 0, 0), 
				      g.Coords(0, d, 0),
				      1.,
				      color(0.,0.99,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(0, 0, 0), 
				      g.Coords(0, 0, d),
				      1.,
				      color(0.,0.,0.99)));



  /*
  board.Add( new ns_whiteboard::line( g.Coords(0, y, 0), 
				      g.Coords(0, y, d),
				      1.,
				      color(0.,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(0, y, d), 
				      g.Coords(d, y+30, d),
				      1.,
				      color(0.,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(d, y+30, d), 
				      g.Coords(d, y+30, 0),
				      1.,
				      color(0.,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(d, y+30, 0), 
				      g.Coords(0, y, 0),
				      1.,
				      color(0.,0.,0.)));
  */
  /*
  y = 100.;
  board.Add( new ns_whiteboard::line( g.Coords(0, y, 0), 
				      g.Coords(0, y, d),
				      1.,
				      color(0.,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(0, y, d), 
				      g.Coords(d, y+30, d),
				      1.,
				      color(0.,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(d, y+30, d), 
				      g.Coords(d, y+30, 0),
				      1.,
				      color(0.,0.,0.)));

  board.Add( new ns_whiteboard::line( g.Coords(d, y+30, 0), 
				      g.Coords(0, y, 0),
				      1.,
				      color(0.,0.,0.)));

  */
  /*
  for (int i=0; i<10; i++) {
    double x = i * 20;
    double y = i * 20;
    double z = 0;
    board.Add( new ns_whiteboard::rect( g.Coords(x, y, z), 
					g.Coords(x+12, y, z+12),
					color(0.,0.,0.)));
  
					}*/


  ofstream out(o.c_str());
  
  ns_whiteboard::ps_display display(out, x_max + 2 * x_offset, y_max + 2 * y_offset);
  board.DisplayOn(&display);
 
}

int main( int argc, char** argv )
{
 
  commandArg<string> aStringO("-o","outfile (post-script)");
 
  
  commandLineParser P(argc,argv);
  P.SetDescription("Animation example");

  P.registerArg(aStringO);

  P.parse();

  string o = P.GetStringValueFor(aStringO);
   
  
  //OneFrame(o, 0, 0);

  int i;
  int k = 1000;
  double x = 0.;
  double y = 0.;

  for (i=0; i<150; i++) {
     char name[256];
    sprintf(name, "%s%d.ps", o.c_str(), k);
    
    x = 0.05*i;
    y = 0.03*i;
    
    OneFrame(name, x, y);
    k++;
  }
 
  return 0;
}
