#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>
#include "visual/Whiteboard.h"
#include "visual/Color.h"

#include "visual/Axes.h"
#include <iostream>
#include "base/RandomStuff.h"

#include <math.h>

int main( int argc, char** argv )
{

  
  string outName = "saguaro.ps";

  //comment. ???
  //FlatFileParser parser;
  
  //parser.Open(fileName);
 

  double x, y;
  double offset = 20.;

  x = offset;
  y = offset;
  double dot = 3.;
  double x_max = 0;
  double y_max = 0;
  
  ns_whiteboard::whiteboard board;

  
  int i, j;
  double rad = 200.;
  for (j=0; j<360; j+= 5) {
    double r, g, b;
    r = g = b = 0.5;

    double phi = 2 * 3.14156*((double)j)/360.;
    x = offset + rad + rad*sin(phi);
    y = offset + rad + rad*cos(phi);
    //cout << x << " " << y << " " << phi << endl;

 
    double size = 2.;

    
    /*if (j > 120 && j < 240) {
      double a = j-180;
      if (a < 0)
	a = -a;
      size += 8. - a/10.;
      }*/

    if (j > 20 && j < 100) {
      double a = j-60;
      if (a < 0)
	a = -a;
      size += 0.9*(8. - a/10.);
      r = 0.6;
      g = b = 0.;
    }

    /*if (j > 300 && j < 360) {
      double a = j-330;
      if (a < 0)
	a = -a;
      size += 0.4*(8. - a/10.);
      b = 0.6;
      r = g = 0.;
      }*/
    color col(r, g, b); 

    if (x > x_max)
      x_max = x;
    if (y > y_max)
      y_max = y;

    board.Add( new ns_whiteboard::arc( ns_whiteboard::xy_coords(x, y), 
				       size, 0., 360, 2.0,
				       col) );
    
    
  }

  ofstream out(outName.c_str());
  
  ns_whiteboard::ps_display display(out, x_max+2*offset, y_max + 2*offset);
  board.DisplayOn(&display);
 

  return 0;
}
