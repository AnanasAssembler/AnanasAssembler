#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>
#include "visual/Whiteboard.h"
#include "visual/Color.h"

#include "visual/Axes.h"
#include <iostream>
#include "base/RandomStuff.h"

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
  for (j=0; j<400; j++) {
    for (i=0; i<20; i++) {
      x = offset + j*dot;
      y = offset + i*dot;

      if (y > y_max)
	y_max = y;
      if (x > x_max)
	x_max = x;
      
 
      double f = RandomFloat(1.);

      double d = 0.9999;

      double thresh = 0.35 + 0.02*((double)i);

      if (f < thresh) {
	d = 0.5;
      }
      double r, g, b;
      r = d;
      g = d;
      b = d;

      if ((j >110 && j < 140) || (j > 283 && j < 329)) {
	f = RandomFloat(1.);
	if (d < 0.8) {
	  r = 0.6;
	  g = 0.;
	  b = 0.;
	}
	if (f > 0.3) {
	  if (i < 10) {
	    r = 0.6;
	    g = 0.;
	    b = 0.;
	  } else {
	    r = 0.9999;
	    g = 0.9999;
	    b = 0.9999;
	  }
	}
	
      }
      
      if (j >200 && j < 230) {
	f = RandomFloat(1.);
	if (d < 0.8) {
	  r = 0.;
	  g = 0.;
	  b = 0.6;
	}
	if (f > 0.3) {
	  if (i % 2 == 0) {
	    r = 0.;
	    g = 0.;
	    b = 0.6;	    	    
	  } else {
	    r = 0.9999;
	    g = 0.9999;
	    b = 0.9999;	    
	  }
	}	
      }
	
	
      color col(r, g, b); 
      board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x, y + dot), 
					  ns_whiteboard::xy_coords(x + dot, y),
					  col) );

    }
  }

  ofstream out(outName.c_str());
  
  ns_whiteboard::ps_display display(out, x_max+2*offset, y_max + 2*offset);
  board.DisplayOn(&display);
 

  return 0;
}
