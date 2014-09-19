#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>
#include "visual/Whiteboard.h"
#include "visual/Color.h"

#include "visual/Axes.h"
#include <iostream>


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> outCmmd("-o","output ps file");
  commandArg<bool> normCmmd("-norm","normalize", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Plots heat maps for gene families.");
  P.registerArg(fileCmmd);
  P.registerArg(outCmmd);
  P.registerArg(normCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  bool bNorm = P.GetBoolValueFor(normCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i;

  double x, y;
  double offset = 20.;

  x = offset;
  y = offset;
  double dot = 10.;
  double x_max = dot * 10. + 90.;
  
  ns_whiteboard::whiteboard board;


  while (parser.ParseLine()) {
    if (parser.GetItemCount() > 6) {
      const string & name = parser.AsString(parser.GetItemCount()-1);
      x = offset;

      double mul = 1.;
      if (bNorm) {
	double sum = 0;
	for (i=7; i<parser.GetItemCount()-1; i++) {
	  double d = parser.AsFloat(i);
	  sum += log(1.+d);   
	}
	if (sum > 0.01)
	  mul = 10. / sum;
      }
   
      for (i=7; i<parser.GetItemCount()-1; i++) {
	double d = parser.AsFloat(i);
	d = log(1.+d);      
	
	d *= mul;

	d /= 6.;
	if (d > 0.99)
	  d = 0.99;
	d = 1. - d;

	color col(d, d, d);
        board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x, y + dot), 
                                            ns_whiteboard::xy_coords(x + dot, y),
                                            col) );
	x += dot;
      }



      board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x + dot, y),
					  name, black, 8., "Times-Roman", 0, true));
      y += dot;
      
    } else {
      y += dot;
    }
  }



  ofstream out(outName.c_str());
  
  ns_whiteboard::ps_display display(out, x_max, y + offset);
  board.DisplayOn(&display);
 

  return 0;
}
