#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>
#include "visual/Whiteboard.h"
#include "visual/Color.h"

#include "visual/Axes.h"
#include <iostream>


string Cut(const string & s) {
  char tmp[256];
  const char * p = s.c_str();
  strcpy(tmp, &p[4]);
  for (int i=0; i<strlen(tmp); i++) {
    if (tmp[i] == 'F')
      tmp[i] = 'M';
  }
    

  string out = tmp;
  return out;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd1("-i1","input file 1");
  commandArg<string> fileCmmd2("-i2","input file 2");
  commandArg<string> outCmmd("-o","output ps file");
  commandArg<bool> normCmmd("-norm","normalize", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Scatterplot with colors.");
  P.registerArg(fileCmmd1);
  P.registerArg(fileCmmd2);
  P.registerArg(outCmmd);
  P.registerArg(normCmmd);
  
  P.parse();
  
  string fileName1 = P.GetStringValueFor(fileCmmd1);
  string fileName2 = P.GetStringValueFor(fileCmmd2);
  string outName = P.GetStringValueFor(outCmmd);
  bool bNorm = P.GetBoolValueFor(normCmmd);
  

  //comment. ???
  FlatFileParser parser1, parser2;
  
  parser1.Open(fileName1);
  parser2.Open(fileName2);
  int i;

  double x, y;
  double offset = 20.;

  x = offset;
  y = offset;
  double dot = 3.;
  double x_max = 0.;
  double y_max = 0.;
  
  ns_whiteboard::whiteboard board;

  double scale_x = 500.;
  double scale_y = 500;


  while (parser1.ParseLine()) {
    parser2.ParseLine();
    if (parser1.GetItemCount() == 0)
      continue;
    if (parser1.AsString(3) == "-nan")
      continue;


    double x = offset + scale_x * parser1.AsFloat(3);
    double y = offset + scale_y * parser2.AsFloat(3);

    if (x > x_max)
      x_max = x;
    if (y > y_max)
      y_max = y;
 
    double r = 0.;
    double g = 0.;
    double b = 0.;
    bool same = false;
    if (Cut(parser1.AsString(0)) == Cut(parser1.AsString(1))) {
      same = true;
      if (x < y) {
	cout << parser1.Line() << endl;
 	cout << parser2.Line() << endl;
	cout << endl;
      }
   }
    bool ok = false;
    if (strstr(parser1.AsString(0).c_str(), "ts") != NULL
	|| strstr(parser1.AsString(1).c_str(), "ts") != NULL) {
      ok = true;
      if (same) {
	r = 0.9;
      } else {
	r = 0.8;
	g = 0.7;
	b = 0.0;
      }
    }
    if (!ok && same) {
      ok = true;
      if (strstr(parser1.AsString(0).c_str(), "kd") != NULL) {
	g = 0.9;
	//b = 0.5;
      }
      if (strstr(parser1.AsString(0).c_str(), "lv") != NULL) {
	b = 0.9;
	g = 0.3;
	r = 0.3;
      }
      if (strstr(parser1.AsString(0).c_str(), "ht") != NULL) {
	b = 0.5;
	g = 0.5;
	r = 0.9;
      }
    } 
    if (!ok && !same) {
      r = g = b = 0.4;
 
      if (!(strstr(parser1.AsString(0).c_str(), "br") != NULL && 
	    strstr(parser1.AsString(1).c_str(), "cb") != NULL) 
	  && !(strstr(parser1.AsString(0).c_str(), "cb") != NULL && 
	       strstr(parser1.AsString(1).c_str(), "br") != NULL)) {
	r = g = b = 0.7;
      }
    }

    color col(r, g, b);
    //board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x, y + dot), 
    //					ns_whiteboard::xy_coords(x + dot, y),
    //					col) );
    board.Add( new ns_whiteboard::arc( ns_whiteboard::xy_coords(x, y), 
				       2, 0., 360, 1.0,
				       col) );
  }
  x_max = scale_x * 0.8;
  y_max = scale_y * 0.8;
  
  double tt = 0.2 * scale_x;

  board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(offset + tt, offset + tt), 
				      ns_whiteboard::xy_coords(x_max+offset, offset + y_max),
				      1, color(0,0,0)) );
  
  board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(offset + tt, offset + tt), 
				      ns_whiteboard::xy_coords(x_max+offset, offset + tt),
				      1, color(0,0,0)) );
  board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(offset + tt, offset + tt), 
				      ns_whiteboard::xy_coords(offset +tt, y_max+offset),
				      1, color(0,0,0)) );

  for (i=2; i<=8; i++) {
    double ss = (i * scale_x)/10.;
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(offset + tt, offset + ss), 
					ns_whiteboard::xy_coords(offset + tt - 10, offset + ss),
					1, color(0,0,0)) );
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(offset + ss, offset + tt), 
					ns_whiteboard::xy_coords(offset + ss, offset + tt - 10),
					1, color(0,0,0)) );

  }


  ofstream out(outName.c_str());
  


  ns_whiteboard::ps_display display(out, x_max + 3*offset, y_max + 3*offset);
  board.DisplayOn(&display);
 

  return 0;
}
