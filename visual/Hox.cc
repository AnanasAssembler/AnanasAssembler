#include "visual/Whiteboard.h"
#include "base/FileParser.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "visual/Compounds.h"

#include <iostream>


void DrawOne(const string & o, int index, const svec< svec < double > > & val, const svec<string> & name)
{
  double x_offset = 20;
  double y_offset = 20;


  int i, j;
  ns_whiteboard::whiteboard board;

 
  double x_max = 550. + 2* x_offset;
  double y_max = 220.;

    
  Box b;

  double x = 0.;
  double y = 0.;
  board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(0, 0), 
				      ns_whiteboard::xy_coords(x_max + 2*x_offset, y_max + 2*y_offset),
				      color(0.99, 0.99, 0.99)) );
  
  for (i=0; i<name.isize(); i++) {
    if (name[i] == "delim") {
      x = 0.;
      y += 30.;
      continue;
    }
    if (name[i] == "empty") {
      x += 45.;
      continue;
    }
     

   const svec<double> & vv = val[i];
    double mean = 0.;
    for (j=0; j<vv.isize(); j++)
      mean += vv[j];
    mean /= (double)vv.isize();

    double v = (vv[index] - mean) / 5.;
    
    color cc = GradientMult(v, color(0., 0., 0.99), color(0.99, 0., 0.), color(0.99, 0.99, 0.99));
    
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y_offset + y + 18), 
					ns_whiteboard::xy_coords(x + x_offset + 36, y + y_offset),
					cc) );


    b.Draw(board,
	   ns_whiteboard::xy_coords(x + x_offset, y_offset + y + 18), 
	   ns_whiteboard::xy_coords(x + x_offset + 36, y + y_offset));
    
    board.Add( new ns_whiteboard::text(  ns_whiteboard::xy_coords(x + x_offset + 3, y_offset + y + 5),
					name[i], black, 14., "Times-Roman", 0, true));
    
    x += 45.;
   }
  

  ofstream out(o.c_str());
  
  ns_whiteboard::ps_display display(out, x_max + 2 * x_offset, y_max + 2 * y_offset);
  board.DisplayOn(&display);
 
}

void Smooth(svec<double> & d, int mul) 
{
  int i, j;
  
  for (j=0; j<500; j++) {
    svec<double> tmp;
    tmp.resize(d.isize(), 0.);
    //cout << "Iter: " << j << endl;
    for (i=0; i<d.isize(); i++) {
      if (i == 0 || i+1 == d.isize() || i % mul == 0) {
	tmp[i] = d[i];
	//cout << i << " " << d[i] << " " << tmp[i] << endl;
	continue;
      }
      tmp[i] = (d[i] + d[i-1] + d[i+1])/3.;
      //cout << i << " " << d[i] << " " << tmp[i] << endl;
    }
    d = tmp;
  }
}

int main( int argc, char** argv )
{
 
  commandArg<string> iStringO("-i","input file");
  commandArg<string> aStringO("-o","outfile (post-script)");
 
  
  commandLineParser P(argc,argv);
  P.SetDescription("Color scale example");

  P.registerArg(iStringO);
  P.registerArg(aStringO);

  P.parse();

  string o = P.GetStringValueFor(aStringO);
  string inFile = P.GetStringValueFor(iStringO);
   

  FlatFileParser parser;
  
  int i, j;

  parser.Open(inFile);
  svec<string> name;
  svec< svec < double > > val;
  int mul = 20;
  int num = 0;
  while(parser.ParseLine()) {
    if (parser.GetItemCount() < 1) {
      name.push_back("delim");
      svec<double> dd;
      val.push_back(dd);
      continue;
    }
    if (parser.GetItemCount() == 1 && parser.AsString(0) == "-") {
      name.push_back("empty");
      svec<double> dd;
      val.push_back(dd);
      continue;
    }
    
    //string n = parser.AsString(18) + " ";
    //n += parser.AsString(19) + " ";
    //n += parser.AsString(20);
    string n = parser.AsString(20);
    cout << n << endl;
    name.push_back(n);
    svec<double> d;
    for (i=1; i<10; i++) {
      d.push_back(parser.AsFloat(i));
      for (j=0; j<mul-1; j++)
	d.push_back(0.);
    }
    
    Smooth(d, mul);
    num = d.isize();
    val.push_back(d);
  }
 
  int k = 1000;
  for (i=0; i<num; i++) {
    char out[256];
    sprintf(out, "%s%d.ps", o.c_str(), k);
  
    DrawOne(out, i, val, name);
    k++;
  }



  return 0;
}
