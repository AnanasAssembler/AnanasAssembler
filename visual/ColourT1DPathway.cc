#include "visual/Whiteboard.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "visual/Axes.h"

#include <iostream>


class Value
{
public:
  Value() {
    m_a = 0.;
    m_b = 0.;
    m_c = 0.;
    m_d = 0.;
    m_e = 0.;
    m_f = 0.;
    m_g = 0.;
    m_h = 0.;
    m_i = 0.; 
  }

  void operator += (const Value & v) {

    m_a += v.m_a;
    m_b += v.m_b;
    m_c += v.m_c;
    m_d += v.m_d;
    m_e += v.m_e;
    m_f += v.m_f;
    m_g += v.m_g;
    m_h += v.m_h;
    m_i += v.m_i; 

  }

  void Set(const string & gene, double a, double b, double d, double e, double f, 
	   double g, double h, double i, double c) {
    /*m_a = a * 1.15839;
    m_b = b * 1.14821;
    m_c = c;
    m_d = d * 1.67326;
    m_e = e * 0.98413;
    m_f = f * 0.845145;
    m_g = g * 2.26896;
    m_h = h * 2.34796;
    m_i = i * 1.28959;*/
    //cout << gene << " " << 
    m_a = a;
    m_b = b;
    m_c = c;
    m_d = d;
    m_e = e;
    m_f = f;
    m_g = g;
    m_h = h;
    m_i = i; 


    m_gene = gene;
  }
  double V1() const {return m_a;}
  double VNorm() const {return m_c;}
  double V2() const {return m_b;}
  double V3() const {return m_d;}
  double V4() const {return m_e;}
  double V5() const {return m_f;}
  double V6() const {return m_g;}
  double V7() const {return m_h;}
  double V8() const {return m_i;}

  double Diff() const {
    if (m_a < 0.5 && m_c < 0.5)
      return 0;
    return m_a - m_c;
  }

  const string & Gene() {return m_gene;}

  void Print() {
    cout << VNorm() << "\t" << V1() << "\t" << V2() << "\t" << V3() << "\t" << V4() << "\t";
    cout << V5() << "\t" << V6() << "\t" << V7() << "\t" << V8() << endl;
   
  }

private:
  double m_a;
  double m_b;
  double m_c;
  double m_d;
  double m_e;
  double m_f;
  double m_g;
  double m_h;
  double m_i;
  string m_gene;
};


void Load(svec<Value> & out) {
  FlatFileParser parser;

  parser.Open("../ugc_572_expr_rpkm_update.txt");
  

  //parser.ParseLine();
  while(parser.ParseLine()) {
    if (parser.GetItemCount() < 1)
      continue;


      
    Value v;
    //v.Set(parser.AsString(1), parser.AsFloat(9), parser.AsFloat(2), parser.AsFloat(6), parser.AsFloat(8),
    //	  (parser.AsFloat(10) + parser.AsFloat(11) + parser.AsFloat(12))/3);
    
    v.Set(parser.AsString(1), parser.AsFloat(9), parser.AsFloat(3), parser.AsFloat(4), parser.AsFloat(8), 
	  parser.AsFloat(2), parser.AsFloat(7), parser.AsFloat(6), parser.AsFloat(5),
	  (parser.AsFloat(10) + parser.AsFloat(11) + parser.AsFloat(12))/3);
    out.push_back(v);
  } 
 
}


int main( int argc, char** argv )
{
  //No call to RunTime() made in order to allow clean ctrl-C exit.
 
  commandArg<string> aStringI1("-i","Pathway file");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<bool> bwCmd("-bw","grey-scale only", 0);
  commandArg<bool> vCmd("-v","print out raw values", false);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Heat map plotter");

  P.registerArg(aStringI1);
  P.registerArg(aStringO);
  P.registerArg(bwCmd);
  P.registerArg(vCmd);

  P.parse();

  string in = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(aStringO);
  bool bBW = P.GetBoolValueFor(bwCmd);
  bool bV = P.GetBoolValueFor(vCmd);
  
  double x_offset = 20;
  double y_offset = 20;

  svec<Value> all;

  Load(all);

  int i, j;
  ns_whiteboard::whiteboard board;

  FlatFileParser parser;
  
  parser.Open(in);
  

  double x_max = 0;
  double y_max = 0;

  double space = 50;
  double dot = 8;

  double x = 0;
  double y = 0;
  double div = 12.;

  double floor = 2.0;

  //parser.ParseLine();
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    
    Value v;

    int count = 0;
    for (int i=0; i<all.isize(); i++) {
      if (all[i].Gene() == parser.AsString(0)) {
	v += all[i];
	count++;
	//break;
      }
    }
    //if (count > 1) {
    //  cout << parser.AsString(0) << "\t" << count << endl;
    //}
     
    y += 10.;
    if (y > y_max)
      y_max = y;


    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset, y_offset + y + dot/4),
					parser.AsString(0), black, 6., "Times-Roman", 0, true));

    x = x_max = 20;
    x += 15;

    double r = 0.;
    double g = 0.;
    double b = 0.;

    double v1, v2;
    double diff;


    //v.Print();

    v1 = v.V1() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << "\t";
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }

    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;
    //----------------------------------------------------------------------------------

    v1 = v.V2() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << "\t";
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;


   //----------------------------------------------------------------------------------

    v1 = v.V3() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << "\t";
     if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;

    //----------------------------------------------------------------------------------

    v1 = v.V4() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << "\t";
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;

   //----------------------------------------------------------------------------------

    v1 = v.V5() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << endl;
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;

   //----------------------------------------------------------------------------------

    v1 = v.V6() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << endl;
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
   if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;
  
   //----------------------------------------------------------------------------------

    v1 = v.V7() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << endl;
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;
   //----------------------------------------------------------------------------------

    v1 = v.V8() + floor;
    v2 = v.VNorm() + floor;
 
    diff = v1 - v2;
    //cout << diff << endl;
    if (diff > 0) {
      diff = v1/v2 - 1.;
      r = diff / div;
      if (r > 0.999)
	r = 0.999;
      g = 1. - r;
      b = 1. - r;
      r = 0.999;
    } else {
      diff = v2/v1 - 1.;
      b = diff / div;
      if (b > 0.999)
	b = 0.999;
      g = 1. - b;
      r = 1. - b;
      b = 0.999;
    }	
    //color black_c(r, g, b);
	//cout << "r=" << r << endl;
	//cout << "Adding at " << localX << endl;
    if (v1 - floor < 0.01) {
      r = g = b = 0.7;
    }
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset, y + y_offset + dot), 
					ns_whiteboard::xy_coords(x + x_offset + dot, y + y_offset),
					color(r, g, b)) );

    x += dot + 5;
  } 


  for (i=-40; i<=40; i++) {
    double diff;
    double r, g, b;
    
    if (i > 0) {
      diff = (double)i/10.;
      r = diff / div;
      if (r > 0.99)
	r = 0.99;
      g = 1. - r;
      b = 1. - r;
      r = 0.99;
    } else {
      diff = (double)-i/10.;
      b = diff / div;
      if (b > 0.99)
	b = 0.99;
      g = 1. - b;
      r = 1. - b;
      b = 0.99;
    }	
    color black_c(r, g, b);
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + i + 80, y_offset/2 + dot), 
					ns_whiteboard::xy_coords(x_offset + 1. + i + 80, y_offset/2),
					black_c) );
    if (i > x_max)
      x_max = i;
  }
  //cout << "xmax=" << x_max << "  ymax=" << y_max << endl;
  ofstream out(o.c_str());
  
  ns_whiteboard::ps_display display(out, x_max + 4 * x_offset + space, y_max + 2 * y_offset + space);
  board.DisplayOn(&display);
 

  return 0;
}
