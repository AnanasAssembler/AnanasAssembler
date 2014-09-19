
#define FORCE_DEBUG
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
  double VAll() const {
    return (m_a + m_b + m_d + m_e + m_f + m_g + m_h + m_i) / 8.;
  }
  double VAllM() const {
    return (m_a + m_b + m_d + m_f + m_g + m_h + m_i) / 7.;
  }
  double VNormDivid() const {
    return (m_b + m_d + m_e + m_g + m_h + m_i) / 6.;
  }

  double Diff() const {
    if (m_a < 0.5 && m_c < 0.5)
      return 0;
    return m_a - m_c;
  }

  const string & Gene() {return m_gene;}


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

bool Split(string & name, string & gene, const string & n) {
  int i;
  char tmp[256];
  strcpy(tmp, n.c_str());
  for (i=0; i<strlen(tmp); i++) {
    if (tmp[i] >= 'A' && tmp[i] <= 'Z') {
      char tmp1[256];
      strcpy(tmp1, &tmp[i]);
      tmp[i] = 0;
      name = tmp;
      gene = tmp1;
    }
  }
  return true;
}


class Block
{
public:
  Block() {
    m_x = 0.;
    m_y = 0.;
  }

  bool Set(const string & name, const string & gene, color col) {
    int i;
    if (m_name != "" && m_name != name)
      return false;
    m_name = name;
    m_gene.push_back(gene);
    m_col.push_back(col);
    return true;
  }

  const string & Name() const {return m_name;}
  int GeneCount() const {return m_gene.isize();}
  const string & Gene(int i) const {return m_gene[i];}
  color Color(int i) const {return m_col[i];}

  void SetCoords(double x, double y) {
    m_x = x;
    m_y = y;
  }

  double X() const {return m_x;}
  double Y() const {return m_y;}

private:
  svec<string> m_gene;
  string m_name;
  svec<color> m_col; 
  double m_x;
  double m_y;

};

void LoadPS(svec<Block> & blocks, string in) {

  FlatFileParser parser;
  int i;
  parser.Open(in);
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.GetItemCount() == 2 && parser.AsString(1) == "newpath" && (parser.AsString(0))[0] == '(') {
      char tmp[256];
      const char * p = parser.AsString(0).c_str();
      strcpy(tmp, &p[1]);
      tmp[strlen(tmp)-1] = 0;
      parser.ParseLine();

      for (i=0; i<blocks.isize(); i++) {
	if (blocks[i].Name() == tmp) {
	  if (blocks[i].X() == 0 && blocks[i].Y() == 0) {
	    blocks[i].SetCoords(parser.AsFloat(0), parser.AsFloat(1));
	    cout << "Set block " << blocks[i].Name() << endl;
	  } else {
	    cout << "Copy block " << blocks[i].Name() << endl;
	    blocks.push_back(blocks[i]);
	    blocks[blocks.isize()-1].SetCoords(parser.AsFloat(0), parser.AsFloat(1));
	  }
	  break;
	}
      }
    }
  }
  cout << "SVN Loaded!!!" << endl;
}


//    transform="translate(276.86441,17.199153)"><text
//         id="text4594"
//         transform="matrix(1,0,0,-1,21,368)"><tspan
//           id="tspan4596"
//           sodipodi:role="line"
//           y="0"
//           x="0 3.0580001 8.5579996"
//           style="font-size:11px;font-variant:normal;font-weight:normal;writing-mode:lr-tb;fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none;font-family:Times New Roman;-inkscape-font-specification:Times-Roman">ins</tspan></text>


bool GetText(string & text, const string & s) 
{
  const char * p = s.c_str();
  if (strstr(p, "style") == NULL) {
    return false;
  }
  
  StringParser parser;
  parser.SetLine(s, ">");
  char tmp[1024];
  strcpy(tmp, parser.AsString(1).c_str());
  for (int i=0; i<strlen(tmp); i++) {
    if (tmp[i] == '<') {
      tmp[i] = 0;
      break;
    }
  }
  text = tmp;
  return true;
}

bool GetXY(double & x, double & y, const string & s)
{
  //x = y = 0.;
  //cout << "Enter GetXY" << endl;
  const char * p = s.c_str();
  if (strstr(p, "transform") == NULL) {
    //cout << "Exit GetXY" << endl;
   return false;
  }
  char tmp[2048];
  int i;
  for (i=0; i<strlen(p); i++) {
    if (p[i] == '(') {
      strcpy(tmp, &p[i+1]);
      break;
    }      
  }
  for (i=0; i<strlen(tmp); i++) {
    if (tmp[i] == ')') {
      tmp[i] = 0;
      break;
    }          
  }
  StringParser parser;
  parser.SetLine(tmp, ",");
  x += parser.AsFloat(parser.GetItemCount()-2);

  double yc = 1.;
  if (parser.GetItemCount() == 6) {
    cout << "CORRECTING: " << endl;
    yc = -parser.AsFloat(3);
  }

  y += yc*parser.AsFloat(parser.GetItemCount()-1);
  cout << "PARSING " << s << "\t" << x << "\t" << y << endl;
  //cout << "Exit GetXY" << endl;
  return true;
}

void LoadSVG(svec<Block> & blocks, string in) {

  FlatFileParser parser;
  int i;
  parser.Open(in);
  double x = 0.;
  double y = 0.;
   while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (GetXY(x, y, parser.Line())) {
      x = y = 0.;
      cout << "real!" << endl;
      GetXY(x, y, parser.Line());
    }
    if (strstr(parser.Line().c_str(), "<text") != NULL) {
      string name;
      if (GetXY(x, y, parser.Line())) {
	x = y = 0.;
	//cout << "real!" << endl;
	GetXY(x, y, parser.Line());	
      }
      while(!GetText(name, parser.Line())) {
	parser.ParseLine();
	GetXY(x, y, parser.Line());	
      }
      cout << "Found: " << name << "\t" << x << "\t" << y << endl;
      for (i=0; i<blocks.isize(); i++) {
	if (blocks[i].Name() == name) {
	  if (blocks[i].X() == 0 && blocks[i].Y() == 0) {
	    blocks[i].SetCoords(x, y);
	    cout << "Set block " << blocks[i].Name() << endl;
	  } else {
	    cout << "Copy block " << blocks[i].Name() << endl;
	    blocks.push_back(blocks[i]);
	    blocks[blocks.isize()-1].SetCoords(x, y);
	  }
	  break;
	}
      }
      x = 0.;
      y = 0.;

    }
  }
   cout << "SVG Loaded!" << endl; 
   
}


int main( int argc, char** argv )
{
  //No call to RunTime() made in order to allow clean ctrl-C exit.
 
  commandArg<string> aStringI1("-i","Pathway file");
  commandArg<string> aStringP("-s","SVG format template (from inkscape)", "");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<bool> bwCmd("-bw","grey-scale only", 0);
  commandArg<int> iCmd("-index","sample index", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Heat map plotter");

  P.registerArg(aStringI1);
  P.registerArg(aStringP);
  P.registerArg(aStringO);
  P.registerArg(bwCmd);
  P.registerArg(iCmd);

  P.parse();

  string in = P.GetStringValueFor(aStringI1);
  string post = P.GetStringValueFor(aStringP);
  string o = P.GetStringValueFor(aStringO);
  bool bBW = P.GetBoolValueFor(bwCmd);
  int index = P.GetIntValueFor(iCmd);
  
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
  double dotX = 48;
  double dotY = 16;

  double x = 0;
  double y = 0;

  //for individuals
  // double div = 9.;
  // For all
  double div = 5.;

  double floor = 2.0;

  double xCoordScale = 60.;
  double yCoordScale = 20.;

 

  //parser.ParseLine();
  svec<Block> blocks;

  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    
    //for (j=0; j<parser.GetItemCount(); j++) {
    bool bb = false;
    
    string name, gene;
    //Split(name, gene, parser.AsString(j));
    name = parser.AsString(1);
    gene = parser.AsString(0);
    
    Value v;
    for (i=0; i<all.isize(); i++) {
      if (all[i].Gene() == gene) {
	v += all[i];
	//break;
      }
    }
    
    x = 0;
    
    double r = 0.;
    double g = 0.;
    double b = 0.;
    
    double v1, v2;
    double diff;
    
    
    //v.Print();
    
    v1 = v.VAll() + floor;
    //v1 = v.VAllM() + floor;

    switch(index) {
    case 0:
      v1 = v.V1();
      break;
    case 1:
      v1 = v.V2();
      break;
    case 2:
      v1 = v.V3();
      break;
    case 3:
      v1 = v.V4();
      break;
    case 4:
      v1 = v.V5();
      break;
    case 5:
      v1 = v.V6();
      break;
    case 6:
      v1 = v.V7();
      break;
    case 7:
      v1 = v.V8();
      break;
    case 8:
      v1 = v.VAll();
      break;
    case 9:
      v1 = v.VAllM();
      break;

    };

    v1 += floor;

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

    if (v1 - floor < 0.01 && v2 - floor < 0.01 ) {
      r = g = b = 0.7;
    }


    color col(r, g, b);
    
    
    for (i=0; i<blocks.isize(); i++) {
      if (blocks[i].Set(name, gene, col)) {
	bb = true;
	break;
      }
    }
    if (!bb) {
      Block tmp;
      tmp.Set(name, gene, col);
      blocks.push_back(tmp);
    }
  }
  
  cout << "Loaded " << blocks.isize() << " blocks." << endl;

  if (post != "") {
    LoadSVG(blocks, post);
  }

  for (j=0; j<blocks.isize(); j++) {
    x = 0.;
    y += yCoordScale;

    if (y + dotY > y_max)
      y_max = y + dotY;

    //cout << x << " " << y << endl;

    const Block & bl = blocks[j];

    if (bl.X() > 0. || bl.Y() > 0.) {
      x = bl.X();
      y = bl.Y();
    }

    double x1 = x;
    

    double lDot = dotX/bl.GeneCount();
    for (i=0; i<bl.GeneCount(); i++) {
      cout << bl.Name() << "\t" << bl.Gene(i) << endl;
      board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x + x_offset + lDot*i, y + y_offset + dotY), 
					  ns_whiteboard::xy_coords(x + x_offset + lDot*(i+1), y + y_offset),
					  bl.Color(i)) );
    }
    x += dotX;
    //----------------------------------------------------------------------------------
    //    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset, y_offset + y + dot/4),
    //                                  parser.AsString(0), black, 6., "Times-Roman", 0, true));

    double fontsize = 14.;
    double shadowX = 0.05;
    double shadowY = 0.05;
    int mm = 10;
    /*for (int xx = 1; xx <= mm; xx++) {
      for (int yy = -mm; yy <0; yy++) {
      
	board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2. + shadowX*(double)xx, y_offset + y + dotY/2 - 3 + shadowY*(double)yy),
					    bl.Name(), color(0.99,0.99,0.99), fontsize, "Times-Roman", 0, true));
      }
      }*/

    /*board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2. + shadowX, y_offset + y + dotY/2 - 3 - shadowY),
					bl.Name(), color(0.99,0.99,0.99), fontsize, "Times-Roman", 0, true));
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2. - shadowX, y_offset + y + dotY/2 - 3 + shadowY),
					bl.Name(), color(0.99,0.99,0.99), fontsize, "Times-Roman", 0, true));
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2. - shadowX, y_offset + y + dotY/2 - 3 - shadowY),
					bl.Name(), color(0.99,0.99,0.99), fontsize, "Times-Roman", 0, true));
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2. + shadowX, y_offset + y + dotY/2 - 3 + shadowY),
    bl.Name(), color(0.99,0.99,0.99), fontsize, "Times-Roman", 0, true));*/

    // board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2. + 2*shadowX, y_offset + y + dotY/2 - 3 - 2*shadowY),
    //					bl.Name(), color(0.99,0.99,0.99), fontsize, "Times-Roman", 0, true));

    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x1 + 2., y_offset + y + dotY/2 - 3),
					bl.Name(), black, fontsize, "Times-Roman", 0, true));

    

    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x1 + x_offset - 0.5, y + y_offset-0.5), 
					ns_whiteboard::xy_coords(x + x_offset + 0.5, y + y_offset-0.5), 1.0,
					color(0, 0, 0)) );
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x1 + x_offset -0.5, y + y_offset + dotY + 0.5), 
					ns_whiteboard::xy_coords(x + x_offset + 0.5, y + y_offset + dotY + 0.5), 1.0,
					color(0, 0, 0)) );
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x + x_offset + 0.5, y + y_offset + dotY + 0.5), 
					ns_whiteboard::xy_coords(x + x_offset + 0.5, y + y_offset - 0.5), 1.0,
					color(0, 0, 0)) );
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x1 + x_offset - 0.5, y + y_offset + dotY + 0.5), 
					ns_whiteboard::xy_coords(x1 + x_offset - 0.5, y + y_offset - 0.5), 1.0,
					color(0, 0, 0)) );


    if (x > x_max)
      x_max = x;
    //cout << x << "\t" << x_max << endl;
  } 

  /*
  for (i=-10*div; i<=10*div; i++) {
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
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + i + 10*div, y_offset/2 + dotY), 
					ns_whiteboard::xy_coords(x_offset + 1. + i + 10*div, y_offset/2),
					black_c) );
    if (i > x_max)
      x_max = i;
      }*/

  //cout << "xmax=" << x_max << "  ymax=" << y_max << endl;
  ofstream out(o.c_str());
  cout << "MAX: " <<  x_max + 4 * x_offset + space << "\t" <<  y_max + 2 * y_offset + space << endl;
  ns_whiteboard::ps_display display(out, x_max + 4 * x_offset + space, y_max + 2 * y_offset + space);
  board.DisplayOn(&display);
 

  return 0;
}
