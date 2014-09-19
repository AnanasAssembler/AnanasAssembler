#define FORCE_DEBUG

#include "visual/Whiteboard.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "visual/Axes.h"

#include <iostream>


string Convert(const string & s) 
{
  //if (strlen(s.c_str()) == 12)
  //continue;
  string o;
  int i;
  for (i=0; i<7; i++) {
    o += s[i];
  }
  for (i=strlen(s.c_str()); i<12; i++) {
    o += '0';
  }
  
  for (i=7; i<strlen(s.c_str()); i++) {
    o += s[i];
  }
  // cout << o << " " << s << endl;
  return o;
}

int Pos(const string & s) 
{
  string o;
  int i;
  for (i=7; i<strlen(s.c_str()); i++) {
    o += s[i];
  }
   
  return atol(o.c_str());
}

bool Comp(const string & a, const string & b) 
{
  for (int i=0; i<strlen(a.c_str()); i++) {
    if (a[i] != b[i])
      return false;
    if (a[i] == '_')
      return true;
    
  }
  return true;
}

int Diff(int a, int b)
{
  int c = a - b;
  if (c < 0)
    c = -c;
  return c;
}

int main( int argc, char** argv )
{
  //No call to RunTime() made in order to allow clean ctrl-C exit.
 
  commandArg<string> aStringI1("-i","ERV file");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<bool> bwCmd("-n","do not filter", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Dot plot plotter for ERVs");

  P.registerArg(aStringI1);
  P.registerArg(aStringO);
  P.registerArg(bwCmd);

  P.parse();

  string in = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(aStringO);
  bool bBW = P.GetBoolValueFor(bwCmd);
  
  double x_offset = 20;
  double y_offset = 20;


  int i, j;
  ns_whiteboard::whiteboard board;

  FlatFileParser parser;
  
  parser.Open(in);
  

  double x_max = 0;
  double y_max = 0;

  double space = 50;
  double dot = 5;

  double x = 0;
  double y = 0;

  svec<string> one, two;

  svec<int> rank;
  svec<int> synt;
  //  svec<double> line;
  int lastPos = 0;
  string last;

  svec<string> inverse;

  while(parser.ParseLine()) {
    if (parser.GetItemCount() < 1)
      continue;

    string a = Convert(parser.AsString(0));
    string b = Convert(parser.AsString(1));

    string in = a + "," + b;
    inverse.push_back(in);

    //cout << parser.AsString(0) << " -> " << a << endl;
    int pos = Pos(parser.AsString(1));
    //cout << "Testing " << a << " " << last << endl;
    if (!bBW && a == last) {
      int diff = pos - lastPos;
      if (diff < 0)
	diff = -diff;
      //cout << pos << " " << lastPos << " " << b << " diff " << diff << endl;
      
      if (diff > 2) {
	continue;
      }
    }
    //cout << pos << endl;
    lastPos = pos;

    //if (!bBW)
    //cout << a << "\t" << b << endl;

    one.push_back(a);
    two.push_back(b);
    rank.push_back(pos);
    //if (pos == 5484) {
    //cout << "DEBUG " << a << " " << b << endl; 
    //}

    last = a;
  }

  Sort(inverse);
  
  synt.resize(rank.isize(), 0);
  for (i=1; i<rank.isize()-1; i++) {

    string in = two[i] + "," + one[i];
    if (BinSearch(inverse, in) < 0)
      continue;
    
    //if (one[i] == "gorGor_00034" && two[i]== "homSap_04916") {
    //  cout << "DEBUG: " << rank[i] << " " << rank[i-1] << " " << rank[i+1] << endl;
    //}

    //cout << rank[i] << endl;
    if (rank[i-1] < rank[i] && rank[i] < rank[i+1]
	|| rank[i-1] > rank[i] && rank[i] > rank[i+1]) {
      int max = 30;
      if (Diff(rank[i-1], rank[i]) < max && Diff(rank[i+1], rank[i]) < max) {      
	synt[i-1] = 1;
	synt[i] = 1;
	synt[i+1] = 1;
      }
    }
  }
  
  for (i=1; i<rank.isize()-1; i++) {
   
    //if (one[i] == "gorGor_00034" && two[i]== "homSap_04916") {
    //  cout << "DEBUG: " << rank[i] << " " << rank[i-1] << " " << rank[i+1] << endl;
    //}

    //cout << rank[i] << endl;
    if (rank[i-1] < rank[i] && rank[i] < rank[i+1]
	|| rank[i-1] > rank[i] && rank[i] > rank[i+1]) {
      if (synt[i-1] >= 1 && synt[i] >= 1 && synt[i+1] >= 1) {
	//synt[i-1] = 2;
	synt[i] = 2;
	//synt[i+1] = 2;
      }
    }
  }

  svec<string> first, second;
  svec<string> uni;

  for (i=0; i<one.isize(); i++) {
    if (synt[i] > 0) {
      cout << one[i] << "\t" << two[i] << endl;
      first.push_back(one[i]);
      second.push_back(two[i]);
      uni.push_back(one[i]);
      uni.push_back(two[i]);
    }
  }

  



  UniqueSort(uni);

  double scale = 0.05;
  dot = 0.5;

  svec<int> hits;
  hits.resize(uni.isize(), 0);

  //cout << "Assigning hits" << endl;
  for (i=0; i<first.isize(); i++) {
    int index = BinSearch(uni, first[i]); 
    hits[index]++;
  }
  //cout << "done." << endl;

  for (i=0; i<first.isize(); i++) {
    int index1 = BinSearch(uni, first[i]); 
    int index2 = BinSearch(uni, second[i]); 

    //cout << "i=" << i << " " << index1 << " " << index2;
    //cout << " " << first[i] << " " << second[i] << endl;
    

    //if (hits[index1] > 30 || hits[index2] > 30)
    //continue;

    x = scale * (double)index1;
    y = scale * (double)index2;

    if (x > x_max)
      x_max = x;
    if (y > y_max)
      y_max = y;

    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + x, y_offset + y + dot), 
					ns_whiteboard::xy_coords(x_offset + x + dot, y_offset + y),
					color(0.99, 0., 0.)));
  }
  for (i=1; i<uni.isize(); i++) {
    
    if (!Comp(uni[i], uni[i-1])) {
      
      double z  = scale * (double)i;      

      board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x_offset, y_offset + z), 
					  ns_whiteboard::xy_coords(x_offset + x_max, y_offset + z),
					  1.0, color(0., 0., 0.)));
      board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x_offset + z, y_offset), 
					  ns_whiteboard::xy_coords(x_offset + z, y_offset + y_max),
					  1.0, color(0., 0., 0.)));
    }
  }

  //cout << "xmax=" << x_max << "  ymax=" << y_max << endl;
  ofstream out(o.c_str());
  
  ns_whiteboard::ps_display display(out, x_max + 2 * x_offset + space, y_max + 2 * y_offset + space);
  board.DisplayOn(&display);
 

  return 0;
}
