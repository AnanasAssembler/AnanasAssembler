
#define FORCE_DEBUG
#include "visual/Whiteboard.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "visual/Axes.h"

#include <iostream>

class Read
{
public:
  Read() {}
  Read(int start, const string & name, const string & partner, const string & seq, const string & ori) {
    m_pos = start;
    m_name = name;
    m_partner = partner;
    m_seq = seq;
    m_ori = ori;
  }
  int Pos() const {return m_pos;}
  const string Name() const {return m_name;}
  const string Partner() const {return m_partner;}
  const string Seq() const {return m_seq;}
  const string Ori() {return m_ori;}

private:
  int m_pos;
  string m_name;
  string m_partner;
  string m_seq;
  string m_ori;
};

int main( int argc, char** argv )
{
  //No call to RunTime() made in order to allow clean ctrl-C exit.
 
  commandArg<string> aStringI1("-i","Layout file");
  commandArg<string> oStringI1("-o","Postscript file");
  commandArg<int> iCmd("-index","contig index", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Displays a transcript/contig");

  P.registerArg(aStringI1);
  P.registerArg(oStringI1);
  P.registerArg(iCmd);

  P.parse();

  string in = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(oStringI1);
  int index = P.GetIntValueFor(iCmd);
  
  double x_offset = 20;
  double y_offset = 20;

  int i, j;
  ns_whiteboard::whiteboard board;

  FlatFileParser parser;
  
  parser.Open(in);
  

  double x_max = 0;
  double y_max = 0;


  double x = 0;
  double y = 0;

  int k = 0;

  svec<Read> reads;

  bool bDo = false;
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "</CONTIG>") {
      if (bDo)
	break;
    }
  
    if (parser.AsString(0) == "<CONTIG>") {
      if (k == index) 
	bDo = true;
      k++;
      continue;
    }
    if (bDo) {
      reads.push_back(Read(parser.AsInt(2), 
			   parser.AsString(5), 
			   parser.AsString(6), 
			   parser.AsString(7), 
			   parser.AsString(1)));
    }
  }


  svec<int> used;
  svec<int> part;
  used.resize(reads.isize());
  part.resize(reads.isize());
  int nn = 0;
  int fw = 0;
  int rc = 0;
  int correct = 0;
  for (i=0; i<reads.isize(); i++) {
    for (j=0; j<reads.isize(); j++) {
      if (reads[i].Name() == reads[j].Partner()) {
	part[i] = 1;
	part[j] = 1;
	nn++;	
	if (reads[i].Ori() == "1")
	  fw++;
	if (reads[i].Ori() == "-1")
	  rc++;
	if (reads[i].Ori() != reads[j].Ori())
	  correct++;
      }
    }
  }
  cout << "Total: " << reads.isize();
  cout << " partnered: " << nn << " forward: " << fw << " reverse: " << rc << " correct: " << correct << endl;

  int d = reads.isize();

 

  double scale = 4.8;
  
  while (d > 0) {
   
    int last = -1;
    x = 0.;
    d = reads.isize();

    for (i=0; i<reads.isize(); i++) {
      if (used[i] > 0) {
	d--;
	continue;
      }
      // if (d == 2)
      //cout << i << " " << reads[i].Pos() << " " << last << endl;

      if (reads[i].Pos() > last) {

	last = reads[i].Pos() + strlen(reads[i].Seq().c_str()) + 2;
	double x = scale * reads[i].Pos();
	double x1 = x + 8 * strlen(reads[i].Seq().c_str());
	if (x1 > x_max)
	  x_max = x1;


	double r, g, b;
	r = g = b = 0.1;
	if (part[i] > 0) {
	  if (reads[i].Ori() == "1") {
	    r = 0.90;
	    g = 0.;
	    b = 0.;
	  }
	  if (reads[i].Ori() == "-1") {
	    r = 0.;
	    g = 0.;
	    b = 092.;
	  }
	}


	board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x, y_offset + y),
					    reads[i].Seq(), color(r, g, b), 8., "Courier", 0, true));
	used[i]++;
	d--;
      }
    }
    y += 9.;
    //cout << d << endl;
  } 

  y_max = y;
 
  ofstream out(o.c_str());
  //cout << "MAX: " <<  x_max + 4 * x_offset + space << "\t" <<  y_max + 2 * y_offset + space << endl;
  ns_whiteboard::ps_display display(out, x_max + x_offset, y_max + y_offset);
  board.DisplayOn(&display);
 

  return 0;
}
