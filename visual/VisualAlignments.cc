#include "visual/Whiteboard.h"

#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "src/SequenceMatch.h"
#include "visual/Axes.h"
#include "visual/Color.h"
#include "visual/VisualAlignments.h"
#include "base/RandomStuff.h"
#include "src/Smorgas/ProtNameCluster.h"

void Limit(double & d)
{
  if (d > 0.999)
    d = 0.999;
  if (d < 0)
    d = 0.;
}

TaxaColor::TaxaColor() {
  m_use.push_back(0);
  m_use.push_back(1);
  m_use.push_back(4);
  m_use.push_back(6);
  m_use.push_back(9);
  m_use.push_back(10);
  m_use.push_back(11);
  m_use.push_back(12);
  m_use.push_back(16);
  m_use.push_back(18);
  m_use.push_back(19);
  m_use.push_back(20);
  m_use.push_back(24);
  m_use.push_back(26);
  m_baseCol = color(0., 0., 0);;
}

void TaxaColor::Add(const string & taxa) {
  if (m_taxa.isize() >= m_use.isize()) {
    m_col.push_back(color(0.2, 0.2, 0.2));
  } else {
    m_col.push_back(MakeUpColor(m_use[m_taxa.isize()]));
  }
  m_taxa.push_back(taxa);
}

void TaxaColor::SetTo(const string & name)
{
  string tax;
  bool b = GetTaxonomy(tax, name);
  //cout << "Name: " << name << endl;
  //cout << "Taxonomy: " << tax << endl;
  if (!b) {
    m_baseCol = color(0.2, 0.2, 0.2);
    return;
  }
  
  StringParser p;
  p.SetLine(tax);
  m_curr = "Unclassified";
  int i, j;
  for (j=p.GetItemCount()-1; j>=0; j--) {
    //cout << j << endl;
    const string & lookup = p.AsString(j);
    for (i=0; i<m_taxa.isize(); i++) {
      //cout << lookup << " " << m_taxa[i] << endl; 
      if (m_taxa[i] == lookup) { 
	m_baseCol = m_col[i];
	//cout << "Found: " << lookup << endl;
	m_curr = lookup;
	return;
      }
    }
  }

  m_baseCol = color(0.2, 0.2, 0.2);
}


color TaxaColor::GetColor(double shade) {
  color c = m_baseCol;
  double r = 1. - (1. - c.R())*shade;
  double g = 1. - (1. - c.G())*shade;
  double b = 1. - (1. - c.B())*shade;    
  Limit(r);
  Limit(g);
  Limit(b);

  return color(r, g, b);
}
  

color TaxaColor::Get(const string & s) {
  int i;
  for (i=0; i<m_taxa.isize(); i++) {
    if (m_taxa[i] == s) 
      return m_col[i];
  }
  return color(0., 0., 0.);
}



void Parse(string & spec, const string & in) 
{
  StringParser a, b;
  a.SetLine(in, "_[");
  //a.SetLine(in, "Tax=");
  if (a.GetItemCount() < 2)
    return;
  //cout << a.AsString(0) << endl;
  string c = a.AsString(1);
  b.SetLine(c, "_");
  spec = b.AsString(0);
}

void Species::Read(FlatFileParser & parser) {
  //cout << parser.Line() << endl;
  m_name = parser.AsString(1);
  int i = 3;
  while (strstr(parser.AsString(i).c_str(), ";") == NULL && i < parser.GetItemCount()) {
    //cout << parser.AsString(i) << endl;
    i++;
  }
  if (i >= parser.GetItemCount()-1) {
    m_name = "";
    return;
  }
  m_tax = "";
  for (; i<parser.GetItemCount(); i++) {
    m_tax += parser.AsString(i);
    m_tax += " ";
  }
  //m_tax += parser.AsString(i+1);
  //if (parser.GetItemCount() > 2) {
  //  m_tax += " ";
  //  m_tax += parser.AsString(i+2);    
  //}
}

void TaxaColor::ReadTaxonomy(const string & file)
{
  FlatFileParser parser1;
  parser1.Open(file);
  
  while (parser1.ParseLine()) {
    if (parser1.GetItemCount() == 0)
      continue;
    Species s;
    if (parser1.GetItemCount() > 9) {
      s.Read(parser1);
      m_spec.push_back(s);
    }
  }
}

bool TaxaColor::GetTaxonomy(string & tax, const string & name)
{
  string s;
  Parse(s, name);
  
  for (int i=0; i<m_spec.isize(); i++) {
    if (m_spec[i].Name() == s) {
      tax = m_spec[i].Tax();
      return true;
    }
  }
  return false;
}


class OneAlign
{
public:
  OneAlign() {
    m_start = 0;
    m_end = 0;
  }

  OneAlign(int start, const string & seq, const string & ref, const string & rName) {
    m_ref = ref;
    m_seq = seq;
    m_start = start;
    m_name = rName;
    m_end = 0;
    ComputeIdent();
  }
  
  int Start() const {return m_start;}
  int End() const {return m_end;}
  const string & Name() const {return m_name;}
  int Len() const {return m_ident.isize();}
  double Val(int i) const {return m_ident[i];}

  
private:
  
  void ComputeIdent() {
    int i, j;
    int n = strlen(m_seq.c_str());
    m_end = m_start;
    int win = 3;
    svec<double> tmp;
    tmp.resize(n, 0.);
    int k = 0;
    //for (i=0; i<n; i++) {
    for (i=0; i<n; i++) {
      if (m_seq[i] != '-' /*&& m_ref[i] != '-'*/) {
        if (m_seq[i] == m_ref[i]) {
	  tmp[k] = 1.;
	  m_end++;
	}
	k++;
      }
    }
    m_ident.resize(k+1, 0.);
    for (i=win; i<m_ident.isize()-win; i++) {
      for (j=i-win; j<=i+win; j++) {
	m_ident[i] += tmp[j]/(2*(double)win+1.); 
      } 
    }
  }
  

  string m_ref;
  string m_seq;
  string m_name;
  svec<double> m_ident;
  int m_start;
  int m_end;
};
 
class SeqAndHits
{
public:
  SeqAndHits(const string & n, int size) {
    m_name = n;
    m_size = size;
  }

  void Add(const OneAlign & a) {
    m_aligns.push_back(a);
  }

  int Num() const {return m_aligns.isize();}
  const string & Name() const {return m_name;}
  const OneAlign & Align(int i) const {return m_aligns[i];} 
  int Size() const {return m_size;}

private:
  string m_name;
  int m_size;
  svec<OneAlign> m_aligns;
};

int Find(const svec<SeqAndHits> & h, const string & name) {
  for (int i=0; i<h.isize(); i++) {
    if (h[i].Name() == name)
      return i;
  }
  return -1;
}

int FindAdd(svec<SeqAndHits> & h, const string & name, int len) {
  //cout << "Size: " << len << endl;

  for (int i=0; i<h.isize(); i++) {
    if (h[i].Name() == name && h[i].Size() == len) 
      return i;
  }
  //cout << "Not found." << endl;
  h.push_back(SeqAndHits(name, len));
  return h.isize()-1;
}


string Shorten(const string &s) {
  //char tmp[4096];
  int i;
  int k = 0;
  int n = 0;
  for (i=0; i<s.size(); i++) {
    if (s[i] == '|') {
      n++;
      if (n == 4) {
	k = i;
 	break;
      }
    }
  }
  string out = &s[k+1];
  for (i=0; i<out.size(); i++) {
    if (out[i] == '_')
      out[i] = ' ';
  }
  return out;
}

class SpeciesLegend
{
public:
  SpeciesLegend() {
    m_bSpecify = true;
  }
  
  void Add(const string & name, const string & taxa, color & c) {
    int i;
    string species;
    Parse(species, name);

    if (m_bSpecify) {
      for (i=0; i<m_species.isize(); i++) {
	if (m_species[i] == species) {
	  c = m_color[i];
	  return;
	}
      }
      //cout << "Add species: " << species << endl;
      m_species.push_back(species);
      m_taxa.push_back(taxa);     
      SpecifyColor(c);
      m_color.push_back(c);
       //SpecifyColor(m_color[m_color.isize()-1]);
    } else {
      for (i=0; i<m_taxa.isize(); i++) {
	if (m_taxa[i] == taxa)
	  c = m_color[i];
	  return;
      }
      m_species.push_back(taxa);
      m_taxa.push_back(taxa);
      m_color.push_back(c);   
    }
  }

  int Num() const {return m_species.isize();}
  const string & Name(int i) const {return m_species[i];}
  const color & Color(int i) const {return m_color[i];}
  string Desc(int i) const {
    
    string d;
    if (m_bSpecify) {
      d = m_species[i];
      d += " (" + m_taxa[i] + ")";
    } else {
      d = m_taxa[i];
    }
    return d;
  }


private:
  void SpecifyColor(color & c) {
    double v = 0.4;
   
    double r = c.R() + RandomFloat(v)-v/2.;
    double g = c.G() + RandomFloat(v)-v/2.;
    double b = c.B() + RandomFloat(v)-v/2.;

    //cout << RandomFloat(v)-v/2. << endl;

    Limit(r);
    Limit(g);
    Limit(b);
    //cout << "SPECIFY: " << c.R() << " " << c.G() << " " << c.B() << " --> " << r << " " << g << " " << b << endl;
    c = color(r, g, b);
   }
    
  bool m_bSpecify;
  svec<string> m_species;
  svec<string> m_taxa;
  svec<color> m_color;
  
};


int PrintAlignment(TaxaColor & tax, const string & inFile, const string & outFile, bool bJpeg, const string & filter)
{

  string i1 = inFile;
  string o = outFile;

  int i, j;
  
  SpeciesLegend legend;

  ns_whiteboard::whiteboard board;

  int x_offset = 20;
  int y_offset = 20;

  double x_max = 0;
  double y_max = 0;

  FlatFileParser parser;

  parser.Open(i1);
  
  string refSeq;
  string qSeq;
  int start = -1;

  ProteinNameClusterer cluster;

  svec<SeqAndHits> hits;
  int num = 0;
  int len = 0;
  int size = 0;
  //parser.ParseLine();
  double maxSize = 0.;
  double maxSizeTemp = 0.;
  bool bGot = false;

  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    //NOTE: This is backwards, it should be the query!!!
    if (parser.AsString(0) == "Target" &&  parser.AsString(2) == "size:") {
      size = parser.AsInt(3);
      if (size > maxSizeTemp)
	maxSizeTemp = size;
      start = -1;
    }

    //BACKWARDS!!
    if (parser.AsString(0) == "Sbjct:") {
      qSeq += parser.AsString(2);
      if (start < 0)
	start = parser.AsInt(1);
    }
    //BACKWARDS!!
    if (parser.AsString(0) == "Query:") {
      refSeq += parser.AsString(2);
    }
    if (parser.AsString(0) == "Summary") {
      if (filter != "" && filter != parser.AsString(1)) {
	maxSizeTemp = 0;
	continue;
      }
      maxSize = maxSizeTemp;
      string name;
      string ref;
      for (i=1; i<parser.GetItemCount(); i++) {
	if (parser.AsString(i) == "vs.")
	  break;
	if (name != "")
	  name += " ";
	name += parser.AsString(i);
      }
      i++;
      for (; i<parser.GetItemCount(); i++) {
	if (parser.AsString(i) == "score:")
	  break;
	if (ref != "")
	  ref += " ";
	ref += parser.AsString(i);
      }
      
  
      int index = FindAdd(hits, name, size);
      //cout << "Query " << index << " size " << size << endl;
      //cout << "   Align " << strlen(refSeq.c_str())
      size = -1;
      hits[index].Add(OneAlign(start, qSeq, refSeq, ref));
      num++;

      refSeq = "";
      qSeq = "";
    }
  }

  cout << "Max size: " << maxSize << endl;

  //double xscale = 3.0;
  double xscale = 1200./maxSize;
  double yscale = 3.;
  double x, y;
  double w = 5.*yscale;
  double ww = 5.*xscale;

  double total = 0.;
  for (i=0; i<hits.isize(); i++) {
    total += 10.;
    const SeqAndHits & s = hits[i];
    total += 6*s.Num();
  }

  y = y_offset;
  double fs = 20.;

  for (i=0; i<hits.isize(); i++) {
    const SeqAndHits & s = hits[i];
    int size = s.Size();
    x = xscale*(double)size;
    // Draw query
    y += 10.*xscale;
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset, y), 
                                        ns_whiteboard::xy_coords(x_offset + x, y+w),
                                        color(0.1,0.1,0.1)) );
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x + 5, y),
                                        s.Name(), color(0,0,0), fs, "Times-Roman", 0, true));
    //x += xscale*14.*strlen(s.Name().c_str());
    x += 600;
    if (x > x_max)
      x_max = x;

    y += 10.*yscale;

   
    for (j=0; j<s.Num(); j++) {
      const OneAlign & a = s.Align(j);
      double x1 = xscale*(double)a.Start();
      double x2 = 0;
      tax.SetTo(a.Name());
      color base = tax.GetColor(1);
      legend.Add(a.Name(), tax.GetCurrent(), base);
      tax.SetBaseColor(base);
      cluster.Add(a.Name());
      

      for (int z=0; z<a.Len(); z++) {
	double v = a.Val(z);
	x2 = x1 + xscale*z;
	//double r = 1. - v;
	//double g = 1. - v;
	//double b = 1.;
	//cout << v << endl;
	color toUse = tax.GetColor(v);

	//toUse = color(0.5, 0.5, 0.5);

	board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + x2, y), 
					    ns_whiteboard::xy_coords(x_offset + x2 + xscale, y+w),
					    toUse) );
      }

      board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + x2 + 5, y),
					  Shorten(a.Name()), color(0,0,0), fs, "Times-Roman", 0, true));

      y += 7.*yscale;
      if (y > y_max)
	y_max = y;
    }
    
  }
  y += 7.*yscale;
  double barwidth = 30;
  double xx = 0.;
  double over = 340.;

  if (hits.isize() == 0) {
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + 30, y),
					"No hits found.", color(0,0,0), fs+4., "Arial", 0, true));
    y += 7.*yscale;
    x_max = 500;
    if (y > y_max)
      y_max = y;
  }


  for (i=0; i<legend.Num(); i++) {
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + xx, y), 
					ns_whiteboard::xy_coords(x_offset + xx + barwidth + yscale, y+w),
					legend.Color(i)) );
    
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + xx + barwidth + 8, y),
					legend.Desc(i), color(0,0,0), fs+2., "Times-Roman", 0, true));
    xx += over;
    if (xx + over + x_offset > x_max) {
      xx = 0;
      y += 7.*yscale;
    }
    if (y > y_max)
      y_max = y;
  }

  y += 14.*yscale;
  xx = 0.;
  cluster.Cluster();
  for (i=0; i<cluster.Num(); i++) {
    stringstream ss;
    ss << cluster.Class(i) << ": " << cluster.Count(i);
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + xx + 2 * barwidth + 8, y),
					ss.str(), color(0,0,0), fs+4., "Arial", 0, true));
 
    y += 8.*yscale;
    
    if (y > y_max)
      y_max = y;

  }

  board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + xx + 2 * barwidth + 8, y),
				      "Summary (similar proteins):", color(0,0,0), fs+4., "Arial", 0, true));
  y += 7.*yscale;
  
  if (y > y_max)
    y_max = y;


  ofstream out (o.c_str());
  ns_whiteboard::ps_display display(out, x_max + x_offset, y_max + x_offset);
  board.DisplayOn(&display);
 
  if (bJpeg) {
    char tmp[2048];
    strcpy(tmp, o.c_str());
    for (i=strlen(tmp); i>=0; i--) {
      if (tmp[i] == '.') {
	tmp[i] = 0;
	break;
      }
    }
    string jpeg = tmp;
    
    //jpeg += ".jpeg";
    jpeg += ".png";
    string cmmd = "convert ";
    cmmd += o;
    cmmd += " ";
    cmmd += jpeg;
    int r = system(cmmd.c_str());
  }

  //board.DeletePointers();
 
  return 0;
}
