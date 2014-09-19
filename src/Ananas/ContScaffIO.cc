#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Ananas/ContScaffIO.h"
#include "base/FileParser.h"


void ContigScaffoldIO::Read(Assembled & assembled, const string &file)
{

  FlatFileParser parser;
  
  parser.Open(file);

  Scaffold scaffold;
  Contig contig;

  bool bHaveScaff = false;

  int offset = 0;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "<SCAFFOLD>") {
      bHaveScaff = true;
      scaffold.SetName(parser.AsString(1));
      continue;
    }

    if (parser.AsString(0) == "</SCAFFOLD>") {
      assembled.push_back(scaffold);
      scaffold.clear();
      continue;
    }

    if (parser.AsString(0) == "</CONTIG>") {
      scaffold.push_back(contig, offset, 1);
      // One scaffold per contig if we don't have explicit info
      if (!bHaveScaff) {
	assembled.push_back(scaffold);
	scaffold.clear();
      }
      contig.clear();
      continue;
    }

    if (parser.AsString(0) == "<CONTIG>") {
      contig.clear();
      contig.SetName(parser.AsString(1));
      if (parser.GetItemCount() >= 3)
	offset = parser.AsInt(2);
      if (parser.GetItemCount() >= 4)
	contig.SetOri(parser.AsInt(3));
      continue;
    }

    int read = parser.AsInt(0);
    int ori = parser.AsInt(1);
    int start = parser.AsInt(2);
    int stop = parser.AsInt(4);
    const string & name = parser.AsString(5);
    contig.push_back(ReadPlacement(read , ori, start, stop, name));
  }  
}


void ContigScaffoldIO::Write(const Assembled & assembled, const string &file)
{
  FILE * pOut = fopen(file.c_str(), "w");
  if (pOut == NULL) {
    cout << "ERROR: Could not open file " << file << " for writing." << endl;
  }

  int count = 0;

  int i, j, k, l;
  for (l=0; l<assembled.isize(); l++) {
    const Scaffold & s = assembled[l];
    if (s.isize() == 0)
      continue;
    string name = s.Name();
    if (name == "" ) {
      //name = "<unknown>";
      char tmp[1024];
      sprintf(tmp, ">Scaffold_%5d", count);
      for (unsigned int x = 0; x<strlen(tmp); x++) {
	if (tmp[x] == ' ')
	  tmp[x] = '0';
      }
      name = tmp;
    }
    count++;
    fprintf(pOut, "<SCAFFOLD>\t%s\n", name.c_str());
    for (i=0; i<s.isize(); i++) {
      const Contig & c = s[i];
      string nameC = c.Name();
      if (nameC == "")
	nameC = "<unknown>";
      fprintf(pOut, "<CONTIG>\t%s\t%d\t%d\n", nameC.c_str(), s.Offset(i), c.Ori());
      
      for (j=0; j<c.isize(); j++) {
	const ReadPlacement & r = c[j];	
	fprintf(pOut, "%d\t%d\t%d - %d\t%s\n", r.Read(), r.Ori(), r.Start(), r.Stop(), r.Name().c_str());
	
      }
      fprintf(pOut, "</CONTIG>\t%s\n", nameC.c_str());
    }
    fprintf(pOut, "</SCAFFOLD>\t%s\n", name.c_str());
  }

  fclose(pOut);
}
  
