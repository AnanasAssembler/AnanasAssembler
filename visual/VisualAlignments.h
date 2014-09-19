#ifndef VISUALALIGNEMNTS_H
#define VISUALALIGNEMNTS_H

#include <string>
#include "visual/Color.h"
#include "base/SVector.h"
#include "base/FileParser.h"

class Species
{
public:
  Species() {}

  void Read(FlatFileParser & parser);


  const string & Name() {return m_name;}
  const string & Tax()  {return m_tax;}

private:
  string m_name;
  string m_tax;
};


class TaxaColor
{
public:
  TaxaColor();

  void ReadTaxonomy(const string & file);

  void Add(const string & taxa);
  color GetColor(double shade);
  void SetTo(const string & name);
  const string & GetCurrent() const {return m_curr;}
  void SetBaseColor(color b) {m_baseCol = b;}
private:
  color Get(const string & s);
  bool GetTaxonomy(string & tax, const string & name);
 
  svec<string> m_taxa;
  svec<color> m_col;
  svec<int> m_use;
  svec<Species> m_spec;
  color m_baseCol;
  string m_curr;

};


int PrintAlignment(TaxaColor & c, const string & inFile, const string & outFile, bool bJpeg = true, const string & filter = "");


#endif

