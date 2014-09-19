#include "visual/VisualAlignments.h"
#include "base/CommandLineParser.h"


int main( int argc, char** argv )
{
 
  commandArg<string> aStringI1("-i","Smörgås output file");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<bool> bJPEG("-j","make jpeg file (using convert)", false);
  //commandArg<string> taxString("t", "taxonomy file", "/references/taxonomy/taxonomy.txt");
  commandArg<string> taxString("-t", "taxonomy file", "data/taxonomy.txt");

  //commandArg<double> dotSize("-d","dot size", 1.);
  //commandArg<double> aScale("-s","scale", 60000.);
  //commandArg<int> cTarget("-t","target id", -1);
  //commandArg<bool> bF("-f","forward only", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Displays alignments.");

  P.registerArg(aStringI1);
  P.registerArg(aStringO);
  P.registerArg(bJPEG);
  P.registerArg(taxString);
  //P.registerArg(dotSize);
  //P.registerArg(aScale);
  //P.registerArg(cTarget);
  //P.registerArg(bF);

  P.parse();

  string i1 = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(aStringO);
  string taxFile = P.GetStringValueFor(taxString);
  //double dd = P.GetDoubleValueFor(dotSize);
  //int targetID = P.GetIntValueFor(cTarget);
  bool bJpeg = P.GetBoolValueFor(bJPEG);

  TaxaColor m_tax;
  m_tax.Add("Metazoa;");
    m_tax.Add("Alveolata;");
    m_tax.Add("Amoebozoa;");
    m_tax.Add("Bacteroidetes;");
    m_tax.Add("Chordata;");
    m_tax.Add("Cyanobacteria;");
    m_tax.Add("Diplomonadida;");
    m_tax.Add("Euglenozoa;");
    m_tax.Add("Euryarchaeota;");
    m_tax.Add("Viridiplantae;");
    m_tax.Add("Firmicutes;");
    m_tax.Add("Fungi;");
    m_tax.Add("Fusobacteria;");
    m_tax.Add("Proteobacteria;");
    m_tax.Add("Tenericutes;");
    m_tax.Add("Actinobacteria;");
    m_tax.Add("Chlorobi;");
    m_tax.Add("Chloroflexi;");
 
  m_tax.ReadTaxonomy(taxFile);

  PrintAlignment(m_tax, i1, o, bJpeg);
 
  return 0;
}
