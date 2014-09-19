#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"


class OneAlign
{
public:
  OneAlign() {
    m_pVal = 1.;
    m_meaningful = 1;
    m_score = 0.;
  }

  void Add(const string & line) {    
    m_data.push_back(line);
  }
  void Set(const string & name, const string & ref, double pval, double score) {
    m_name = name;
    m_pVal = pval;
    m_score = score;
    m_meaningful = 0;
    if (Meaningless(ref, "unnamed_protein_product"))
      m_meaningful = 1;
    if (Meaningless(ref, "uncharacterized_protein"))
      m_meaningful = 1;
    if (Meaningless(ref, "_unknown_"))
      m_meaningful = 1;
    if (Meaningless(ref, "hypothetical_protein"))
      m_meaningful = 1;
    if (Meaningless(ref, "predicted_protein"))
      m_meaningful = 1;
    m_ref = ref;
  }


  void Print() const {
    for (int i=0; i<m_data.isize(); i++) 
      cout << m_data[i] << endl;
    cout << endl;
    //cout << m_name << "\t" << m_pVal << "\t" << m_meaningful << endl;
  }
  

  bool operator < (const OneAlign & a) const {
    if (m_name != a.m_name)
      return (m_name < a.m_name);
    if (m_meaningful != a.m_meaningful)
      return (m_meaningful < a.m_meaningful);
    //return (m_pVal < a.m_pVal);
    return (-m_score < -a.m_score);
  }

  const string & Name() const {return m_name;}
  const string & Ref() const {return m_ref;}
  double PVal() const {return m_pVal;}
  double Score() const {return m_score;}

  void Clear() {
    m_data.clear();
    m_pVal = 1.;
    m_meaningful = 1;
    m_name = "";
  }
  
  void SetPlainName(const string & n) {
    m_plain = n;
  }
  const string & PlainName() const {
    return m_plain;
  }


private:
  bool Meaningless(const string & name, const string & phrase) {
    if (strstr(name.c_str(), phrase.c_str()) != NULL)
      return true;
    return false;
  }


  svec<string> m_data;
  double m_pVal;
  string m_name;
  string m_ref;
  int m_meaningful;
  string m_plain;
  double m_score;
};



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> fCmmd("-filter","filter for only those sequences", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Sorts the protein aligner parser.");
  P.registerArg(fileCmmd);
  P.registerArg(fCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string filter = ">" + P.GetStringValueFor(fCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  svec<OneAlign> all;
  OneAlign tmp;
  int i;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() > 4) {
      if (parser.AsString(0) == "Sequence" && 
	  parser.AsString(parser.GetItemCount()-1) == "more.")
	continue;
      if (parser.AsString(0) == "Summary") {
	string name;
	for (i=1; i<parser.GetItemCount(); i++) {
	  if (parser.AsString(i) == "vs.")
	    break;
	  name += parser.AsString(i);
	}
	double score = (parser.AsFloat(i+9)-parser.AsFloat(i+8))*parser.AsFloat(i+10);
	//cout << parser.Line() << " SCORE: " << score << endl;
	tmp.Add(parser.Line());	
	tmp.Set(name, parser.AsString(i+1), parser.AsFloat(i+3), score);
	tmp.SetPlainName(parser.AsString(1));
	all.push_back(tmp);
	tmp.Clear();
	continue;
      }
    }
    tmp.Add(parser.Line());
  }

  Sort(all);
  int max_per = 20;
  int k = 0;
  string last;

  if (all.isize() == 0) {
    cout << "<br>No hits found." << endl;
    cout << "No hits found." << endl;
    cout << "No hits found." << endl;
    cout << "There are no results to display.<br><br>" << endl;
  }

  svec<string> best, bestRef;
  svec<double> bestP;

  for (i=0; i<all.isize(); i++) {
    if (all[i].Name() == "")
      continue;

    string localName;
    StringParser local;
    local.SetLine(all[i].Name(), "frame:");
    localName = local.AsString(0);

    if (localName == last) {
      //if (all[i].PVal() < bestP[bestP.isize()-1]) {
      //bestP[bestP.isize()-1] = all[i].PVal();
      if (all[i].Score() > bestP[bestP.isize()-1]) {
	bestP[bestP.isize()-1] = all[i].Score();
	bestRef[bestP.isize()-1] = all[i].Ref();
      }
    } else {
      best.push_back(localName);
      bestRef.push_back(all[i].Ref());
      bestP.push_back(all[i].Score());
      //cout << localName << endl;
      last = localName;
    }
  }

  for (i=0; i<best.isize(); i++) {
    cout << "Best hit for " << best[i] << " " << bestRef[i] << " " << bestP[i] << endl;
    if (filter == ">")
      filter = best[i];
  }
  //cout << "Filter: " << filter << endl;
  last = "";
  for (i=0; i<all.isize(); i++) {
    if (all[i].Name() == "")
      continue;
    // Filter here!!
    if (filter != "" && all[i].PlainName() != filter)
      continue;

    if (all[i].Name() == last) {
      k++;
    } else {
     
      k = 0;
    }
    last = all[i].Name();
    if (k <max_per) {
      all[i].Print();
    }
  }


  return 0;
}
