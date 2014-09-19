//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include <string>

#include "base/CommandLineParser.h"
#include "src/Saguaro/HMMDecode.h"
#include "src/Saguaro/HMMDistance.h"
#include "base/FileParser.h"

#include "src/Saguaro/SOM.h"
#include "src/Saguaro/HMMDistance.h"
#include "base/RandomStuff.h"


class HMMMatrixScores : public HMMStateFrameScoreSource
{
 public:
  HMMMatrixScores(HMMTrees * pTrees) {
    m_pTrees = pTrees;
  } 
  virtual ~HMMMatrixScores() {}

  void SetFeature(const HMMFeature & f) {
    m_feature = f;
    m_pTrees->ClearCache();
  }

  virtual double GetScore(const string & stateName) const {
    //out << "Getting score..." << endl;
    return m_pTrees->GetScore(stateName, m_feature);
  }

private:
  HMMTrees * m_pTrees;
  HMMFeature m_feature;
};


bool F2F(svec<double> & feat, svec<int> & valid, const svec<double> & f)
{  
  int k = f.isize();
  if (f[f.isize()-1] < 0.)
    return false;
     
  feat.resize(k, 0.);
  
  valid.resize(k, 1);
  for (int x=0; x<f.isize(); x++) {
    
    if (f[x] < 0.) {
      valid[x] = 0;
    } else {
      valid[x] = 1;
      feat[x] = f[x];
    }
  }
  return true;
}



 
int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-c","Configuration file that specifies the structure of the HMM");
  commandArg<string> bStringCmmd("-i","Feature vector", "");
  commandArg<string> lStringCmmd("-l","Feature vector list file", "");
  commandArg<string> cStringCmmd("-m","Distance matrices");
  commandArg<string> singleCmmd("-single","Name of the cactus (for one only)", "");
  commandArg<int> minCmmd("-min","Minimum length", 3);
  commandLineParser P(argc,argv);
  P.SetDescription("Generic and configurable HMM training.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(lStringCmmd);
  P.registerArg(cStringCmmd);
  P.registerArg(singleCmmd);
  P.registerArg(minCmmd);
  P.parse();

  string config = P.GetStringValueFor(aStringCmmd);
  string in = P.GetStringValueFor(bStringCmmd);
  string lst = P.GetStringValueFor(lStringCmmd);
  string matrix = P.GetStringValueFor(cStringCmmd);
  string single = P.GetStringValueFor(singleCmmd);
  int minlen = P.GetIntValueFor(minCmmd);
 
  if (lst == "" && in == "") {
    cout << "You MUST specify either -i or -l!" << endl;
    return -1;
  }

  int i, j, k;


  HMMFeatureVector features;

  if (in != "") {
    cout << "Reading features..." << endl;
    features.Read(in);
    cout << "done!" << endl;
  } else {
    FlatFileParser parserIn;
    parserIn.Open(lst);

    while (parserIn.ParseLine()) {
      if (parserIn.GetItemCount() == 0)
	continue;
      cout << "Reading features from " << parserIn.AsString(0) << endl;
      features.MergeRead(parserIn.AsString(0));
      
    } 
  }



  HMMTrees model;
  HMMMatrixScores scores(&model);

  cout << "Reading models." << endl;
 
  model.AddRead(matrix);


  int wordCount;

  cout << "Dynprog'ing..." << endl;
  HMMUpdate update;
  Read(update, config);
  wordCount = update.GetWordCount();
  cout << "Number of frames: " << features.isize() << endl;
  for (i=0; i<features.isize(); i++) {
    scores.SetFeature(features[i]);
    update.Update(scores);
    if (i > 0 && i % 10000 == 0) {
      cout << "Processed: " << i << " (" << 100. * (double)i/(double)features.isize() << " %)"  << endl;
    }
  }

  cout << "REPORTING Traceback and Update" << endl;  
  vec<HMMTracebackNode> out;
  update.TraceBack(out);

  int lastScore = 0;
  for (i=0; i<out.isize(); i++) {  
    int len = features[out[i].To()].GetPosition() - features[out[i].From()+1].GetPosition();
    if (len < minlen)
      continue;

    int frames = out[i].To() - out[i].From() - 1;
    cout << out[i].Name() << "\t" << features[out[i].To()].Name() << ": " << features[out[i].From()+1].GetPosition();
    cout << " - " << features[out[i].To()].GetPosition() << "\tlength: " << len;
    cout << "\t(frames " << out[i].From()+1 << "-" << out[i].To() << " l=" << frames << ") \tscore=" << (out[i].Score() - lastScore)/double(out[i].To()-out[i].From());
    cout << endl;
      

    if (single == "" || single == out[i].Name()) { 
      model.ResetSameDiff(); 
      for (j=out[i].From()+1; j<=out[i].To(); j++) {
	model.Update(out[i].Name(), features[j]);
      }
      //cout << "Model " << endl;
      model.PrettyPrintOne(out[i].Name());      
    }
  }

  return 0;

}
  
