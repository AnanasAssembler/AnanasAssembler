#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>

#include "base/CommandLineParser.h"
#include "src/Saguaro/HMMDecode.h"
#include "base/FileParser.h"


int main(int argc,char** argv)
{

  
  commandArg<string> aStringCmmd("-c","Configuration file that specifies the structure og the HMM");
  commandArg<string> bStringCmmd("-i","Stream of scores per frame for each state (by name)");
  commandLineParser P(argc,argv);
  P.SetDescription("Generic and configurable HMM decoder/classifier.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);

  P.parse();

  string config = P.GetStringValueFor(aStringCmmd);
  string in = P.GetStringValueFor(bStringCmmd);
 


  HMMUpdate update;


  Read(update, config);

  int i;



  cout << "Feeding frames..." << endl;
  HMMStateFrameScores a;

  FlatFileParser parserIn;
  
  parserIn.Open(in);


  svec<string> labels;
  svec<int> frameID;

  parserIn.ParseLine();

  for (i=0; i<parserIn.GetItemCount(); i++) {
    labels.push_back(parserIn.AsString(i));
  }
    
  int k = 0;
  while (parserIn.ParseLine()) {
    if (parserIn.GetItemCount() == 0)
      continue;

    a.Reset();

 
    for (i=1; i<parserIn.GetItemCount(); i++) {
      double val = parserIn.AsFloat(i);
      if (val < 0.)
	val = -val;
      //cout << k << "\tSetting " << labels[i] << "\t" << val << endl; 
      a.Set(labels[i], val);
    }

    if (k >= frameID.isize())
      frameID.resize(k + 2000000);
    frameID[k] = parserIn.AsInt(0);

    k++;

    update.Update(a);
    
    if (k % 10000 == 0)
      cout << "Processed: " << k << endl;
  }


  cout << "=================== Traceback ===========================" << endl;

  svec<HMMTracebackNode> out;

  update.TraceBack(out);


  for (i=0; i<out.isize(); i++) {
    int len = frameID[out[i].To()] - frameID[out[i].From()];
    int frames = out[i].To() - out[i].From();
    cout << out[i].Name() << "\t" << in << ": " << frameID[out[i].From()] << " - " << frameID[out[i].To()] << "\tlength: " << len;
    cout << "\t(frames " << out[i].From() << "-" << out[i].To() << " l=" << frames << ") \tscore=" << out[i].Score();
    cout << endl;
  }


  return 0;

}
  
