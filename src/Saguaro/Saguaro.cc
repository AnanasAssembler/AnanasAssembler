#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/Saguaro/HMMDecode.h"
#include "src/Saguaro/HMMDistance.h"



class CactusLine
{
public:
  CactusLine() {
  }

  double & operator [] (int i) {return m_data[i];}
  void resize(int i) {m_data.resize(i);}
  int isize() const {return m_data.isize();}

private:
  svec<double> m_data;
};



class Cactus
{
public:
  Cactus(int size) {
    m_data.resize(size);
    for (int i=0; i<size; i++)
      m_data[i].resize(size);
  }

  void Read(FlatFileParser & parser, const string & name, int size) {
    int i, j;
    for (i=0; i<size; i++)  {
      parser.ParseLine();
      CactusLine & l = m_data[i];
      for (j=1; j<parser.GetItemCount(); j++) {
	l[j-1] = parser.AsFloat(j);
	//if (l[j-1] > 1.0)
	//l[j-1] = 1.0;
	if (l[j-1] < 0.0)
	  l[j-1] = 0.0;
	  
      }
    }
  }

  void WriteLine(FILE * p, int i) {
     CactusLine & l = m_data[i];
     for (int j=0; j<l.isize(); j++) {
       fprintf(p, " %f", l[j]);
     }
     fprintf(p, "\n");
  }

  void Add(double v)  {
    int i, j;
    for (j=0; j<m_data.isize(); j++) {
      CactusLine & l = m_data[j];
      for (i=0; i<l.isize(); i++) {
	l[i] += v;
	if (l[i] < 0.01)
	  l[i] = 0.01;
      }
    }

  }
  void Scale(double v)  {
    int i, j;
    for (j=0; j<m_data.isize(); j++) {
      CactusLine & l = m_data[j];
      for (i=0; i<l.isize(); i++) {
	l[i] *= v;
	if (l[i] < 0.001)
	  l[i] = 0.001;
      }
    }

  }
  
private:
  string m_name;
  svec<CactusLine> m_data;

};



class HMMControl
{

public:
  HMMControl(const string & outDir) {
    m_outDir = outDir;
    Exec("mkdir " + outDir, true);
    m_numStates = 3;
    m_config = m_outDir + "/saguaro.config";
    m_trees = m_outDir + "/saguaro.cactus";
    m_log = m_outDir + "/HMMTrain.out";
    m_garbage = m_outDir + "/saguaro.garbage";
    
    string tmp = "rm " + m_garbage + ".vec";
    Exec(tmp, true);

    m_run = 0;
    m_exec = "./";
    m_trans = 150.;
    m_threads = 1;
  }

  void SetTransition(double d) {
    m_trans = d;
  }

  void SetThreads(int nt) {
    m_threads = nt;
  }

  void SetExec(const string &s) {
    m_exec = s;
  }

  void SetupNames(const string & in, const string & lst) {
    HMMFeatureVector features;

    if (in != "") {
      cout << "Reading features..." << endl;
      features.Read(in);
      m_feat = in;
      cout << "done!" << endl;
    } else {
      FlatFileParser parserIn;
      parserIn.Open(lst);
      m_list = lst;
      while (parserIn.ParseLine()) {
	if (parserIn.GetItemCount() == 0)
	  continue;
	cout << "Reading features from " << parserIn.AsString(0) << endl;
	features.MergeRead(parserIn.AsString(0));
      }
    } 
    int i;
    m_names.resize(features.GetNameCount());
    for (i=0; i<features.GetNameCount(); i++) {
      m_names[i] = features.GetName(i);
      cout << "Adding individual " << m_names[i] << endl;
    }
  }


  void UpdateConfigs() {
    FlatFileParser parser;
    parser.Open(m_trees);

    svec<Cactus> cacti;

    while (parser.ParseLine()) {
      if (parser.GetItemCount() == 0)
	continue;
      string name = parser.AsString(0);
      parser.ParseLine();
      Cactus tmp(m_names.isize());
      tmp.Read(parser, name, m_names.isize());
      cacti.push_back(tmp);      
    }
    cout << "# of active cacti: " << cacti.isize() << endl;
    int half = cacti.isize();

    /*
    if (cacti.isize() == 1) {
      string cp = "cp " + m_trees + " " + m_garbage;
      Exec(cp);
      }*/


    // if (cacti.isize() == 1) {
    //Cactus tmp = cacti[0];
    //tmp.Scale(0.95);
    //cacti.push_back(tmp);
    //cacti[0].Scale(1.05);
    // } else {
    
      FlatFileParser parserGarb;
      parserGarb.Open(m_garbage);
      
      cout << "Reading garbage..." << endl;
      while (parserGarb.ParseLine()) {
	if (parserGarb.GetItemCount() == 0)
	  continue;
	string name = parserGarb.AsString(0);
	parserGarb.ParseLine();
	Cactus tmp(m_names.isize());
	tmp.Read(parserGarb, name, m_names.isize());
	//tmp.Scale(0.6);
	cacti.push_back(tmp);      
      }
      //}

    int i, j;
    /*
    double up = 1.2;
    double down = 0.5;
    if (cacti.isize() > 4) {
      down = 0.8;
      up = 1.;
    }
    

    for (i=0; i<half; i++) {
      Cactus tmp = cacti[i];
      cacti[i].Scale(down);
      tmp.Scale(up);
      cacti.push_back(tmp);
      }*/

    FILE * pOut = fopen(m_trees.c_str(), "w");

    for (i=0; i<cacti.isize(); i++) {
      fprintf(pOut, "cactus%d\n", i);
      for (j=0; j<m_names.isize(); j++) {
	fprintf(pOut, "%s ", m_names[j].c_str());
      }
      fprintf(pOut, "\n");
      
      for (j=0; j<m_names.isize(); j++) {
	fprintf(pOut, "%s ", m_names[j].c_str());
	cacti[i].WriteLine(pOut, j);
      }      
    }

    fclose(pOut);


    FILE * pConfig = fopen(m_config.c_str(), "w");
    for (j=0; j<cacti.isize(); j++) {
      fprintf(pConfig, "cactus%d", j);
      for (i=0; i<m_numStates; i++)
	fprintf(pConfig, "\tcactus%d", j);
      
      fprintf(pConfig, "\n");
    }

    fprintf(pConfig, "\nall_transitions %f\n", m_trans);

    fclose(pConfig);
  }



  void WriteInitialFiles() {

    int i, j;

    FILE * pConfig = fopen(m_config.c_str(), "w");
    fprintf(pConfig, "cactus0");
    for (i=0; i<m_numStates; i++)
      fprintf(pConfig, "\tcactus0");
      
    fprintf(pConfig, "\n");
    fclose(pConfig);

    
    FILE * pCactus = fopen(m_trees.c_str(), "w");
    fprintf(pCactus, "cactus0\n");
    for (i=0; i<m_names.isize(); i++) {
      fprintf(pCactus, "%s ", m_names[i].c_str());
    }
    fprintf(pCactus, "\n");

    for (i=0; i<m_names.isize(); i++) {
      fprintf(pCactus, "%s ", m_names[i].c_str());
      for (j=0; j<m_names.isize(); j++) {
	fprintf(pCactus, "0.2 ");
      }
      fprintf(pCactus, "\n");
    }
 
    fclose(pCactus);

    string cp = "cp " + m_trees + " " + m_garbage;
    Exec(cp, true);


  }


  void Run(int iter) {
    string cmd = m_exec + "HMMTrain ";
    cmd += m_neurons;
    if (m_list != "") {
      cmd += "-l " + m_list;
    } else {
      cmd += "-i " + m_feat;
    }

    cmd += " -c " + m_config;
    cmd += " -m " + m_trees;
    cmd += " -o " + m_trees;
    cmd += " -g " + m_garbage;


    char tmpp[256];
    sprintf(tmpp, " -n %d ", m_threads);
    cmd += tmpp;


    char tmp[256];
    sprintf(tmp, " -iter %d ", iter);
    cmd += tmp;
    cmd += " > " + m_log;

    sprintf(tmp, ".%d ", m_run);
    cmd += tmp;
    m_run++;
    Exec(cmd);
    cout << "Done running." << endl;
  }
  
  void SetRunNum(int i) {
    m_run = i;
  }

  void SetNumNeurons(int k) {
    char n[1024];
    sprintf(n, " -neurons %d ", k);
    m_neurons = n;
    cout << m_neurons << endl;
  }

  void RunLocalTrees() {
    string cmd = m_exec + "LocalTrees ";
    if (m_list != "") {
      cmd += "-l " + m_list;
    } else {
      cmd += "-i " + m_feat;
    }
    cmd += " -c " + m_config;
    cmd += " -m " + m_trees;
    cmd += " > " + m_outDir + "/LocalTrees.out";
    Exec(cmd);
  }

private:
  void Exec(const string & cmd, bool bIgnore = false) {
    cout << "Running " << cmd << endl;
    int ret = system(cmd.c_str());
    if (ret != 0 && bIgnore == false) {
      cout << "ERROR: HMMTrain failed with error code " << ret << endl;
      exit(-1);
    }
      
  }

  int m_run;
  int m_threads;
  int m_numStates;
  svec<string> m_names;
  string m_outDir;
  string m_config;
  string m_trees;
  string m_feat;
  string m_list;
  string m_log;
  string m_garbage;

  string m_exec;
  string m_neurons;


  double m_trans;
};


int main( int argc, char** argv )
{

  //cout << argv[0] << endl;

  char tmpPath[512];
  strcpy(tmpPath, argv[0]);
  tmpPath[strlen(tmpPath)-7] = 0;
  cout << "Path to executable: " << tmpPath << endl;
  

  //commandArg<string> fileCmmd("-i","feature input file");
  commandArg<string> bStringCmmd("-f","Feature vector", "");
  commandArg<string> lStringCmmd("-l","Feature vector list file", "");
  commandArg<string> oStringCmmd("-o","output directory");
  commandArg<int> iterCmmd("-cycle","iterations per cycle", 2);
  commandArg<int> allCmmd("-iter","iterations with split", 40);
  commandArg<double> transCmmd("-t","transition penalty", 150.);
  commandArg<int> resumeCmmd("-resume","resume w/ iteration #", 0);
  commandArg<int> numCmmd("-n","number of threads", 1);
  commandArg<int> neuronCmmd("-neurons","number of neurons in the SOM", 800);
  commandLineParser P(argc,argv);
  P.SetDescription("Smoothly, automatically and generically uncover the ancestry of related organisms.");
  //P.registerArg(fileCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(lStringCmmd);
  P.registerArg(oStringCmmd);
  P.registerArg(iterCmmd);
  P.registerArg(allCmmd);
  P.registerArg(transCmmd);
  P.registerArg(resumeCmmd);
  P.registerArg(numCmmd);
  P.registerArg(neuronCmmd);
 
  P.parse();
  
  //string fileName = P.GetStringValueFor(fileCmmd);
  string in = P.GetStringValueFor(bStringCmmd);
  string lst = P.GetStringValueFor(lStringCmmd);
  string outdir = P.GetStringValueFor(oStringCmmd);
  int numIter = P.GetIntValueFor(iterCmmd);
  int allIter = P.GetIntValueFor(allCmmd);
  double transPen = P.GetDoubleValueFor(transCmmd);
  int resume = P.GetIntValueFor(resumeCmmd);
  int numThreads = P.GetIntValueFor(numCmmd);
  int neurons = P.GetIntValueFor(neuronCmmd);
  
  HMMControl control(outdir);
  control.SetExec(tmpPath);
  control.SetTransition(transPen);
  control.SetThreads(numThreads);
  control.SetNumNeurons(neurons);
  control.SetupNames(in, lst);

  if (resume == 0) {
    control.WriteInitialFiles();
    control.Run(1);
  } else {
    control.SetRunNum(resume);
    cout << "WARNING: Skipping initial file write due to the -resume option" << endl;
  }


  //control.Run(numIter);

  int i;
  for (i=resume; i<allIter; i++) {
    cout << "Iteration # " << i << endl;
    control.UpdateConfigs();
    control.Run(numIter);
  }

  cout << "Computing local phylogienies." << endl;
  control.RunLocalTrees();
  cout << "All done!" << endl;

  return 0;
}
