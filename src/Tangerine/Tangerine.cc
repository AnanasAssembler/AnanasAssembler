
#include <string>

#include "base/CommandLineParser.h"

int main(int argc,char** argv)
{
  int i;
  const char * pExec = argv[0];
  cout << "Executing " << pExec << endl;

  char exec_dir[1024];

  strcpy(exec_dir, pExec);
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      break;
    }
  }

  
  commandArg<string> qCmmd("-q", "protein query");
  commandArg<string> kCmmd("-k", "default k-mer size", "4");
  commandArg<int> cpuCmmd("-cpu", "number of CPUs to use", 16);
  commandArg<string> dbCmmd("-db", "database", "/references/databases/uniref90/uniref90.fasta");
  commandArg<string> outputCmmd("-o", "output folder", "");

  commandLineParser P(argc,argv);
  P.SetDescription("Runs SatsumaProtein in parallel.");
  P.registerArg(qCmmd);
  P.registerArg(kCmmd);
  P.registerArg(cpuCmmd);
  P.registerArg(dbCmmd);
  P.registerArg(outputCmmd);

  P.parse();

  string query = P.GetStringValueFor(qCmmd);
  string k = P.GetStringValueFor(kCmmd);
  string db = P.GetStringValueFor(dbCmmd);
  string outputFolder = P.GetStringValueFor(outputCmmd);

  int n= P.GetIntValueFor(cpuCmmd);

  for (i=0; i<n; i++) {
    string cmmd = exec_dir;
    cmmd += "/runSatsumaProt -t ";
    cmmd += db;
    char tmp[128];
    sprintf(tmp, ".%d", i);

    cmmd += tmp;
    cmmd += " -q ";
    cmmd += query;
    cmmd += " > ";
    cmmd += outputFolder;
    cmmd += "/runSatsumaProt";
    cmmd += tmp;
    cmmd += " &";
    cout << cmmd << endl;
    system (cmmd.c_str());
  }


}
