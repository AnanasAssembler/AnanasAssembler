
#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "ryggrad/src/base/ThreadHandler.h"
#include "ryggrad/src/util/SysTime.h"
#include "Version.h"

int Run(const string & exec, const string & cmmd, bool bIgnoreFailure = false)
{
    string c = exec + cmmd;
    cout << GetTimeStatic();
    cout << " Executing: " << c << endl;
    int ret = system(c.c_str());
    if (ret != 0 && !bIgnoreFailure) {
      cout << GetTimeStatic() << " ERROR running " << c << endl;
      cout << "Terminated abnormally." << endl;
      exit(-1);
    }

    cout << GetTimeStatic() << " Completed " << endl;
    return ret;
}

bool Exists(const string & file) 
{
    FILE * p = fopen(file.c_str(), "r");
    if ( p == NULL)
        return false;
    fclose(p);
    return true;
}

class LayoutThread : public IOneThread
{
public:
    LayoutThread(const string & exec,
                 const string & command,
                 int i,
                 int total)
    {
        m_index = i;
        m_total = total;
        char tmp[1024];
        sprintf(tmp, " -num %d -of %d ", i, total);
        m_cmmd = command;
        m_cmmd += tmp;
        m_exec = exec;
    }


protected:

    virtual bool OnDie() {
        //cout << "Killed!!" << endl;
        return true;
    }

    virtual bool OnDo(const string & msg) {
        //int startRead = atol(msg);
        cout << "Starting layout guided w/ " << msg << endl;
        cout << "Spawning: " << m_cmmd << endl;
        //int ret = system(m_cmmd.c_str());
        int ret = Run(m_exec, m_cmmd);
        cout << "Process " << m_index << " is done, return code: " << ret << endl;

        //cout << "Done!" << endl;
        return true;
    }

    virtual bool OnInitialize(const string & msg) {
        //cout << "Initializing!" << endl;
        return true;
    }
private:
    string m_exec;
    string m_cmmd;
    int m_index;
    int m_total;
};



string Number(int i)
{
    char tmp[256];
    sprintf(tmp, "%d", i);
    string f = tmp;
    return f;
}
string NumberFloat(double d)
{
    char tmp[256];
    sprintf(tmp, "%f", d);
    string f = tmp;
    return f;
}


void PrintLogo()
{
    cout << endl;
    cout << "\\\\\\||///" << endl;
    cout << " \\\\||//" << endl;
    cout << "  \\||/" << endl;
    cout << " .<><><>." << endl;
    cout << ".<><><><>." << endl;
    cout << "'<><><><>'" << endl;
    cout << " '<><><>'" << endl;
    cout << endl;
}




int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input fasta file");
    commandArg<string> outCmmd("-o","output directory", "ananas_out");
    commandArg<double> minCmmd("-m","minimum overlap identity", 0.98);
    commandArg<double> minGroupCmmd("-mg","minimum identity for grouping", 0.99);
    commandArg<int> cpuCmmd("-n","number of CPU cores", 1);
    commandArg<int> cpuCmmd2("-n2","number of CPU cores for isoform enumeration", 1);
    commandArg<string> dirCmmd("-dir","direction of pairs: fr fowards each other, ff same direction, na unpaired");
    commandArg<int> cpuLapCmmd("-no","number of processes for overlap finding", 1);
    commandArg<int> bandCmmd("-b","bandwidth of alignments (maximum indel size)", 0);
    commandArg<int> mlCmmd("-ml","minimum overlap (for alignments)", 25);
    commandArg<int> stepCmmd("-s","step size (for alignments)", 30);
    commandArg<string> ssCmmd("-strand","strand specificity (0=no 1=yes)", "0");
    commandArg<int> contigSizeCmmd("-minContigLen","minimum length of a single-contig scaffold to report", 200);
    commandArg<int> libSizeCmmd("-libSize","Maximum library size", 500);
    //commandArg<bool> filtCmmd("-group","group identical reads (recommended for large data sets)", false);
    commandArg<int>    overlapIterCmmd("-overlapIter", "Number of iterations with increasing liniency for overlap computation", 2);
    commandArg<string> readGroupFileCmmd("-readGroupFile","read groupin information file if available","");
    commandArg<string> prefixCmmd("-prefix","The prefix to add to all generated contig names", "Sample1");
    commandArg<bool> redunCmmd("-rr","Remove redundant transcripts", false);
    commandLineParser P(argc,argv);
    P.SetDescription("Assembles reads from overlaps.");
    P.registerArg(fileCmmd);
    P.registerArg(outCmmd);
    P.registerArg(minCmmd);
    P.registerArg(minGroupCmmd);
    P.registerArg(dirCmmd);
    P.registerArg(bandCmmd);
    P.registerArg(ssCmmd);
    P.registerArg(libSizeCmmd);
    P.registerArg(contigSizeCmmd);
    P.registerArg(mlCmmd);
    P.registerArg(stepCmmd);
    P.registerArg(cpuCmmd);
    P.registerArg(cpuCmmd2);
    P.registerArg(cpuLapCmmd);
    P.registerArg(overlapIterCmmd);
    //P.registerArg(filtCmmd);
    P.registerArg(readGroupFileCmmd);
    P.registerArg(prefixCmmd);
    P.registerArg(redunCmmd);
  
    P.parse();
  
    string readsFileName = P.GetStringValueFor(fileCmmd);
    string outName = P.GetStringValueFor(outCmmd);
    string dir = P.GetStringValueFor(dirCmmd);
    string ss = P.GetStringValueFor(ssCmmd);
    int minContig = P.GetIntValueFor(contigSizeCmmd);
    int libSize = P.GetIntValueFor(libSizeCmmd);
    double mI = P.GetDoubleValueFor(minCmmd);
    double minGroupI = P.GetDoubleValueFor(minGroupCmmd);
    int cpu = P.GetIntValueFor(cpuCmmd);
    int cpu2 = P.GetIntValueFor(cpuCmmd2);
    int cpuLap = P.GetIntValueFor(cpuLapCmmd);
    int step = P.GetIntValueFor(stepCmmd);
    int bandwidth = P.GetIntValueFor(bandCmmd);
    int minoverlap = P.GetIntValueFor(mlCmmd);
    int    overlapIter     = P.GetIntValueFor(overlapIterCmmd); 
    //bool bGroup = P.GetBoolValueFor(filtCmmd);
    string readGroupFile   = P.GetStringValueFor(readGroupFileCmmd);
    string prefix = P.GetStringValueFor(prefixCmmd);
    bool bRemoveRedundant = P.GetBoolValueFor(redunCmmd);

    PrintLogo();
    cout << "Welcome to package " << ananas_software << " " << ananas_version << endl;
    cout << endl;

 


    if (dir != "fr" && dir != "ff" && dir != "na") {
        cout << "ERROR: Don't understand -dir " << dir << endl;
        return 0;
    }
    bool bUnpaired = false;
    if (dir == "na") {
        //dir = "fr";
        bUnpaired = true;
    }

    int i;
  
    const char * pExec = argv[0];
    cout << "Executing " << pExec << endl;

    char exec_dir[8192];
    strcpy(exec_dir, pExec);
    bool bSlash = false;
    for (i = strlen(exec_dir)-1; i>=0; i--) {
        if (exec_dir[i] == '/') {
            exec_dir[i+1] = 0;
	    bSlash = true;
            break;
        }
    }
    if (!bSlash)
        strcpy(exec_dir, "");


    Run("mkdir ", outName, true);
    outName += "/";
    string cmmd;

    /*if (bGroup) {
      cmmd = "SimplicityFilter "
      Run(exec_dir, cmmd);
      cmmd = "BuildReadGroups -i "
      Run(exec_dir, cmmd);
      }*/

    ////////////////////////////////////////////////////////////////////
    ///// findOverlaps ////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////


    string pairSzFile = outName + "/pairSz.tmp";
    string lapFile    = outName + "/allOverlaps.out";
    string groupFile  = outName + "/consensusReads.out";

    if (bUnpaired)
      //cmmd = "findOverlaps -I 0.98 -b 30 -B 2 -O 75 -s 1 -i " + readsFileName;
      //cmmd = "findOverlaps -I 0.98 -b 30 -B 0 -O 75 -s 0 -i " + readsFileName;
      cmmd = "findOverlaps -S 25 -I " + NumberFloat(mI) + " -b " + Number(step);
    else
      cmmd = "findOverlaps -I " + NumberFloat(mI) + " -b " + Number(step);

    cmmd += " -d " + NumberFloat(minGroupI) + " -B " + Number(bandwidth) +  " -O " + Number(minoverlap) + " -s " 
             + ss + " -i " + readsFileName + " -t " + pairSzFile + " -T " + Number(cpu) + " -g " 
             + readGroupFile + " -C " + groupFile + " -o " + lapFile + " -overlapIter " + Number(overlapIter);

  if   (Exists(lapFile)) {
      cout << "Overlaps exist, skipping." << endl;
      cout << "Remove " << lapFile << " to re-compute read overlaps." << endl;
    } else {
      Run(exec_dir, cmmd);
    }



    ////////////////////////////////////////////////////////////////////
    ////  Layout //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////



    string layoutFile = outName + "contigs.layout.all_done";
  
    if (Exists(layoutFile)) {
        cout << "Layout files exis, skipping." << endl;
        cout << "Remove " << layoutFile << " to re-compute initial stage layout of contigs" << endl;

    } else {
        char cpuNum[256];
        sprintf(cpuNum, "%d", cpu);
    
    
        cmmd = "Layout -i " + pairSzFile;
        cmmd += " -n ";
        cmmd += cpuNum;
        cmmd += " -l " + outName + "/allOverlaps.out";
        cmmd += " -o " + outName + "/contigs.layout";
        cmmd += " -g " + groupFile;
        cmmd += " -dir " + dir;
        cmmd += " -libSize " + Number(libSize);
        cmmd += " -prefix " + prefix;
    
        Run(exec_dir, cmmd);
    }

    if (cpu > 1) {
        layoutFile = outName + "contigs.layout";
        if (!Exists(layoutFile)) {
            // Note: this is extremely STUPID and dangerous!!!
            cmmd = "cat " + outName + "/contigs.layout.* > " + outName + "/contigs.layout";
            cout << "Concatenating layout files: " << cmmd << endl;
            int ret = system(cmmd.c_str());
        }
    }


    ///////////////////////////////////////////////////////////
    ///// Scaffolder /////////////////////////////////////////
    //////////////////////////////////////////////////////////

    cmmd = "Scaffolder -f " + pairSzFile;
    cmmd += " -l " + outName + "/allOverlaps.out";
    cmmd += " -o " + outName + "/scaffolds.layout";
    cmmd += " -i " + outName + "/contigs.layout";
    cmmd += " -g " + groupFile;
    cmmd += " -dir " + dir;

    Run(exec_dir, cmmd);



    ///////////////////////////////////////////////////////////
    //// LayoutGuided ////////////////////////////////////////
    /////////////////////////////////////////////////////////

    cmmd = "LayoutGuided -i " + pairSzFile;
    cmmd += " -l " + outName + "/allOverlaps.out";
    cmmd += " -o " + outName + "/contigs_altsplic_raw.layout";
    //  cmmd += " -f " + outName + "/contigs_altsplic.fasta";
    cmmd += " -s " + outName + "/scaffolds.layout";
    cmmd += " -g " + groupFile;
    cmmd += " -dir " + dir;
    cmmd += " -libSize " + Number(libSize);
    cmmd += " -prefix " + prefix;

    if (cpu2 == 1) {
        Run(exec_dir, cmmd);
    } else {
        ThreadHandler th;
        for (i=0; i<cpu2; i++) {
            th.AddThread(new LayoutThread(exec_dir, cmmd, i, cpu2), "");
        }
        for (i=0; i<cpu2; i++) {
            th.Feed(i, "run");
        }
        while (!th.AllDone()) {
            usleep(10000);
        }
        cout << "LayoutGuided is done!!!" << endl;
        // Note: this is extremely STUPID and dangerous!!!
        cmmd = "cat " + outName + "/contigs_altsplic_raw.layout.* > " + outName + "/contigs_altsplic_raw.layout";
        cout << "Concatenating layout files: " << cmmd << endl;
        int ret = system(cmmd.c_str());    
    }

    //////////////////////////////////////////////////////////
    // RunIsoEM         /////////////////////////////////////
    ////////////////////////////////////////////////////////

  
    cmmd = "RunIsoEM " ;
    cmmd += " -i " + outName + "/contigs_altsplic_raw.layout";
    cmmd += " -o " + outName + "/contigs_altsplic.layout";
    Run(exec_dir, cmmd);


    //////////////////////////////////////////////////////////
    // GenAssemblyFasta /////////////////////////////////////
    ////////////////////////////////////////////////////////

    string partitionsOutName = outName + "partitions";
#if defined(FORCE_DEBUG)
    Run("mkdir ", partitionsOutName, true);
#endif
    partitionsOutName += "/";
  
    cmmd = "GenAssemblyFasta " ;
    cmmd += " -r " + readsFileName;
    cmmd += " -i " + outName + "/contigs_altsplic.layout";
    cmmd += " -c " + outName + "/consensusReads.out";
    cmmd += " -o " + outName + "/final.fa";
    cmmd += " -minContig " + Number(minContig);
    cmmd += " -prefix " + prefix;
    cmmd += " -readsOutDir " + partitionsOutName;
    Run(exec_dir, cmmd);


    if (bRemoveRedundant) {
      cmmd = "RemoveRedundantSeqs " ;
      cmmd += " -i " + outName + "/final.fa";
      cmmd += " -o " + outName + "/final.fa";
      Run(exec_dir, cmmd);
    }


    cout << GetTimeStatic() << " All DONE! " << endl;

    return 0;
}