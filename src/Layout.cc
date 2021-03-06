#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "SearchOverlaps.h"
#include "GlobUsage.h"
#include "ryggrad/src/base/ThreadHandler.h"

class SearchThread : public IOneThread
{
public:
    SearchThread(const string & layoutName, 
                 int index,
                 string prefix,
                 GlobalUsageHandler * pGlob, 
                 const ConsensOverlapUnit * pReads, 
                 int startRead,
                 const string & dir,
                 int libSize,
                 bool bEx, bool pRest)
    {
        m_startRead = startRead;
        m_pReads = pReads;
        m_search.SetExhaustive(bEx);
        m_search.SetPairRestrict(pRest);
        m_search.SetOutput(layoutName);
        m_search.SetGlobalUsage(pGlob);
        m_search.SetIndex(index);
        m_search.SetPrefix(prefix);
        m_search.SetDir(dir);
        m_search.SetLibSize(libSize);
        m_index = index;
    }


protected:

    virtual bool OnDie() {
        //cout << "Killed!!" << endl;
        return true;
    }

    virtual bool OnDo(const string & msg) {
        m_search.DoSearchAll(*m_pReads, m_pReads->GetNumReads(), m_startRead);
        return true;
    }

    virtual bool OnInitialize(const string & msg) {
        return true;
    }
private:
    Search m_search;
    const ConsensOverlapUnit * m_pReads;
    int m_startRead;
    int m_index;
};



int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input read pair/size info file");
    commandArg<string> lapCmmd("-l","input read overlap file");
    commandArg<string> consCmmd("-g","input read consensus group file");
    commandArg<string> layoutCmmd("-o","output layout file", "contigs.layout");
    commandArg<int> libSizeCmmd("-libSize","Maximum library size", 500);
    commandArg<string> dirCmmd("-dir","direction of pairs (fr or ff)");
    commandArg<int> cpuCmmd("-n","number of CPUs/cores to use", 1);
    commandArg<bool> exCmmd("-e","Exhaustive search (report top-n)", false);
    commandArg<bool> pRestCmmd("-pairRestrict","Restrict contigs with pair-end support", true);
    commandArg<bool> fCmmd("-force","Force using the thread handler", false);
    commandArg<string> prefixCmmd("-prefix","The prefix to add to all generated contig names", "Sample1");
    commandLineParser P(argc,argv);
    P.SetDescription("Assembles COUnit from overlaps.");
    P.registerArg(fileCmmd);
    P.registerArg(lapCmmd);
    P.registerArg(consCmmd);
    P.registerArg(layoutCmmd);
    P.registerArg(libSizeCmmd);
    P.registerArg(dirCmmd);
    P.registerArg(cpuCmmd);
    P.registerArg(exCmmd);
    P.registerArg(pRestCmmd);
    P.registerArg(fCmmd);
    P.registerArg(prefixCmmd);
  
    P.parse();
  
    string pairSzFileName = P.GetStringValueFor(fileCmmd);
    string lapName = P.GetStringValueFor(lapCmmd);
    string consName = P.GetStringValueFor(consCmmd);
    string layoutName = P.GetStringValueFor(layoutCmmd);
    int libSize = P.GetIntValueFor(libSizeCmmd);
    bool bEx = P.GetBoolValueFor(exCmmd);
    bool pRest = P.GetBoolValueFor(pRestCmmd);
    bool bForce = P.GetBoolValueFor(fCmmd);
    int cpu = P.GetIntValueFor(cpuCmmd);
    int n = cpu;
    string dir = P.GetStringValueFor(dirCmmd);
    string prefix = P.GetStringValueFor(prefixCmmd);

    ConsensOverlapUnit COUnit(pairSzFileName, consName, lapName);
 
    if (n == 0) {
        cout << "ERROR: -n was set to 0. We need at least one core. Re-setting to 1." << endl;
        n = 1;
    }
    

    // Special-case single CPU here (for extreme compatibility)
    if (n == 1 && !bForce) {
        Search search;
        search.SetExhaustive(bEx);
        search.SetPairRestrict(pRest);
        search.SetDir(dir);
        search.SetOutput(layoutName);
        search.SetLibSize(libSize);
        search.DoSearchAll(COUnit, COUnit.GetNumReads());
    } else {
        int i;
    
        GlobalUsageHandler glob;
        glob.SetSize(COUnit.GetNumReads());
    
        int step = COUnit.GetNumReads() / n / 2;
        ThreadHandler th;
        int startRead = 0;
    
        for (i=0; i<n; i++) {
            char tmp[256];
            sprintf(tmp, ".%d", i);
            string layout = layoutName;
            layout += tmp;
      
            th.AddThread(new SearchThread(layout, 
                                          i,
                                          prefix,
                                          &glob, 
                                          &COUnit, 
                                          startRead,
                                          dir,
                                          libSize,
                                          bEx, pRest), "init");
      
            startRead += step;
        }
    
        for (i=0; i<n; i++) {
            th.Feed(i, "do_it");
        }
    
        while (!th.AllDone()) {
            usleep(10000);
        }
    }

    cout << "Layout is done!! Exiting normal." << endl;

    string alldone = layoutName + ".all_done";
    cout << "Signing off: " << alldone << endl;
    FILE * pFinal = fopen(alldone.c_str(), "w");
    fprintf(pFinal, "Layout is done. Exited normally.\n");
    fclose(pFinal);

    return 0;
}
