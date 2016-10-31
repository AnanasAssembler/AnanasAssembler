#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ryggrad/src/base/CommandLineParser.h"
#include "ContScaff.h"
#include "ContScaffIO.h"
#include "ConsensOverlapUnit.h"


class MergeItem
{
public:
    MergeItem() {
        m_count = 0;
        m_scaff1 = -1;
        m_scaff2 = -1;
        m_dir = 0;
    }
    MergeItem(int s1, int s2, int count = 0, int dir = 0) {
        m_count = count;   
        if (s1 < s2) {
            m_scaff1 = s1;
            m_scaff2 = s2;
        } else {
            m_scaff1 = s2;
            m_scaff2 = s1;
        }
        m_dir = dir;
    }

    bool operator < (const MergeItem & s) const {
        if (m_count != s.m_count) {
            return m_count < s.m_count;
        }
        if (m_scaff1 != s.m_scaff1) {
            return m_scaff1 < s.m_scaff1;
        }
        return m_scaff2 < s.m_scaff2;
    }

    int Dir() const {return m_dir;}
    int Count() const {return m_count;}
    int Scaff1() const {return m_scaff1;}
    int Scaff2() const {return m_scaff2;}

private:
    int m_count;
    int m_scaff1;
    int m_scaff2;
    int m_dir;
};

class MergeQueue
{
public:
    MergeQueue() {}

    void Add(int s1, int s2, int dir) {
        m_queue.push_back(MergeItem(s1, s2, 0, dir));
    }

    void Compact();
 
    int isize() const {return m_queue.isize();}
    const MergeItem & operator [] (int i) const {return m_queue[i];}

private:
    svec<MergeItem> m_queue;
};



class Scaffolder
{
public:
    Scaffolder() {
        m_pairs = 1;
    }
    void SetPairing(const string &s) {
        if (s == "fr")
            m_pairs = -1;
        if (s == "ff")
            m_pairs = 1;
        if (s == "na")
            m_pairs = 0;    
    }

    void Join(Assembled & assembly, const ConsensOverlapUnit & COUnit);
    void RemoveIdentical(Assembled & assembled, const ConsensOverlapUnit & COUnit);


private:
    void BuildScaffLocTable(const Assembled & assembly, const ConsensOverlapUnit & COUnit);
    void Merge(Assembled & assembly, const ConsensOverlapUnit & COUnit, int s1, int s2, int dir);
    bool CheckOverlaps(const Scaffold & s1, const Scaffold & s2, const ConsensOverlapUnit & COUnit);

    svec<int> m_scaffLoc;
    svec<int> m_scaffLocDir;
    int m_pairs;
};

void MergeQueue::Compact()
{
    int i;
    Sort(m_queue);
    svec<MergeItem> tmp;

    int l1 = -1; 
    int l2 = -1;
    int np = 0;
    int nm = 0;

    for (i=0; i<m_queue.isize(); i++) {
        if (m_queue[i].Scaff1() != l1 || m_queue[i].Scaff2() != l2) {
            if (np > 1 || nm > 1) { // Stupid heuristics....
                int n = np;
                int dir = 1;
                if (nm > np) {
                    dir = -1;
                    n = nm;
                }
                tmp.push_back(MergeItem(l1, l2, n, dir));
            }
            np = 0;
            nm = 0;
            l1 = m_queue[i].Scaff1();
            l2 = m_queue[i].Scaff2();
        }
        if (m_queue[i].Dir() == -1)
            nm++;
        else
            np++;
    }
  
    if (np > 1 || nm > 1) { // Stupid heuristics....
        int n = np;
        int dir = 1;
        if (nm > np) {
            dir = -1;
            n = nm;
        }
        tmp.push_back(MergeItem(l1, l2, n, dir));
    }

    Sort(tmp);

    m_queue.clear();
    m_queue.resize(tmp.isize());
    int k = 0;
    for (i=tmp.isize()-1; i>=0; i--) {
        m_queue[i] = tmp[k];
        k++;
    }
}

void Scaffolder::Join(Assembled & assembly, const ConsensOverlapUnit & COUnit)
{
    BuildScaffLocTable(assembly, COUnit);
    MergeQueue q;
    int i, j;

    for (i=0; i<COUnit.GetNumReads(); i++) {
        int numPartner = COUnit.getNumOfPartners(i);
        //if (partner < 0)
        //continue;
        for (j=0; j<numPartner; j++) {
            int partner = COUnit.getPartner(i, j);
            int s1 = m_scaffLoc[i];
            int s2 = m_scaffLoc[partner];
            int d1 = m_scaffLocDir[i];
            int d2 = m_scaffLocDir[partner];

            int dir = d1 * d2 * m_pairs;
      
            if (s1 != -1 && s2 != -1) {
                if (s1 != s2) {
                    q.Add(s1, s2, dir);
                }
            }
        }
    }
    q.Compact();
  
    for (i=0; i<q.isize(); i++) {
        const MergeItem & m = q[i];

        int count = m.Count();
        int r1 = assembly[m.Scaff1()].NumReads();
        int r2 = assembly[m.Scaff2()].NumReads();
        int h1 = assembly[m.Scaff1()].Length();
        int h2 = assembly[m.Scaff2()].Length();

        double total = r1;
        double len = h1;
        if (r2 < r1) {
            total = r2;
            len = h2;
        }
    
        if (total == 0) {
            continue;
        }

        int expect = (int)(150. * total / len + 0.5); 
//        cout << m.Scaff1() << " vs. " << m.Scaff2() << " c=" << count << " e=" << expect << " r1=" << r1 << " r2=" << r2 << endl;
        if (count >= expect) {
//            cout << "Merging." << endl;
            Merge(assembly, COUnit, m.Scaff1(), m.Scaff2(), m.Dir());
        } else {
//            cout << "Discard." << endl;
        }

        //    if (CheckOverlaps(assembly[m.Scaff1()], assembly[m.Scaff2()], COUnit))
        //Merge(assembly, COUnit, m.Scaff1(), m.Scaff2(), m.Dir());
    }

}

bool Scaffolder::CheckOverlaps(const Scaffold & s1, const Scaffold & s2, const ConsensOverlapUnit & COUnit)
{
    int i, j, k, l;
    int total = 0;
    int good = 0;
    for (i=0; i<s2.isize(); i++) {
        const Contig & c2 = s2[i];
        for (j=0; j<c2.isize(); j++) {
            const ReadPlacement & r2 = c2[j]; 
            total++;
            int id2 = r2.Read();
            int hasLap = 0;
            for (k=0; k<s1.isize(); k++) {
                const Contig & c1 = s1[k];
                for (l=0; l<c1.isize(); l++) {
                    const ReadPlacement & r1 = c1[l];
                    if (COUnit.HasLap(r1.Read(), id2))
                        hasLap++;
                } 
            }
            if (hasLap > 0) {
                good++;
            }
        }
    }
    if (total == 0)
        return false;
    double frac = (double)good/(double)total;
//    cout << "Fraction: " << frac << ", " << good << " out of " << total << endl;
    if (frac >= 0.1)
        //if (good > 0)
        return true;
    else
        return false;
}

void Scaffolder::Merge(Assembled & assembly, const ConsensOverlapUnit & COUnit, int s1, int s2, int dir)
{
    if (assembly[s1].isize() == 0 || assembly[s2].isize() == 0) {
        return;
    }

//    cout << "Merging scaffolds " << s1 << " <- " << s2 << endl;
    if (dir == -1)
        assembly[s2].Reverse();
    assembly[s1].push_back(assembly[s2], 0, 1);
    assembly[s2].clear();
}

void Scaffolder::BuildScaffLocTable(const Assembled & assembled, const ConsensOverlapUnit & COUnit)
{

    m_scaffLoc.clear();
    m_scaffLoc.resize(COUnit.GetNumReads(), -1);
    m_scaffLocDir.clear();
    m_scaffLocDir.resize(COUnit.GetNumReads(), 0);

    int i, j, k, l;
    for (l=0; l<assembled.isize(); l++) {
        const Scaffold & s = assembled[l];
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                const ReadPlacement & r = c[j]; 
                m_scaffLoc[r.Read()] = l;
                m_scaffLocDir[r.Read()] = r.Ori();	
            }
        }
    }
}

void Scaffolder::RemoveIdentical(Assembled & assembled, const ConsensOverlapUnit & COUnit)
{
    //cout << "Sorting assembly... " << endl;
    //assembled.Sort();
    //cout << "Done, start merging." << endl;

    m_scaffLoc.clear();
    m_scaffLoc.resize(COUnit.GetNumReads(), -1);

    int i, j, k, l;
    for (l=0; l<assembled.isize(); l++) {
        const Scaffold & s = assembled[l];
        // Check for identical scaffolds
        svec<int> other;
        int total = 0;
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                const ReadPlacement & r = c[j]; 
                total++;
                if (m_scaffLoc[r.Read()] == -1) { 
                    m_scaffLoc[r.Read()] = l;
                } else {
                    other.push_back(m_scaffLoc[r.Read()]);
                }
            }
            // Let's see if we have duplicates.
        }
        int n = other.isize();
        //cout << "n=" << n << endl;
        UniqueSort(other);
        //cout << "size=" << other.isize() << endl;
        if (other.isize() == 1 && total - n < 2) {
            //cout << "Removing scaffold " << l << " as duplicate of scaffold " << other[0] << endl; 
            assembled[l].clear();
        }
    }
}


int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input scaffold/contig/layout file");
    commandArg<string> lapCmmd("-l","input read overlap file");
    commandArg<string> fastaCmmd("-f","read fasta file"); // Do we really need this???
    commandArg<string> consCmmd("-g","input read consensus group file");
    commandArg<string> outCmmd("-o","output file");
    commandArg<string> pairCmmd("-dir","pairing (fr, ff, or na)");

    //commandArg<string> layoutCmmd("-o","output layout file", "contigs.layout");
    //commandArg<double> minCmmd("-m","minimum overlap identity", 0.99);
    //commandArg<bool> exCmmd("-e","Exhaustive search (report top-n)", false);
  
    commandLineParser P(argc,argv);
    P.SetDescription("Joins contigs into scaffolds.");
    P.registerArg(fileCmmd);
    P.registerArg(lapCmmd);
    P.registerArg(consCmmd);
    P.registerArg(fastaCmmd);
    P.registerArg(outCmmd);
    P.registerArg(pairCmmd);
    //P.registerArg(layoutCmmd);
    //P.registerArg(minCmmd);
    //P.registerArg(exCmmd);
  
    P.parse();
  
    string fileName = P.GetStringValueFor(fileCmmd);
    string lapName = P.GetStringValueFor(lapCmmd);
    string fastaName = P.GetStringValueFor(fastaCmmd);
    string outName = P.GetStringValueFor(outCmmd);
    string pairing = P.GetStringValueFor(pairCmmd);
    string consName = P.GetStringValueFor(consCmmd); 
 
    Assembled assembly;

    ContigScaffoldIO io;

    io.Read(assembly, fileName);

    ConsensOverlapUnit laps(fastaName, consName, "");

    Scaffolder scaffolder;
    scaffolder.SetPairing(pairing);
    scaffolder.RemoveIdentical(assembly, laps);
  
    int i, j;
  
    for (i=0; i<5; i++) {
        cout << "\rScaffolding, iteration " << i << "  " << flush;
        scaffolder.Join(assembly, laps);
        string tmpName = outName + ".iter.";
        char tmp[64];
        sprintf(tmp, "%d", i);
        tmpName += tmp;
        io.Write(assembly, tmpName);
    }

    io.Write(assembly, outName);

    return 0;
}

