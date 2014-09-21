#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"

void ConnectUp(svec<int> & partner, const vecDNAVector & seq)
{
    int cc = 0;
    for (int i=0; i<seq.isize(); i++) {
        char name[2048];
        strcpy(name, seq.Name(i).c_str());
        int n = strlen(name);
        bool b = false;
        if (name[n-2] != '/')
            continue;
        if (name[n-1] == '1') {
            name[n-1] = '2';
            b = true;
        } else {
            if (name[n-1] == '2') {
                name[n-1] = '1';
                b = true;
            }
        }
        if (!b)
            continue;
        int index = seq.NameIndex(name);
        if (index < 0)
            continue;
        //cout << "Connecting partners: " << name << " " << seq.Name(i) << endl;
        partner[i] = index;
        cc++;
    }
    cout << "Found partner for " << cc << " reads, out of " << seq.isize() << endl;
}


int Mis(const DNAVector & a, const DNAVector & b)
{
    if (a.isize() != b.isize())
        return a.isize();
    int m = 0;
    for (int i=0; i<a.isize(); i++) {
        if (a[i] != b[i])
            m++;
    }
    return m;
}


int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input fasta file");
    commandArg<string> outCmmd("-o","output fasta file");
    commandArg<int> nCmmd("-n","minimum redundancy", 5);
    commandArg<int> mCmmd("-m","maximum mismatch", 1);
    commandArg<bool> invCmmd("-inverse","output redundant reads", false);
    commandLineParser P(argc,argv);
    P.SetDescription("Groups reads by identity.");
    P.registerArg(fileCmmd);
    P.registerArg(outCmmd);
    P.registerArg(nCmmd);
    P.registerArg(mCmmd);
    P.registerArg(invCmmd);
  
    P.parse();
  
    string fileName = P.GetStringValueFor(fileCmmd);
    string outName = P.GetStringValueFor(outCmmd);
    int max = P.GetIntValueFor(nCmmd);
    int mis = P.GetIntValueFor(mCmmd);
    bool bInv = P.GetBoolValueFor(invCmmd);

    vecDNAVector dnaR;
    cout << "Reading input fasta." << endl;
    dnaR.Read(fileName);
    svec<int> check;
    check.resize(dnaR.isize(), -1);
    ConnectUp(check, dnaR);

    cout << "Sorting reads." << endl;
    dnaR.Sort();

    vecDNAVector dna;
    dna.resize(dnaR.isize());
    for (int x=0; x<dna.isize(); x++) {
        dna[x] = dnaR[x];
        dna.SetName(x, dnaR.Name(x));
    }

    check.clear();
    check.resize(dna.isize(), -1);
    ConnectUp(check, dna);
    //dna.Write(outName);
    //return 0;

    int i, j;
  
    svec<int> readClass;
    readClass.resize(dna.isize(), -1);
    svec<int> partner;
    partner.resize(dna.isize(), -1);
    svec<int> counts;

    int k = 0;
    readClass[0] = k;
    counts.push_back(1);
    for (i=1; i<dna.isize(); i++) {
        if (Mis(dna[i-1], dna[i]) > mis) {
            k++;
            counts.push_back(0);
        }
        readClass[i] = k;
        counts[k]++;
    }
  
    cout << "Reads: " << dna.isize() << " classes: " << k << endl;

    cout << "Finding partners" << endl;
    ConnectUp(partner, dna);
  

    svec<int> good;
    good.resize(dna.isize(), 0);

    //  int keep = 0;

    for (i=0; i<dna.isize(); i++) {
        int c = readClass[i];
        if (counts[c] < max) {
            good[i] = 1;
            //keep++;
            continue;
        }
        int maxPartner = 0;
        int index = -1;

        int from = i;

        for (; i<dna.isize(); i++) {
            if (readClass[i] != c) {
                i--;
                break;
            }
            int p = partner[i];
            if (p < 0)
                continue;
            int cp = readClass[p];
            if (counts[cp] > maxPartner) {
                maxPartner = counts[cp];
                index = i;
            } 
        }
        if (index == -1) {
            good[from] = 1;
        } else {
            //cout << "Found good pair." << endl;
            good[index] = 1;
            good[partner[index]] = 1;
        }
    
    }

    if (bInv) {
        for (i=0; i<good.isize(); i++) {
            if (good[i] == 1)
                good[i] = 0;
            else
                good[i] = 1;
        }
    }


    int keep = 0;
    for (i=0; i<good.isize(); i++) {
        int p = partner[i];
        //if (good[i])
        if (good[i] && (p < 0 || good[p])) {
            keep++;
        }
    }
    cout << "Will keep " << keep << " reads." << endl;

    vecDNAVector out;
    out.resize(keep);
    k = 0;
    for (i=0; i<good.isize(); i++) {
        //if (good[i]) {
        int p = partner[i];
        if (good[i] && (p < 0 || good[p])) {
            out[k] = dna[i];
            out.SetName(k, dna.Name(i));
            k++;
        }
    }

    cout << "Almost done. " << endl;
    check.clear();
    check.resize(out.isize(), -1);
    ConnectUp(check, out);

    out.Write(outName);
  
    return 0;
}
