#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"

class Kmer
{
public:
    Kmer() {
        m_pos = -1;
        m_mul = -1;
    }
    Kmer(const DNAVector & d, int from, int size, int pos_base) {
        m_mul = -1;
        Set(d, from, size, pos_base);
    }

    void Set(const DNAVector & d, int from, int size, int pos_base) {
        d.Substring(m_kmer, from, size);
        m_pos = pos_base + from;
    }
 
    void IncMult() {
        m_mul++;
    }
 
    void SetMult(int m) {
        m_mul = m;
    }
    void SetPos(int i) {
        m_pos = i;
    }

    const string & GetKmer() const {return m_kmer;}
    int Pos() const {return m_pos;}
    int Mult() const {return m_mul;}

    bool operator < (const Kmer & k) const {
        //return m_kmer < k.m_kmer;
        if (m_kmer != k.m_kmer)
            return m_kmer < k.m_kmer;   
        return m_pos < k.m_pos;
    }
    bool operator != (const Kmer & k) const {
        return !operator == (k);
    }
    bool operator == (const Kmer & k) const {
        int diff = m_pos - k.m_pos;
        if (diff < 0)
            diff = -diff;
        // HARD CODED!!!!
        if (m_kmer == k.m_kmer && diff < 10)
            return true;
        return false;
    }

    void Extend(char l) {
        int i;
        for (i=1; i<m_kmer.size(); i++) {
            m_kmer[i-1] = m_kmer[i];
        }
        m_kmer[m_kmer.size()-1] = l;
    }

private:
    string m_kmer;
    int m_pos;
    int m_mul;
};


class KmerAssembly
{
public:
    KmerAssembly() {
        m_size = 24;
    }

    void Clear() {
        m_kmers.clear();
    }

    void Add(const DNAVector & d, int base_pos) {
        int i;
        for (i=0; i<=d.isize()-m_size; i++) {
            m_kmers.push_back(Kmer(d, i, m_size, base_pos));
        }
    }

    void Build(DNAVector & out) {
        Sort(m_kmers);

        int start = -1;
   
        int startMul = 0;
        int i, j;

        svec<Kmer> comp;
    
        Kmer tmp;
        for (i=0; i<m_kmers.isize(); i++) {
            tmp = m_kmers[i];
            tmp.SetMult(1);
            for (j=i+1; j<m_kmers.isize(); j++) {
                if (m_kmers[j] != m_kmers[i]) {
                    comp.push_back(tmp);
                    //cout << "Adding: " << tmp.GetKmer() << " " << tmp.Mult() << endl;
                    if (tmp.Pos() == 0 && tmp.Mult() > startMul) {
                        startMul = tmp.Mult();
                        start = comp.isize()-1;
                    } 
                    break;
                } else {
                    tmp.IncMult();
                    //if (m_kmers[j].Pos() < tmp.Pos()) {
                    //cout << "WERID!!" << endl;
                    //tmp.SetPos(m_kmers[j].Pos());
                    //}
                }
            }
            i = j-1;
        }

        // Last element (stupid!!!)
        comp.push_back(tmp);
        //cout << "Adding: " << tmp.GetKmer() << " " << tmp.Mult() << endl;
        if (tmp.Pos() == 0 && tmp.Mult() > startMul) {
            startMul = tmp.Mult();
            start = comp.isize()-1;
        } 


        // Now we have unique kmers w/ counts.
        comp.push_back(tmp);
        Sort(comp);
        svec<int> used;
        used.resize(comp.isize(), 0);
        if (start == -1) {
            cout << "ERROR: Start " << start << endl;
            ///for (i=0; i<comp.isize()
        }
        m_kmers = comp; // Stupid...


        string result = comp[start].GetKmer();
        int best = -1;
        int lastPos = 0;
        do {
            const Kmer & k = comp[start];
            //cout << "Kmer: " << k.GetKmer() << endl;
      
            svec<Kmer> ext;
            ext.resize(4);
            for (i=0; i<ext.isize(); i++) {
                ext[i] = k;
                ext[i].SetPos(0); // Get the first one
                ext[i].Extend(NucLetter(i));
            }

            best = -1;
            int bestBase = -1;
            int bestCount = 0;
            for (i=0; i<ext.isize(); i++) {
                //cout << "Lookup " << ext[i].GetKmer() << endl;
                int index = BinSearchFuzzy(comp, ext[i]);
                index = FindBest(comp, ext[i], index, lastPos);
                if (index < 0)
                    continue;
                //cout << "Index=" << index << " pos " <<  comp[index].Pos() << endl;
	
                if (used[index] == 0 && comp[index].Mult() > bestCount) {
                    bestCount = comp[index].Mult();
                    best = index;
                    bestBase = i;
                }
            }
            if (bestBase != -1) {
                //cout << "Best ext: " << NucLetter(bestBase) << endl;
                if (best == -1)
                    cout << "ERROR!" << endl;
                used[best] = 1;
                result += NucLetter(bestBase);
                lastPos = comp[best].Pos();
                //cout << result << endl;
            }
            start = best;
        } while (best != -1);
        out.SetFromBases(result);
    }
    

private:
    int FindBest(const svec<Kmer> & comp, Kmer & k, int index, int pos) {
        if (index < 0)
            return index;
        int i;
        int min = m_size;
        int bestIndex = -1;
        for (i=index; i<comp.isize(); i++) {
            if (comp[i].GetKmer() != k.GetKmer())
                break;
            int p = comp[i].Pos();
            int diff = p - pos - 1;
            //cout << "Compare " << pos << " " << p << " " << diff << " i=" << i << endl;
            if (diff < 0)
                diff = -diff * 2;
            if (diff < min) {
                min = diff;
                bestIndex = i;
            }
        }
        return bestIndex;
    }


    svec<Kmer> m_kmers;
    int m_size;
};




int main( int argc, char** argv )
{

    commandArg<string> fileCmmd("-i","input file (layout)");
    commandLineParser P(argc,argv);
    P.SetDescription("Rebuilds consensus sequence for ananas contigs.");
    P.registerArg(fileCmmd);
  
    P.parse();
  
    string fileName = P.GetStringValueFor(fileCmmd);
    int i;

    //comment. ???
    FlatFileParser parser;
  
    parser.Open(fileName);
    string name;

    KmerAssembly kmers;

    vecDNAVector out;

    while (parser.ParseLine()) {
        if (parser.GetItemCount() == 0)
            continue;
        if (parser.AsString(0) == "<CONTIG>") {
            name = parser.AsString(1);
            continue;
        }
        if (parser.AsString(0) == "</CONTIG>") {
            // Build consensus
            DNAVector dd; 
            kmers.Build(dd);
            cout << name << endl;
            for (i=0; i<dd.isize(); i++) 
                cout << dd [i];
            cout << endl;
            kmers.Clear();
            continue;
        }
        DNAVector d;
        d.SetFromBases(parser.AsString(7));
        kmers.Add(d, parser.AsInt(2));
    }
    return 0;
}
