#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "ContScaffIO.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/base/StringUtil.h"


void ContigScaffoldIO::Read(Assembled & assembled, const string &file)
{

    FlatFileParser parser;
  
    parser.Open(file);

    Scaffold scaffold;
    Contig contig;

    bool bHaveScaff = false;

    int offset = 0;

    while (parser.ParseLine()) {
        if (parser.GetItemCount() == 0)
            continue;
        if (parser.AsString(0) == "<SCAFFOLD>") {
            bHaveScaff = true;
            scaffold.SetName(parser.AsString(1));
            continue;
        }

        if (parser.AsString(0) == "</SCAFFOLD>") {
            assembled.push_back(scaffold);
            scaffold.clear();
            continue;
        }

        if (parser.AsString(0) == "</CONTIG>") {
            scaffold.push_back(contig, offset, 1);
            // One scaffold per contig if we don't have explicit info
            if (!bHaveScaff) {
                assembled.push_back(scaffold);
                scaffold.clear();
            }
            contig.clear();
            continue;
        }

        if (parser.AsString(0) == "<CONTIG>") {
            contig.clear();
            contig.SetName(parser.AsString(1));
            if (parser.GetItemCount() >= 3)
                offset = parser.AsInt(2);
            if (parser.GetItemCount() >= 4)
                contig.SetOri(parser.AsInt(3));
            continue;
        }

        if (parser.AsString(0) == "<CONTIG_READCOUNT>") {
            contig.SetNumReads(parser.AsInt(2));
            continue;
        }

        if (parser.AsString(0) == "<CONTIG_PAIRCOUNT>") {
            contig.SetNumPairs(parser.AsInt(2));
            continue;
        }

        if(parser.GetItemCount()>6) {
            int read = parser.AsInt(0);
            int ori = parser.AsInt(1);
            int start = parser.AsInt(2);
            int stop = parser.AsInt(4);
            int pair = parser.AsInt(5);
            int pairOrient = parser.AsInt(6);
            contig.push_back(ReadPlacement(read , ori, start, stop, pair, pairOrient, ""));
        }
    }  
}


void ContigScaffoldIO::Write(const Assembled & assembled, const string &file)
{
    FILE * pOut = fopen(file.c_str(), "w");
    if (pOut == NULL) {
        cout << "ERROR: Could not open file " << file << " for writing." << endl;
    }

    int count = 0;

    int i, j, k, l;
    for (l=0; l<assembled.isize(); l++) {
        const Scaffold & s = assembled[l];
        if (s.isize() == 0)
            continue;
        string name = s.Name();
        if (name == "" ) {
            //name = "<unknown>";
            char tmp[1024];
            sprintf(tmp, ">Scaffold_%5d", count);
            for (unsigned int x = 0; x<strlen(tmp); x++) {
                if (tmp[x] == ' ')
                    tmp[x] = '0';
            }
            name = tmp;
        }
        count++;
        fprintf(pOut, "<SCAFFOLD>\t%s\n", name.c_str());
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            string nameC = c.Name();
            if (nameC == "")
                nameC = "<unknown>";
            fprintf(pOut, "<CONTIG>\t%s\t%d\t%d\n", nameC.c_str(), s.Offset(i), c.Ori());
      
            for (j=0; j<c.isize(); j++) {
                const ReadPlacement & r = c[j];	
                fprintf(pOut, "%d\t%d\t%d - %d\t%d\t%d\n", r.Read(), r.Ori(), r.Start(), r.Stop(), r.Pair(), r.PairOrient());
            }

            fprintf(pOut, "<CONTIG_READCOUNT> %s %d </CONTIG_READCOUNT>\n", nameC.c_str(), c.NumReads());
            fprintf(pOut, "<CONTIG_PAIRCOUNT> %s %d </CONTIG_PAIRCOUNT>\n", nameC.c_str(), c.NumPairs());
            fprintf(pOut, "</CONTIG>\t%s\n", nameC.c_str());
        }
        fprintf(pOut, "<SCAFFOLD_READCOUNT> %s %d </SCAFFOLD_READCOUNT>\n", name.c_str(), s.NumUniqReads());
        fprintf(pOut, "<SCAFFOLD_PAIRCOUNT> %s %d </SCAFFOLD_PAIRCOUNT>\n", name.c_str(), s.NumUniqPairs());
        fprintf(pOut, "</SCAFFOLD>\t%s\n", name.c_str());
    }

    fclose(pOut);
}

void ContigScaffoldIO::WriteAssembledRawReads(const Assembled & assembled,  const ConsensOverlapUnit& COUnit, const string &outFile) {
    int i, j, k, l;
    FILE * pOutCurr = fopen(outFile.c_str(), "w");
    if (pOutCurr == NULL) {
        cout << "ERROR: Could not open file " << outFile << " for writing." << endl;
    }
    for (l=0; l<assembled.isize(); l++) {
        const Scaffold & s = assembled[l];
        if (s.isize() == 0)
            continue;
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                int readIdx = c[j].Read();	
                const svec<int>&  rawReadIdxs = COUnit.getConsReads().getConsMembers(readIdx);
                for(int rr=0; rr<rawReadIdxs.isize(); rr++) {
                    fprintf(pOutCurr, "%s\n", COUnit.getRawReadName(rawReadIdxs[rr]).c_str());
                }
            }
        }
    }
    fclose(pOutCurr);
}

 
void ContigScaffoldIO::WriteScaffoldReads(const Assembled & assembled,  const ConsensOverlapUnit& COUnit, const string &outDir) {
    int i, j, k, l;
    for (l=0; l<assembled.isize(); l++) {
        const Scaffold & s = assembled[l];
        if (s.isize() == 0)
            continue;
        FILE * pOutCurr = fopen((outDir+Stringify(l)+".reads").c_str(), "w");
        if (pOutCurr == NULL) {
            cout << "ERROR: Could not open file " << outDir+Stringify(l)+".reads" << " for writing." << endl;
            break;
        }
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            for (j=0; j<c.isize(); j++) {
                int readIdx = c[j].Read();	
                const DNAVector& dna = COUnit.getConsReadDNA(readIdx);
                fprintf(pOutCurr, "%s\n%s\n", dna.Name().c_str(), dna.AsString().c_str());
            }
        }
        fclose(pOutCurr);
    }
    cout<<"Total number of scaffold produced: " << l << endl;
}

void ContigScaffoldIO::WriteReadCountSummary(const Assembled & assembled, const string &file)
{
    FILE * pOut = fopen(file.c_str(), "w");
    if (pOut == NULL) {
        cout << "ERROR: Could not open file " << file << " for writing." << endl;
    }

    int i, l;
    for (l=0; l<assembled.isize(); l++) {
        const Scaffold & s = assembled[l];
        if (s.isize() == 0)
            continue;
        string name = s.Name();
        for (i=0; i<s.isize(); i++) {
            const Contig & c = s[i];
            string nameC = c.Name();
            if (nameC == "") nameC = "<unknown>";
            fprintf(pOut, "%s %d %d \n", nameC.c_str(), c.NumReads(), c.NumPairs());
        }
    }
    fclose(pOut);
}
