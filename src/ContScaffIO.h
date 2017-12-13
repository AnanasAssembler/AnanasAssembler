#ifndef CONTSCAFFIO_H_
#define CONTSCAFFIO_H_

#include "ContScaff.h"
#include "ConsensOverlapUnit.h"


//Forward declaration
class ConsensOverlapUnit;

class ContigScaffoldIO
{

public:
    ContigScaffoldIO() {}

    void Read(Assembled & assembled, const string &file);
    void Write(const Assembled & assembled, const string &file);
    void WriteScaffoldReads(const Assembled & assembled, const ConsensOverlapUnit& COUnit, const string &outDir); 
    void WriteAssembledRawReads(const Assembled & assembled,  const ConsensOverlapUnit& COUnit, const string &outFile); 
    void WriteReadCountSummary(const Assembled & assembled, const string &file);
};


#endif //CONTSCAFFIO_H_
