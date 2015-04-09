#ifndef CONTSCAFFIO_H_
#define CONTSCAFFIO_H_

#include "src/Ananas/ContScaff.h"
#include "src/Ananas/ConsensOverlapUnit.h"


class ContigScaffoldIO
{

 public:
  ContigScaffoldIO() {}

  void Read(Assembled & assembled, const string &file);
  void Write(const Assembled & assembled, const string &file);
  void WriteScaffoldReads(const Assembled & assembled,  
                          const ConsensOverlapUnit& COUnit, const string &outDir); 
  

};


#endif //CONTSCAFFIO_H_
