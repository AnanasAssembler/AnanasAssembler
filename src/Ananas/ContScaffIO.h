#ifndef CONTSCAFFIO_H_
#define CONTSCAFFIO_H_

#include "src/Ananas/ContScaff.h"


class ContigScaffoldIO
{

 public:
  ContigScaffoldIO() {}

  void Read(Assembled & assembled, const string &file);
  void Write(const Assembled & assembled, const string &file);
  

};


#endif //CONTSCAFFIO_H_
