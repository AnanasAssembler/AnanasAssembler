#ifndef _REDUNDANCEREMOVAL__H_
#define _REDUNDANCEREMOVAL_H_

#include <string>
#include <vector>
#include <set>
#include "cola/src/cola/Cola.h"
#include "cola/src/fastAlign/AlignmentParams.h"

class RedundanceRemoval {
public:
  void removeRedundants(const AlignmentParams& aParams, const string& fastaFileName, const string& layoutFileName, int numThreads); 
private:
  void handleAlignments(svec<AlignmentInfo> currAlignments, set<string>& contigsToRemove, float identRatio); 
  void removeItems(const string& fastaFileName, const string& layoutFileName, set<string>& contigsToRemove); 
};

#endif // _REDUNDANCEREMOVAL_H_
