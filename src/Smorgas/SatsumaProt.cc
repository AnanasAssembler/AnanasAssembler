#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "extern/logger/log.h"
#include "src/Smorgas/SatsumaProt.h"

void SatsumaProt::init(string& refFile, string& dbName, string& queryFile) {
  trans.SetSize(params.protein_K);
  core.SetNumTables(1);
  core.SetTranslator(&trans);

  FILE_LOG(logINFO)<< "Reading ref & query...";
  readTarget(refFile, dbName);
  if (queryFile != "")
    readQuery(queryFile);

  if (params.bExhaust) {
    if (params.n_blocks == 0) {
      allhits.resize(ref.isize());    
      for (int i=0; i<ref.isize(); i++) {      
        allhits[i] = i; //.Set(i, 0);
      }
    } else {
      int size = (ref.isize() + params.n_blocks)/params.n_blocks;
      allhits.resize(size);    
      for (int i=size*params.block; i<size*(params.block+1); i++) {
        if (i >= ref.isize())
          break;
	allhits[i-size*params.block] = i; //.Set(i, 0);
      }
    }
  } else {
    contrib.resize(ref.isize(), 0);
  }

  if(params.prefilterType==1) { //Need to use quals for this filter
    FILE_LOG(logDEBUG2) << "Resetting quals." ;
    ResetQuals(ref);
  } else if(params.prefilterType==3) {
    FILE_LOG(logDEBUG2) << "Resetting quals." ;
    ResetQuals(ref);
    filter1Passed.resize(ref.isize(), false);
  }
}

 
void SatsumaProt::readTarget(string& refFile, string& dbName) { 
  if (refFile != "") {
    ref.Read(refFile);
  } else {
    FILE_LOG(logINFO) << "Reading database " << dbName;
    ref.ReadV(dbName);
  }
  if (!params.bExhaust)
    core.AddData(ref, true, 1);
  FILE_LOG(logINFO) << "done!";
}

void SatsumaProt::setQuery(const vecDNAVector & d)
{
  queryStream.SetSource(d);
}

void SatsumaProt::readQuery(string& queryFile) {
  queryStream.Open(queryFile);
}

void SatsumaProt::alignAll(ostream & o) { 
  int totalHits = 0;
  string name;
  //for (int i=0; i<query.isize() ; i++) {
  //DNAVector& d = query[i];
  //name = query.NameClean(i);
  DNAVector d;
  int i = -1;
  while (queryStream.GetNext(d, name)) {
    i++;
    if (!params.bExhaust) {
      switch(params.prefilterType) {
      case 1:
          filterTopHits_1(d);
          break;
      case 2:
          filterTopHits_2(d);
          break;
      case 3:
        filterTopHits_3(d);
        break;
      }
    } else { FILE_LOG(logDEBUG1)<<"Exhaustive search. "; }
    FILE_LOG(logDEBUG1) << endl << "Start XC alignments. " <<  GetTimeStatic(params.bQuiet);
    if (params.bSame) {
      allhits.clear();
      allhits.push_back(Hit(i, 10));
    }
    int failedAttempts = 0; //Keep track of failed attempts & stop if reach given limit
    bool hit = false;
    for (int x=0; x<allhits.isize() && failedAttempts<params.allowFailHits; x++) {
      if (!params.bExhaust && x >= params.limit) {
	FILE_LOG(logDEBUG1) << "Sequence " << name << ". "
                            << "Displaying the first " << params.limit 
                            << ", there are " << allhits.isize()-params.limit << " more.";
        break;
      }
      DNAVector & t = ref[allhits[x].Contig()];

      // Determine the size.
      int size = xc.Size(t.isize(), d.isize());
      svec<float> signal;

      CCSignalProtein target_sig, query_sig;      
      target_sig.SetSequence(d, size);
      query_sig.SetSequence(t, size);
      xc.CrossCorrelate(signal, target_sig, query_sig);

      SignalFilter filter;
      filter.Do(signal);

      XCDynProg dp;
      MultiProtein tt(d);
      MultiProtein dd(t);

      double pVal = dp.Align(tt, dd, filter, params.ungapped_cutoff, o);
      double eVal = pVal * ref.size();
	
      if (eVal< params.eValThresh) {
        hit = true;
        o << "Summary " << name << " vs. " << ref.NameClean(allhits[x].Contig()) << " score: " << pVal;
        o << " q-coords: " << dp.GetQueryStart() << " " << dp.GetQueryStop();
        o << " t-coords: " << dp.GetTargetStart() << " " << dp.GetTargetStop();
        o << " " << dp.GetIdentity() << " " << dp.GetPVal();
        o << endl << endl << endl;
      } else  {
        failedAttempts++;
      }
    }
    if(hit) {totalHits++; }
  }
  FILE_LOG(logINFO)<<"********************  "<<totalHits<<"     out of "<<i+1;    
}

void SatsumaProt::filterTopHits_1(DNAVector d) {
  int maxContrib = 0;
  int allCount = 0;
  allhits.clear();
  allhits.reserve(1000);

  int hitCount = 0;

    svec<int> hitstats;
    hitstats.resize(ref.isize(), 0);
 
  bool bFull = true;
  FILE_LOG(logDEBUG1) << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
  FILE_LOG(logDEBUG1) << "Looking up k-mers for " << "  size=" << d.size() << " " <<  GetTimeStatic(params.bQuiet);

  for (int j=0; j<=d.size()-params.protein_K; j++) {
    const svec<KmerAlignCoreRecord> & matches = core.GetMatchesDirectly(d, j);
    if (matches.isize() == 0)
      bFull = false;
    for (int y=0; y<matches.isize(); y++) {
/*TODO
if (params.bSelf) {
  if (i == matches[y].GetContig()) {
    continue;
  }
}
*/
      int diff = matches[y].GetPosition() - j; 
      if (hitCount >= hit_ref.isize()) {
        hit_ref.resize(hitCount + 1000000);
        hit_diff.resize(hitCount + 1000000);
      }
      hit_ref[hitCount] = matches[y].GetContig();
      hit_diff[hitCount] = diff;
      hitCount++;
      hitstats[matches[y].GetContig()]++;
      int max = SetQual(ref[matches[y].GetContig()], diff);
      if (max > contrib[matches[y].GetContig()])
        contrib[matches[y].GetContig()] = max;
      if (max > maxContrib)
	maxContrib = max;
      }
    }
    FILE_LOG(logDEBUG2) << "Generating stats: ";
    int min2 = 0;
    int min3 = 0;
    int min4 = 0;
    for (int j=0; j<hitstats.isize(); j++) {
      if (hitstats[j] >= 2)
        min2++;
      if (hitstats[j] >= 3)
        min3++;
      if (hitstats[j] >= 4)
        min4++;
    }
    FILE_LOG(logDEBUG2) << "Min2=" << min2 << " Min3=" << min3 << " Min4=" << min4;
    FILE_LOG(logDEBUG2) << "Got k-mers, count. " <<  GetTimeStatic(params.bQuiet);
    for (int j=0; j<contrib.isize(); j++) {
      if (contrib[j] >= max(maxContrib/2,4) || contrib[j] >= 8) //TODO Parmeterize
        allhits.push_back(Hit(j, contrib[j]));
       contrib[j] = 0;
    }
    FILE_LOG(logDEBUG2) << "Reset the quals. Total count=" << hitCount << "  " << GetTimeStatic(params.bQuiet);
    for (int j=0; j<hitCount; j++) {
      ResetOneQual(ref[hit_ref[j]], hit_diff[j]);
    }
    hitCount = 0;
    FILE_LOG(logDEBUG2) << "Done. " <<  GetTimeStatic(params.bQuiet);
    Sort(allhits);
}

void SatsumaProt::filterTopHits_2(DNAVector d) {
  allhits.clear();
  allhits.reserve(1000);
  FILE_LOG(logDEBUG2) << "Looking up k-mers for " << "  size=" << d.size() 
                      << " " <<  GetTimeStatic(params.bQuiet);
  const int lim1 = d.size() - params.protein_K;
  for (int j=0; j<lim1-params.filterWSlide; j+=params.filterWSlide) {
    const svec<KmerAlignCoreRecord> & matches = core.GetMatchesDirectly(d, j);
    const int lim2 = matches.isize();
    for (int y=0; y<lim2; y++) {
      contrib[matches[y].GetContig()]++;
    }
  }
  FILE_LOG(logDEBUG2) << "Got k-mers, count. " <<  GetTimeStatic(params.bQuiet);
  const int lim3 = contrib.isize();
  double threshold = 0.1/params.filterWSlide;
  for (int j=0; j<lim3; j++) {
    if ((float)contrib[j]/ref[j].size() >= threshold)
      allhits.push_back(Hit(j, contrib[j]));
    contrib[j] = 0;
  }
  Sort(allhits);
}


void SatsumaProt::filterTopHits_3(DNAVector d) {
  allhits.clear();
  allhits.reserve(1000);
  FILE_LOG(logDEBUG2) << "Looking up k-mers for " << "  size=" << d.size() << " " <<  GetTimeStatic(params.bQuiet);
  const int lim1 = d.size() - params.protein_K;
  for (int j=0; j<lim1-params.filterWSlide; j+=params.filterWSlide) {
    const svec<KmerAlignCoreRecord> & matches = core.GetMatchesDirectly(d, j);
    const int lim2 = matches.isize();
    for (int y=0; y<lim2; y++) {
      contrib[matches[y].GetContig()]++;
    }
  }
  FILE_LOG(logDEBUG2) << "Got k-mers, count. " <<  GetTimeStatic(params.bQuiet);
  const int lim3 = contrib.isize();
  double threshold = 0.05/params.filterWSlide;
  for (int j=0; j<lim3; j++) {
    if ((float)contrib[j]/ref[j].size() >= threshold)
      filter1Passed.at(j) = true;
    contrib[j] = 0;
  }
  int maxContrib = 0;
  int allCount = 0;
  allhits.clear();
  allhits.reserve(1000);
 
  int hitCount = 0;

  for (int j=0; j<=lim1; j++) {
    const svec<KmerAlignCoreRecord> & matches = core.GetMatchesDirectly(d, j);
    for (int y=0; y<matches.isize(); y++) {
      if(filter1Passed.at(matches[y].GetContig())) {
        int diff = matches[y].GetPosition() - j; 
        if (hitCount >= hit_ref.isize()) {
          hit_ref.resize(hitCount + 1000000);
	  hit_diff.resize(hitCount + 1000000);
	}
        hit_ref[hitCount] = matches[y].GetContig();
        hit_diff[hitCount] = diff;
        hitCount++;
        int max = SetQual(ref[matches[y].GetContig()], diff);
        if (max > contrib[matches[y].GetContig()])
	  contrib[matches[y].GetContig()] = max;
	if (max > maxContrib)
	  maxContrib = max;
      }
    }
  }
  for (int j=0; j<contrib.isize(); j++) {
    if (contrib[j] >= max(maxContrib/2,4) || contrib[j] >= 8)
      allhits.push_back(Hit(j, contrib[j]));
    contrib[j] = 0;
    filter1Passed.at(j) = false;
//TODO combine in pair to reduce access time
  }
  FILE_LOG(logDEBUG2) << "Reset the quals. Total count=" << hitCount << "  " << GetTimeStatic(params.bQuiet);
  for (int j=0; j<hitCount; j++) {
    ResetOneQual(ref[hit_ref[j]], hit_diff[j]);
  }
  hitCount = 0;
  FILE_LOG(logINFO) << "Done. " <<  GetTimeStatic(params.bQuiet);     
  Sort(allhits);
}

