#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include<fstream>
#include "src/Cola/NSaligner.h"
#include "src/Cola/SWGAaligner.h"
#include "src/Cola/NSGAaligner.h"
#include "src/Cola/NOIAligner.h"
#include "src/Satsuma/RefineSatsuma.h"

//=====================================================================
void RefineSatsuma::alignAll() {
	AlignmentBlock last, curr, lastGap;
	while (curr.parse(satsumaParser, false)) {
		if (last.merge(curr)) {
			if(outputMode==0) {
				realignFout << "Merged blocks." << endl << last.getTargetChrom() << " : " 
			        	    << last.getTargetStart() << " - " << last.getTargetStop()-1;
			}
			lastGap = curr;
			continue;
		}
		if(last.getTargetChrom() == "" || last.getTargetStart() < 0 || last.getQueryStart() < 0) {
			last = curr;
			lastGap = curr;
			continue;
		}
		if(outputMode==0) { realignFout << "Re-aligning Satsuma block: " << endl; }
		alignBlock(last); 
		processGap(lastGap, curr);
		last = curr;
		lastGap = curr;
	}
	alignBlock(last);
}

void RefineSatsuma::alignBlock(const AlignmentBlock& aBlock) {
	printRealignInfo(aBlock);
	DNAVector q, t;
	q.SetToSubOf(query(aBlock.getQueryChrom()), aBlock.getQueryStart(),
		 aBlock.getQueryStop()-aBlock.getQueryStart());
	t.SetToSubOf(target(aBlock.getTargetChrom()), aBlock.getTargetStart(),
		 aBlock.getTargetStop()-aBlock.getTargetStart());
	if (aBlock.isReversed()) { q.ReverseComplement(); }

	alignSeqs(t, q, realignFout);
}

void RefineSatsuma::alignSeqs(const DNAVector& t, const DNAVector& q, ostream& sout) {
	Cola cola1 = Cola();
	cola1.createAlignment(t, q, aligner).print(outputMode, 1.0, sout, screenWidth); 
}

void RefineSatsuma::processGap(const AlignmentBlock& lastGap,
					    const AlignmentBlock& curr) {
	if (!curr.isCompatible(lastGap)) { return; } 
	int tGapStart = lastGap.getTargetStop();
	int tGapStop  = curr.getTargetStart();
	int qGapStart = lastGap.getQueryStop();
	int qGapStop  = curr.getQueryStart();
	if (curr.isReversed()) {
	    qGapStart = curr.getQueryStop();
	    qGapStop  = lastGap.getQueryStart(); 
	}
	switch(outputMode) {
		case 0:	
		gapFout << "Aligning gap: " << lastGap.getTargetStop() 
			<< " " << curr.getTargetStart() << ", " << qGapStart << " " << qGapStop << endl
			<< lastGap.getTargetChrom() << " " << curr.getTargetChrom() 
			<< " " << lastGap.getQueryChrom() << " " << curr.getQueryChrom() 
			<< " " << lastGap.getOrient() << " " << curr.getOrient() 
			<< tGapStop - tGapStart << " - "
			<< qGapStop - qGapStart << endl;
			break;
		case 1:
		gapFout << lastGap.getTargetChrom() << "," << tGapStart << "," 
			<< tGapStop << "," << tGapStop - tGapStart << "," 
			<< curr.getQueryChrom() << "," << qGapStart << "," 
			<< qGapStop << "," << qGapStop - qGapStart << ",";
	}

	// Check length of gap not to exceed certain limit
	if((tGapStop - tGapStart < MIN_GAP_SIZE || tGapStop - tGapStart > MAX_GAP_SIZE) 
	    || (qGapStop - qGapStart < MIN_GAP_SIZE || qGapStop - qGapStart > MAX_GAP_SIZE)) {
		gapFout << endl; //End line for this block
		return; // Gap will not be processed
	}
	DNAVector q, t;
	q.SetToSubOf(query(lastGap.getQueryChrom()), qGapStart, qGapStop - qGapStart);
	t.SetToSubOf(target(lastGap.getTargetChrom()), tGapStart, tGapStop - tGapStart);
	if (lastGap.isReversed()) { q.ReverseComplement(); }
	alignSeqs(t, q, gapFout);
}

void RefineSatsuma::printRealignInfo(const AlignmentBlock& aBlock) {
	switch(outputMode) {
		case 0:
		realignFout << "Alignment (full): " << aBlock.getTargetChrom() << " : " 
			    << aBlock.getTargetStart() << " - " << aBlock.getTargetStop()-1
	         	    << " vs. " << aBlock.getQueryChrom() << " : " 
	            	    << aBlock.getQueryStart() << " - " << aBlock.getQueryStop()-1 
	            	    << " " << aBlock.getOrient() << endl
	            	    << aBlock.getTargetStop() - aBlock.getTargetStart() << " - "
	            	    << aBlock.getQueryStop() - aBlock.getQueryStart() << endl;
		case 1: 
		realignFout << aBlock.getTargetChrom() << "," << aBlock.getTargetStart()
			    << "," << aBlock.getTargetStop() 
			    << "," << aBlock.getTargetStop() - aBlock.getTargetStart() 
			    << "," << aBlock.getQueryChrom()
			    << "," << aBlock.getQueryStart() << "," << aBlock.getQueryStop()
			    << "," << aBlock.getQueryStop() - aBlock.getQueryStart()
			    << "," << aBlock.getOrient() << ",";
	}
}
			
void RefineSatsuma::output(const string& str) {
	cout<<str;
}
