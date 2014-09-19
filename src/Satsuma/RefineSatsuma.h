#ifndef _REFINESATSUMA_H_
#define _REFINESATSUMA_H_

#include "src/DNAVector.h"
#include "src/AlignmentBlock.h"
#include "src/Cola/Cola.h"

// Maximum and minimum gap size allowed for alignment 
#define MAX_GAP_SIZE 15000 
#define MIN_GAP_SIZE 1 

//==================================================================
/**
 * This class facilitates taking satsuma alignment coordinates and 
 * the sequences that were used in the satsuma alignment and to postprocess
 * these with a more optimal alignment tool such as Cola
 * The class is templated over the refining alignment class which is assumed
 * to have an "align" function associated with it. 
 **/
class RefineSatsuma 
{
public:
	RefineSatsuma(const vecDNAVector& t, const vecDNAVector& q, 
	               const string& satsumaFile, 
                       ofstream& reFout, ofstream& gFout, int oMode
	             ):target(t), query(q), realignFout(reFout),
		       gapFout(gFout), outputMode(oMode), screenWidth(170),
		       aligner() 
	{
		satsumaParser.Open(satsumaFile);
	}
	
        /** 
         *  Set the Type of aligner to use and also the relevant parameters 
         *  If this function is not called the class will use its defaults,
         *  which is the NSGA aligner and its default parameters.
         */
        void setAligner(AlignerParams algn) { aligner = algn; }

	/** Realign the aligned regions and also attempt to align gaps 
	 * @param outputMode: 
	 */
	void alignAll(); 
private:
	/** 
	 * maps the start and end coordinates of an alignment block on to the 
	 * origniating query/target sequences and realigns those regions
	 */ 
	void alignBlock(const AlignmentBlock& aBlock); 

	/**
	 * Processes a gap by aligning it if it is between the min-max limits
	 */
	void processGap(const AlignmentBlock& lastGap, const AlignmentBlock& curr); 

	/** 
	 * Create an instance of the alignerType class and align the given sequences 
	 */
	void alignSeqs(const DNAVector& t, const DNAVector& q, ostream& sout=cout); 
	
	void printRealignInfo(const AlignmentBlock& aBlock); 
	void printGapInfo(const AlignmentBlock& lastGap, const AlignmentBlock& curr); 

	/** handle output of alignments and other info */
	// TODO handle output instead of scattered couts
	void output(const string& str);

	const vecDNAVector& target;   /// target genome contained in chromosomes
	const vecDNAVector& query;    /// query genome contained in chromosomes
	FlatFileParser satsumaParser; /// parser object containing the satsuma alignment coords
	ofstream& realignFout;        /// The outputstream for the realignments
	ofstream& gapFout;            /// The outputstream for the gap alignments
	int outputMode;               /// 0-full alignments and info 1- only info needed for extracting stats 
	int screenWidth;              /// Width to print alignments horizontally
	AlignerParams aligner;        /// Object containg aligner type and necessary parameters
};

#endif //_REFINESATSUMA_H_
