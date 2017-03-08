cmake_minimum_required(VERSION 3.3)
project(ananas)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -O3 -w")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

# include directory in find path where all dependency modules exist
include_directories(./)


# ananas binaries
set(SOURCE_FILES_ANANAS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/util/mutil.cc ryggrad/src/util/SysTime.cc src/Ananas.cc) 
set(SOURCE_FILES_BIN2ASCIIOVERLAPS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/util/mutil.cc src/Bin2AsciiOverlaps.cc src/ReadOverlap.cc) 
set(SOURCE_FILES_BUILDREADGROUPS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc ryggrad/src/util/SysTime.cc src/BuildReadGroups.cc) 
set(SOURCE_FILES_CLUSTERCONTIGS  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/ClusterContigs.cc src/ContigClusterUnit.cc src/ReadOverlap.cc src/Reads.cc) 
set(SOURCE_FILES_CONVERTREADCOUNTSTEMP  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/ContScaffIO.cc src/ConvertReadCounts_temp.cc src/ReadOverlap.cc src/Reads.cc) 
set(SOURCE_FILES_ENUMERATESEQS  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/aligns/KmerAlignCore.cc ryggrad/src/aligns/KmerDynProg.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/ContScaffIO.cc src/EnumerateSeqs.cc src/ReadOverlap.cc src/Reads.cc src/SequenceGraph.cc) 
set(SOURCE_FILES_FINDALTERNATIVESPLICING  ryggrad/src/aligns/KmerAlignCore.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/FindAlternativeSplicing.cc) 
set(SOURCE_FILES_GENASSEMBLYFASTA  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/aligns/KmerAlignCore.cc ryggrad/src/aligns/KmerDynProg.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/Consensus.cc src/ContScaffIO.cc src/GenAssemblyFasta.cc src/ReadOverlap.cc src/Reads.cc) 
set(SOURCE_FILES_KMERKONSENSUS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/KmerKonsensus.cc) 
set(SOURCE_FILES_LAYOUT  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/aligns/KmerAlignCore.cc ryggrad/src/aligns/KmerDynProg.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/ConsensOverlapUnit.cc src/Consensus.cc src/GlobUsage.cc src/Layout.cc src/ReadOverlap.cc src/Reads.cc src/SearchOverlaps.cc) 
set(SOURCE_FILES_LAYOUTGUIDED  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/Consensus.cc src/ContScaffIO.cc src/GlobUsage.cc src/LayoutGuided.cc src/ReadOverlap.cc src/Reads.cc src/SearchOverlaps.cc) 
set(SOURCE_FILES_OVERLAPSTAT  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/Consensus.cc src/OverlapStat.cc src/ReadOverlap.cc src/Reads.cc) 
set(SOURCE_FILES_REMOVEREDUNDANTSEQS  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/RemoveRedundantSeqs) 
set(SOURCE_FILES_RUNISOEM  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/ContScaffIO.cc src/ReadOverlap.cc src/Reads.cc src/RunIsoEM.cc) 
set(SOURCE_FILES_SCAFFOLDER  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/ContScaffIO.cc src/ReadOverlap.cc src/Reads.cc src/Scaffolder.cc) 
set(SOURCE_FILES_TESTCONSPAIRING  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/ContScaffIO.cc src/ReadOverlap.cc src/Reads.cc src/TestConsPairing.cc) 
set(SOURCE_FILES_FINDOVERLAPS  cola/src/cola/AlignmentCola.cc cola/src/cola/Cola.cc cola/src/cola/EditGraph.cc cola/src/cola/NSaligner.cc cola/src/cola/NSGAaligner.cc cola/src/cola/SWGAaligner.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc ryggrad/src/base/StringUtil.cc ryggrad/src/base/ThreadHandler.cc ryggrad/src/general/Alignment.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc  src/ConsensOverlapUnit.cc src/ContScaffIO.cc src/ReadOverlap.cc src/Reads.cc src/FindOverlaps.cc) 

add_executable(Ananas                  ${SOURCE_FILES_ANANAS})
add_executable(Bin2AsciiOverlaps       ${SOURCE_FILES_BIN2ASCIIOVERLAPS})
add_executable(BuildReadGroups         ${SOURCE_FILES_BUILDREADGROUPS})
add_executable(ClusterContigs          ${SOURCE_FILES_CLUSTERCONTIGS})
add_executable(ConvertReadCountsTemp   ${SOURCE_FILES_CONVERTREADCOUNTSTEMP})
add_executable(EnumerateSeqs           ${SOURCE_FILES_ENUMERATESEQS})
add_executable(FindAlternativeSplicing ${SOURCE_FILES_FINDALTERNATIVESPLICING})
add_executable(GenAssemblyFasta        ${SOURCE_FILES_GENASSEMBLYFASTA})
add_executable(KmerKonsensus           ${SOURCE_FILES_KMERKONSENSUS})
add_executable(Layout                  ${SOURCE_FILES_LAYOUT})
add_executable(LayoutGuided            ${SOURCE_FILES_LAYOUTGUIDED})
add_executable(OverlapStat             ${SOURCE_FILES_OVERLAPSTAT})
add_executable(RemoveRedundantSeqs     ${SOURCE_FILES_REMOVEREDUNDANTSEQS})
add_executable(RunIsoEM                ${SOURCE_FILES_RUNISOEM})
add_executable(Scaffolder              ${SOURCE_FILES_SCAFFOLDER})
add_executable(TestConsPairing         ${SOURCE_FILES_TESTCONSPAIRING})
add_executable(FindOverlaps            ${SOURCE_FILES_FINDOVERLAPS})

# ananas analysis binaries
set(SOURCE_FILES_ADDREADCOUNTS src/analysis/AddReadCounts.cc ryggrad/src/base/FileParser.cc ryggrad/src/util/mutil.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/StringUtil.cc) 
set(SOURCE_FILES_GENERATESAMPLEGENEMATRIX src/analysis/GenerateSampleGeneMatrix.cc src/analysis/SampleGeneMatrix.cc ryggrad/src/base/FileParser.cc ryggrad/src/util/mutil.cc ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/StringUtil.cc) 

add_executable(AddReadCounts              ${SOURCE_FILES_ADDREADCOUNTS})
add_executable(GenerateSampleGeneMatrix   ${SOURCE_FILES_GENERATESAMPLEGENEMATRIX})
