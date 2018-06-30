all: iMapSplice-phased iMapSplice-unphased getSNPmer-phased getSNPmer-unphased buildWholeGenome build2ndLevelIndex buildSNPmerIndex

OPTFLAGS = -O0 -g -rdynamic -std=c++0x
OUTPUTDIR = ./

CFLAGS += $(OPTFLAGS) -fopenmp

MapSplice:
	g++  $(CFLAGS) -o $(OUTPUTDIR)bin/MapSplice src/align/src/alignAll.cpp

iMapSplice:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -o $(OUTPUTDIR)bin/iMapSplice src/align/src/alignAll.cpp

iMapSplice-phased:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -D TWO_HAP -o $(OUTPUTDIR)bin/iMapSplice-phased src/align/src/alignAll.cpp

iMapSplice-unphased:
	g++ $(CFLAGS) -D PERSONALIZED_CHR_SEQ -D VARY_SNP_MER -o $(OUTPUTDIR)bin/iMapSplice-unphased src/align/src/alignAll.cpp

getSNPmer-phased:
	g++ $(CFLAGS) -o $(OUTPUTDIR)bin/getSNPmer-phased src/align/src/otherProjects/incorporateGenomicVariants/utils/generateSyntheticSNPgenomSeq.cpp

getSNPmer-unphased:
	g++ $(CFLAGS) -o $(OUTPUTDIR)bin/getSNPmer-unphased src/align/src/otherProjects/incorporateGenomicVariants/utils/generateSyntheticSNPgenomSeq_allComb_sampling.cpp

buildWholeGenome:
	g++  $(CFLAGS) -o $(OUTPUTDIR)bin/buildWholeGenome src/buildIndex/src/buildGlobalAndLocalIndex/buildWholeGenomeIndex.cpp

build2ndLevelIndex:
	g++  $(CFLAGS) -o $(OUTPUTDIR)bin/build2ndLevelIndex src/buildIndex/src/buildGlobalAndLocalIndex/buildAll2ndLevelIndex.cpp

buildSNPmerIndex:
	g++  $(CFLAGS) -o $(OUTPUTDIR)bin/buildSNPmerIndex src/buildIndex/src/buildGlobalAndLocalIndex/buildWholeGenomeIndex_mergedFa_noPreIndex.cpp