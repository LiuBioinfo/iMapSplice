// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>

#include "option_info.h"
#include "read_block_test.h"
#include "bwtmap_info.h"
#include "DoubleAnchorScore.h"
#include "sbndm.h"
#include "otherFunc.h"
#include "index_info.h"
#include "annotation_info.h"

using namespace std;

int main(int argc, char**argv)
{
	if(argc < 3)
	{
		cout << "execution <inputIndexFilePath> <inputAnnotationFile> <outputSpliceJunctionFile>" << endl;
		exit(1);
	}
	string indexStr = argv[1];
	string inputAnnotationFile = argv[2];
	string outputSpliceJunctionFile = argv[3];

	string logFile = outputSpliceJunctionFile + ".log";
	ofstream log_ofs(logFile.c_str());

	log_ofs << "start to load index file ..." << endl;

	string chrom_bit_file = indexStr; chrom_bit_file.append("/_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("/_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	log_ofs << "index: " << indexStr << endl;
	///////////////////////////////////////
 
	log_ofs << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	log_ofs << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	log_ofs << "start to load every chromosome" << endl;
	//chromStr[0] = 
	(indexInfo->chromStr).push_back((indexInfo->chromString).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
	(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[0]+1));
	
	for(int tmp = 1; tmp < indexInfo->chromNum; tmp++)
	{
		//chromStr[tmp] = 
		
		(indexInfo->chromStr).push_back((indexInfo->chromString).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));	
		(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));
	}	

	log_ofs << "finish loading all index ..." << endl;

	ofstream outputSJ_ofs(outputSpliceJunctionFile.c_str());
	ifstream annotation_ifs(inputAnnotationFile.c_str());

	AnnotatedSJ_Info* annotationSJ_info = new AnnotatedSJ_Info(indexInfo);
	annotationSJ_info->extractSJfromAnnotation(annotation_ifs, indexInfo);
	annotationSJ_info->outputSJ(outputSJ_ofs, indexInfo);

	delete annotationSJ_info;

	return 0;
}