// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdlib.h>
#include <string>
#include <string.h>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>

#include "read_block_test.h"
#include "index_info.h"
#include "transcript_set.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		// inputIndexFolder -- real genome index folder; 
		// inputAnnotationFile -- annotation file;
		// outputRefFolder -- transcript seq output folder;
		// outputTranscriptInfoFile -- transcript information file;
		// format -- GAF / GFF/ GTF / BEER (...transcripts...txt)
		cout << "executable <inputIndexFolder> <inputAnnotationFile> <outputRefFolder> <outputTranscriptInfoFile> <Format>" << endl;
		exit(1);
	}

	string indexStr = argv[1];
	string inputAnnotationFileStr = argv[2];
	string outputRefFolderStr = argv[3];
	string outputTranscriptInfoFileStr = argv[4];
	string formatStr = argv[5];

	if((formatStr == "GAF")||(formatStr == "GFF")||(formatStr == "GTF")
		||(formatStr == "BEER"))
	{}
	else
	{
		cout << "annotation format error !" << endl;
		exit(1);
	}


	ifstream annotation_ifs(inputAnnotationFileStr.c_str());
	ofstream transcriptInfo_ofs(outputTranscriptInfoFileStr.c_str());


	cout << "start to load index file ..." << endl;

	string chrom_bit_file = indexStr; chrom_bit_file.append("/_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("/_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);

	cout << "index: " << indexStr << endl;
	cout << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	cout << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	cout << "start to load every chromosome" << endl;
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

	cout << "finish loading all index ..." << endl;

	Transcript_Set* transcript_set = new Transcript_Set();
	cout << "start to extractTranscript from annotation" << endl;
	transcript_set->extractTranscript(annotation_ifs, indexInfo, formatStr);
	cout << "transcript #: " << transcript_set->returnTranscriptNum() << endl;		
	cout << "start to output transcript info" << endl;
	transcript_set->outputTranscriptInfo(transcriptInfo_ofs);
	cout << "start to outputRefSeq" << endl;
	transcript_set->outputRefSeq(outputRefFolderStr, indexInfo, 50);

	return 0;
}