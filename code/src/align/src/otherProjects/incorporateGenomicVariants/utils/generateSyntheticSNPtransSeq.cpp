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
#include <hash_map>
#include <map>
#include <set>
#include <sstream>

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "../general/SNPinTranscript_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputAnnotationPath ";
		cout << "inputFormattedSNPfile outputSyntheticSNPtransSeqFile syntheticSeqLength" << endl;
		exit(1);
	}
	string syntheticSeqLengthStr = argv[5];
	int syntheticSeqLength = atoi(syntheticSeqLengthStr.c_str());

	cout << "creating results folder ...." << endl;
	string outputFilePath = argv[4];

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	////////////////// transcript info loading .... //////////
	string transcript_file_path = argv[2];
	cout << "start to load transcriptInfo " << endl;
	Transcript_Set* transcriptInfo = new Transcript_Set();
	ifstream transcript_file_ifs(transcript_file_path.c_str());
	string transcript_type_GAF = "GAF";
	transcriptInfo->extractTranscript(transcript_file_ifs, indexInfo, transcript_type_GAF);

	///////////////////////  load SNPs //////////////////////////////////
	string inputFormattedSNPfile = argv[3];
	SNPhash_Info tmpSNPhashInfo;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo.generateSNPhash_formattedSNPfile(inputFormattedSNPfile, indexInfo);

	////////////////  insert SNPs into transcripts //////////////////////
	cout << "start to insert SNPs 2 transcripts " << endl;
	SNPinTranscript_Info SNPinTranscriptInfo;
	SNPinTranscriptInfo.addSNP2transcript(transcriptInfo, tmpSNPhashInfo, indexInfo);

	//////////////  update chromStr with loaded SNPs //////////////
	cout << "start to update chromStr with loaded SNPs" << endl;
	tmpSNPhashInfo.updateChromStrWithSNP(indexInfo);		

	cout << "start to outputSyntheticTransSeq with SNP inserted... " << endl;
	SNPinTranscriptInfo.outputSimulatedTransReadSeqAroundSNP(
		outputFilePath, transcriptInfo, indexInfo, syntheticSeqLength, false);

	cout << "all jobs done ..." << endl;
	transcriptInfo->memoryFree();
	delete transcriptInfo;
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}