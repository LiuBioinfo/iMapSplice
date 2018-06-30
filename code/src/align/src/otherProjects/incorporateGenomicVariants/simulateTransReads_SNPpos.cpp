// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//**************************************************************************************//
//step 1: take SNP files and indexes (chrom only) as input, 
//        and initiate the indexInfo, and get the new chrStr for each chromosome
//step 2: go through all the transcripts, and scan each transcript sequence, 
//        if it contains a SNP pos, get the simulated read from the new chrStr
//        be careful when generating the cigar strings.
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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
#include "./general/SNPinTranscript_info.h"

//#define SIMULATED_READ_LENGTH 51

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputAnnotationPath ";
		cout << "inputFormattedSNPfile outputSimulatedReadFolderPath" << endl;
		exit(1);
	}
	vector<int> readLengthVec;
	readLengthVec.push_back(21);
	readLengthVec.push_back(31);
	readLengthVec.push_back(41);
	readLengthVec.push_back(51);
	readLengthVec.push_back(61);
	readLengthVec.push_back(71);
	readLengthVec.push_back(81);
	readLengthVec.push_back(91);		
	readLengthVec.push_back(101);

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());

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

	cout << "start to outputSimulatedTransReadSeq with original bases... " << endl;
	bool outputFastaReadWithOriOrSNPbaseBool = true;
	for(int tmp = 0; tmp < readLengthVec.size(); tmp++)
	{
		int tmpReadLength = readLengthVec[tmp];
		string tmpReadLengthStr = int_to_str(tmpReadLength);
		string tmpOutputSimulatedReadFilePath = outputFolderStr + "simulatedRead_" + tmpReadLengthStr + "_ori" + ".fa";
		SNPinTranscriptInfo.outputSimulatedTransReadSeqAroundSNP(
			tmpOutputSimulatedReadFilePath, transcriptInfo, indexInfo, tmpReadLength, outputFastaReadWithOriOrSNPbaseBool);
	}
	outputFastaReadWithOriOrSNPbaseBool = false;
	cout << "start to outputSimulatedTransReadSeq with SNP inserted... " << endl; 
	for(int tmp = 0; tmp < readLengthVec.size(); tmp++)
	{
		int tmpReadLength = readLengthVec[tmp];
		string tmpReadLengthStr = int_to_str(tmpReadLength);
		string tmpOutputSimulatedReadFilePath = outputFolderStr + "simulatedRead_" + tmpReadLengthStr + "_SNP" + ".fa";
		SNPinTranscriptInfo.outputSimulatedTransReadSeqAroundSNP(
			tmpOutputSimulatedReadFilePath, transcriptInfo, indexInfo, tmpReadLength, outputFastaReadWithOriOrSNPbaseBool);
	}
	cout << "all jobs done ..." << endl;
	log_ofs.close();
	transcriptInfo->memoryFree();
	delete transcriptInfo;
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}