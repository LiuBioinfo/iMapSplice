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

//#include "option_info.h"
#include "../extractSJfromAnnotation/read_block_test.h"
#include "../extractSJfromAnnotation/index_info.h"
#include "../extractSJfromAnnotation/normal_transcript_info.h"
#include "../extractSJfromAnnotation/transcript_set.h"
#include "../extractSJfromAnnotation/annotation_info.h"

using namespace std;

/*
void filterOutInvalidChrGTF(ifstream& ori_gtf_ifs, ofstream& invalidChrGTF_ofs,
	ofstream& validChrGTF_ofs, Index_Info* indexInfo)
{
	while(!(ori_gtf_ifs.eof()))
	{
		string tmpLineStr;
		getline(ori_gtf_ifs, tmpLineStr);
		int	tabLoc = tmpLineStr.find("\t", 0);			
		string tmpChrStr = tmpLineStr.substr(0, tabLoc);

		int tmpChrInt = indexInfo->convertStringToInt(tmpChrStr);

		if(tmpChrInt == -1)
			invalidChrGTF_ofs << tmpLineStr << endl;
		else
			validChrGTF_ofs << tmpLineStr << endl;
	}
}

void filterValidChrNonExonGTF(ifstream& input_ifs, ofstream& nonExon_ofs, ofstream& exon_ofs)
{
	while(!(input_ifs.eof()))
	{
		string tmpLineStr;
		getline(input_ifs, tmpLineStr);
		vector<string> fieldStrVec;
		int searchStartLoc = 0;
		int tabLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			tabLoc = tmpLineStr.find("\t", searchStartLoc);			
			string tmpFieldStr = tmpLineStr.substr(searchStartLoc, tabLoc - 1 - searchStartLoc + 1);
			fieldStrVec.push_back(tmpFieldStr);
			searchStartLoc = tabLoc + 1;
		}
		if(fieldStrVec[2] == "exon")
			exon_ofs << tmpLineStr << endl;
		else
			nonExon_ofs << tmpLineStr << endl;
	}
}*/



int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "executable indexFolderPath inputBeersTranscriptInfoFile outputFolder annotation_name"<< endl;
		exit(1);
	}

	string indexFolderPath = argv[1];
	string inputBeersTranscriptInfoFile = argv[2];
	string outputFolder = argv[3];
	string annotation_name = argv[4];
	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";

	string invalidBeersTranscriptInfoFile = outputFolder + "invalidBeersTranscriptInfo.txt";
	string outputGTFfileStr = outputFolder + "output.gtf";
	string logFileStr = outputFolder + "BeersTranscriptInfo2gtf.log";

	//ofstream invalidBeersTranscriptInfo_ofs(invalidBeersTranscriptInfoFile.c_str());
	//ofstream outputGTF_ofs(outputGTFfileStr.c_str());
	ofstream log_ofs(logFileStr.c_str());

	indexFolderPath += "/";
	string indexParameterFile = indexFolderPath + "_parameter";
	ifstream indexParameter_ifs(indexParameterFile.c_str());

	log_ofs << "index: " << indexFolderPath << endl;
	log_ofs << "start to load index information ...." << endl;
	cout << "start to load index information ...." << endl;
	Index_Info* indexInfo = new Index_Info(indexParameter_ifs);
	log_ofs << "finish loading index information ...." << endl;
	cout << "finish loading index information ...." << endl;

	log_ofs << "start to generate transcriptSetInfo ..." << endl;
	cout << "start to generate transcriptSetInfo ..." << endl;
	Transcript_Set* transcriptSetInfo = new Transcript_Set();
	transcriptSetInfo->extractTranscript_BEER_outputInvalidTranscriptInfo(indexInfo,
		inputBeersTranscriptInfoFile, invalidBeersTranscriptInfoFile, log_ofs);
	log_ofs << "end of generating transcriptSetInfo ..." << endl;
	cout << "end of generating transcriptSetInfo ..." << endl;

	log_ofs << "start to output as GTF " << endl;
	cout << "start to output as GTF " << endl;
	ofstream GTF_ofs(outputGTFfileStr.c_str());
	transcriptSetInfo->output_GTF(GTF_ofs, annotation_name);

	GTF_ofs.close();
	indexParameter_ifs.close();
	log_ofs.close();
	delete indexInfo;
	return 0;
}