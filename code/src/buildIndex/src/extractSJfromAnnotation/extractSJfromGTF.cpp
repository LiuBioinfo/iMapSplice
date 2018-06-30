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
#include "read_block_test.h"
#include "index_info.h"
#include "normal_transcript_info.h"
#include "transcript_set.h"
#include "annotation_info.h"

using namespace std;

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
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "executable indexFolderPath inputGTFpath outputFolder" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputGTFpath = argv[2];
	string outputFolder = argv[3];

	indexFolderPath += "/";
	string indexParameterFile = indexFolderPath + "_parameter";
	ifstream indexParameter_ifs(indexParameterFile.c_str());


	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";
	string outputLog = outputFolder + "process.log";
	ofstream log_ofs(outputLog.c_str());

	string invalidChrGTF = outputFolder + "invalidChr.gtf";
	string validChrGTF = outputFolder + "validChr.gtf";
	string validChrExonGTF = outputFolder + "validChr.exon.gtf";
	string validChrNonExonGTF = outputFolder + "validChr.nonExon.gtf";

	string outputSJ = outputFolder + "output.SJ";

	log_ofs << "index: " << indexFolderPath << endl;
	log_ofs << "start to load index information ...." << endl;
	cout << "start to load index information ...." << endl;
	Index_Info* indexInfo = new Index_Info(indexParameter_ifs);
	log_ofs << "finish loading index information ...." << endl;
	cout << "finish loading index information ...." << endl;

	log_ofs << "start to filter out invalidChrGTF from original gtf file ...." << endl;
	cout << "start to filter out invalidChrGTF from original gtf file ...." << endl;
	ifstream ori_gtf_ifs(inputGTFpath.c_str());	
	ofstream invalidChrGTF_ofs(invalidChrGTF.c_str());
	ofstream validChrGTF_ofs(validChrGTF.c_str());	
	//ofstream SJ_ofs(outputSJ.c_str());
	filterOutInvalidChrGTF(ori_gtf_ifs, invalidChrGTF_ofs, validChrGTF_ofs, indexInfo);
	ori_gtf_ifs.close();
	invalidChrGTF_ofs.close();
	validChrGTF_ofs.close();
	log_ofs << "finish filtering out invalidChrGTF from original gtf file ...." << endl;
	cout << "finish filtering out invalidChrGTF from original gtf file ...." << endl;

	log_ofs << "start to filter out validChr non-exon gtf file ...." << endl;
	cout << "start to filter out validChr non-exon gtf file ...." << endl;
	ifstream validChrGTF_ifs(validChrGTF.c_str());
	ofstream validChrExonGTF_ofs(validChrExonGTF.c_str());
	ofstream validChrNonExonGTF_ofs(validChrNonExonGTF.c_str());
	filterValidChrNonExonGTF(validChrGTF_ifs, validChrNonExonGTF_ofs, validChrExonGTF_ofs);
	validChrGTF_ifs.close();
	validChrExonGTF_ofs.close();
	validChrNonExonGTF_ofs.close();
	log_ofs << "finish filtering out validChr non-exon gtf file ...." << endl;
	cout << "finish filtering out validChr non-exon gtf file ...." << endl;

	log_ofs << "start to generate transcriptSetInfo ..." << endl;
	cout << "start to generate transcriptSetInfo ..." << endl;
	ifstream gtf_ifs(validChrExonGTF.c_str());
	Transcript_Set* transcriptSetInfo = new Transcript_Set();
	string transcriptType = "GTF";
	transcriptSetInfo->extractTranscript(gtf_ifs, indexInfo, transcriptType);
	gtf_ifs.close();
	log_ofs << "finish generating transcriptSetInfo ..." << endl;
	cout << "finish generating transcriptSetInfo ..." << endl;
	log_ofs << "transcriptSetSize: " << transcriptSetInfo->returnTranscriptNum() << endl;;
	cout << "transcriptSetSize: " << transcriptSetInfo->returnTranscriptNum() << endl;;
	log_ofs << "start to insert SJs of transcriptSetInfo into annotationInfo ..." << endl;
	cout << "start to insert SJs of transcriptSetInfo into annotationInfo ..." << endl;
	Annotation_Info* annotationInfo = new Annotation_Info();
	annotationInfo->initiateHash_vec(indexInfo);
	annotationInfo->extractSJfromTranscriptSet(transcriptSetInfo, transcriptType, indexInfo);
	log_ofs << "finish inserting SJs of transcriptSetInfo into annotationInfo ..." << endl;
	cout << "finish inserting SJs of transcriptSetInfo into annotationInfo ..." << endl;
	log_ofs << "start to output SJs" << endl;
	cout << "start to output SJs" << endl;
	ofstream SJ_ofs(outputSJ.c_str());
	annotationInfo->outputSJ(SJ_ofs, indexInfo);
	SJ_ofs.close();
	log_ofs << "finish outputing SJs" << endl;
	cout << "finish outputing SJs" << endl;
	log_ofs.close();
	delete annotationInfo;
	transcriptSetInfo->memoryFree();
	delete transcriptSetInfo;
	delete indexInfo;
	return 0;
}