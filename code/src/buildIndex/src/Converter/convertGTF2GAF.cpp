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
		//cout << "tmpChrInt: " << tmpChrInt << endl;
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
	if(argc != 5)
	{
		cout << "executable indexFolderPath inputGTFfile outputFolder chromNameInGTFwithPrefixChrBool" << endl;
		exit(1);
	}

	string indexFolderPath = argv[1];
	string inputGTFpath = argv[2];
	string outputFolder = argv[3];

	bool chromNameInGTFwithPrefixChrBool = true;
	string chromNameInGTFwithPrefixChrBoolStr = argv[4];
	if((chromNameInGTFwithPrefixChrBoolStr == "True")||(chromNameInGTFwithPrefixChrBoolStr == "true")||(chromNameInGTFwithPrefixChrBoolStr == "True"))
		chromNameInGTFwithPrefixChrBool = true;
	else if((chromNameInGTFwithPrefixChrBoolStr == "FALSE")||(chromNameInGTFwithPrefixChrBoolStr == "false")||(chromNameInGTFwithPrefixChrBoolStr == "False"))
		chromNameInGTFwithPrefixChrBool = false;
	else
	{
		cout << "chromNameInGTFwithPrefixChrBoolStr: " << chromNameInGTFwithPrefixChrBoolStr << endl;
		exit(1);
	}

	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";

	string gtfFile_withPrefixChr;
	if(!chromNameInGTFwithPrefixChrBool)
	{
		cout << "start to add chr prefix to ori gtf file" << endl;
		string inputGTFpath_ori = inputGTFpath;
		string outputGTF_withPrefixChr = outputFolder + "chrPrefixAdded2oriGTF.gtf";
		ifstream oriGTF_ifs(inputGTFpath_ori.c_str());
		ofstream GTF_withPrefixChr_ofs(outputGTF_withPrefixChr.c_str());
		while(!oriGTF_ifs.eof())
		{
			string gtfStr;
			getline(oriGTF_ifs, gtfStr);
			string GTF_withPrefixChr_str = "chr" + gtfStr;
			GTF_withPrefixChr_ofs << GTF_withPrefixChr_str << endl;
		}
		oriGTF_ifs.close();
		GTF_withPrefixChr_ofs.close();
		gtfFile_withPrefixChr = outputGTF_withPrefixChr;
		cout << "finish to add chr prefix to ori gtf file" << endl;
	}
	else
		gtfFile_withPrefixChr = inputGTFpath;
	cout << "start to refine GTF with transcriptID and gene name" << endl;
	string gtfFile_withPrefixChr_transcriptID_geneName_only = outputFolder + "chrPrefixAdded2oriGTF_transcriptID_geneName.gtf";
	string gtfFile_withPrefixChr_withOut_transcriptID_geneName_only = outputFolder + "chrPrefixAdded2oriGTF_withOut_transcriptID_geneName.gtf";
	ofstream refinedGTF_withPrefixChr_transcriptID_geneName_ofs(gtfFile_withPrefixChr_transcriptID_geneName_only.c_str());
	ofstream refinedGTF_withPrefixChr_withOut_transcriptID_geneName_ofs(gtfFile_withPrefixChr_withOut_transcriptID_geneName_only.c_str());
	ifstream gtfFile_withPrefixChr_ifs(gtfFile_withPrefixChr.c_str());
	while(!gtfFile_withPrefixChr_ifs.eof())
	{
		string tmpGtfStr;
		getline(gtfFile_withPrefixChr_ifs,tmpGtfStr);
		int transcriptIDlocInStr = tmpGtfStr.find("transcript_id", 0);
		int geneNameLocInStr = tmpGtfStr.find("gene_name",0);
		int geneTypeLocInStr = tmpGtfStr.find("gene_biotype",0);
		if((transcriptIDlocInStr == string::npos)
			||(geneNameLocInStr == string::npos)
			||(geneTypeLocInStr == string::npos))
		{
			refinedGTF_withPrefixChr_withOut_transcriptID_geneName_ofs << tmpGtfStr << endl;
			continue;
		}
		int transcriptIDlocInStr_semicolon = tmpGtfStr.find(";", transcriptIDlocInStr + 1);
		int geneNameLocInStr_semicolon = tmpGtfStr.find(";", geneNameLocInStr + 1);
		int geneTypeLocInStr_semicolon = tmpGtfStr.find(";", geneTypeLocInStr + 1);
		string transcriptIDstr = tmpGtfStr.substr(transcriptIDlocInStr, 
			transcriptIDlocInStr_semicolon - transcriptIDlocInStr + 1);
		string geneNameStr = tmpGtfStr.substr(geneNameLocInStr,
			geneNameLocInStr_semicolon - geneNameLocInStr + 1);
		string geneTypeStr = tmpGtfStr.substr(geneTypeLocInStr,
			geneTypeLocInStr_semicolon - geneTypeLocInStr + 1);
		int startLoc = 0;
		for(int tmp = 0; tmp < 8; tmp++)
		{
			int tabLoc = tmpGtfStr.find("\t",startLoc);
			string tmpGtfStr_field = tmpGtfStr.substr(startLoc, tabLoc - 1 - startLoc + 1);
			refinedGTF_withPrefixChr_transcriptID_geneName_ofs << tmpGtfStr_field << "\t";
			startLoc = tabLoc + 1;
		}
		refinedGTF_withPrefixChr_transcriptID_geneName_ofs 
			<< transcriptIDstr << " " 
			<< geneNameStr << " " << geneTypeStr << endl;
	}
	gtfFile_withPrefixChr_ifs.close();
	refinedGTF_withPrefixChr_transcriptID_geneName_ofs.close();
	refinedGTF_withPrefixChr_withOut_transcriptID_geneName_ofs.close();
	inputGTFpath = gtfFile_withPrefixChr_transcriptID_geneName_only;

	cout << "finish to refine GTF with transcriptID and gene name" << endl;

	string invalidChrGTF = outputFolder + "invalidChr.gtf";
	string validChrGTF = outputFolder + "validChr.gtf";
	string validChrExonGTF = outputFolder + "validChr.exon.gtf";
	string validChrNonExonGTF = outputFolder + "validChr.nonExon.gtf";

	string outputGAFfileStr = outputFolder + "output.gaf";
	string outputGAFfileStr_standardID = outputGAFfileStr + ".reissuedID.gaf";
	string logFileStr = outputGAFfileStr + ".gtf2gaf.log";

	ofstream GAF_standardID_ofs(outputGAFfileStr_standardID.c_str());
	ofstream GAF_ofs(outputGAFfileStr.c_str());
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

	log_ofs << "start to output transcriptInfo in GAF format " << endl;
	cout << "start to output transcriptInfo in GAF format" << endl;
	transcriptSetInfo->output_GAF(GAF_ofs);
	transcriptSetInfo->output_GAF_standardID(GAF_standardID_ofs);
	log_ofs << "finish to output transcriptInfo in GAF format " << endl;
	cout << "finish to output transcriptInfo in GAF format" << endl;

	log_ofs << "start to output reissuedTranscriptID_geneName_oriTranscriptName_map" << endl; 
	cout << "start to output reissuedTranscriptID_geneName_oriTranscriptName_map" << endl; 
	string reissuedTranscriptID_geneName_oriTranscriptName_map_path = outputFolder + "reissuedTranscriptID_geneName_oriTranscriptName_map.txt";
	transcriptSetInfo->output_reissuedTranscriptID_geneName_oriTranscriptName(reissuedTranscriptID_geneName_oriTranscriptName_map_path);
	log_ofs << "finish to output reissuedTranscriptID_geneName_oriTranscriptName_map" << endl; 
	cout << "finish to output reissuedTranscriptID_geneName_oriTranscriptName_map" << endl; 

	GAF_standardID_ofs.close();
	GAF_ofs.close();
	log_ofs.close();
	//delete annotationInfo;
	transcriptSetInfo->memoryFree();
	delete transcriptSetInfo;
	
	//remove(invalidChrGTF.c_str());
	//remove(validChrGTF.c_str());
	//remove(validChrExonGTF.c_str());
	//remove(validChrNonExonGTF.c_str());

	return 0;
}