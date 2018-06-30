// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: SJinfoTable from compareSJfromMultiSamples output
// output: unannotated SJ num, unannotated SJ supportingReadNum distribution in each region along the refGenome
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

#include "../general/read_block_test.h"
#include "../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable inputIndexInfoPath inputCompareSJfromMultiSampleResultsPath outputFolder sampleName_1 sampleName_2 ...." << endl;
		exit(1);
	}

	int binwidth = 3000000;

	string inputIndexInfoPath = argv[1];
	string inputCompareSJfromMultiSampleResultsPath = argv[2];
	int sampleNum = argc - 4;
	vector<string> sampleNameVec;
	for(int tmpSample = 0; tmpSample < sampleNum; tmpSample++)
	{
		string tmpSampleName = argv[4+tmpSample];
		sampleNameVec.push_back(tmpSampleName);
	}
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	cout << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;
	int chromNum = indexInfo->returnChromNum();

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_log = outputFolderStr + "log.txt";
	ofstream log_ofs(output_log.c_str());
	string output_annotatedSJ_num = outputFolderStr + "annotatedSJ_num_distribution.txt";
	string output_unannotatedSJ_num = outputFolderStr + "unannotatedSJ_num_distribution.txt";
	string output_annotatedSJ_supNum = outputFolderStr + "annotatedSJ_supNum_distribution.txt";
	string output_unannotatedSJ_supNum = outputFolderStr + "unannotatedSJ_supNum_distribution.txt";
	ofstream annotatedSJ_num_distribution_ofs(output_annotatedSJ_num.c_str());
	ofstream annotatedSJ_supNum_distribution_ofs(output_annotatedSJ_supNum.c_str());
	ofstream unannotatedSJ_num_distribution_ofs(output_unannotatedSJ_num.c_str());
	ofstream unannotatedSJ_supNum_distribution_ofs(output_unannotatedSJ_supNum.c_str());

	// initiate SJnum, SJsupNum vectors 
	cout << "start to initiate SJnum, SJsupNum vectors ...." << endl;
	log_ofs << "start to initiate SJnum, SJsupNum vectors ...." << endl; 
	vector< vector< vector<int> > > annotatedSJnumVecVecVec;
	vector< vector< vector<int> > > unannotatedSJnumVecVecVec;
	vector< vector< vector<int> > > annotatedSJsupNumVecVecVec;
	vector< vector< vector<int> > > unannotatedSJsupNumVecVecVec;

	for(int tmpSample = 0; tmpSample < sampleNum; tmpSample ++)
	{
		vector< vector<int> > tmpAnnotatedSJnumVecVec;
		vector< vector<int> > tmpUnannotatedSJnumVecVec;
		vector< vector<int> > tmpAnnotatedSJsupNumVecVec; 
		vector< vector<int> > tmpUnannotatedSJsupNumVecVec;
		for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
		{
			int tmpChrLength = indexInfo->returnChromLength(tmpChr);
			int tmpChr_binNum = tmpChrLength/binwidth + 1;
			vector<int> tmpAnnotatedSJnumVec;
			vector<int> tmpUnannotatedSJnumVec;
			vector<int> tmpAnnotatedSJsupNumVec;
			vector<int> tmpUnannotatedSJsupNumVec;
			for(int tmpBin = 0; tmpBin < tmpChr_binNum; tmpBin++)
			{
				tmpAnnotatedSJnumVec.push_back(0);
				tmpUnannotatedSJnumVec.push_back(0);
				tmpAnnotatedSJsupNumVec.push_back(0);
				tmpUnannotatedSJsupNumVec.push_back(0);
			}
			tmpAnnotatedSJnumVecVec.push_back(tmpAnnotatedSJnumVec);
			tmpUnannotatedSJnumVecVec.push_back(tmpUnannotatedSJnumVec);
			tmpAnnotatedSJsupNumVecVec.push_back(tmpAnnotatedSJsupNumVec); 
			tmpUnannotatedSJsupNumVecVec.push_back(tmpUnannotatedSJsupNumVec);						
		}
		annotatedSJnumVecVecVec.push_back(tmpAnnotatedSJnumVecVec);
		unannotatedSJnumVecVecVec.push_back(tmpUnannotatedSJnumVecVec);
		annotatedSJsupNumVecVecVec.push_back(tmpAnnotatedSJsupNumVecVec);
		unannotatedSJsupNumVecVecVec.push_back(tmpUnannotatedSJsupNumVecVec);
	}

	cout << "start to fill numbers in those vectors" << endl;
	log_ofs << "start to fill numbers in those vectors" << endl; 
	ifstream SJinfo_ifs(inputCompareSJfromMultiSampleResultsPath.c_str());
	int fieldNumInSJinfoTable = sampleNum + 4;
	while(!SJinfo_ifs.eof())
	{
		string tmpSJinfoStr;
		getline(SJinfo_ifs, tmpSJinfoStr);
		if((SJinfo_ifs.eof())||(tmpSJinfoStr == ""))
			break;
		vector<string> tmpSJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < fieldNumInSJinfoTable - 1; tmp++)
		{
			int tabLoc = tmpSJinfoStr.find("\t", startLoc);
			string tmpSJfieldStr = tmpSJinfoStr.substr(startLoc, tabLoc - startLoc);
			tmpSJfieldVec.push_back(tmpSJfieldStr);
			startLoc = tabLoc + 1;
		}
		tmpSJfieldVec.push_back(tmpSJinfoStr.substr(startLoc));
		string tmpSJ_chrNameStr = tmpSJfieldVec[0];
		string tmpSJ_chrPosStr_start = tmpSJfieldVec[1];
		string tmpSJ_chrPosStr_end = tmpSJfieldVec[2];
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrNameStr);
		int tmpSJ_chrPosInt_start = atoi(tmpSJ_chrPosStr_start.c_str());
		int tmpSJ_chrPosInt_end = atoi(tmpSJ_chrPosStr_end.c_str());
		int tmpSJ_chrPosInt_mid = (tmpSJ_chrPosInt_start + tmpSJ_chrPosInt_end)/2;
		int tmpSJ_chrPosInt_mid_binIndex = tmpSJ_chrPosInt_mid/binwidth;
		//cout << endl << "tmpSJ_chrNameStr: " << tmpSJ_chrNameStr << endl;
		//cout << "tmpSJ_chrNameInt: " << tmpSJ_chrNameInt << endl;
		//cout << "tmpSJ_chrPosInt_start: " << tmpSJ_chrPosInt_start << endl;
		//cout << "tmpSJ_chrPosInt_end: " << tmpSJ_chrPosInt_end << endl;
		string tmpSJ_annotatedOrNot = tmpSJfieldVec[3];
		//cout << "tmpSJ_annotatedOrNot: " << tmpSJ_annotatedOrNot << endl;
		bool annotatedOrNot_bool;
		if(tmpSJ_annotatedOrNot == "Y")
			annotatedOrNot_bool = true;
		else
			annotatedOrNot_bool = false;
		//vector<int> tmpSupNumForEachSampleVec;
		for(int tmpSample = 0; tmpSample < sampleNum; tmpSample ++)
		{
			//cout << "tmpSample: " << tmpSample << endl;
			string tmpSampleSupNumStr = tmpSJfieldVec[tmpSample + 4];
			int tmpSampleSupNumInt = atoi(tmpSampleSupNumStr.c_str());
			//cout << "tmpSampleSupNumInt: " << tmpSampleSupNumInt << endl;

			//tmpSupNumForEachSampleVec.push_back(tmpSampleSupNumInt);
			if(tmpSampleSupNumInt > 0)
			{
				if(annotatedOrNot_bool)
				{
					((annotatedSJnumVecVecVec[tmpSample])[tmpSJ_chrNameInt])[tmpSJ_chrPosInt_mid_binIndex] ++;
					((annotatedSJsupNumVecVecVec[tmpSample])[tmpSJ_chrNameInt])[tmpSJ_chrPosInt_mid_binIndex] += tmpSampleSupNumInt;
				}
				else
				{
					((unannotatedSJnumVecVecVec[tmpSample])[tmpSJ_chrNameInt])[tmpSJ_chrPosInt_mid_binIndex] ++;
					((unannotatedSJsupNumVecVecVec[tmpSample])[tmpSJ_chrNameInt])[tmpSJ_chrPosInt_mid_binIndex] += tmpSampleSupNumInt;
				}
			}
		}
	}
	SJinfo_ifs.close();

	cout << "start to output SJnum SJsupNum distribution along refGenome" << endl;
	log_ofs << "start to output SJnum SJsupNum distribution along refGenome" << endl;

	int accumulativeBinNum = 0;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr ++)
	{
		int tmpChrLength = indexInfo->returnChromLength(tmpChr);
		int tmpChr_binNum = tmpChrLength/binwidth;
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		for(int tmpBin = 0; tmpBin < tmpChr_binNum; tmpBin++)
		{
			int tmpAccumulativeBinNum = accumulativeBinNum + tmpBin + 1;
			int tmpBin_startPos = tmpBin * binwidth;
			int tmpBin_endPos = (tmpBin + 1) * binwidth - 1;
			for(int tmpSample = 0; tmpSample < sampleNum; tmpSample ++)
			{
				string tmpSampleName = sampleNameVec[tmpSample];
				int tmpAnnotatedSJnum 
					= ((annotatedSJnumVecVecVec[tmpSample])[tmpChr])[tmpBin];
				annotatedSJ_num_distribution_ofs 
					<< tmpAccumulativeBinNum << "\t" << tmpChrName << "\t"
					<< tmpBin_startPos << "\t" << tmpBin_endPos << "\t"
					<< tmpAnnotatedSJnum << "\t" 
					<< tmpSampleName << endl;
				int tmpAnnotatedSJsupNum 
					= ((annotatedSJsupNumVecVecVec[tmpSample])[tmpChr])[tmpBin];
				annotatedSJ_supNum_distribution_ofs 
					<< tmpAccumulativeBinNum << "\t" << tmpChrName << "\t"
					<< tmpBin_startPos << "\t" << tmpBin_endPos << "\t"
					<< tmpAnnotatedSJsupNum << "\t" 
					<< tmpSampleName << endl;
				int tmpUnannotatedSJnum 
					= ((unannotatedSJnumVecVecVec[tmpSample])[tmpChr])[tmpBin];
				unannotatedSJ_num_distribution_ofs 
					<< tmpAccumulativeBinNum << "\t" << tmpChrName << "\t"
					<< tmpBin_startPos << "\t" << tmpBin_endPos << "\t"
					<< tmpUnannotatedSJnum << "\t" 
					<< tmpSampleName << endl;				
				int tmpUnannotatedSJsupNum 
					= ((unannotatedSJsupNumVecVecVec[tmpSample])[tmpChr])[tmpBin];
				unannotatedSJ_supNum_distribution_ofs
					<< tmpAccumulativeBinNum << "\t" << tmpChrName << "\t"
					<< tmpBin_startPos << "\t" << tmpBin_endPos << "\t"
					<< tmpUnannotatedSJsupNum << "\t" 
					<< tmpSampleName << endl;	
			}
		}
		accumulativeBinNum += tmpChr_binNum;
	}
	return 0;
}