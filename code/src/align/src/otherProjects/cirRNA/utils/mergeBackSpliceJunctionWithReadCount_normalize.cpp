// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input:  backSplice counts files from different groups or samples;
// output: read counts for merged backSplices (for some backSplices, read counts may be 0 for some samples)

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
#include "../../../general/index_info.h"
using namespace std;

typedef map<int, int> SJendPos2supNumMap;
typedef map<int, SJendPos2supNumMap > SJchrPosPairMap;

void generate_sampleSJchrPosPairMapVec(
	vector<SJchrPosPairMap>& sampleSJchrPosPairMapVec,
	int tmpBackSplice_chrNameInt,
	int tmpBackSplice_chrPosInt_start,
	int tmpBackSplice_chrPosInt_end,
	int tmpBackSplice_supportNumInt)
{
	SJchrPosPairMap::iterator tmpSJchrPosPairMapIter
		= sampleSJchrPosPairMapVec[tmpBackSplice_chrNameInt].find(
			tmpBackSplice_chrPosInt_start);
	if(tmpSJchrPosPairMapIter 
		== sampleSJchrPosPairMapVec[tmpBackSplice_chrNameInt].end())// SJstartPos notFound, addSJ
	{
		SJendPos2supNumMap tmpSJendPos2supNumMap;
		tmpSJendPos2supNumMap.insert(pair<int,int>(tmpBackSplice_chrPosInt_end, 
			tmpBackSplice_supportNumInt));
		sampleSJchrPosPairMapVec[tmpBackSplice_chrNameInt].insert(
			pair<int,SJendPos2supNumMap>(
				tmpBackSplice_chrPosInt_start, tmpSJendPos2supNumMap));		
	}	
	else// SJstartPos found, check SJendPos
	{
		SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter
			= (tmpSJchrPosPairMapIter->second).find(tmpBackSplice_chrPosInt_end);		
		if(tmpSJendPos2supNumMapIter 
			== (tmpSJchrPosPairMapIter->second).end()) // SJendPos notFound, addSJ
			(tmpSJchrPosPairMapIter->second).insert(pair<int,int>(
				tmpBackSplice_chrPosInt_end, tmpBackSplice_supportNumInt));
		else
		{
			cout << "error !, all SJs should not exist ......" << endl;
			exit(1);
		}
	}
}

void generate_mergedBackSpliceSitePairSetVec_sampleSJchrPosPairMapVec(
	string& sampleBackSpliceReadCountFilePath,
	vector< set< pair<int,int> > >& backSpliceSitePairSetVec,
	vector<SJchrPosPairMap>& sampleSJchrPosPairMapVec,
	Index_Info* indexInfo)
{
	ifstream sampleBackSpliceReadCount_ifs(sampleBackSpliceReadCountFilePath.c_str());
	while(!sampleBackSpliceReadCount_ifs.eof())
	{
		string tmpBackSpliceSiteReadCountStr;
		getline(sampleBackSpliceReadCount_ifs, tmpBackSpliceSiteReadCountStr);
		if((sampleBackSpliceReadCount_ifs.eof())||(tmpBackSpliceSiteReadCountStr == ""))
			break;
		//cout << "tmpBackSpliceSiteReadCountStr: " << tmpBackSpliceSiteReadCountStr << endl;
		vector<string> tmpBackSpliceSiteReadCountFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = tmpBackSpliceSiteReadCountStr.find("\t", startLoc);
			string tmpBackSpliceSiteReadCountFieldStr = tmpBackSpliceSiteReadCountStr.substr(
				startLoc, tabLoc - startLoc);
			tmpBackSpliceSiteReadCountFieldVec.push_back(tmpBackSpliceSiteReadCountFieldStr);
			//cout << "tmpBackSpliceSiteReadCountFieldStr: " << tmpBackSpliceSiteReadCountFieldStr << endl;
			startLoc = tabLoc + 1;
		}
		int lastTabLoc = tmpBackSpliceSiteReadCountStr.find("\t", startLoc);
		if(lastTabLoc == string::npos)
			tmpBackSpliceSiteReadCountFieldVec.push_back(tmpBackSpliceSiteReadCountStr.substr(startLoc));
		else
			tmpBackSpliceSiteReadCountFieldVec.push_back(tmpBackSpliceSiteReadCountStr.substr(
				startLoc, lastTabLoc - startLoc));
		string tmpBackSplice_chrNameStr = tmpBackSpliceSiteReadCountFieldVec[0];
		//cout << "tmpBackSplice_chrNameStr: " << tmpBackSplice_chrNameStr << endl;
		string tmpBackSplice_chrPosStr_start = tmpBackSpliceSiteReadCountFieldVec[1];
		//cout << "tmpBackSplice_chrPosStr_start: " << tmpBackSplice_chrPosStr_start << endl;
		string tmpBackSplice_chrPosStr_end = tmpBackSpliceSiteReadCountFieldVec[2];
		//cout << "tmpBackSplice_chrPosStr_end: " << tmpBackSplice_chrPosStr_end << endl;
		string tmpBackSplice_juncName = tmpBackSpliceSiteReadCountFieldVec[3];
		string tmpBackSplice_supportNumStr = tmpBackSpliceSiteReadCountFieldVec[4];
		int tmpBackSplice_chrNameInt = indexInfo->convertStringToInt(tmpBackSplice_chrNameStr);
		int tmpBackSplice_chrPosInt_start = atoi(tmpBackSplice_chrPosStr_start.c_str());
		int tmpBackSplice_chrPosInt_end = atoi(tmpBackSplice_chrPosStr_end.c_str());
		int tmpBackSplice_supportNumInt = atoi(tmpBackSplice_supportNumStr.c_str());
		// cout << "tmpBackSplice_chrNameInt: " << tmpBackSplice_chrNameInt << endl;
		// cout << "tmpBackSplice_chrPosInt_start: " << tmpBackSplice_chrPosInt_start << endl;
		// cout << "tmpBackSplice_chrPosInt_end: " << tmpBackSplice_chrPosInt_end << endl;
		// cout << "tmpBackSplice_supportNumInt: " << tmpBackSplice_supportNumInt << endl;
		backSpliceSitePairSetVec[tmpBackSplice_chrNameInt].insert(pair<int,int>(
			tmpBackSplice_chrPosInt_start, tmpBackSplice_chrPosInt_end));
		generate_sampleSJchrPosPairMapVec(sampleSJchrPosPairMapVec,
			tmpBackSplice_chrNameInt, tmpBackSplice_chrPosInt_start,
			tmpBackSplice_chrPosInt_end, tmpBackSplice_supportNumInt);
	}
	sampleBackSpliceReadCount_ifs.close();
}

void generate_mergedBackSpliceSitePairSetVec_sampleSJchrPosPairMapVecVec(
	vector<string>& sampleBackSpliceReadCountFilePathVec,
	vector< set< pair<int,int> > >& backSpliceSitePairSetVec,
	vector< vector<SJchrPosPairMap> >& sampleSJchrPosPairMapVecVec,	
	Index_Info* indexInfo)
{
	for(int tmp = 0; tmp < sampleBackSpliceReadCountFilePathVec.size(); tmp++)
	{
		//cout << "tmp: " << tmp << endl;
		//cout << "file name: " << endl << sampleBackSpliceReadCountFilePathVec[tmp] << endl;
		generate_mergedBackSpliceSitePairSetVec_sampleSJchrPosPairMapVec(
			sampleBackSpliceReadCountFilePathVec[tmp],
			backSpliceSitePairSetVec,
			sampleSJchrPosPairMapVecVec[tmp], indexInfo);
	}
}

int searchAndReturnBackSpliceSupportNum(
	vector<SJchrPosPairMap>& sampleSJchrPosPairMapVec,
	int tmpChrNameInt, int tmpChrPos_start, int tmpChrPos_end)
{
	SJchrPosPairMap::iterator tmpIter = sampleSJchrPosPairMapVec[tmpChrNameInt].find(tmpChrPos_start);
	if(tmpIter == sampleSJchrPosPairMapVec[tmpChrNameInt].end()) // start pos not found
		return 0;
	else
	{	
		SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter = (tmpIter->second).find(tmpChrPos_end);
		if(tmpSJendPos2supNumMapIter == (tmpIter->second).end()) // end pos not found
			return 0;
		else
			return (tmpSJendPos2supNumMapIter->second);
	}
}

void getDataNameReadNumFromStr(string& tmpDataName, int& tmpReadNum, string& tmpReadNumStr)
{
	int tabLoc = tmpReadNumStr.find("\t", 0);
	tmpDataName = tmpReadNumStr.substr(0, tabLoc);
	string tmpReadNumFieldStr = tmpReadNumStr.substr(tabLoc + 1);
	tmpReadNum = atoi(tmpReadNumFieldStr.c_str());
}

void generateReadNumVecFromFile(string& readNumFilePath, 
	vector<int>& readNumVec, vector<double>& normalizeFactorVec)
{
	ifstream readNum_ifs(readNumFilePath.c_str());
	while(!readNum_ifs.eof())
	{
		string tmpReadNumStr;
		getline(readNum_ifs, tmpReadNumStr);
		if(tmpReadNumStr == "")
			break;
		string tmpDataName;
		int tmpReadNum;
		getDataNameReadNumFromStr(tmpDataName, tmpReadNum, tmpReadNumStr);
		readNumVec.push_back(tmpReadNum);
	}
	readNum_ifs.eof();
	normalizeFactorVec.push_back(1.0);
	for(int tmp = 1; tmp < readNumVec.size(); tmp++)
	{
		double tmpNormalizeFactor = ((double)readNumVec[0])/((double)readNumVec[tmp]);
		normalizeFactorVec.push_back(tmpNormalizeFactor);
	}
}

int main(int argc, char** argv)
{
	if(argc < 8)
	{
		cout << "Executable inputIndexFolderPath outputFolderPath inputReadNumFile name_1 backSplice_readCount_1 ";
		cout << " name_2 backSplice_readCount_2 (name_3 backSplice_readCount_3 ...)" << endl;
		exit(1);
	}
	// cout << "start to initiate indexInfo" << endl;
	// string indexFolderPath = argv[1];
	// string indexStr = indexFolderPath;
	// indexStr += "/";
	// string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	// ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	// string indexParameterFileStr = indexFolderPath + "/_parameter";
	// ifstream parameter_ifs(indexParameterFileStr.c_str());
	// cout << "initiate indexInfo" << endl;
	// Index_Info* indexInfo = new Index_Info(parameter_ifs);
	// int chromNum = indexInfo->returnChromNum();

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();


	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputMergedBackSpliceCoveragePath = outputFolderStr + "mergedBackSplice.readCount";
	ofstream mergedBackSpliceCoverage_ofs(outputMergedBackSpliceCoveragePath.c_str());
	string outputMergedBackSpliceCoveragePath_normalize = outputFolderStr + "mergedBackSplice_normalize.readCount";
	ofstream mergedBackSpliceCoverage_normalize_ofs(outputMergedBackSpliceCoveragePath_normalize.c_str());

	int sampleSize = (argc - 4)/2;
	if(sampleSize * 2 + 4 != argc)
	{
		cout << "incorrect parameter number" << endl;
		exit(1);
	}

	cout << "start to extract read num in each data set" << endl;
	vector<int> readNumVec;
	vector<double> normalizeFactorVec;
	string readNumFilePath = argv[3];
	generateReadNumVecFromFile(readNumFilePath, readNumVec, normalizeFactorVec);

	cout << "start to generate sampleNameVec, sampleBackSpliceReadCountFilePathVec" << endl;
	vector<string> sampleNameVec;
	vector<string> sampleBackSpliceReadCountFilePathVec;
	vector< vector<SJchrPosPairMap> > sampleSJchrPosPairMapVecVec;
	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample ++)
	{
		string tmpSampleNameStr = argv[4+tmpSample*2];
		string tmpSampleBackSpliceReadCountFilePath = argv[5+tmpSample*2];
		sampleNameVec.push_back(tmpSampleNameStr);
		sampleBackSpliceReadCountFilePathVec.push_back(tmpSampleBackSpliceReadCountFilePath);
		vector<SJchrPosPairMap> tmpSJchrPosPairMapVec;
		for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
		{
			SJchrPosPairMap tmpSJchrPosPairMap;
			tmpSJchrPosPairMapVec.push_back(tmpSJchrPosPairMap);
		}
		sampleSJchrPosPairMapVecVec.push_back(tmpSJchrPosPairMapVec);
	}
	cout << "start to initiate backSpliceSitePairSetVec_merged" << endl;
	vector< set< pair<int,int> > > backSpliceSitePairSetVec_merged;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		set< pair<int,int> > tmpBackSpliceSitePairSet;
		backSpliceSitePairSetVec_merged.push_back(tmpBackSpliceSitePairSet);
	}
	cout << "start to generate backSpliceSitePairSetVec_merged" << endl;
	generate_mergedBackSpliceSitePairSetVec_sampleSJchrPosPairMapVecVec(
		sampleBackSpliceReadCountFilePathVec, backSpliceSitePairSetVec_merged,
		sampleSJchrPosPairMapVecVec, indexInfo);
	cout << "start to output backSplice read count info ..." << endl;
	for(int tmp = 0; tmp < backSpliceSitePairSetVec_merged.size(); tmp++)
	{
		for(set< pair<int,int> >::iterator tmpSetIter = backSpliceSitePairSetVec_merged[tmp].begin();
			tmpSetIter != backSpliceSitePairSetVec_merged[tmp].end(); tmpSetIter ++)
		{
			int tmpBackSplice_chrNameInt = tmp;
			int tmpBackSplice_chrPos_start = (*tmpSetIter).first;
			int tmpBackSplice_chrPos_end = (*tmpSetIter).second;
			string tmpBackSplice_chrNameStr = indexInfo->returnChrNameStr(tmpBackSplice_chrNameInt);
			//cout << "tmpBackSplice_chrNameStr: " << tmpBackSplice_chrNameStr << endl;
			//cout << "tmpBackSplice_chrPos_start: " << tmpBackSplice_chrPos_start << endl;
			//cout << "tmpBackSplice_chrPos_end: " << tmpBackSplice_chrPos_end << endl;
			string tmpBackSplice_flankString = indexInfo->returnBackSpliceFlankString(tmpBackSplice_chrNameInt,
				tmpBackSplice_chrPos_start, tmpBackSplice_chrPos_end);
			//cout << "tmpBackSplice_flankString: " << tmpBackSplice_flankString << endl;
			mergedBackSpliceCoverage_ofs << tmpBackSplice_chrNameStr << "\t"
				<< tmpBackSplice_chrPos_start << "\t" 
				<< tmpBackSplice_chrPos_end << "\t" << tmpBackSplice_flankString;
			mergedBackSpliceCoverage_normalize_ofs << tmpBackSplice_chrNameStr << "\t"
				<< tmpBackSplice_chrPos_start << "\t" 
				<< tmpBackSplice_chrPos_end << "\t" << tmpBackSplice_flankString;
			int tmpTotalBackSpliceNum = 0;
			for(int tmpSample = 0; tmpSample < sampleSize; tmpSample++)
			{
				//cout << "tmpSample: " << tmpSample << endl;
				int tmpBackSplice_supportNum = searchAndReturnBackSpliceSupportNum(
					sampleSJchrPosPairMapVecVec[tmpSample],
					tmpBackSplice_chrNameInt, tmpBackSplice_chrPos_start, tmpBackSplice_chrPos_end);
				double tmpBackSplice_supportNum_normalized 
					= ((double)tmpBackSplice_supportNum)*(normalizeFactorVec[tmpSample]); 
				mergedBackSpliceCoverage_ofs << "\t" << tmpBackSplice_supportNum;
				mergedBackSpliceCoverage_normalize_ofs << "\t" << tmpBackSplice_supportNum_normalized;
				tmpTotalBackSpliceNum += tmpBackSplice_supportNum;
			}
			mergedBackSpliceCoverage_ofs << "\t" << tmpTotalBackSpliceNum << endl;
			mergedBackSpliceCoverage_normalize_ofs << "\t" << tmpTotalBackSpliceNum << endl;
		}
	}

	delete indexInfo;
	mergedBackSpliceCoverage_ofs.close();
	mergedBackSpliceCoverage_normalize_ofs.close();
	return 0;
}