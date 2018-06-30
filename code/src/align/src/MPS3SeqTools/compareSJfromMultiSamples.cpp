// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
// input: annotated SJ file, SJ from each sample
// output: a table of overall info for all SJs
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

typedef map<int, int> SJendPos2supNumMap;
typedef map<int, SJendPos2supNumMap > SJchrPosSupNumMap;

typedef set<int> SJendPosSet;
typedef map<int, SJendPosSet> SJchrPosMap;

bool returnSJannotatedOrNot(
	vector<SJchrPosMap>& annotatedSJchrPosMapVec,
	int tmpChrNameInt, int tmpChrPos_start, int tmpChrPos_end)
{
	SJchrPosMap::iterator tmpSJchrPosMapIter
		= annotatedSJchrPosMapVec[tmpChrNameInt].find(tmpChrPos_start);
	if(tmpSJchrPosMapIter == annotatedSJchrPosMapVec[tmpChrNameInt].end())
		return false;
	else
	{
		SJendPosSet::iterator tmpSJendPosSetIter
			= (tmpSJchrPosMapIter->second).find(tmpChrPos_end);
		if(tmpSJendPosSetIter == (tmpSJchrPosMapIter->second).end())
			return false;
		else
			return true;
	}
}

int searchAndReturnSJsupportNum(
	vector<SJchrPosSupNumMap>& sampleSJchrPosSupNumMapVec,
	int tmpChrNameInt, int tmpChrPos_start, int tmpChrPos_end)
{
	SJchrPosSupNumMap::iterator tmpIter 
		= sampleSJchrPosSupNumMapVec[tmpChrNameInt].find(tmpChrPos_start);
	if(tmpIter == sampleSJchrPosSupNumMapVec[tmpChrNameInt].end()) // start pos not found
		return 0;
	else
	{	
		SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter 
			= (tmpIter->second).find(tmpChrPos_end);
		if(tmpSJendPos2supNumMapIter == (tmpIter->second).end()) // end pos not found
			return 0;
		else
			return (tmpSJendPos2supNumMapIter->second);
	}
}

void outputSJchrPosMapVec(
	vector<SJchrPosMap>& tmpSJchrPosMapVec,
	string& outputPath, Index_Info* indexInfo)
{
	ofstream tmp_ofs(outputPath.c_str());
	for(int tmpChr = 0; tmpChr < tmpSJchrPosMapVec.size(); tmpChr++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		for(SJchrPosMap::iterator tmpSJchrPosMapIter 
			= tmpSJchrPosMapVec[tmpChr].begin();
			tmpSJchrPosMapIter != tmpSJchrPosMapVec[tmpChr].end();
			tmpSJchrPosMapIter ++)
		{
			int tmpSJ_chrPosInt_start = tmpSJchrPosMapIter->first;
			for(SJendPosSet::iterator tmpSJendPosSetIter
				= (tmpSJchrPosMapIter->second).begin();
				tmpSJendPosSetIter != (tmpSJchrPosMapIter->second).end();
				tmpSJendPosSetIter ++)
			{
				int tmpSJ_chrPosInt_end = (*tmpSJendPosSetIter);
				tmp_ofs << tmpChrName << "\t"
					<< tmpSJ_chrPosInt_start << "\t"
					<< tmpSJ_chrPosInt_end << endl;
			}
		}
	}
	tmp_ofs.close();
}

void outputSJchrPosSupNumMapVec(
	vector<SJchrPosSupNumMap>& tmpSJchrPosSupNumMapVec,
	string& outputPath, Index_Info* indexInfo)
{
	ofstream tmp_ofs(outputPath.c_str());
	int tmpJuncNum = 0;
	for(int tmpChr = 0; tmpChr < tmpSJchrPosSupNumMapVec.size();
		tmpChr ++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		for(SJchrPosSupNumMap::iterator tmpSJchrPosSupNumMapIter
			= tmpSJchrPosSupNumMapVec[tmpChr].begin();
			tmpSJchrPosSupNumMapIter != tmpSJchrPosSupNumMapVec[tmpChr].end();
			tmpSJchrPosSupNumMapIter ++)
		{
			int tmpSJ_chrPosInt_start = tmpSJchrPosSupNumMapIter->first;
			for(SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter
				= (tmpSJchrPosSupNumMapIter->second).begin();
				tmpSJendPos2supNumMapIter != (tmpSJchrPosSupNumMapIter->second).end();
				tmpSJendPos2supNumMapIter ++)
			{
				int tmpSJ_chrPosInt_end = tmpSJendPos2supNumMapIter->first;
				int tmpSJ_supportNumInt = tmpSJendPos2supNumMapIter->second;
				tmpJuncNum ++;
				tmp_ofs << tmpChrName << "\t" << tmpSJ_chrPosInt_start << "\t"
					<< tmpSJ_chrPosInt_end << "\tJUNC_" << tmpJuncNum << "\t"
					<< tmpSJ_supportNumInt << endl;
			}
		}
	}
	tmp_ofs.close();
}

void insert2SJchrPosMapVec(
	vector<SJchrPosMap>& targetSJchrPosMapVec,
	int tmpSJ_chrNameInt, int tmpSJ_chrPosInt_start,
	int tmpSJ_chrPosInt_end)
{
	SJchrPosMap::iterator tmpSJchrPosMapIter
		= targetSJchrPosMapVec[tmpSJ_chrNameInt].find(
			tmpSJ_chrPosInt_start);
	if(tmpSJchrPosMapIter 
		== targetSJchrPosMapVec[tmpSJ_chrNameInt].end())
	{
		SJendPosSet tmpSJendPosSet;
		tmpSJendPosSet.insert(tmpSJ_chrPosInt_end);
		targetSJchrPosMapVec[tmpSJ_chrNameInt].insert(
			pair<int, set<int> >(tmpSJ_chrPosInt_start,
				tmpSJendPosSet));
	}
	else
	{
		SJendPosSet::iterator tmpSJendPosSetIter
			= (tmpSJchrPosMapIter->second).find(tmpSJ_chrPosInt_end);
		if(tmpSJendPosSetIter == (tmpSJchrPosMapIter->second).end())
		{
			(tmpSJchrPosMapIter->second).insert(tmpSJ_chrPosInt_end);
		}
		else
		{}
	}
}

void insert2sampleSJchrPosSupNumMapVec(
	vector<SJchrPosSupNumMap>& sampleSJchrPosSupNumMapVec,
	int tmpSJ_chrNameInt,
	int tmpSJ_chrPosInt_start,
	int tmpSJ_chrPosInt_end,
	int tmpSJ_supportNumInt)
{
	SJchrPosSupNumMap::iterator tmpSJchrPosSupNumMapIter
		= sampleSJchrPosSupNumMapVec[tmpSJ_chrNameInt].find(
			tmpSJ_chrPosInt_start);
	if(tmpSJchrPosSupNumMapIter
		== sampleSJchrPosSupNumMapVec[tmpSJ_chrNameInt].end())// SJstartPos notFound, addSJ
	{
		SJendPos2supNumMap tmpSJendPos2supNumMap;
		tmpSJendPos2supNumMap.insert(pair<int,int>(tmpSJ_chrPosInt_end, 
			tmpSJ_supportNumInt));
		sampleSJchrPosSupNumMapVec[tmpSJ_chrNameInt].insert(
			pair<int,SJendPos2supNumMap>(
				tmpSJ_chrPosInt_start, tmpSJendPos2supNumMap));		
	}	
	else// SJstartPos found, check SJendPos
	{
		SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter
			= (tmpSJchrPosSupNumMapIter->second).find(tmpSJ_chrPosInt_end);		
		if(tmpSJendPos2supNumMapIter 
			== (tmpSJchrPosSupNumMapIter->second).end()) // SJendPos notFound, addSJ
			(tmpSJchrPosSupNumMapIter->second).insert(pair<int,int>(
				tmpSJ_chrPosInt_end, tmpSJ_supportNumInt));
		else
		{
			cout << "error !, no SJ should exist ......" << endl;
			exit(1);
		}
	}
}

void mergeSJchrPosSupNumMapVecVec2SJchrPosMapVec(
	vector<SJchrPosMap>& targetSJchrPosMapVec,
	vector< vector<SJchrPosSupNumMap> >& sampleSJchrPosSupNumMapVecVec)
{
	int chromNum = targetSJchrPosMapVec.size();
	for(int tmpSample = 0; tmpSample < sampleSJchrPosSupNumMapVecVec.size();
		tmpSample ++)
	{
		for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
		{
			for(SJchrPosSupNumMap::iterator tmpSJchrPosSupNumMapIter 
				= ((sampleSJchrPosSupNumMapVecVec[tmpSample])[tmpChr]).begin();
				tmpSJchrPosSupNumMapIter != ((sampleSJchrPosSupNumMapVecVec[tmpSample])[tmpChr]).end();
				tmpSJchrPosSupNumMapIter ++)
			{	
				int tmpSJ_chrPosInt_start = tmpSJchrPosSupNumMapIter->first;
				for(SJendPos2supNumMap::iterator tmpSJendPos2supNumMapIter
					= (tmpSJchrPosSupNumMapIter->second).begin();
					tmpSJendPos2supNumMapIter != (tmpSJchrPosSupNumMapIter->second).end();
					tmpSJendPos2supNumMapIter ++)
				{	
					int tmpSJ_chrPosInt_end = tmpSJendPos2supNumMapIter->first;
					insert2SJchrPosMapVec(targetSJchrPosMapVec,
						tmpChr, tmpSJ_chrPosInt_start, tmpSJ_chrPosInt_end);
				}
			}
		}
	}
}

void generateAnnotatedSJchrPosMapVec(
	vector<SJchrPosMap>& annotatedSJchrPosMapVec,
	string& annotatedSJpath, Index_Info* indexInfo)
{
	ifstream annotatedSJ_ifs(annotatedSJpath.c_str());
	while(!annotatedSJ_ifs.eof())
	{
		string tmpSJstr;
		getline(annotatedSJ_ifs, tmpSJstr);
		if((annotatedSJ_ifs.eof())||(tmpSJstr == ""))
			break;		
		if(tmpSJstr.substr(0,3) != "chr")
			continue;
		vector<string> tmpSJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = tmpSJstr.find("\t", startLoc);
			string tmpSJfieldStr = tmpSJstr.substr(startLoc, tabLoc - startLoc);
			tmpSJfieldVec.push_back(tmpSJfieldStr);
			startLoc = tabLoc + 1;
		}
		if((tmpSJstr.substr(startLoc)).find("\t") == string::npos)
			tmpSJfieldVec.push_back(tmpSJstr.substr(startLoc));
		else
		{
			int tabLoc = (tmpSJstr.substr(startLoc)).find("\t", 0);
			tmpSJfieldVec.push_back((tmpSJstr.substr(startLoc)).substr(0, tabLoc));
		}

		string tmpSJ_chrNameStr = tmpSJfieldVec[0];
		string tmpSJ_chrPosStr_start = tmpSJfieldVec[1];
		string tmpSJ_chrPosStr_end = tmpSJfieldVec[2];
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrNameStr);
		int tmpSJ_chrPosInt_start = atoi(tmpSJ_chrPosStr_start.c_str());
		int tmpSJ_chrPosInt_end = atoi(tmpSJ_chrPosStr_end.c_str());
		insert2SJchrPosMapVec(annotatedSJchrPosMapVec, tmpSJ_chrNameInt, 
			tmpSJ_chrPosInt_start, tmpSJ_chrPosInt_end);
	}
	annotatedSJ_ifs.close();
}

void generateSampleSJchrPosSupNumMapVec(
	vector<SJchrPosSupNumMap>& sampleSJchrPosSupNumMapVec,
	string& sampleSJpath, Index_Info* indexInfo)
{
	ifstream sampleSJ_ifs(sampleSJpath.c_str());
	while(!sampleSJ_ifs.eof())
	{
		string tmpSJstr;
		getline(sampleSJ_ifs, tmpSJstr);
		//cout << "tmpSJstr: " << tmpSJstr << endl;
		if((sampleSJ_ifs.eof())||(tmpSJstr == ""))
			break;
		//cout << "tmpSJstr: " << tmpSJstr << endl;
		if(tmpSJstr.substr(0,3) != "chr")
			continue;
		vector<string> tmpSJfieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 4; tmp++)
		{
			int tabLoc = tmpSJstr.find("\t", startLoc);
			string tmpSJfieldStr = tmpSJstr.substr(startLoc, tabLoc - startLoc);
			tmpSJfieldVec.push_back(tmpSJfieldStr);
			startLoc = tabLoc + 1;
		}
		if((tmpSJstr.substr(startLoc)).find("\t") == string::npos)
			tmpSJfieldVec.push_back(tmpSJstr.substr(startLoc));
		else
		{
			int tabLoc = (tmpSJstr.substr(startLoc)).find("\t", 0);
			tmpSJfieldVec.push_back((tmpSJstr.substr(startLoc)).substr(0, tabLoc));
		}

		string tmpSJ_chrNameStr = tmpSJfieldVec[0];
		string tmpSJ_chrPosStr_start = tmpSJfieldVec[1];
		string tmpSJ_chrPosStr_end = tmpSJfieldVec[2];
		string tmpSJ_supportNumStr = tmpSJfieldVec[4];
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrNameStr);
		int tmpSJ_chrPosInt_start = atoi(tmpSJ_chrPosStr_start.c_str());
		int tmpSJ_chrPosInt_end = atoi(tmpSJ_chrPosStr_end.c_str());
		int tmpSJ_supportNumInt = atoi(tmpSJ_supportNumStr.c_str());
		insert2sampleSJchrPosSupNumMapVec(
			sampleSJchrPosSupNumMapVec, tmpSJ_chrNameInt,
			tmpSJ_chrPosInt_start, tmpSJ_chrPosInt_end, tmpSJ_supportNumInt);
	}
	sampleSJ_ifs.close();
}


int main(int argc, char** argv)
{
	if(argc < 8)
	{
		cout << "Executable inputIndexFolderPath inputAnnotatedSJ outputFolderPath ";
		cout << "name_1 SJ_1 name_2 SJ_2 (name_3 SJ_3 ....." << endl;
		exit(1);
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
	int chromNum = indexInfo->returnChromNum();

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	int sampleSize = (argc - 4)/2;
	if(sampleSize * 2 + 4 != argc)
	{
		cout << "incorrect parameter number" << endl;
		exit(1);
	}
	cout << "start to initiate sampleNameVec, sampleSJpathVec";
	cout << " sampleSJchrPosSupNumMapVecVec" << endl;
	vector<string> sampleNameVec;
	vector<string> sampleSJpathVec;
	vector< vector<SJchrPosSupNumMap> > sampleSJchrPosSupNumMapVecVec;
	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample ++)
	{
		string tmpSampleNameStr = argv[4+tmpSample*2];
		string tmpSampleSJpath = argv[5+tmpSample*2];
		sampleNameVec.push_back(tmpSampleNameStr);
		sampleSJpathVec.push_back(tmpSampleSJpath);
		vector<SJchrPosSupNumMap> tmpSJchrPosSupNumMapVec;
		for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
		{
			SJchrPosSupNumMap tmpSJchrPosSupNumMap;
			tmpSJchrPosSupNumMapVec.push_back(tmpSJchrPosSupNumMap);
		}
		sampleSJchrPosSupNumMapVecVec.push_back(tmpSJchrPosSupNumMapVec);
	}
	cout << "start to initiate mergedSJchrPosMapVec, annotatedSJchrPosMapVec" << endl;
	vector<SJchrPosMap> mergedSJchrPosMapVec;
	vector<SJchrPosMap> annotatedSJchrPosMapVec;	
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		SJchrPosMap tmpSJchrPosMap_merged;
		mergedSJchrPosMapVec.push_back(tmpSJchrPosMap_merged);
		SJchrPosMap tmpSJchrPosMap_annotated;
		annotatedSJchrPosMapVec.push_back(tmpSJchrPosMap_annotated);
	}
	cout << "start to generate annotatedSJ chrPosMapVec" << endl;
	string inputAnnotatedSJpath = argv[2];
	string outputAnnotatedSJpath = outputFolderStr + "annotatedSJ.junc";
	generateAnnotatedSJchrPosMapVec(annotatedSJchrPosMapVec,
		inputAnnotatedSJpath, indexInfo);

	cout << "start to output annotatedSJ chrPosMapVec" << endl;
	outputSJchrPosMapVec(annotatedSJchrPosMapVec,
		outputAnnotatedSJpath, indexInfo);

	cout << "start to generate sampleSJchrPosSupNumMapVecVec and output each ..." << endl;
	for(int tmpSample = 0; tmpSample < sampleSize; tmpSample ++)
	{
		cout << "tmpSampleSJpath: " << sampleSJpathVec[tmpSample] << endl;
		generateSampleSJchrPosSupNumMapVec(
			sampleSJchrPosSupNumMapVecVec[tmpSample],
			sampleSJpathVec[tmpSample], indexInfo);
		string tmpSJoutputPath = outputFolderStr + sampleNameVec[tmpSample] + ".junc";
		outputSJchrPosSupNumMapVec(
			sampleSJchrPosSupNumMapVecVec[tmpSample], tmpSJoutputPath, indexInfo);
	}

	cout << "start to mergeSJchrPosSupNumMapVecVec2SJchrPosMapVec " << endl;
	mergeSJchrPosSupNumMapVecVec2SJchrPosMapVec(
		mergedSJchrPosMapVec, sampleSJchrPosSupNumMapVecVec);

	cout << "start to output mergedSJchrPosMapVec " << endl;
	string outputMergedSJpath = outputFolderStr + "mergedSJ.junc";
	outputSJchrPosMapVec(mergedSJchrPosMapVec, outputMergedSJpath, indexInfo);

	cout << "start to creat and output annotatedOrNot info and supNum info for all juncs" << endl;
	string outputSJ_annotatedOrNot_supportNumInEachSample_path
		= outputFolderStr + "SJinfoTable.txt";
	ofstream SJinfoTable_ofs(outputSJ_annotatedOrNot_supportNumInEachSample_path.c_str());
	for(int tmpChr = 0; tmpChr < mergedSJchrPosMapVec.size(); tmpChr++)
	{
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		for(SJchrPosMap::iterator tmpSJchrPosMapIter 
			= mergedSJchrPosMapVec[tmpChr].begin();
			tmpSJchrPosMapIter != mergedSJchrPosMapVec[tmpChr].end();
			tmpSJchrPosMapIter ++)
		{
			int tmpSJ_chrPosInt_start = tmpSJchrPosMapIter->first;
			for(SJendPosSet::iterator tmpSJendPosSetIter
				= (tmpSJchrPosMapIter->second).begin();
				tmpSJendPosSetIter != (tmpSJchrPosMapIter->second).end();
				tmpSJendPosSetIter ++)
			{
				int tmpSJ_chrPosInt_end = (*tmpSJendPosSetIter);
				SJinfoTable_ofs << tmpChrName << "\t"
					<< tmpSJ_chrPosInt_start << "\t"
					<< tmpSJ_chrPosInt_end << "\t";
				bool annotatedOrNot_bool = returnSJannotatedOrNot(
					annotatedSJchrPosMapVec, tmpChr,
					tmpSJ_chrPosInt_start, tmpSJ_chrPosInt_end);
				if(annotatedOrNot_bool)
					SJinfoTable_ofs << "Y";
				else
					SJinfoTable_ofs << "N";
				for(int tmpSample = 0; tmpSample < sampleSize; tmpSample ++)
				{
					int tmpSJ_supportNum_inTmpSample 
						= searchAndReturnSJsupportNum(
							sampleSJchrPosSupNumMapVecVec[tmpSample],
							tmpChr, tmpSJ_chrPosInt_start, tmpSJ_chrPosInt_end);
					SJinfoTable_ofs << "\t" << tmpSJ_supportNum_inTmpSample;
				}
				SJinfoTable_ofs << endl;
			}
		}
	}	

	return 0;
}