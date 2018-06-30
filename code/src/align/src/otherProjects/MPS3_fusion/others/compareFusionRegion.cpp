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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

using namespace std;

typedef set< pair<int,int> > FusionRegionPairSetType;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable fusionRegionFile_1 fusionRegionFile_2 outputFolder" << endl;
		exit(1);
	}

	string outputFolderStr = argv[3];
	cout << "creat folders and files ......" << endl;	
	string mkdirOutputCommand = "mkdir -p " + outputFolderStr;
	system(mkdirOutputCommand.c_str());
	string outputFilePrefix = outputFolderStr + "/";
	string output_fusionRegion_inBothFile = outputFilePrefix + "inBoth.genomicRegionIndex";
	ofstream fusionRegion_inBothFile_ofs(output_fusionRegion_inBothFile.c_str());
	string output_fusionRegion_inFile1_notInFile2 = outputFilePrefix + "inFile1_notInFile2.genomicRegionIndex";
	ofstream fusionRegion_inFile1_notInFile2_ofs(output_fusionRegion_inFile1_notInFile2.c_str());
	string output_fusionRegion_inFile2_notInFile1 = outputFilePrefix + "inFile2_notInFile1.genomicRegionIndex";
	ofstream fusionRegion_inFile2_notInFile1_ofs(output_fusionRegion_inFile2_notInFile1.c_str());

	string fusionRegionFile_1 = argv[1];
	string fusionRegionFile_2 = argv[2];

	ifstream fusionRegionFile_ifs_1(fusionRegionFile_1.c_str());
	ifstream fusionRegionFile_ifs_2(fusionRegionFile_2.c_str());

	cout << "generate fusionRegionPairSet_1 ......" << endl;	
	FusionRegionPairSetType fusionRegionPairSet_1;
	while(!fusionRegionFile_ifs_1.eof())
	{
		string tmpFusionRegionStr;
		getline(fusionRegionFile_ifs_1, tmpFusionRegionStr);
		if(tmpFusionRegionStr == "")
			break;
		vector<string> fusionRegionStrFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = tmpFusionRegionStr.find("\t", startLoc);
			fusionRegionStrFieldVec.push_back(tmpFusionRegionStr.substr(startLoc, tabLoc-startLoc));
			startLoc = tabLoc + 1;
		}
		int genomicRegion_index_1 = atoi(fusionRegionStrFieldVec[0].c_str());
		int genomicRegion_index_2 = atoi(fusionRegionStrFieldVec[1].c_str());
		fusionRegionPairSet_1.insert(pair<int,int>(genomicRegion_index_1, genomicRegion_index_2));
	}
	fusionRegionFile_ifs_1.close();
	cout << "fusionRegionPairSet_1.size(): " << fusionRegionPairSet_1.size() << endl;

	cout << "generate fusionRegionPairSet_2 ......" << endl;
	FusionRegionPairSetType fusionRegionPairSet_2;
	while(!fusionRegionFile_ifs_2.eof())
	{
		string tmpFusionRegionStr;
		getline(fusionRegionFile_ifs_2, tmpFusionRegionStr);
		if(tmpFusionRegionStr == "")
			break;
		vector<string> fusionRegionStrFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = tmpFusionRegionStr.find("\t", startLoc);
			fusionRegionStrFieldVec.push_back(tmpFusionRegionStr.substr(startLoc, tabLoc-startLoc));
			startLoc = tabLoc + 1;
		}
		int genomicRegion_index_1 = atoi(fusionRegionStrFieldVec[0].c_str());
		int genomicRegion_index_2 = atoi(fusionRegionStrFieldVec[1].c_str());
		fusionRegionPairSet_2.insert(pair<int,int>(genomicRegion_index_1, genomicRegion_index_2));		
	}
	fusionRegionFile_ifs_2.close();
	cout << "fusionRegionPairSet_2.size(): " << fusionRegionPairSet_2.size() << endl;

	ifstream fusionRegionFile_ifs_1_reCheck(fusionRegionFile_1.c_str());
	ifstream fusionRegionFile_ifs_2_reCheck(fusionRegionFile_2.c_str());
	while(!fusionRegionFile_ifs_1_reCheck.eof())
	{
		string tmpFusionRegionStr;
		getline(fusionRegionFile_ifs_1_reCheck, tmpFusionRegionStr);
		if(tmpFusionRegionStr == "")
			break;
		vector<string> fusionRegionStrFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = tmpFusionRegionStr.find("\t", startLoc);
			fusionRegionStrFieldVec.push_back(tmpFusionRegionStr.substr(startLoc, tabLoc-startLoc));
			startLoc = tabLoc + 1;
		}
		int genomicRegion_index_1 = atoi(fusionRegionStrFieldVec[0].c_str());
		int genomicRegion_index_2 = atoi(fusionRegionStrFieldVec[1].c_str());
		FusionRegionPairSetType::iterator tmpIter_2 = fusionRegionPairSet_2.find(pair<int,int>(
			genomicRegion_index_1, genomicRegion_index_2));
		if(tmpIter_2 == fusionRegionPairSet_2.end())
		{
			//cout << "notFound: " << genomicRegion_index_1 << "\t" << genomicRegion_index_2 << endl;
			fusionRegion_inFile1_notInFile2_ofs << tmpFusionRegionStr << endl;
		}
		else
		{
			//cout << "found: " << genomicRegion_index_1 << "\t" << genomicRegion_index_2 << endl;
			fusionRegion_inBothFile_ofs << tmpFusionRegionStr << endl;
		}
	}
	fusionRegionFile_ifs_1_reCheck.close();

	while(!fusionRegionFile_ifs_2_reCheck.eof())
	{
		string tmpFusionRegionStr;
		getline(fusionRegionFile_ifs_2_reCheck, tmpFusionRegionStr);
		if(tmpFusionRegionStr == "")
			break;
		vector<string> fusionRegionStrFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 2; tmp++)
		{
			int tabLoc = tmpFusionRegionStr.find("\t", startLoc);
			fusionRegionStrFieldVec.push_back(tmpFusionRegionStr.substr(startLoc, tabLoc-startLoc));
			startLoc = tabLoc + 1;
		}
		int genomicRegion_index_1 = atoi(fusionRegionStrFieldVec[0].c_str());
		int genomicRegion_index_2 = atoi(fusionRegionStrFieldVec[1].c_str());
		FusionRegionPairSetType::iterator tmpIter_1 = fusionRegionPairSet_1.find(pair<int,int>(
			genomicRegion_index_1, genomicRegion_index_2));
		if(tmpIter_1 == fusionRegionPairSet_1.end())
		{
			fusionRegion_inFile2_notInFile1_ofs << tmpFusionRegionStr << endl;
		}
		else
		{
			//fusionRegion_inBothFile_ofs << tmpFusionRegionStr << endl;
		}
	}
	fusionRegionFile_ifs_2_reCheck.close();



	return 0;
}