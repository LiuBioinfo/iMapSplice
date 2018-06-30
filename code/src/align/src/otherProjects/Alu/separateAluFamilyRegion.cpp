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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"

using namespace std;

void extractAluRegionInfoFromStr(string& tmpAluRegionStr,
	int& startPos, int& endPos, string& tmpAluFamilyName, string& tmpChrName)
{
	int tabLoc_1 = tmpAluRegionStr.find("\t");
	int tabLoc_2 = tmpAluRegionStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpAluRegionStr.find("\t", tabLoc_2 + 1);	
	int tabLoc_4 = tmpAluRegionStr.find("\t", tabLoc_3 + 1);
	//cout << "tabLoc_3: " << tabLoc_3 << " tabLoc_4: " << tabLoc_4 << endl;
	//cout << "string::npos: " << string::npos << endl;
	string startPosStr = tmpAluRegionStr.substr(0, tabLoc_1);
	startPos = atoi(startPosStr.c_str()) + 1;
	string endPosStr = tmpAluRegionStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	endPos = atoi(endPosStr.c_str());
	tmpAluFamilyName = tmpAluRegionStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	if(tabLoc_4 == string::npos)
	{
		//cout << "tabLoc_4 == string::npos " << endl;
		//cout << "length of tmpAluRegionStr: " << tmpAluRegionStr.length() << endl;
		//cout << "last base: " << tmpAluRegionStr.substr(tmpAluRegionStr.length()-1, 1) << endl;
		tmpChrName = tmpAluRegionStr.substr(tabLoc_3 + 1, tmpAluRegionStr.length() - 1 - tabLoc_3 - 1);
	}
	else
		tmpChrName = tmpAluRegionStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);	
}

void add2aluFamilyNameVecAndAluRegionVecVec(
	vector <string>& aluFamilyNameVec,
	vector < vector< pair< int, pair<int,int> > > >& aluRegionVecVec,
	int tmpChrNameInt, int startPos, int endPos, string& tmpAluFamilyName)
{
	int matchedIndexInVec = -1;
	for(int tmp = 0; tmp < aluFamilyNameVec.size(); tmp++)
	{
		string tmpNameInVec = aluFamilyNameVec[tmp];
		if(tmpNameInVec == tmpAluFamilyName)
		{
			matchedIndexInVec = tmp;
			break;
		}
	}
	if(matchedIndexInVec < 0)
	{
		aluFamilyNameVec.push_back(tmpAluFamilyName);
		vector< pair< int, pair<int,int> > > tmpRegionVec;
		tmpRegionVec.push_back(pair< int, pair<int,int> >(tmpChrNameInt, 
			pair<int,int>(startPos, endPos)));
		aluRegionVecVec.push_back(tmpRegionVec);
	}
	else
		aluRegionVecVec[matchedIndexInVec].push_back(pair< int, pair<int,int> >(
			tmpChrNameInt, pair<int,int>(startPos, endPos)));
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolderPath inputAluRegionInfoFile outputFolder" << endl;
		exit(1);
	}
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());	

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	//parameter_ifs.close();
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	parameter_ifs.close();
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;
	//exit(1);
	string inputAluRegionInfoFile = argv[2];
	vector <string> aluFamilyNameVec;
	vector < vector< pair< int, pair<int,int> > > > aluRegionVecVec;
	ifstream aluRegion_ifs(inputAluRegionInfoFile.c_str());
	while(!aluRegion_ifs.eof())
	{
		string tmpAluRegionStr;
		getline(aluRegion_ifs, tmpAluRegionStr);
		if(tmpAluRegionStr == "")
			break;
		string tmpChrName, tmpAluFamilyName;
		int tmpStartPos, tmpEndPos;
		extractAluRegionInfoFromStr(tmpAluRegionStr, tmpStartPos, 
			tmpEndPos, tmpAluFamilyName, tmpChrName);
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			cout << "chrName: " << tmpChrName << endl;
			cout << "length of chrNameStr: " << tmpChrName.length() << endl;
			cout << "invalid chrName:" << endl << tmpAluRegionStr << endl;
			log_ofs << "invalid chrName:" << endl << tmpAluRegionStr << endl;
		}
		else
			add2aluFamilyNameVecAndAluRegionVecVec(aluFamilyNameVec, aluRegionVecVec, 
				tmpChrNameInt, tmpStartPos, tmpEndPos, tmpAluFamilyName);
	}
	aluRegion_ifs.close();
	
	int aluFamilyNum = aluFamilyNameVec.size();
	cout << endl << "detected Alu family #: " << aluFamilyNum << endl << endl;
	log_ofs << endl << "detected Alu family #: " << aluFamilyNum << endl << endl;

	string AluFamilyRegionFile_Ygroup = outputFolderStr + "AluYgroup_region.txt";
	ofstream AluFamilyRegion_Ygroup_ofs(AluFamilyRegionFile_Ygroup.c_str());
	string AluFamilyRegionFile_merged = outputFolderStr + "AluMerged_region.txt";
	ofstream AluFamilyRegionFile_merged_ofs(AluFamilyRegionFile_merged.c_str());
	for(int tmp = 0; tmp < aluFamilyNum; tmp++)
	{
		string tmpAluFamilyName = aluFamilyNameVec[tmp];
		bool belong2aluYgroup_bool = (tmpAluFamilyName.substr(0,4) == "AluY");
		string tmpAluFamilyRegionFile = outputFolderStr + tmpAluFamilyName + "_region.txt";
		ofstream tmpAluFamilyRegion_ofs(tmpAluFamilyRegionFile.c_str());
		int tmpRegionNum = (aluRegionVecVec[tmp]).size();
		cout << "alu family [" << tmp+1 << "]:\t" << tmpAluFamilyName 
			<< "\tregion #:\t" << tmpRegionNum << endl;
		log_ofs << "alu family [" << tmp+1 << "]:\t" << tmpAluFamilyName 
			<< "\tregion #:\t" << tmpRegionNum << endl;
		for(int tmpRegion = 0; tmpRegion < tmpRegionNum; tmpRegion++)
		{
			int tmpRegion_chrNameInt = ((aluRegionVecVec[tmp])[tmpRegion]).first;
			string tmpChrName = indexInfo->returnChrNameStr(tmpRegion_chrNameInt);
			int tmpStartPos = (((aluRegionVecVec[tmp])[tmpRegion]).second).first;
			int tmpEndPos = (((aluRegionVecVec[tmp])[tmpRegion]).second).second;
			tmpAluFamilyRegion_ofs << tmpChrName << "\t" << tmpStartPos << "\t"
				<< tmpEndPos << "\t" << tmpEndPos - tmpStartPos + 1 << endl;
			AluFamilyRegionFile_merged_ofs << tmpChrName << "\t" << tmpStartPos << "\t"
				<< tmpEndPos << "\t" << tmpEndPos - tmpStartPos + 1 << "\t" 
				<< tmpAluFamilyName << endl; 
			if(belong2aluYgroup_bool)
			{
				AluFamilyRegion_Ygroup_ofs << tmpChrName << "\t" << tmpStartPos << "\t"
					<< tmpEndPos << "\t" << tmpEndPos - tmpStartPos + 1 << endl;				
			}
		}
		tmpAluFamilyRegion_ofs.close();
	}
	AluFamilyRegionFile_merged_ofs.close();
	AluFamilyRegion_Ygroup_ofs.close();
	delete indexInfo;	
	log_ofs.close();
	return 0;
}