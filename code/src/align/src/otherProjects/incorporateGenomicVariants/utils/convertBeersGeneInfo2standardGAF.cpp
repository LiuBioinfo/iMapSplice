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

using namespace std;

bool sortExonPosAsceOrderBool(int i, int j) { return (i < j);}
bool sortExonPosDescOrderBool(int i, int j) { return (i > j);}

void updatedBeersGeneInfoExonStr2standardGAF(string& tmpStrandStr, int tmpExonNum,
	string& tmpExonPosVecStr, string& tmpExonPosVecStr_updated)
{
	bool positiveStrandBool = (tmpStrandStr == "+");
	vector<string> exonPosStrVec;
	if(tmpExonNum == 1)
		exonPosStrVec.push_back(tmpExonPosVecStr);
	else
	{	
		int startLoc = 0;
		for(int tmp = 0; tmp < tmpExonNum-1; tmp++)
		{
			int tabLoc = tmpExonPosVecStr.find(",", startLoc);
			string tmpExonPosStr = tmpExonPosVecStr.substr(startLoc, tabLoc - startLoc);
			exonPosStrVec.push_back(tmpExonPosStr);
			startLoc = tabLoc + 1;
		}
		exonPosStrVec.push_back(tmpExonPosVecStr.substr(startLoc));
	}
	vector<int> exonPosVec;	
	for(int tmp = 0; tmp < exonPosStrVec.size(); tmp++)
	{
		string tmpExonPosStr = exonPosStrVec[tmp];
		int tmpExonPos = atoi(tmpExonPosStr.c_str());
		exonPosVec.push_back(tmpExonPos);
	}
	if(positiveStrandBool)
		sort(exonPosVec.begin(), exonPosVec.end(), sortExonPosAsceOrderBool);
	else
		sort(exonPosVec.begin(), exonPosVec.end(), sortExonPosDescOrderBool);
	tmpExonPosVecStr_updated = "";
	for(int tmp = 0; tmp < exonPosVec.size(); tmp++)
	{
		int tmpExonPos = exonPosVec[tmp];
		string tmpExonPosStr = int_to_str(tmpExonPos);
		tmpExonPosVecStr_updated += tmpExonPosStr;
		tmpExonPosVecStr_updated += ",";
	}
}


int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexPath inputBeersGeneInfoFile outputStandardGAF" << endl;
		exit(1);
	}

	cout << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string inputBeersGeneInfoFile = argv[2];
	string outputStandardGAF = argv[3];
	string outputStandardGAF_invalidChr = outputStandardGAF + ".invalidChr";
	ifstream Beers_ifs(inputBeersGeneInfoFile.c_str()); 
	ofstream gaf_ofs(outputStandardGAF.c_str());
	ofstream gaf_invalid_ofs(outputStandardGAF_invalidChr.c_str());
	while(!Beers_ifs.eof())
	{
		string tmpStr;
		getline(Beers_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 7; tmp++)
		{
			int tabLoc = tmpStr.find('\t', startLoc);
			if(tabLoc == string::npos)
				return false;
			string tmpFieldStr = tmpStr.substr(startLoc, tabLoc - startLoc);
			tmpFieldVec.push_back(tmpFieldStr);
			startLoc = tabLoc + 1;
		}
		int nextTabLoc = tmpStr.find("\t", startLoc);
		if(nextTabLoc == string::npos)
			tmpFieldVec.push_back(tmpStr.substr(startLoc));
		else
		{
			cout << "not the last field ......" << endl;
			exit(1);
		}
		string tmpChrNameStr = tmpFieldVec[0];
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
		if(tmpChrNameInt < 0)
		{
			cout << "tmpChrNameStr: " << tmpChrNameStr << endl;
			cout << "tmpChrNameInt: " << tmpChrNameInt << endl;
			gaf_invalid_ofs << tmpStr << endl;
			continue;
		}
		string tmpStrandStr = tmpFieldVec[1];
		string tmpTranscriptLeftPosStr = tmpFieldVec[2];
		string tmpTranscriptRightPosStr = tmpFieldVec[3];
		string tmpExonNumStr = tmpFieldVec[4];
		int tmpExonNum = atoi(tmpExonNumStr.c_str());
		string tmpExonStartPosStr = tmpFieldVec[5];
		string tmpExonEndPosStr = tmpFieldVec[6];
		string tmpTranscriptNameStr = tmpFieldVec[7];

		// generate tmpExonStartPosStr_updated, tmpExonEndPosStr_updated
		string tmpExonStartPosStr_updated, tmpExonEndPosStr_updated;
		updatedBeersGeneInfoExonStr2standardGAF(tmpStrandStr, tmpExonNum,
			tmpExonStartPosStr, tmpExonStartPosStr_updated);
		updatedBeersGeneInfoExonStr2standardGAF(tmpStrandStr, tmpExonNum,
			tmpExonEndPosStr, tmpExonEndPosStr_updated);		
		// output
		if(tmpStrandStr == "+")
		{
			gaf_ofs << tmpTranscriptNameStr << "\t" << tmpChrNameStr << "\t+\t"
				<< tmpTranscriptLeftPosStr << "\t" << tmpTranscriptRightPosStr << "\t"
				<< tmpTranscriptLeftPosStr << "\t" << tmpTranscriptRightPosStr << "\t"
				<< tmpExonNumStr << "\t" << tmpExonStartPosStr_updated << "\t"
				<< tmpExonEndPosStr_updated << "\t" << tmpTranscriptNameStr << "\t" << tmpTranscriptNameStr << endl;
		}
		else if(tmpStrandStr == "-")
		{
			gaf_ofs << tmpTranscriptNameStr << "\t" << tmpChrNameStr << "\t-\t"
				<< tmpTranscriptRightPosStr << "\t" << tmpTranscriptLeftPosStr << "\t"
				<< tmpTranscriptRightPosStr << "\t" << tmpTranscriptLeftPosStr << "\t"
				<< tmpExonNumStr << "\t" << tmpExonStartPosStr_updated << "\t"
				<< tmpExonEndPosStr_updated << "\t" << tmpTranscriptNameStr << "\t" << tmpTranscriptNameStr << endl;
		}
		else
		{
			cout << "error in tmpStrandStr ..." << endl;
			exit(1);
		}
	}
	gaf_ofs.close();
	Beers_ifs.close();
	return 0;
}