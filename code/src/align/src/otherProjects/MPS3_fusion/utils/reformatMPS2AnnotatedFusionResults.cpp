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

#include "../../../general/index_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable inputMPS2AnnotatedFusionResults outputReformatedFusionfile" << endl;
		exit(1);
	}	

	string inputFilePath = argv[1];
	string outputFilePath = argv[2];
	ifstream MPS2fusion_ifs(inputFilePath.c_str());
	ofstream reformatedFusion_ofs(outputFilePath.c_str());
	while(!MPS2fusion_ifs.eof())
	{
		string tmpStr;
		getline(MPS2fusion_ifs, tmpStr);
		if((MPS2fusion_ifs.eof())||(tmpStr == ""))
			break;		
		vector<string> fusionFieldVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			if(tabLoc != string::npos)
			{	
				string tmpFieldStr = tmpStr.substr(startLoc, tabLoc-startLoc);
				fusionFieldVec.push_back(tmpFieldStr);
				startLoc = tabLoc + 1;
			}
			else
				break;
		}
		fusionFieldVec.push_back(tmpStr.substr(startLoc));
		string tmpChrInfoStr = fusionFieldVec[0];
		string tmpStartPosStr = fusionFieldVec[1];
		string tmpEndPosStr = fusionFieldVec[2];
		string tmpJuncIdStr = fusionFieldVec[3];
		string tmpSupNumStr = fusionFieldVec[4];
		string tmpStrandStr = fusionFieldVec[5];
		string tmpFlankStringStr = fusionFieldVec[12];
		string tmpGeneNameStr_1 = fusionFieldVec[fusionFieldVec.size()-3];
		string tmpGeneNameStr_2 = fusionFieldVec[fusionFieldVec.size()-2];
		int tlideLoc = tmpChrInfoStr.find("~");
		string tmpChrName_1 = tmpChrInfoStr.substr(0, tlideLoc);
		string tmpChrName_2 = tmpChrInfoStr.substr(tlideLoc+1);
		string tmpStrand_1 = tmpStrandStr.substr(0,1);
		string tmpStrand_2 = tmpStrandStr.substr(1,1);
		reformatedFusion_ofs << tmpChrName_1 << "\t" << tmpChrName_2 << "\t"
			<< tmpStartPosStr << "\t" << tmpEndPosStr << "\t"
			<< tmpStrand_1 << "\t" << tmpStrand_2 << "\t" << tmpFlankStringStr << "\t"
			<< tmpSupNumStr << "\t" << tmpGeneNameStr_1 << "\t" << tmpGeneNameStr_2 << endl;
	}
	MPS2fusion_ifs.close();
	reformatedFusion_ofs.close();
	return 0;
}