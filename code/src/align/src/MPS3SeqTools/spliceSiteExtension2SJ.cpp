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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"
#include "../general/alignmentToJunc_supportNum.h"
#include "../general/fixSingleAnchorNWDP_info.h"

using namespace std;

void outputSpliceSiteExtensionSJ(string& tmpJuncStr, 
	ofstream& outputSJ_lowPenalty_ofs,
	ofstream& outputSJ_highPenalty_ofs)
{
	vector<string> juncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 10; tmp++)
	{
		int tabLoc = tmpJuncStr.find("\t", startLoc);
		string tmpJuncField = tmpJuncStr.substr(startLoc, tabLoc-startLoc);
		juncFieldVec.push_back(tmpJuncField);
		startLoc = tabLoc + 1;
	}
	string chrNameStr = juncFieldVec[0];
	string donerEndPosStr = juncFieldVec[1];
	string acceptorStartPosStr = juncFieldVec[2];
	string juncSupNumStr = juncFieldVec[3];
	string donerAnchorLengthStr = juncFieldVec[4];
	string acceptorAnchorLengthStr = juncFieldVec[5];
	string flankStringCaseStr = juncFieldVec[6];
	string donerAnchorNWDPcigarSting = juncFieldVec[7];
	string acceptorAnchorNWDPcigarSting = juncFieldVec[8];
	string donerAnchorNWDPpenaltyStr = juncFieldVec[9];
	string acceptorAnchorNWDPpenaltyStr = tmpJuncStr.substr(startLoc);

	int donerAnchorLength = atoi(donerAnchorLengthStr.c_str());
	int acceptorAnchorLength = atoi(acceptorAnchorLengthStr.c_str());

	int donerAnchorNWDPpenalty = atoi(donerAnchorNWDPpenaltyStr.c_str());
	int acceptorAnchorNWDPpenalty = atoi(acceptorAnchorNWDPpenaltyStr.c_str());

	bool donerMatchBool = false;
	bool acceptorMatchBool = false;

	if((1 + donerAnchorLength/8) >= donerAnchorNWDPpenalty)
		donerMatchBool = true;
	if((1 + acceptorAnchorLength/8) >= acceptorAnchorNWDPpenalty)
		acceptorMatchBool = true;

	string junctionStr;
	junctionStr = chrNameStr + "\t" + donerEndPosStr + "\t"
		+ acceptorStartPosStr + "\tJUNC\t" + juncSupNumStr + "\t*\t*\t*\t*\t*\t"
		+ donerAnchorLengthStr + "," + acceptorAnchorLengthStr + ",\t*\t*\t" + flankStringCaseStr
		+ "\t*\t*\t*\t*\t*\t*";  
	if(donerMatchBool || acceptorMatchBool)
		outputSJ_lowPenalty_ofs << junctionStr << endl;
	else
		outputSJ_highPenalty_ofs << junctionStr << endl;
}



int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Executable <inputSJfile> <outputSJfolder>" << endl;
		exit(0);
	}
	string inputSJpath = argv[1];
	string outputFolder = argv[2];
	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";
	string outputSJ_lowPenalty = outputFolder + "lowPenaltySJ.junc";
	string outputSJ_highPenalty = outputFolder + "highPenaltySJ.junc";

	ifstream inputSJ_ifs(inputSJpath.c_str());
	ofstream outputSJ_lowPenalty_ofs(outputSJ_lowPenalty.c_str());
	ofstream outputSJ_highPenalty_ofs(outputSJ_highPenalty.c_str());	
	while(1)
	{
		if(inputSJ_ifs.eof())
			break;
		string tmpJunctionStr;
		getline(inputSJ_ifs, tmpJunctionStr);
		if(inputSJ_ifs.eof())
			break;
		if(tmpJunctionStr.substr(0,3) != "chr")
			continue;
		outputSpliceSiteExtensionSJ(tmpJunctionStr, 
			outputSJ_lowPenalty_ofs, outputSJ_highPenalty_ofs);
	}
	inputSJ_ifs.close();
	outputSJ_lowPenalty_ofs.close();
	outputSJ_highPenalty_ofs.close();
	return 0;
}