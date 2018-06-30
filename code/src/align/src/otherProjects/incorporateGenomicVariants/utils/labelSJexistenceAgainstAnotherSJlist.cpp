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
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputJuncListCompared2 inputJuncList2compare outputLabeledJuncList" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "end of initiating indexInfo" << endl;

	string inputJuncListCompared2 = argv[2];
	AlignInferJunctionHash_Info* juncHash_compared2 = new AlignInferJunctionHash_Info();
	juncHash_compared2->initiateAlignInferJunctionInfo(chromNum);
	juncHash_compared2->insertJuncFromJuncFile_chrNamePosOnly(inputJuncListCompared2, indexInfo);

	string inputJuncList2compare = argv[3];
	string outputLabeledJuncList = argv[4];
	ifstream juncList2compare_ifs(inputJuncList2compare.c_str());
	ofstream labeledJuncList_ofs(outputLabeledJuncList.c_str());
	while(!juncList2compare_ifs.eof())
	{
		string tmpStr;
		getline(juncList2compare_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpJunc_chrNameStr = tmpStr.substr(0, tabLoc_1);
		string tmpJunc_chrStartPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpJunc_chrEndPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		string tmpJunc_otherStr = tmpStr.substr(tabLoc_3 + 1);
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrNameStr);
		int tmpJunc_startPos = atoi(tmpJunc_chrStartPosStr.c_str());
		int tmpJunc_endPos = atoi(tmpJunc_chrEndPosStr.c_str());
		int tmpJunc_indexInTheJuncHashCompared2 = juncHash_compared2->searchAndReturnAlignInferInfoVecIndex(
				tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
		if(tmpJunc_indexInTheJuncHashCompared2 > 0)
			labeledJuncList_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_chrStartPosStr << "\t"
				<< tmpJunc_chrEndPosStr << "\tY\t" << tmpJunc_otherStr << endl;
		else
			labeledJuncList_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_chrStartPosStr << "\t"
				<< tmpJunc_chrEndPosStr << "\tN\t" << tmpJunc_otherStr << endl;			
	}
	juncList2compare_ifs.close();
	labeledJuncList_ofs.close();
	delete juncHash_compared2;
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}	