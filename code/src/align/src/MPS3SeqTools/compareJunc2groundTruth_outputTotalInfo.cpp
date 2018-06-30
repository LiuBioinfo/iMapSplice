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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/alignInferJunctionHash_info.h"

using namespace std;

void extractChrNamePosFromSJstr(int& tmpSJ_chrNameInt, 
	int& tmpSJ_donerEndPos, int& tmpSJ_acceptorStartPos, 
	string& tmpRegularSJstr, Index_Info* indexInfo)
{
	vector<string> tmpSJfieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = tmpRegularSJstr.find("\t",startLoc);
		string tmpSJfield = tmpRegularSJstr.substr(startLoc, tabLoc - startLoc);
		tmpSJfieldVec.push_back(tmpSJfield);
		startLoc = tabLoc + 1;
	}
	string tmpSJ_chrNameStr = tmpSJfieldVec[0];
	string tmpSJ_donerEndPosStr = tmpSJfieldVec[1];
	string tmpSJ_acceptorStartPosStr = tmpSJfieldVec[2];
	tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrNameStr);
	tmpSJ_donerEndPos = atoi(tmpSJ_donerEndPosStr.c_str());
	tmpSJ_acceptorStartPos = atoi(tmpSJ_acceptorStartPosStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolder <BeersGroundTruthSJ> <inputRegularSJ> <offset> <outputFilePrefix>" << endl;
		exit(1);
	}

	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	string offsetStr = argv[4];
	int offset = atoi(offsetStr.c_str());

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[5];
	string outputDirStr = outputFolderStr + "/";
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	

	string outputFilePrefix = outputDirStr;
	string outputFile_compared2groundTruth = outputFilePrefix + "SJ_compared2groundTruth.txt";
	string outputFile_inGroundTruth = outputFilePrefix + "SJ_compared2groundTruth_correct.txt";
	string outputFile_outOfGroundTruth = outputFilePrefix + "SJ_compared2groundTruth_incorrect.txt";
	string outputFile_log = outputFilePrefix + "log.txt";

	ofstream compared2groundTruth_ofs(outputFile_compared2groundTruth.c_str());
	ofstream correct_ofs(outputFile_inGroundTruth.c_str());
	ofstream incorrect_ofs(outputFile_outOfGroundTruth.c_str());
	ofstream log_ofs(outputFile_log.c_str());

	string inputBeersGroundTruthSJpath = argv[2];
	AlignInferJunctionHash_Info* juncHash_groundTruth = new AlignInferJunctionHash_Info();
	juncHash_groundTruth->initiateAlignInferJunctionInfo(chromNum);
	juncHash_groundTruth->insertJuncFromJuncFile_chrNamePosOnly(inputBeersGroundTruthSJpath, indexInfo);

	int SJnum_correct = 0;
	int SJnum_incorrect = 0;

	string inputRegularSJ_path = argv[3];
	ifstream regularSJ_ifs(inputRegularSJ_path.c_str());
	while(!regularSJ_ifs.eof())
	{
		string tmpRegularSJstr;
		getline(regularSJ_ifs, tmpRegularSJstr);
		if(tmpRegularSJstr == "")
			break;
		int tmpSJ_chrNameInt, tmpSJ_donerEndPos, tmpSJ_acceptorStartPos;
		extractChrNamePosFromSJstr(tmpSJ_chrNameInt, tmpSJ_donerEndPos, 
			tmpSJ_acceptorStartPos, tmpRegularSJstr, indexInfo);
		
		bool foundInGroundTruthJuncHash_bool = false;
		for(int tmp = 0 - offset; tmp <= offset; tmp++)
		{
			for(int tmp2 = 0 - offset; tmp2 <= offset; tmp2++)
			{
				int tmpSJ_donerEndPos_withOffset = tmpSJ_donerEndPos + tmp;
				int tmpSJ_acceptorStartPos_withOffset = tmpSJ_acceptorStartPos + tmp2;
				int tmpSJ_index_inGroundTruthJuncHash
					= juncHash_groundTruth->searchAndReturnAlignInferInfoVecIndex(
						tmpSJ_chrNameInt, tmpSJ_donerEndPos_withOffset, tmpSJ_acceptorStartPos_withOffset);
				if(tmpSJ_index_inGroundTruthJuncHash >= 0)
				{
					foundInGroundTruthJuncHash_bool = true;
					break;
				}				
			}
			if(foundInGroundTruthJuncHash_bool)
				break;
		}

		if(foundInGroundTruthJuncHash_bool)
		{	
			SJnum_correct ++;
			correct_ofs << tmpRegularSJstr << endl;
			compared2groundTruth_ofs << tmpRegularSJstr << "\tTRUE" << endl;
		}
		else
		{	
			SJnum_incorrect ++;
			incorrect_ofs << tmpRegularSJstr << endl;
			compared2groundTruth_ofs << tmpRegularSJstr << "\tFALSE" << endl;
		}
	}
	int SJnum_total = SJnum_correct + SJnum_incorrect;
	log_ofs << "SJnum_total: " << SJnum_total << endl << endl;
	log_ofs << "SJnum_correct: " << SJnum_correct << endl;
	log_ofs << "SJnum_incorrect: " << SJnum_incorrect << endl;
	double SJ_specificity = ((double)SJnum_correct / (double) SJnum_total) * 100;
	log_ofs << endl << "SJ_specificity: " << SJ_specificity << "%" << endl << endl;
	
	parameter_ifs.close();
	compared2groundTruth_ofs.close();
	correct_ofs.close();
	incorrect_ofs.close();
	log_ofs.close();
	delete juncHash_groundTruth;
	delete indexInfo;
	return 0;
}