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

void extractAluRegionInfoFromStr(string& tmpAluRegionStr, string& tmpChrName, 
	int& startPos, int& endPos, int& aluSeqLength, string& tmpAluFamilyName)
{
	int tabLoc_1 = tmpAluRegionStr.find("\t");
	int tabLoc_2 = tmpAluRegionStr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpAluRegionStr.find("\t", tabLoc_2 + 1);	
	int tabLoc_4 = tmpAluRegionStr.find("\t", tabLoc_3 + 1);
	tmpChrName = tmpAluRegionStr.substr(0, tabLoc_1);
	string startPosStr = tmpAluRegionStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	startPos = atoi(startPosStr.c_str());
	string endPosStr = tmpAluRegionStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	endPos = atoi(endPosStr.c_str());
	string aluSeqLengthStr = tmpAluRegionStr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	aluSeqLength = atoi(aluSeqLengthStr.c_str());
	tmpAluFamilyName = tmpAluRegionStr.substr(tabLoc_4 + 1);
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder inputMergedAluRegionFile outputAluSeqFastaFile" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string inputMergedAluRegionFile = argv[2];
	string outputAluSeqFastaFile = argv[3];
	ifstream aluRegion_ifs(inputMergedAluRegionFile.c_str());
	ofstream aluSeq_ofs(outputAluSeqFastaFile.c_str());
	while(!aluRegion_ifs.eof())
	{
		string tmpAluRegionStr;
		getline(aluRegion_ifs, tmpAluRegionStr);
		if(tmpAluRegionStr == "")
			break;
		string tmpChrName, tmpAluFamilyName;
		int tmpStartPos, tmpEndPos, tmpAluSeqLength;
		extractAluRegionInfoFromStr(tmpAluRegionStr, tmpChrName, 
			tmpStartPos, tmpEndPos, tmpAluSeqLength, tmpAluFamilyName);		
		int tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);
		if(tmpChrNameInt < 0)
		{
			cout << "invalid chrom, existing ......" << endl;
			exit(1);
		}
		string tmpAlu_name = tmpAluFamilyName + "_" + tmpChrName + "_" 
			+ int_to_str(tmpStartPos) + "_" + int_to_str(tmpEndPos) 
			+ "_" + int_to_str(tmpAluSeqLength);
		string tmpAlu_seq = indexInfo->returnChromStrSubstr(
			tmpChrNameInt, tmpStartPos, tmpAluSeqLength);
		aluSeq_ofs << ">" << tmpAlu_name << endl << tmpAlu_seq << endl;
	}
	delete indexInfo;
	parameter_ifs.close();
	aluRegion_ifs.close();
	aluSeq_ofs.close();
	return 0;
}