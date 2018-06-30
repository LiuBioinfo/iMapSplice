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
#include "../general/SNPhash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputBeersSNPfile inputSJfilePath outputClosestSNP2SJfilePath" << endl;
		exit(1);
	}
	string inputBeersSNPfile = argv[2];
	string inputSJfilePath = argv[3];
	string outputFilePath = argv[4];

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
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	SNPhash_Info tmpSNPhashInfo;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_BeersSNPfile ..." << endl;
	tmpSNPhashInfo.generateSNPhash_BeersSNPfile(inputBeersSNPfile, indexInfo);
	cout << "start to do SNP search beside splice site " << endl;
	ifstream SJ_ifs(inputSJfilePath.c_str());
	ofstream closestSNP2SJ_ofs(outputFilePath.c_str());
	while(!SJ_ifs.eof())
	{
		string tmpSJstr;
		getline(SJ_ifs, tmpSJstr);
		if(tmpSJstr == "")
			break;
		//cout << "tmpSJstr: " << endl << tmpSJstr << endl;
		int firstTabLoc = tmpSJstr.find("\t");
		int secondTabLoc = tmpSJstr.find("\t", firstTabLoc + 1);
		int thirdTabLoc = tmpSJstr.find("\t", secondTabLoc + 1);
		string tmpSJ_chrName = tmpSJstr.substr(0, firstTabLoc);
		string tmpSJ_startPosStr = tmpSJstr.substr(firstTabLoc + 1, secondTabLoc - firstTabLoc - 1);
		string tmpSJ_endPosStr;
		if(thirdTabLoc != string::npos)
			tmpSJ_endPosStr = tmpSJstr.substr(secondTabLoc + 1, thirdTabLoc - secondTabLoc - 1);
		else
			tmpSJ_endPosStr = tmpSJstr.substr(secondTabLoc + 1);
		//cout << "tmpSJ_chrName: " << tmpSJ_chrName << endl;
		//cout << "tmpSJ_startPosStr: " << tmpSJ_startPosStr << endl;
		//cout << "tmpSJ_endPosStr: " << tmpSJ_endPosStr << endl;
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrName);
		int tmpSJ_startPos = atoi(tmpSJ_startPosStr.c_str());
		int tmpSJ_endPos = atoi(tmpSJ_endPosStr.c_str());
		//cout << "tmpSJ_chrNameInt: " << tmpSJ_chrNameInt << endl;
		//cout << "tmpSJ_startPos: " << tmpSJ_startPos << endl;
		//cout << "tmpSJ_endPos: " << tmpSJ_endPos << endl;
		int tmpRightMostSNPpos_SJdonerEnd = tmpSNPhashInfo.returnRightMostSNPpos(tmpSJ_chrNameInt,
			tmpSJ_startPos - 100, tmpSJ_startPos);
		//cout << "tmpRightMostSNPpos_SJdonerEnd: " << tmpRightMostSNPpos_SJdonerEnd << endl;
		int tmpLeftMostSNPpos_SJacceptorStart = tmpSNPhashInfo.returnLeftMostSNPpos(tmpSJ_chrNameInt,
			tmpSJ_endPos, tmpSJ_endPos + 100);
		//cout << "tmpLeftMostSNPpos_SJacceptorStart: " << tmpLeftMostSNPpos_SJacceptorStart << endl;
		int tmpClosestSNPdistance_left;
		if(tmpRightMostSNPpos_SJdonerEnd > 0)
			tmpClosestSNPdistance_left = tmpSJ_startPos - tmpRightMostSNPpos_SJdonerEnd;
		else
			tmpClosestSNPdistance_left = 9999;
		int tmpClosestSNPdistance_right;
		if(tmpLeftMostSNPpos_SJacceptorStart > 0)
			tmpClosestSNPdistance_right = tmpLeftMostSNPpos_SJacceptorStart - tmpSJ_endPos;
		else
			tmpClosestSNPdistance_right = 9999;
		int tmpClosestSNPdistance = 9999;
		if(tmpClosestSNPdistance_left < tmpClosestSNPdistance)
			tmpClosestSNPdistance = tmpClosestSNPdistance_left;
		if(tmpClosestSNPdistance_right < tmpClosestSNPdistance)
			tmpClosestSNPdistance = tmpClosestSNPdistance_right;
		closestSNP2SJ_ofs << tmpSJstr << "\t" << tmpClosestSNPdistance_left 
			<< "\t" << tmpClosestSNPdistance_right 
			<< "\t" << tmpClosestSNPdistance << endl;
	}
	cout << "all jobs done ..." << endl;
	SJ_ifs.close();
	closestSNP2SJ_ofs.close();
	return 0;
}