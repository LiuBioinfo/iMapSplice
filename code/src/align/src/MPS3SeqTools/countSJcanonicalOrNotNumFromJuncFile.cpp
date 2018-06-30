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

using namespace std;

int main(int argc, char**argv)
{
	if(argc != 3)
	{
		cout << "Executable <inputIndexFolderPath> <InputSJ>" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	//cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	//cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	//cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	//cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	cout << "finish loading chromosomes" << endl;

	string inputJuncFile = argv[2];
	//string outputSpliceSiteNumPath = argv[3];
	ifstream junc_ifs(inputJuncFile.c_str());

	
	int SJcanonicalOrNotTypeNumTotal = 7;
	vector<int> SJcanonicalOrNotNumVec;
	for(int tmp = 0; tmp < SJcanonicalOrNotTypeNumTotal; tmp++)
		SJcanonicalOrNotNumVec.push_back(0);
	//cout << "start to read splice junction file!" << endl;
	int totalJuncNum = 0;
	while(!junc_ifs.eof())
	{
		string tmpJuncStr;
		getline(junc_ifs, tmpJuncStr);
		if(tmpJuncStr == "")
			break;
		int tabLoc_1 = tmpJuncStr.find("\t");
		int tabLoc_2 = tmpJuncStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpJuncStr.find("\t", tabLoc_2 + 1);
		string tmpJunc_chrName = tmpJuncStr.substr(0,tabLoc_1);
		string tmpJunc_startPosStr = tmpJuncStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpJunc_endPosStr = tmpJuncStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
		int tmpJunc_chrNameInt = indexInfo->convertStringToInt(tmpJunc_chrName);
		int tmpJunc_startPos = atoi(tmpJunc_startPosStr.c_str());
		int tmpJunc_endPos = atoi(tmpJunc_endPosStr.c_str());
		string tmpJunc_flankString = indexInfo->returnFlankString(tmpJunc_chrNameInt, 
			tmpJunc_startPos, tmpJunc_endPos);
		int tmpJunc_flankStringCase = indexInfo->returnFlankStringCaseFromFlankString(tmpJunc_flankString);
		SJcanonicalOrNotNumVec[tmpJunc_flankStringCase] ++;
		totalJuncNum ++;
	}
	cout << "***************************" << endl << endl;
	int canonicalSJnum = SJcanonicalOrNotNumVec[5] + SJcanonicalOrNotNumVec[6];
	double canonicalSJnum_perc = ((double)canonicalSJnum / (double)totalJuncNum) * 100;
	int semiCanonicalSJnum = SJcanonicalOrNotNumVec[1] + SJcanonicalOrNotNumVec[2]
		+ SJcanonicalOrNotNumVec[3] + SJcanonicalOrNotNumVec[4];
	double semiCanonicalSJnum_perc = ((double)semiCanonicalSJnum / (double)totalJuncNum) * 100;
	int nonCanonicalSJnum = SJcanonicalOrNotNumVec[0];
	double nonCanonicalSJnum_perc = ((double)nonCanonicalSJnum / (double)totalJuncNum) * 100;
	cout << "total SJ #: " << totalJuncNum << endl;
	cout << "canonical SJ #: " << canonicalSJnum << " -- " << canonicalSJnum_perc << "%" << endl;
	cout << "semi-canonical SJ #: " << semiCanonicalSJnum << " -- " << semiCanonicalSJnum_perc << "%" << endl;
	cout << "non-canonical SJ #: " << nonCanonicalSJnum << " -- " << nonCanonicalSJnum_perc << "%" << endl;
	cout << endl << "***************************" << endl;
	cout << "GTAG: " << SJcanonicalOrNotNumVec[5] << " -- " << ((double)SJcanonicalOrNotNumVec[5]/(double)totalJuncNum)*100 << "%" << endl;
	cout << "CTAC: " << SJcanonicalOrNotNumVec[6] << " -- " << ((double)SJcanonicalOrNotNumVec[6]/(double)totalJuncNum)*100 << "%" << endl;
	cout << endl;
	cout << "ATAC: " << SJcanonicalOrNotNumVec[1] << " -- " << ((double)SJcanonicalOrNotNumVec[1]/(double)totalJuncNum)*100 << "%" << endl;
	cout << "GTAT: " << SJcanonicalOrNotNumVec[2] << " -- " << ((double)SJcanonicalOrNotNumVec[2]/(double)totalJuncNum)*100 << "%" << endl;
	cout << "CTGC: " << SJcanonicalOrNotNumVec[3] << " -- " << ((double)SJcanonicalOrNotNumVec[3]/(double)totalJuncNum)*100 << "%" << endl;
	cout << "GCAG: " << SJcanonicalOrNotNumVec[4] << " -- " << ((double)SJcanonicalOrNotNumVec[4]/(double)totalJuncNum)*100 << "%" << endl;
	cout << endl;
	cout << "others: " << SJcanonicalOrNotNumVec[0] << " -- " << ((double)SJcanonicalOrNotNumVec[0]/(double)totalJuncNum)*100 << "%" << endl;
	cout << endl;
	delete indexInfo; 
	free(chrom);
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	junc_ifs.close();
	//spliceSiteNum_ofs.close();
	return 0;
}	