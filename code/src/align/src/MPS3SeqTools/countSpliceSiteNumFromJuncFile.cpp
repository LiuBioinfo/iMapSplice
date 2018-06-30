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
		cout << "Executable <InputIndexInfo> <InputSJ>" << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	//string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	//char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	//cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	//indexInfo->readGenome(chrom);
	//cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	//cout << "start to load every chromosome" << endl;
	//indexInfo->initiate();	
	//cout << "start to initiate chrNameIndexArray" << endl;
	//indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();
	cout << "finish loading chromosomes" << endl;

	string inputJuncFile = argv[2];
	//string outputSpliceSiteNumPath = argv[3];


	ifstream junc_ifs(inputJuncFile.c_str());
	//ofstream spliceSiteNum_ofs(outputSpliceSiteNumPath.c_str());

	vector< set<int> > donerSiteSetVec;
	vector< set<int> > acceptorSiteSetVec;
	for(int tmp = 0; tmp < chromNum; tmp ++)
	{
		set<int> tmpDonerSiteSet;
		set<int> tmpAcceptorSiteSet;
		donerSiteSetVec.push_back(tmpDonerSiteSet);
		acceptorSiteSetVec.push_back(tmpAcceptorSiteSet);
	}
	cout << "start to read splice junction file!" << endl;
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
		donerSiteSetVec[tmpJunc_chrNameInt].insert(tmpJunc_startPos);
		acceptorSiteSetVec[tmpJunc_chrNameInt].insert(tmpJunc_endPos);
	}

	int donerSiteNum_total = 0;
	int acceptorSiteNum_total = 0;
	for(int tmp = 0; tmp < chromNum; tmp ++)
	{
		int tmpDonerSiteNum = donerSiteSetVec[tmp].size();
		int tmpAcceptorSiteNum = acceptorSiteSetVec[tmp].size();
		donerSiteNum_total += tmpDonerSiteNum;
		acceptorSiteNum_total += tmpAcceptorSiteNum;
	}

	cout << "total number of doner site: " << donerSiteNum_total << endl;
	cout << "total number of acceptor site: " << acceptorSiteNum_total << endl;

	delete indexInfo;
	//free(chrom);
	//chrom_bit_file_ifs.close();
	parameter_ifs.close();
	junc_ifs.close();
	//spliceSiteNum_ofs.close();
	return 0;
}	