// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "../general/read_block_test.h"
#include "../general/index_info.h"

using namespace std;

void extractIntraChromJuncInfoFromStr(
	int& tmpChrNameInt, int& tmpStartPos, int& tmpEndPos,
	string& tmpIntraChromSJstr, Index_Info* indexInfo)
{
	vector<string> tmpIntraChromSJfieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = tmpIntraChromSJstr.find("\t", startLoc);
		string tmpIntraChromSJfield = tmpIntraChromSJstr.substr(startLoc, tabLoc - startLoc);
		tmpIntraChromSJfieldVec.push_back(tmpIntraChromSJfield);
		startLoc = tabLoc + 1;
	}
	string tmpChrNameStr = tmpIntraChromSJfieldVec[0];
	string tmpIntraChromJuncStartPosStr = tmpIntraChromSJfieldVec[1];
	string tmpIntraChromJuncEndPosStr = tmpIntraChromSJfieldVec[2];
	tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
	tmpStartPos = atoi(tmpIntraChromJuncStartPosStr.c_str());
	tmpEndPos = atoi(tmpIntraChromJuncEndPosStr.c_str());
}

bool extractExonChrNamePosFromGTFstr(
	int& tmpExonChrNameInt, int& tmpExonChrStartPos, int& tmpExonChrEndPos,
	string& tmpExonStrand, string& tmpGeneName,
	string& tmpExonStr, Index_Info* indexInfo)
{
	vector<string> tmpExonFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 7; tmp++)
	{
		int tabLoc = tmpExonStr.find("\t", startLoc);
		if(tabLoc == string::npos)
			return false;
		tmpExonFieldVec.push_back(tmpExonStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}
	string tmpExon_chrNameStr = tmpExonFieldVec[0];
	string tmpExon_sourceStr = tmpExonFieldVec[1];
	string tmpExon_featureStr = tmpExonFieldVec[2];
	string tmpExon_startPosStr = tmpExonFieldVec[3];
	string tmpExon_endPosStr = tmpExonFieldVec[4];
	
	tmpExonStrand = tmpExonFieldVec[6];

	int tmpGeneNameLoc = tmpExonStr.find("gene_name");
	if(tmpGeneNameLoc == string::npos)
		return false;
	int tmpGeneNameLoc_nextSemiColon = tmpExonStr.find(";", tmpGeneNameLoc);
	if(tmpGeneNameLoc_nextSemiColon == string::npos)
		return false;
	tmpGeneName = tmpExonStr.substr(tmpGeneNameLoc+11, tmpGeneNameLoc_nextSemiColon - 2 - tmpGeneNameLoc - 11 + 1);

	//if(tmpExon_featureStr != "exon")
	//	return false;
	if(tmpExon_chrNameStr.substr(0,3) == "chr")
	{}
	else
	{
		string tmpExon_chrNameStrOri = tmpExon_chrNameStr;
		tmpExon_chrNameStr = "chr" + tmpExon_chrNameStrOri;
	}
	tmpExonChrNameInt = indexInfo->convertStringToInt(tmpExon_chrNameStr);
	tmpExonChrStartPos = atoi(tmpExon_startPosStr.c_str());
	tmpExonChrEndPos = atoi(tmpExon_endPosStr.c_str());
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexInfoFolderPath inputGTFannotation inputIntraChromJuncPath";
		cout << " outputFusionJuncWithGene offset" << endl;
		exit(1);
	}
	string offsetStr = argv[5];
	int offsetInt = atoi(offsetStr.c_str());
	string inputGTFannotationPath = argv[2];
	string inputIntraChromJuncPath = argv[3];
	//string outputIntraChromJuncWithGeneInfoPath = argv[4];

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	cout << "end of initiating indexInfo" << endl;

	//int tmpFusionJuncNum = 0;
	ifstream intraChromJunc_ifs(inputIntraChromJuncPath.c_str());

	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	string outputIntraChromJuncWithGeneInfoPath = outputFolderStr + "junc_geneInfo.txt";
	string outputIntraChromJuncWithGeneInfoPath_intraGene
		= outputFolderStr + "junc_geneInfo_intraGene.txt";
	string outputInterChromJuncWithGeneInfoPath_interGene
		= outputFolderStr + "junc_geneInfo_interGene.txt";

	ofstream intraChromJuncWithGeneInfo_ofs(outputIntraChromJuncWithGeneInfoPath.c_str());
	ofstream intraChromJuncWithGeneInfo_intraGene_ofs(
		outputIntraChromJuncWithGeneInfoPath_intraGene.c_str());
	ofstream intraChromJuncWithGeneInfo_interGene_ofs(
		outputInterChromJuncWithGeneInfoPath_interGene.c_str());
	while(!intraChromJunc_ifs.eof())
	{
		string tmpIntraChromJuncStr;
		getline(intraChromJunc_ifs, tmpIntraChromJuncStr);
		if((intraChromJunc_ifs.eof())||(tmpIntraChromJuncStr == ""))
			break;
		int tmpChrNameInt;
		int tmpStartPos;
		int tmpEndPos;
		extractIntraChromJuncInfoFromStr(
			tmpChrNameInt, tmpStartPos, tmpEndPos,
			tmpIntraChromJuncStr, indexInfo);
		vector<string> tmpCandiExonGTFstrVec_1;
		set<string> tmpExonGeneNameSet_1;
		vector<string> tmpCandiExonGTFstrVec_2;
		set<string> tmpExonGeneNameSet_2;
		ifstream gtf_ifs(inputGTFannotationPath.c_str());
		while(!gtf_ifs.eof())
		{
			string tmpGTFstr;
			getline(gtf_ifs, tmpGTFstr);
			if((gtf_ifs.eof())||(tmpGTFstr == ""))
				break;
			int tmpExonChrNameInt;
			int tmpExonChrStartPos;
			int tmpExonChrEndPos;
			string tmpExonStrand;
			string tmpExonGeneName;
			bool exonInfoExtractedSuccessFully_bool 
				= extractExonChrNamePosFromGTFstr(
					tmpExonChrNameInt, tmpExonChrStartPos, tmpExonChrEndPos, 
					tmpExonStrand, tmpExonGeneName, tmpGTFstr, indexInfo);
			//cout << "exonInfoExtractedSuccessFully_bool: " << exonInfoExtractedSuccessFully_bool << endl;
			if(exonInfoExtractedSuccessFully_bool)
			{
				if(tmpChrNameInt == tmpExonChrNameInt)
				{
					if((tmpStartPos + offsetInt >= tmpExonChrStartPos)
						&&(tmpStartPos <= tmpExonChrEndPos))
					{
						tmpCandiExonGTFstrVec_1.push_back(tmpGTFstr);
						tmpExonGeneNameSet_1.insert(tmpExonGeneName + ":" + tmpExonStrand);
					}
					if((tmpEndPos + offsetInt >= tmpExonChrStartPos)
						&&(tmpEndPos - offsetInt <= tmpExonChrEndPos))
					{	
						tmpCandiExonGTFstrVec_2.push_back(tmpGTFstr);
						tmpExonGeneNameSet_2.insert(tmpExonGeneName + ":" + tmpExonStrand);
					}
				}
			}
		}

		string tmpIntraChromJuncStr_geneInfo = "";
		tmpIntraChromJuncStr_geneInfo += tmpIntraChromJuncStr;
		tmpIntraChromJuncStr_geneInfo += "\t";
		//intraChromJuncWithGeneInfo_ofs << tmpIntraChromJuncStr << "\t";
		if(tmpExonGeneNameSet_1.size() == 0)
			tmpIntraChromJuncStr_geneInfo += "NULL";
		for(set<string>::iterator tmpSetIter = tmpExonGeneNameSet_1.begin();
			tmpSetIter != tmpExonGeneNameSet_1.end(); tmpSetIter ++)
		{
			tmpIntraChromJuncStr_geneInfo += (*tmpSetIter);
			tmpIntraChromJuncStr_geneInfo += ",";
		}
		tmpIntraChromJuncStr_geneInfo += "\t";
		if(tmpExonGeneNameSet_2.size() == 0)
			tmpIntraChromJuncStr_geneInfo += "NULL";
		for(set<string>::iterator tmpSetIter = tmpExonGeneNameSet_2.begin();
			tmpSetIter != tmpExonGeneNameSet_2.end(); tmpSetIter ++)
		{
			tmpIntraChromJuncStr_geneInfo += (*tmpSetIter);
			tmpIntraChromJuncStr_geneInfo += ",";	
		}
		intraChromJuncWithGeneInfo_ofs << tmpIntraChromJuncStr_geneInfo << endl;

		bool intraGene_bool = false;
		for(set<string>::iterator tmpSetIter = tmpExonGeneNameSet_1.begin();
			tmpSetIter != tmpExonGeneNameSet_1.end(); tmpSetIter ++)
		{
			string tmpGeneName = (*tmpSetIter);
			set<string>::iterator tmpSetIter_2 = tmpExonGeneNameSet_2.find(tmpGeneName);
			if(tmpSetIter_2 != tmpExonGeneNameSet_2.end())
			{
				intraGene_bool = true;
				break;
			}
		}
		if(intraGene_bool)
			intraChromJuncWithGeneInfo_intraGene_ofs << tmpIntraChromJuncStr_geneInfo << endl;
		else
			intraChromJuncWithGeneInfo_interGene_ofs << tmpIntraChromJuncStr_geneInfo << endl;	
		gtf_ifs.close();
	}
	intraChromJuncWithGeneInfo_ofs.close();
	intraChromJuncWithGeneInfo_intraGene_ofs.close();
	intraChromJuncWithGeneInfo_interGene_ofs.close();	
	intraChromJunc_ifs.close();
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}