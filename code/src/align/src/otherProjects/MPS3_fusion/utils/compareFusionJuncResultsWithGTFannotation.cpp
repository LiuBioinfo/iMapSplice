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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;

void extractFusionJuncChrNamePosFromFusionJuncStr(
	int& tmpChrNameInt_1, int& tmpChrNameInt_2, 
	int& tmpChrPosInt_1, int& tmpChrPosInt_2, //int& tmpSupNum,
	string& tmpStrand_1, string& tmpStrand_2,
	string& tmpFusionJuncStr, Index_Info* indexInfo)
{
	vector<string> tmpFusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpFusionJuncStr.find("\t", startLoc);
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}
	//tmpFusionJuncFieldVec.push_back(tmpFusionJuncStr.substr(startLoc));

	string tmpChrNameStr_1 = tmpFusionJuncFieldVec[0];
	string tmpChrNameStr_2 = tmpFusionJuncFieldVec[1];
 	string tmpChrPosStr_1 = tmpFusionJuncFieldVec[2];
	string tmpChrPosStr_2 = tmpFusionJuncFieldVec[3];
	//string tmpSupportNumStr = tmpFusionJuncFieldVec[4];
	tmpStrand_1 = tmpFusionJuncFieldVec[4];
	tmpStrand_2 = tmpFusionJuncFieldVec[5];

	tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpChrNameStr_1);
	tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpChrNameStr_2);
	tmpChrPosInt_1 = atoi(tmpChrPosStr_1.c_str());
	tmpChrPosInt_2 = atoi(tmpChrPosStr_2.c_str());
	//tmpSupNum = atoi(tmpSupportNumStr.c_str());
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
	if((tmpExonStrand != "+")&&(tmpExonStrand != "-"))
		return false;

	int tmpGeneNameLoc = tmpExonStr.find("gene_name");
	if(tmpGeneNameLoc == string::npos)
		return false;
	int tmpGeneNameLoc_nextSemiColon = tmpExonStr.find(";", tmpGeneNameLoc);
	if(tmpGeneNameLoc_nextSemiColon == string::npos)
		return false;
	tmpGeneName = tmpExonStr.substr(tmpGeneNameLoc+11, tmpGeneNameLoc_nextSemiColon - 2 - tmpGeneNameLoc - 11 + 1);

	int tmpGeneTypeLoc = tmpExonStr.find("gene_biotype");
	if(tmpGeneTypeLoc == string::npos)
		return false;
	int tmpGeneTypeLoc_nextSemiColon = tmpExonStr.find(";", tmpGeneTypeLoc);
	if(tmpGeneTypeLoc_nextSemiColon == string::npos)
		return false;
	string tmpGeneType = tmpExonStr.substr(tmpGeneTypeLoc+14, tmpGeneTypeLoc_nextSemiColon - 2 - tmpGeneTypeLoc - 14 + 1);
	//cout << "tmpGeneType: " << tmpGeneType << endl;
	if(tmpGeneType == "protein_coding")
	{}
	else
		return false;

	if(tmpExon_featureStr != "exon")
		return false;
	if(tmpExon_chrNameStr.substr(0,3) == "chr")
	{}
	else
	{
		string tmpExon_chrNameStrOri = tmpExon_chrNameStr;
		tmpExon_chrNameStr = "chr" + tmpExon_chrNameStrOri;
	}
	//cout << "tmpExon_chrNameStr: " << tmpExon_chrNameStr << endl;
	tmpExonChrNameInt = indexInfo->convertStringToInt(tmpExon_chrNameStr);
	if((tmpExonChrNameInt < 0)||(tmpExonChrNameInt >= indexInfo->returnChromNum()))
		return false;
	tmpExonChrStartPos = atoi(tmpExon_startPosStr.c_str());
	tmpExonChrEndPos = atoi(tmpExon_endPosStr.c_str());
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputIndexInfoFolderPath inputGTFannotation";
		cout << " inputFusionJuncPath outputFolderPath offset checkStrandOrNot" << endl;
		exit(1);
	}
	string offsetStr = argv[5];
	int offsetInt = atoi(offsetStr.c_str());
	string inputGTFannotationPath = argv[2];
	string inputFusionJuncPath = argv[3];

	bool checkStrandOrNot_bool;
	string checkStrandOrNotStr = argv[6];
	if((checkStrandOrNotStr == "Y")||(checkStrandOrNotStr == "y")||(checkStrandOrNotStr == "Yes"))
		checkStrandOrNot_bool = true;
	else if((checkStrandOrNotStr == "N")||(checkStrandOrNotStr == "n")||(checkStrandOrNotStr == "No"))
		checkStrandOrNot_bool = false;
	else
	{
		cout << "incorrect settings for checkStrandOrNot, argv[6]: " << argv[6] << endl;
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

	vector< vector< pair<int,int> > > regionVecVec;
	vector< vector<string> > geneNameVecVec;
	vector< vector<string> > strandVecVec;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		vector< pair<int,int> > tmpRegionVec;
		regionVecVec.push_back(tmpRegionVec);
		vector< string > tmpGeneNameVec;
		geneNameVecVec.push_back(tmpGeneNameVec);
		vector<string> tmpStrandVec;
		strandVecVec.push_back(tmpStrandVec);
	}

	ifstream gtf_ifs(inputGTFannotationPath.c_str());
	while(!gtf_ifs.eof())
	{
		string tmpGTFstr;
		getline(gtf_ifs, tmpGTFstr);
		if(tmpGTFstr == "")
			break;
		int tmpExonChrNameInt;
		int tmpExonChrStartPos;
		int tmpExonChrEndPos;
		string tmpExonStrand;
		string tmpExonGeneName;
		//cout << "tmpGTFstr: " << endl << tmpGTFstr << endl;
		bool exonInfoExtractedSuccessFully_bool 
			= extractExonChrNamePosFromGTFstr(
				tmpExonChrNameInt, tmpExonChrStartPos, tmpExonChrEndPos, 
				tmpExonStrand, tmpExonGeneName, tmpGTFstr, indexInfo);
		
		if(exonInfoExtractedSuccessFully_bool)
		{
			// cout << "exonInfoExtractedSuccessFully_bool: " << exonInfoExtractedSuccessFully_bool << endl;
			// cout << "tmpExonChrNameInt: " << tmpExonChrNameInt << endl;
			// cout << "startPos: " << tmpExonChrStartPos << endl;
			// cout << "endPos: " << tmpExonChrEndPos << endl;
			// cout << "tmpExonGeneName: " << tmpExonGeneName << endl;
			// cout << "tmpExonStrand: " << tmpExonStrand << endl;
			regionVecVec[tmpExonChrNameInt].push_back(pair<int,int>(tmpExonChrStartPos, tmpExonChrEndPos));
			geneNameVecVec[tmpExonChrNameInt].push_back(tmpExonGeneName);
			strandVecVec[tmpExonChrNameInt].push_back(tmpExonStrand);
		}
	}
	gtf_ifs.close();

	cout << "chromNum: " << chromNum << endl;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		cout << "tmpChr: " << tmp << endl;
		cout << "tmpRegionVec.size(): " << regionVecVec[tmp].size() << endl;
		cout << "tmpGeneNameVec.size(): " << geneNameVecVec[tmp].size() << endl;
		cout << "tmpStrandVec.size(): " << strandVecVec[tmp].size() << endl;
	}	

	int tmpFusionJuncNum = 0;
	ifstream fusionJunc_ifs(inputFusionJuncPath.c_str());

	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	string outputFusionJuncWithGeneInfoPath 
		= outputFolderStr + "fusionJunc_geneInfo.txt";
	ofstream fusionJuncWithGeneInfo_ofs(outputFusionJuncWithGeneInfoPath.c_str());
	string outputFusionJuncWithGeneInfoPath_intraGene 
		= outputFolderStr + "fusionJunc_geneInfo_intraGene.txt";
	ofstream fusionJuncWithGeneInfo_intraGene_ofs(outputFusionJuncWithGeneInfoPath_intraGene.c_str());
	string geneList_intraGene
		= outputFolderStr + "geneList_intraGene.txt";
	ofstream geneList_intraGene_ofs(geneList_intraGene.c_str());

	string outputFusionJuncWithGeneInfoPath_interGene 
		= outputFolderStr + "fusionJunc_geneInfo_interGene.txt";
	ofstream fusionJuncWithGeneInfo_interGene_ofs(outputFusionJuncWithGeneInfoPath_interGene.c_str());
	string geneList_interGene
		= outputFolderStr + "geneList_interGene.txt";
	ofstream geneList_interGene_ofs(geneList_interGene.c_str());

	string outputFusionJuncWithGeneInfoPath_outOfGene 
		= outputFolderStr + "fusionJunc_geneInfo_outOfGene.txt";
	ofstream fusionJuncWithGeneInfo_outOfGene_ofs(outputFusionJuncWithGeneInfoPath_outOfGene.c_str());
	string geneList_outOfGene;
	ofstream geneList_outOfGene_ofs(geneList_outOfGene.c_str());

	while(!fusionJunc_ifs.eof())
	{
		string tmpFusionJuncStr;
		getline(fusionJunc_ifs, tmpFusionJuncStr);
		if((fusionJunc_ifs.eof())||(tmpFusionJuncStr == ""))
			break;
		int tmpChrNameInt_1;
		int tmpChrNameInt_2;
		int tmpChrPosInt_1;
		int tmpChrPosInt_2;
		string tmpStrand_1;
		string tmpStrand_2;
		//int tmpSupNum;
		extractFusionJuncChrNamePosFromFusionJuncStr(
			tmpChrNameInt_1, tmpChrNameInt_2, 
			tmpChrPosInt_1, tmpChrPosInt_2, //tmpSupNum,
			tmpStrand_1, tmpStrand_2, 
			tmpFusionJuncStr, indexInfo);

		set<string> tmpExonGeneNameSet_1;
		set<string> tmpExonGeneNameSet_2;

		int tmpRegionVecSize_1 = regionVecVec[tmpChrNameInt_1].size();
		for(int tmp = 0; tmp < tmpRegionVecSize_1; tmp++)
		{
			int tmpRegion_startPos = (regionVecVec[tmpChrNameInt_1])[tmp].first;
			int tmpRegion_endPos = (regionVecVec[tmpChrNameInt_1])[tmp].second;
			string tmpRegion_geneName = (geneNameVecVec[tmpChrNameInt_1])[tmp];
			string tmpRegion_strand = (strandVecVec[tmpChrNameInt_1])[tmp];
			if(//(tmpStrand_1 == tmpRegion_strand)&&
				(tmpChrPosInt_1 + offsetInt >= tmpRegion_startPos)
				&&(tmpChrPosInt_1 - offsetInt <= tmpRegion_endPos))
			{
				if(checkStrandOrNot_bool)	
				{
					if(tmpStrand_1 == tmpRegion_strand)					
						tmpExonGeneNameSet_1.insert(tmpRegion_geneName);
				}
				else
					tmpExonGeneNameSet_1.insert(tmpRegion_geneName);
			}
		}

		int tmpRegionVecSize_2 = regionVecVec[tmpChrNameInt_2].size();
		for(int tmp = 0; tmp < tmpRegionVecSize_2; tmp++)
		{
			int tmpRegion_startPos = (regionVecVec[tmpChrNameInt_2])[tmp].first;
			int tmpRegion_endPos = (regionVecVec[tmpChrNameInt_2])[tmp].second;
			string tmpRegion_geneName = (geneNameVecVec[tmpChrNameInt_2])[tmp];
			string tmpRegion_strand = (strandVecVec[tmpChrNameInt_2])[tmp];
			if(//(tmpStrand_2 == tmpRegion_strand)&&
				(tmpChrPosInt_2 + offsetInt >= tmpRegion_startPos)
				&&(tmpChrPosInt_2 - offsetInt <= tmpRegion_endPos))
			{
				if(checkStrandOrNot_bool)
				{
					if(tmpStrand_2 == tmpRegion_strand)
						tmpExonGeneNameSet_2.insert(tmpRegion_geneName);
				}
				else
					tmpExonGeneNameSet_2.insert(tmpRegion_geneName);
			}
		}

		string tmpFusionJuncStr_geneInfo;
		string tmpGeneList_geneInfo;

		tmpFusionJuncStr_geneInfo += tmpFusionJuncStr;
		tmpFusionJuncStr_geneInfo += "\t";
		if(tmpExonGeneNameSet_1.size() == 0)
		{	
			tmpFusionJuncStr_geneInfo += "NULL";
			tmpGeneList_geneInfo += "NULL";
		}
		for(set<string>::iterator tmpSetIter = tmpExonGeneNameSet_1.begin();
			tmpSetIter != tmpExonGeneNameSet_1.end(); tmpSetIter ++)
		{
			tmpFusionJuncStr_geneInfo += (*tmpSetIter);
			tmpFusionJuncStr_geneInfo += ",";
			tmpGeneList_geneInfo += (*tmpSetIter);
			tmpGeneList_geneInfo += ",";
		}
		tmpFusionJuncStr_geneInfo += "\t";
		if(tmpExonGeneNameSet_2.size() == 0)
		{	
			tmpFusionJuncStr_geneInfo += "NULL";
			tmpGeneList_geneInfo += "NULL";
		}
		for(set<string>::iterator tmpSetIter = tmpExonGeneNameSet_2.begin();
			tmpSetIter != tmpExonGeneNameSet_2.end(); tmpSetIter ++)
		{
			tmpFusionJuncStr_geneInfo += (*tmpSetIter);
			tmpFusionJuncStr_geneInfo += ",";
			tmpGeneList_geneInfo += (*tmpSetIter);
			tmpGeneList_geneInfo += ",";
		}
		fusionJuncWithGeneInfo_ofs << tmpFusionJuncStr_geneInfo << endl;

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
		{
			fusionJuncWithGeneInfo_intraGene_ofs << tmpFusionJuncStr_geneInfo << endl;
			geneList_intraGene_ofs << tmpGeneList_geneInfo << endl;
		}
		else if((tmpExonGeneNameSet_1.size() > 0) && (tmpExonGeneNameSet_2.size() > 0))
		{
			fusionJuncWithGeneInfo_interGene_ofs << tmpFusionJuncStr_geneInfo << endl;
			geneList_interGene_ofs << tmpGeneList_geneInfo << endl;
		}
		else
		{
			fusionJuncWithGeneInfo_outOfGene_ofs << tmpFusionJuncStr_geneInfo << endl;
			geneList_outOfGene_ofs << tmpGeneList_geneInfo << endl;
		}
		gtf_ifs.close();
	}
	fusionJuncWithGeneInfo_ofs.close();
	geneList_intraGene_ofs.close();
	geneList_interGene_ofs.close();
	geneList_outOfGene_ofs.close();
	fusionJuncWithGeneInfo_intraGene_ofs.close();
	fusionJuncWithGeneInfo_interGene_ofs.close();
	fusionJuncWithGeneInfo_outOfGene_ofs.close();
	fusionJunc_ifs.close();
	delete indexInfo;
	parameter_ifs.close();
	return 0;
}