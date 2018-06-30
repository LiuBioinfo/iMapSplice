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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"
#include "../../../general/alignInferJunctionHash_info.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputIndexFolder inputTotalJuncFile outputFolder" << endl;
		exit(1);
	}
	string inputJuncFilePath = argv[2];
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();

	string output_backSpliceSite_total = outputFolderStr + "backSpliceSite_total.txt";
	string output_backSpliceSite_alterSite = outputFolderStr + "backSpliceSite_alterSite.txt";
	string output_backSpliceSite_alterBackSiteOnly = outputFolderStr + "backSpliceSite_alterBackSiteOnly.txt";
	string output_backSpliceSite_alterNormalSiteOnly = outputFolderStr + "backSpliceSite_alterNormalSiteOnly.txt";
	ofstream backSpliceSite_total_ofs(output_backSpliceSite_total.c_str());
	ofstream backSpliceSite_alterSite_ofs(output_backSpliceSite_alterSite.c_str());
	ofstream backSpliceSite_alterBackSiteOnly_ofs(output_backSpliceSite_alterBackSiteOnly.c_str());
	ofstream backSpliceSite_alterNormalSiteOnly_ofs(output_backSpliceSite_alterNormalSiteOnly.c_str());

	cout << "start to initaite merged alignInferJunctionHashInfo " << endl;
	log_ofs << "start to initaite merged alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_total = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_total->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo_total->insertJuncFromJuncFile_chrNamePos_supportNum_includeBackSplice(inputJuncFilePath, indexInfo);


	cout << "start to initiate SJhashInfo" << endl;
	log_ofs << "start to initiate SJhashInfo" << endl;	
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	cout << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	log_ofs << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo_total->convert2SJhashInfo(SJ, indexInfo);


	cout << "start to initiate donerSiteSetVec and acceptorSiteSetVec ..." << endl;
	vector< set<int> > donerSiteSetVec;
	vector< set<int> > acceptorSiteSetVec;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpDonerSiteSet;
		set<int> tmpAcceptorSiteSet;
		donerSiteSetVec.push_back(tmpDonerSiteSet);
		acceptorSiteSetVec.push_back(tmpAcceptorSiteSet);
	}

	cout << "start to insert splice sites 2 donerSiteSetVec and acceptorSiteSetVec ..." << endl;
	int alignInferInfoVecSize = alignInferJunctionHashInfo_total->returnAlignInferInfoVecSize();
	for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
	{	
		int tmpSJ_chrNameInt = alignInferJunctionHashInfo_total->returnAlignInferInfo_chrNameInt(tmp);
		//string tmpSJ_chrNameStr = indexInfo->returnChrNameStr(tmpSJ_chrNameInt);
		int tmpSJ_donerEndPos = alignInferJunctionHashInfo_total->returnAlignInferInfo_donerEndPos(tmp);
		int tmpSJ_acceptorStartPos = alignInferJunctionHashInfo_total->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpSJ_supportNum = alignInferJunctionHashInfo_total->returnAlignInferInfo_supportNum(tmp);
		int tmpSJ_juncSize = tmpSJ_acceptorStartPos - tmpSJ_donerEndPos - 1;
		if(tmpSJ_juncSize >= 0)
			continue;
		donerSiteSetVec[tmpSJ_chrNameInt].insert(tmpSJ_donerEndPos);
		acceptorSiteSetVec[tmpSJ_chrNameInt].insert(tmpSJ_acceptorStartPos);
	}			

	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		cout << "tmpChr: " << tmpChr << endl;
		string tmpChrName = indexInfo->returnChrNameStr(tmpChr);
		for(set<int>::iterator tmpIntSetIter = donerSiteSetVec[tmpChr].begin();
			tmpIntSetIter != donerSiteSetVec[tmpChr].end(); tmpIntSetIter ++)
		{
			int tmpDonerSite = *tmpIntSetIter;
			vector<int> tmpAlterAcceptorSiteVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChr, tmpDonerSite, tmpAlterAcceptorSiteVec);
			vector< pair<int,int> > tmpAlterAcceptorSiteSupNumPairVec_backSplice; // (acceptorSite, supNum)
			vector< pair<int,int> > tmpAlterAcceptorSiteSupNumPairVec_normalSplice;
			for(int tmpAlterSiteIndex = 0; tmpAlterSiteIndex < tmpAlterAcceptorSiteVec.size(); tmpAlterSiteIndex ++)
			{
				int tmpAlterAcceptorSite = tmpAlterAcceptorSiteVec[tmpAlterSiteIndex];
				int tmpSJ_supNum = alignInferJunctionHashInfo_total->searchAndReturnAlignInferJuncHashSupNum(tmpChr, tmpDonerSite, tmpAlterAcceptorSite);
				if(tmpDonerSite < tmpAlterAcceptorSite) // normal splice
					tmpAlterAcceptorSiteSupNumPairVec_normalSplice.push_back(pair<int,int> (tmpAlterAcceptorSite, tmpSJ_supNum));
				else // back splice
					tmpAlterAcceptorSiteSupNumPairVec_backSplice.push_back(pair<int,int> (tmpAlterAcceptorSite, tmpSJ_supNum));
			}
			
			int tmpNormalSpliceNum = tmpAlterAcceptorSiteSupNumPairVec_normalSplice.size();
			int tmpBackSpliceNum = tmpAlterAcceptorSiteSupNumPairVec_backSplice.size();
			int tmpNormalSplice_supNumSum = 0;
			for(int tmpNormalSpliceIndex = 0; tmpNormalSpliceIndex < tmpNormalSpliceNum; tmpNormalSpliceIndex++)
				tmpNormalSplice_supNumSum += (tmpAlterAcceptorSiteSupNumPairVec_normalSplice[tmpNormalSpliceIndex]).second;
			int tmpBackSplice_supNumSum = 0;
			for(int tmpBackSpliceIndex = 0; tmpBackSpliceIndex < tmpBackSpliceNum; tmpBackSpliceIndex ++)
				tmpBackSplice_supNumSum += (tmpAlterAcceptorSiteSupNumPairVec_backSplice[tmpBackSpliceIndex]).second;
			
			string tmpSiteStr = tmpChrName + "\t" + int_to_str(tmpDonerSite) + "\tDoner\t" 
				+ int_to_str(tmpBackSpliceNum) + "\t" + int_to_str(tmpBackSplice_supNumSum)
				+ "\t" + int_to_str(tmpNormalSpliceNum) + "\t" + int_to_str(tmpNormalSplice_supNumSum);
			tmpSiteStr += "\t";
			for(int tmpBackSpliceIndex = 0; tmpBackSpliceIndex < tmpBackSpliceNum; tmpBackSpliceIndex ++)
			{
				tmpSiteStr += int_to_str((tmpAlterAcceptorSiteSupNumPairVec_backSplice[tmpBackSpliceIndex]).first);
				tmpSiteStr += ":";
				tmpSiteStr += int_to_str((tmpAlterAcceptorSiteSupNumPairVec_backSplice[tmpBackSpliceIndex]).second);
				tmpSiteStr += ";";
			}
			if(tmpBackSpliceNum == 0)
				tmpSiteStr += "NULL";
			tmpSiteStr += "\t";
			for(int tmpNormalSpliceIndex = 0; tmpNormalSpliceIndex < tmpNormalSpliceNum; tmpNormalSpliceIndex++)
			{
				tmpSiteStr += int_to_str((tmpAlterAcceptorSiteSupNumPairVec_normalSplice[tmpNormalSpliceIndex]).first);
				tmpSiteStr += ":";
				tmpSiteStr += int_to_str((tmpAlterAcceptorSiteSupNumPairVec_normalSplice[tmpNormalSpliceIndex]).second);
				tmpSiteStr += ";";
			}
			if(tmpNormalSpliceNum == 0)
				tmpSiteStr += "NULL";

			if(tmpNormalSpliceNum + tmpBackSpliceNum == 0)
			{
				cout << "error ! tmpNormalSpliceNum + tmpBackSpliceNum == 0" << endl;
				exit(1);
			}
			else
			{
				backSpliceSite_total_ofs << tmpSiteStr << endl;
				if(tmpNormalSpliceNum + tmpBackSpliceNum != 1)
				{
					backSpliceSite_alterSite_ofs << tmpSiteStr << endl;
					if(tmpNormalSpliceNum == 0)
						backSpliceSite_alterBackSiteOnly_ofs << tmpSiteStr << endl;
					if(tmpBackSpliceNum == 0)
						backSpliceSite_alterNormalSiteOnly_ofs << tmpSiteStr << endl;
				}
			}
		}
		for(set<int>::iterator tmpIntSetIter = acceptorSiteSetVec[tmpChr].begin();
			tmpIntSetIter != acceptorSiteSetVec[tmpChr].end(); tmpIntSetIter ++)
		{
			int tmpAcceptorSite = *tmpIntSetIter;
			vector<int> tmpAlterDonerSiteVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChr, tmpAcceptorSite, tmpAlterDonerSiteVec);
			vector< pair<int,int> > tmpAlterDonerSiteSupNumPairVec_backSplice;
			vector< pair<int,int> > tmpAlterDonerSiteSupNumPairVec_normalSplice;
			for(int tmpAlterSiteIndex = 0; tmpAlterSiteIndex < tmpAlterDonerSiteVec.size(); tmpAlterSiteIndex ++)
			{
				int tmpAlterDonerSite = tmpAlterDonerSiteVec[tmpAlterSiteIndex];
				int tmpSJ_supNum = alignInferJunctionHashInfo_total->searchAndReturnAlignInferJuncHashSupNum(tmpChr, tmpAlterDonerSite, tmpAcceptorSite);
				if(tmpAlterDonerSite < tmpAcceptorSite) // normal splice
					tmpAlterDonerSiteSupNumPairVec_normalSplice.push_back(pair<int,int> (tmpAlterDonerSite, tmpSJ_supNum));
				else // back splice
					tmpAlterDonerSiteSupNumPairVec_backSplice.push_back(pair<int,int> (tmpAlterDonerSite, tmpSJ_supNum));
			}
			int tmpNormalSpliceNum = tmpAlterDonerSiteSupNumPairVec_normalSplice.size();
			int tmpBackSpliceNum = tmpAlterDonerSiteSupNumPairVec_backSplice.size();
			int tmpNormalSplice_supNumSum = 0;
			for(int tmpNormalSpliceIndex = 0; tmpNormalSpliceIndex < tmpNormalSpliceNum; tmpNormalSpliceIndex++)
				tmpNormalSplice_supNumSum += (tmpAlterDonerSiteSupNumPairVec_normalSplice[tmpNormalSpliceIndex]).second;
			int tmpBackSplice_supNumSum = 0;
			for(int tmpBackSpliceIndex = 0; tmpBackSpliceIndex < tmpBackSpliceNum; tmpBackSpliceIndex ++)
				tmpBackSplice_supNumSum += (tmpAlterDonerSiteSupNumPairVec_backSplice[tmpBackSpliceIndex]).second;
			string tmpSiteStr = tmpChrName + "\t" + int_to_str(tmpAcceptorSite) + "\tAcceptor\t" 
				+ int_to_str(tmpBackSpliceNum) + "\t" + int_to_str(tmpBackSplice_supNumSum)
				+ "\t" + int_to_str(tmpNormalSpliceNum) + "\t" + int_to_str(tmpNormalSplice_supNumSum);			
			tmpSiteStr += "\t";
			for(int tmpBackSpliceIndex = 0; tmpBackSpliceIndex < tmpBackSpliceNum; tmpBackSpliceIndex ++)
			{
				tmpSiteStr += int_to_str((tmpAlterDonerSiteSupNumPairVec_backSplice[tmpBackSpliceIndex]).first);
				tmpSiteStr += ":";
				tmpSiteStr += int_to_str((tmpAlterDonerSiteSupNumPairVec_backSplice[tmpBackSpliceIndex]).second);
				tmpSiteStr += ";";
			}
			if(tmpBackSpliceNum == 0)
				tmpSiteStr += "NULL";
			tmpSiteStr += "\t";
			for(int tmpNormalSpliceIndex = 0; tmpNormalSpliceIndex < tmpNormalSpliceNum; tmpNormalSpliceIndex++)
			{
				tmpSiteStr += int_to_str((tmpAlterDonerSiteSupNumPairVec_normalSplice[tmpNormalSpliceIndex]).first);
				tmpSiteStr += ":";
				tmpSiteStr += int_to_str((tmpAlterDonerSiteSupNumPairVec_normalSplice[tmpNormalSpliceIndex]).second);
				tmpSiteStr += ";";
			}
			if(tmpNormalSpliceNum == 0)
				tmpSiteStr += "NULL";

			if(tmpNormalSpliceNum + tmpBackSpliceNum == 0)
			{
				cout << "error ! tmpNormalSpliceNum + tmpBackSpliceNum == 0" << endl;
				exit(1);
			}
			else
			{
				backSpliceSite_total_ofs << tmpSiteStr << endl;
				if(tmpNormalSpliceNum + tmpBackSpliceNum != 1)
				{
					backSpliceSite_alterSite_ofs << tmpSiteStr << endl;
					if(tmpNormalSpliceNum == 0)
						backSpliceSite_alterBackSiteOnly_ofs << tmpSiteStr << endl;
					if(tmpBackSpliceNum == 0)
						backSpliceSite_alterNormalSiteOnly_ofs << tmpSiteStr << endl;
				}
			}			
		}
	}

	parameter_ifs.close();
	chrom_bit_file_ifs.close();
	delete alignInferJunctionHashInfo_total;
	delete indexInfo;
	backSpliceSite_alterNormalSiteOnly_ofs.close();
	backSpliceSite_alterBackSiteOnly_ofs.close();
	backSpliceSite_alterSite_ofs.close();
	backSpliceSite_total_ofs.close();
	log_ofs.close();
	return 0;
}