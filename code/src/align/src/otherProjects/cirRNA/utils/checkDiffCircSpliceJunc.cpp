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
	// string output_backSpliceOnly = outputFolderStr + "backSpliceOnly_junc.txt";
	// ofstream backSpliceOnly_ofs(output_backSpliceOnly.c_str());
	// string output_backNormalSpliceBoth = outputFolderStr + "backNormalSpliceBoth_junc.txt";
	// ofstream backNormalSpliceBoth_ofs(output_backNormalSpliceBoth.c_str());
	string output_backSpliceTotal = outputFolderStr + "backSplice_total.txt";
	ofstream backSpliceTotal_ofs(output_backSpliceTotal.c_str());

	string output_backSplice_alterSiteOnBothEnds = outputFolderStr + "backSplice_alterSiteOnBothEnds.txt";
	ofstream backSplice_alterSiteOnBothEnds_ofs(output_backSplice_alterSiteOnBothEnds.c_str());

	string output_backSplice_alterSiteOnOneEnd = outputFolderStr + "backSplice_alterSiteOnOneEnd.txt";
	ofstream backSplice_alterSiteOnOneEnd_ofs(output_backSplice_alterSiteOnOneEnd.c_str());

	string output_backSplice_noAlterSite = outputFolderStr + "backSplice_noAlterSite.txt";
	ofstream backSplice_noAlterSite_ofs(output_backSplice_noAlterSite.c_str());

	string output_backSplice_withAlterSite = outputFolderStr + "backSplice_withAlterSite.txt";
	ofstream backSplice_withAlterSite_ofs(output_backSplice_withAlterSite.c_str());

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	indexInfo->readGenome(chrom);
	indexInfo->initiate();	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	indexInfo->initiateChrNameIndexArray(1000);
	int chromNum = indexInfo->returnChromNum();

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite merged alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite merged alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_total = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo_total->initiateAlignInferJunctionHashInfo(chromNum);
	alignInferJunctionHashInfo_total->insertJuncFromJuncFile_chrNamePos_supportNum_includeBackSplice(
		inputJuncFilePath, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initiate SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initiate SJhashInfo" << endl;	
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo_total->convert2SJhashInfo(SJ, indexInfo);

	int backSpliceNum_total = 0;
	int backSpliceNum_alterSiteOnBothEnds = 0;
	int backSpliceNum_alterSiteOnOneEnd = 0;
	int backSpliceNum_withAlterSite = 0;
	int backSpliceNum_noAlterSite = 0;
	int alignInferInfoVecSize = alignInferJunctionHashInfo_total->returnAlignInferInfoVecSize();
	for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
	{	
		int tmpSJ_chrNameInt = alignInferJunctionHashInfo_total->returnAlignInferInfo_chrNameInt(tmp);
		string tmpSJ_chrNameStr = indexInfo->returnChrNameStr(tmpSJ_chrNameInt);
		int tmpSJ_donerEndPos = alignInferJunctionHashInfo_total->returnAlignInferInfo_donerEndPos(tmp);
		int tmpSJ_acceptorStartPos = alignInferJunctionHashInfo_total->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpSJ_supportNum = alignInferJunctionHashInfo_total->returnAlignInferInfo_supportNum(tmp);
		int tmpSJ_juncSize = tmpSJ_acceptorStartPos - tmpSJ_donerEndPos - 1;
		if(tmpSJ_juncSize >= 0)
			continue;
		string tmpSJinfoStr_withAlterSite = tmpSJ_chrNameStr + "\t" + int_to_str(tmpSJ_donerEndPos) + "\t"
			+ int_to_str(tmpSJ_acceptorStartPos) + "\t" + int_to_str(tmpSJ_supportNum) + "\t";
		// backSpliceTotal_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_donerEndPos << "\t" << tmpSJ_acceptorStartPos << "\t"
		// 		<< tmpSJ_supportNum << "\t";
		vector<int> alterAcceptorSiteVec;
		vector<int> alterDonerSiteVec;
		SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpSJ_chrNameInt, tmpSJ_donerEndPos, alterAcceptorSiteVec);
		SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpSJ_chrNameInt, tmpSJ_acceptorStartPos, alterDonerSiteVec);
		string alterAcceptorStr;
		string alterDonerStr;
		if((alterAcceptorSiteVec.size() == 0)||(alterDonerSiteVec.size() == 0))
		{
			cout << "error in alterSiteVecSize !" << endl;
			cout << "alterAcceptorSiteVec.size(): " << alterAcceptorSiteVec.size() << endl;
			cout << "alterDonerSiteVec.size(): " << alterDonerSiteVec.size() << endl;
			exit(1);
		}
		int tmpSJ_alterAcceptorSiteSupNumSum = 0;
		for(int tmp = 0; tmp < alterAcceptorSiteVec.size(); tmp++)
		{
			int alterAcceptorSite = alterAcceptorSiteVec[tmp];
			if(alterAcceptorSite != tmpSJ_acceptorStartPos)
			{
				int tmpAlterSJ_supNum = alignInferJunctionHashInfo_total->searchAndReturnAlignInferJuncHashSupNum(
					tmpSJ_chrNameInt, tmpSJ_donerEndPos, alterAcceptorSite);
				if(tmpAlterSJ_supNum == 0)
				{
					cout << "error in tmpAlterSJ_supNum: " << tmpAlterSJ_supNum << endl;
					exit(1);
				}
				tmpSJ_alterAcceptorSiteSupNumSum += tmpAlterSJ_supNum;
				alterAcceptorStr += int_to_str(alterAcceptorSite);
				alterAcceptorStr += ":";
				alterAcceptorStr += int_to_str(tmpAlterSJ_supNum);
				alterAcceptorStr += ";";
			}
		}
		int tmpSJ_alterDonerSiteSupNumSum = 0;
		for(int tmp = 0; tmp < alterDonerSiteVec.size(); tmp++)
		{
			int alterDonerSite = alterDonerSiteVec[tmp];
			if(alterDonerSite != tmpSJ_donerEndPos)
			{
				int tmpAlterSJ_supNum = alignInferJunctionHashInfo_total->searchAndReturnAlignInferJuncHashSupNum(
					tmpSJ_chrNameInt, alterDonerSite, tmpSJ_acceptorStartPos);
				if(tmpAlterSJ_supNum == 0)
				{
					cout << "error in tmpAlterSJ_supNum: " << tmpAlterSJ_supNum << endl;
					exit(1);
				}
				tmpSJ_alterDonerSiteSupNumSum += tmpAlterSJ_supNum;
				alterDonerStr += int_to_str(alterDonerSite);
				alterDonerStr += ":";
				alterDonerStr += int_to_str(tmpAlterSJ_supNum);
				alterDonerStr += ";";
			}
		}
		tmpSJinfoStr_withAlterSite = tmpSJinfoStr_withAlterSite + int_to_str(tmpSJ_alterAcceptorSiteSupNumSum) + "\t"
			+ int_to_str(tmpSJ_alterDonerSiteSupNumSum) + "\t";
		
		if(tmpSJ_alterAcceptorSiteSupNumSum == 0)
			tmpSJinfoStr_withAlterSite += "NULL";
		else
			tmpSJinfoStr_withAlterSite += alterAcceptorStr;
		tmpSJinfoStr_withAlterSite += "\t";
		if(tmpSJ_alterDonerSiteSupNumSum == 0)
			tmpSJinfoStr_withAlterSite += "NULL";
		else
			tmpSJinfoStr_withAlterSite += alterDonerStr;
		backSpliceTotal_ofs << tmpSJinfoStr_withAlterSite << endl;
		backSpliceNum_total ++;
		if((tmpSJ_alterAcceptorSiteSupNumSum != 0)&&(tmpSJ_alterDonerSiteSupNumSum != 0))
		{
			backSplice_alterSiteOnBothEnds_ofs << tmpSJinfoStr_withAlterSite << endl;
			backSpliceNum_alterSiteOnBothEnds ++;
			backSplice_withAlterSite_ofs << tmpSJinfoStr_withAlterSite << endl;
			backSpliceNum_withAlterSite ++;
		}
		else if((tmpSJ_alterAcceptorSiteSupNumSum != 0)||(tmpSJ_alterDonerSiteSupNumSum != 0))
		{
			backSplice_alterSiteOnOneEnd_ofs << tmpSJinfoStr_withAlterSite << endl;
			backSpliceNum_alterSiteOnOneEnd ++;
			backSplice_withAlterSite_ofs << tmpSJinfoStr_withAlterSite << endl;
			backSpliceNum_withAlterSite ++;
		}
		else
		{
			backSplice_noAlterSite_ofs << tmpSJinfoStr_withAlterSite << endl;
			backSpliceNum_noAlterSite ++;
		}
	}
	log_ofs << "backSpliceNum_total: " << backSpliceNum_total << endl;
	log_ofs << "    backSpliceNum_withAlterSite: " << backSpliceNum_withAlterSite << endl;
	log_ofs << "        backSpliceNum_alterSiteOnBothEnds: " << backSpliceNum_alterSiteOnBothEnds << endl;
	log_ofs << "        backSpliceNum_alterSiteOnOneEnd: " << backSpliceNum_alterSiteOnOneEnd << endl;
	log_ofs << "    backSpliceNum_noAlterSite: " << backSpliceNum_noAlterSite << endl;
	backSpliceTotal_ofs.close();
	delete alignInferJunctionHashInfo_total;
	backSplice_alterSiteOnBothEnds_ofs.close();
	backSplice_alterSiteOnOneEnd_ofs.close();
	backSplice_noAlterSite_ofs.close();
	backSplice_withAlterSite_ofs.close();
	log_ofs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}