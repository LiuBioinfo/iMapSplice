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

void extractSJchrNamePosSupNumFromSJstr(string& tmpSJstr, string& tmpSJchrName, 
	int& tmpSJstartPos, int& tmpSJendPos, int& tmpSJsupNum)
{
	int tabLoc_1 = tmpSJstr.find("\t");
	int tabLoc_2 = tmpSJstr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSJstr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpSJstr.find("\t", tabLoc_3 + 1);
	if(tabLoc_4 == string::npos)
	{
		cout << "incorrect SJ format" << endl;
		exit(1);
	}
	int tabLoc_5 = tmpSJstr.find("\t", tabLoc_4 + 1);	
	tmpSJchrName = tmpSJstr.substr(0,tabLoc_1);
	string tmpStartPosStr = tmpSJstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpEndPosStr = tmpSJstr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
	string tmpIDstr = tmpSJstr.substr(tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
	string tmpSupNumStr;
	if(tabLoc_5 == string::npos)
		tmpSupNumStr = tmpSJstr.substr(tabLoc_4 + 1);
	else
		tmpSupNumStr = tmpSJstr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);
	tmpSJstartPos = atoi(tmpStartPosStr.c_str());
	tmpSJendPos = atoi(tmpEndPosStr.c_str());
	tmpSJsupNum = atoi(tmpSupNumStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndex inputSdSJ inputPsSJ outputFolder" << endl;
		exit(1);
	}
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());	

	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	//parameter_ifs.close();
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	parameter_ifs.close();
	int chromNum = indexInfo->returnChromNum();
	cout << "chromNum: " << chromNum << endl;	

	cout << "start to initiate Standard and Personal alignInferJunctionHashInfo ...." << endl;
	log_ofs << "start to initiate Standard and Personal alignInferJunctionHashInfo ...." << endl;
	string juncHash_file_sd_str = argv[2];
	string juncHash_file_ps_str = argv[3];
	AlignInferJunctionHash_Info* juncHash_sd = new AlignInferJunctionHash_Info();
	AlignInferJunctionHash_Info* juncHash_ps = new AlignInferJunctionHash_Info();
	juncHash_sd->initiateAlignInferJunctionInfo(chromNum);
	juncHash_ps->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to read 2 juncfiles ...." << endl;
	log_ofs << "start to read 2 juncfiles ...." << endl;
	juncHash_sd->insertJuncFromJuncFile_chrNamePos_supportNum(juncHash_file_sd_str, indexInfo);
	juncHash_ps->insertJuncFromJuncFile_chrNamePos_supportNum(juncHash_file_ps_str, indexInfo);

	cout << "start to compare junctions" << endl;
	log_ofs << "start to compare junctions" << endl;	
	string compare_path = outputFolderStr + "/compare.txt";
	ofstream compare_ofs(compare_path.c_str());
	string sharedSJ_path = outputFolderStr + "/shared_SJ.txt";
	ofstream sharedSJ_ofs(sharedSJ_path.c_str());
	string sharedSJ_Sd_path = outputFolderStr + "/shared_SJ_Sd.txt";
	ofstream sharedSJ_Sd_ofs(sharedSJ_Sd_path.c_str());
	string sharedSJ_Ps_path = outputFolderStr + "/shared_SJ_Ps.txt";
	ofstream sharedSJ_Ps_ofs(sharedSJ_Ps_path.c_str());
	string SdOnlySJ_path = outputFolderStr + "/Sd_only_SJ.txt";
	ofstream SdOnlySJ_ofs(SdOnlySJ_path.c_str());
	string PsOnlySJ_path = outputFolderStr + "/Ps_only_SJ.txt";
	ofstream PsOnlySJ_ofs(PsOnlySJ_path.c_str());

	ifstream SdSJ_ifs(juncHash_file_sd_str.c_str());
	ifstream PsSJ_ifs(juncHash_file_ps_str.c_str());

	int SdSJfoundInPs_num = 0;
	int PsSJfoundInSd_num = 0;
	int SdSJnotFoundInPs_num = 0;
	int PsSJnotFoundInSd_num = 0;

	cout << "start to compare Sd SJ to Ps Sj" << endl;
	log_ofs << "start to compare Sd SJ to Ps Sj" << endl;
	while(!SdSJ_ifs.eof())
	{
		string tmpSdSJstr;
		getline(SdSJ_ifs, tmpSdSJstr);
		if(tmpSdSJstr == "")
			break;
		string tmpSdSJchrName;
		int tmpSdSJstartPos, tmpSdSJendPos, tmpSdSJsupNum;
		extractSJchrNamePosSupNumFromSJstr(tmpSdSJstr, tmpSdSJchrName, 
			tmpSdSJstartPos, tmpSdSJendPos, tmpSdSJsupNum);
		int tmpPsSJsupNum;
		bool foundInPsJuncHashOrNotBool 
			= juncHash_ps->searchAndReturnSupNumInAlignInferJuncHash(
				tmpPsSJsupNum, tmpSdSJchrName, tmpSdSJstartPos, tmpSdSJendPos, indexInfo);
		if(foundInPsJuncHashOrNotBool)
		{
			SdSJfoundInPs_num ++;
			sharedSJ_ofs << tmpSdSJchrName << "\t" << tmpSdSJstartPos << "\t"
			 	<< tmpSdSJendPos << "\t" << tmpSdSJsupNum << "\t" << tmpPsSJsupNum << endl;
			sharedSJ_Sd_ofs << tmpSdSJstr << endl; 
		}
		else
		{
			SdSJnotFoundInPs_num ++;
			SdOnlySJ_ofs << tmpSdSJstr << endl;
		}
	}

	cout << "start to compare Ps SJ to Sd Sj" << endl;
	log_ofs << "start to compare Ps SJ to Sd Sj" << endl;
	while(!PsSJ_ifs.eof())
	{
		string tmpPsSJstr;
		getline(PsSJ_ifs, tmpPsSJstr);
		if(tmpPsSJstr == "")
			break;
		string tmpPsSJchrName;
		int tmpPsSJstartPos, tmpPsSJendPos, tmpPsSJsupNum;
		extractSJchrNamePosSupNumFromSJstr(tmpPsSJstr, tmpPsSJchrName, 
			tmpPsSJstartPos, tmpPsSJendPos, tmpPsSJsupNum);
		int tmpSdSJsupNum;
		bool foundInSdJuncHashOrNotBool 
			= juncHash_sd->searchAndReturnSupNumInAlignInferJuncHash(
				tmpSdSJsupNum, tmpPsSJchrName, tmpPsSJstartPos, tmpPsSJendPos, indexInfo);
		if(foundInSdJuncHashOrNotBool)
		{
			PsSJfoundInSd_num ++;
			sharedSJ_Ps_ofs << tmpPsSJstr << endl; 
		}
		else
		{
			PsSJnotFoundInSd_num ++;
			PsOnlySJ_ofs << tmpPsSJstr << endl;
		}
	}
	cout << "end of all the comparison" << endl;
	log_ofs << "end of all the comparison" << endl;

	cout << endl << "SdSJfoundInPs_num: " << SdSJfoundInPs_num << endl;
	cout << "SdSJnotFoundInPs_num: " << SdSJnotFoundInPs_num << endl;
	cout << endl << "PsSJfoundInSd_num: " << PsSJfoundInSd_num << endl;
	cout << "PsSJnotFoundInSd_num: " << PsSJnotFoundInSd_num << endl;
	log_ofs << endl << "SdSJfoundInPs_num: " << SdSJfoundInPs_num << endl;
	log_ofs << "SdSJnotFoundInPs_num: " << SdSJnotFoundInPs_num << endl;
	log_ofs << endl << "PsSJfoundInSd_num: " << PsSJfoundInSd_num << endl;
	log_ofs << "PsSJnotFoundInSd_num: " << PsSJnotFoundInSd_num << endl;

	double sharedSJpercInSdSJ 
		= (double)SdSJfoundInPs_num * 100 / (double)(SdSJfoundInPs_num + SdSJnotFoundInPs_num);
	double uniqueSJpercInSdSJ = 100 - sharedSJpercInSdSJ;
	double sharedSJpercInPsSJ 
		= (double)PsSJfoundInSd_num * 100 / (double)(PsSJfoundInSd_num + PsSJnotFoundInSd_num);
	double uniqueSJpercInPsSJ = 100 - sharedSJpercInPsSJ;
	compare_ofs << endl << "SdSJ_num: " << SdSJfoundInPs_num + SdSJnotFoundInPs_num << endl;
	compare_ofs << "shared_Standard_SJ_num: " << SdSJfoundInPs_num << " -- " << sharedSJpercInSdSJ << "%" << endl;
	compare_ofs << "unique_Standard_SJ_num: " << SdSJnotFoundInPs_num << " -- " << uniqueSJpercInSdSJ << "%" << endl;
	compare_ofs << endl << "PsSJ_num: " << PsSJfoundInSd_num + PsSJnotFoundInSd_num << endl;
	compare_ofs << "shared_Personal_SJ_num: " << PsSJfoundInSd_num << " -- " << sharedSJpercInPsSJ << "%" << endl;
	compare_ofs << "unique_Personal_SJ_num: " << PsSJnotFoundInSd_num << " -- " << uniqueSJpercInPsSJ << "%" << endl;

	SdSJ_ifs.close();
	PsSJ_ifs.close();
	PsOnlySJ_ofs.close();
	SdOnlySJ_ofs.close();
	sharedSJ_ofs.close();
	sharedSJ_Sd_ofs.close();
	sharedSJ_Ps_ofs.close();
	compare_ofs.close();
	delete juncHash_sd;
	delete juncHash_ps;
	delete indexInfo;
	log_ofs.close();
	parameter_ifs.close();
	return 0;
}