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

int main(int argc, char** argv)
{
	if(argc <= 7)
	{
		cout << "Executable <InputIndexFolderPath> <SJsizeMin> <SJmaxMin> ";
		cout << " <outputFolder> SJ_1 SJ_2 (SJ_3 ... SJ_n)" << endl;
		exit(1);
	}
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "SJsizeMin: " << argv[2] << endl;
	log_ofs << "SJmaxMin: " << argv[3] << endl;
	log_ofs << "outputFolder: " << argv[4] << endl;

	string juncExistence = outputFolderStr + "juncExistence.txt";
	ofstream juncExistence_ofs(juncExistence.c_str());
	string invalidJuncFile = outputFolderStr + "invalidJunc.txt";
	ofstream invalidJunc_ofs(invalidJuncFile.c_str());
	vector<string> juncFileVec;
	for(int tmp = 5; tmp <= argc - 1; tmp++)
	{
		string tmpJuncFile = argv[tmp];
		juncFileVec.push_back(tmpJuncFile);
		log_ofs << "totalJuncFile: " << tmpJuncFile << endl;
	}
	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	cout << "assigning parameters ......" << endl;
	log_ofs << "assigning parameters ......" << endl;
	string min_intron_size_str = argv[2];
	int min_intron_size_int = atoi(min_intron_size_str.c_str());
	log_ofs << "min_intron_size_int: " << min_intron_size_int << endl;
	string max_intron_size_str = argv[3];
	int max_intron_size_int = atoi(max_intron_size_str.c_str());	
	log_ofs << "max_intron_size_int: " << max_intron_size_int << endl;

	cout << "start to initiate 2 alignInferJunctionHashInfo_total ...." << endl;
	log_ofs << "start to initiate 2 alignInferJunctionHashInfo ...." << endl;
	AlignInferJunctionHash_Info* juncHash_total = new AlignInferJunctionHash_Info();
	juncHash_total->initiateAlignInferJunctionInfo(chromNum);
	juncHash_total->insertJuncFromJuncFileVec_chrNamePosOnly(juncFileVec, indexInfo);

	vector<AlignInferJunctionHash_Info*> juncHashVec;
	AlignInferJunctionHash_Info* tmpJuncHash_1st = new AlignInferJunctionHash_Info();
	tmpJuncHash_1st->initiateAlignInferJunctionInfo(chromNum);
	tmpJuncHash_1st->insertJuncFromJuncFile_chrNamePos_supportNum(juncFileVec[0], indexInfo);
	juncHashVec.push_back(tmpJuncHash_1st);
	for(int tmp = 1; tmp < juncFileVec.size(); tmp++)
	{
		string tmpJuncFile = juncFileVec[tmp];
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
		tmpJuncHash->insertJuncFromJuncFile_chrNamePosOnly(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
	}

	int totalJuncNum = juncHash_total->returnAlignInferInfoVecSize();
	int totalJuncNum_valid = 0;
	int totalJuncNum_invalid = 0;
	for(int tmp = 0; tmp < totalJuncNum; tmp++)
	{
		int tmpJunc_chrNameInt = juncHash_total->returnAlignInferInfo_chrNameInt(tmp);
		int tmpJunc_startPos = juncHash_total->returnAlignInferInfo_donerEndPos(tmp);
		int tmpJunc_endPos = juncHash_total->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpJunc_size = tmpJunc_endPos - tmpJunc_startPos - 1;
		if((tmpJunc_size > max_intron_size_int)||(tmpJunc_size < min_intron_size_int))
		{
			invalidJunc_ofs << indexInfo->returnChrNameStr(tmpJunc_chrNameInt) << "\t"
				<< tmpJunc_startPos << "\t" << tmpJunc_endPos << endl;
			totalJuncNum_invalid ++;
			continue;
		}
		totalJuncNum_valid ++;

		int tmpJunc_supportNum = 0;
		bool search_bool = tmpJuncHash_1st->searchAndReturnSupNumInAlignInferJuncHash(
			tmpJunc_supportNum, tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
		if(!search_bool)
			tmpJunc_supportNum = 0;
		juncExistence_ofs << indexInfo->returnChrNameStr(tmpJunc_chrNameInt) << "\t"
			<< tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" << tmpJunc_supportNum;
		for(int tmpJuncHash = 0; tmpJuncHash < juncHashVec.size(); tmpJuncHash ++)
		{
			bool SJexist_tmp_bool = juncHashVec[tmpJuncHash]->SJexistInAlignInferJuncHash(
				tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
			if(SJexist_tmp_bool)
				juncExistence_ofs << "\tY";
			else
				juncExistence_ofs << "\tN"; 
		}
		juncExistence_ofs << endl;;
	}
	cout << "totalJuncNum: " << totalJuncNum << endl;
	cout << "totalJuncNum_valid: " << totalJuncNum_valid << endl;
	cout << "totalJuncNum_invalid: " << totalJuncNum_invalid << endl;

	cout << "All jobs done ! " << endl;
	log_ofs << "All jobs done ! " << endl;	
	delete indexInfo;
	delete juncHash_total;
	for(int tmp = 0; tmp < juncHashVec.size(); tmp++)
		delete juncHashVec[tmp];
	juncExistence_ofs.close();
	//compare_ofs.close();
	parameter_ifs.close();
	log_ofs.close();
}	