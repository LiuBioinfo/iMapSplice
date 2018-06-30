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
	if(argc <= 5)
	{
		cout << "Executable inputIndexFolder outputFolder populationSize SJ_1 SJ_2 (SJ_3 ... SJ_n)" << endl;
		exit(1);
	}
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "outputFolder: " << argv[2] << endl;
	log_ofs << "populationSize: " << argv[3] << endl;
	int SJfileNum = argc - 4;
	log_ofs << "SJfileNum: " << SJfileNum << endl;
	string populationSizeStr = argv[3];
	int populationSize = atoi(populationSizeStr.c_str());
	if(populationSize != SJfileNum)
	{
		cout << "populationSize != SJfileNum " << endl;
		cout << "populationSize: " << populationSize << endl;
		cout << "SJfileNum: " << SJfileNum << endl;
		exit(1);
	}
	cout << "loading indexInfo parameters ......" << endl;
	log_ofs << "loading indexInfo parameters ......" << endl;
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

	cout << "loading SJ file paths......" << endl;
	log_ofs << "loading SJ file paths......" << endl;
	vector<string> SJfilePathVec;
	for(int tmp = 0; tmp < populationSize; tmp++)
	{
		string tmpSJfile = argv[4 + tmp];
		log_ofs << "SJfile[" << tmp + 1 << "]" << tmpSJfile << endl;
		SJfilePathVec.push_back(tmpSJfile);
	}

	cout << "generating juncHash_merged ......" << endl;
	log_ofs << "generating juncHash_merged ......" << endl;	
	AlignInferJunctionHash_Info* juncHash_merged = new AlignInferJunctionHash_Info();
	juncHash_merged->initiateAlignInferJunctionInfo(chromNum);
	juncHash_merged->insertJuncFromJuncFileVec_chrNamePosOnly(SJfilePathVec, indexInfo);

	cout << "generating juncHashVec ......" << endl;
	log_ofs << "generating juncHashVec ......" << endl;	
	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 0; tmp < populationSize; tmp++)
	{
		string tmpJuncFile = SJfilePathVec[tmp];
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum_flankStringChange(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
	}
	int totalJuncNum = juncHash_merged->returnAlignInferInfoVecSize();
	cout << "totalJuncNum: " << totalJuncNum << endl;
	log_ofs << "totalJuncNum: " << totalJuncNum << endl;
	string junc_merged_file = outputFolderStr + "/mergedJunc.txt";
	ofstream junc_merged_ofs(junc_merged_file.c_str());
	for(int tmp = 0; tmp < totalJuncNum; tmp++)
	{	
		int tmpJunc_chrNameInt = juncHash_merged->returnAlignInferInfo_chrNameInt(tmp);
		string tmpJunc_chrNameStr = indexInfo->returnChrNameStr(tmpJunc_chrNameInt);
		int tmpJunc_startPos = juncHash_merged->returnAlignInferInfo_donerEndPos(tmp);
		int tmpJunc_endPos = juncHash_merged->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpJunc_size = tmpJunc_endPos - tmpJunc_startPos - 1;
		if((tmpJunc_size < 50)||(tmpJunc_size > 300000))
			continue;
		string tmpJunc_flankString_ori = juncHash_merged->returnAlignInferInfo_flankString(tmp, indexInfo);
		junc_merged_ofs << tmpJunc_chrNameStr << "\t" << tmpJunc_startPos << "\t" << tmpJunc_endPos << "\t" << tmpJunc_flankString_ori;
		int tmpJunc_totalSupNum = 0;
		int tmpJunc_existNumInSJvec = 0;
		string tmpJuncSupNumInEachSampleStr = "";
		for(int tmpJuncHash = 0; tmpJuncHash < juncHashVec.size(); tmpJuncHash ++)
		{
			int tmpSJsupNum = juncHashVec[tmpJuncHash]->searchAndReturnAlignInferJuncHashSupNum(
				tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
			string tmpSJflankStringChanged;
			if(tmpSJsupNum > 0)
			{
				tmpJunc_existNumInSJvec ++;
				int tmpSJindex = juncHashVec[tmpJuncHash]->searchAndReturnAlignInferInfoVecIndex(tmpJunc_chrNameInt, tmpJunc_startPos, tmpJunc_endPos);
				if(tmpSJindex < 0)
				{
					cout << "error ! SJ does not exist in tmpAlignInferJuncHash" << endl;
					exit(1);
				}
				tmpSJflankStringChanged = juncHashVec[tmpJuncHash]->returnAlignInferInfo_flankStringChanged(tmpSJindex);
			}
			else
				tmpSJflankStringChanged = "NULL";
			tmpJunc_totalSupNum += tmpSJsupNum;
			tmpJuncSupNumInEachSampleStr += "\t";
			tmpJuncSupNumInEachSampleStr += tmpSJflankStringChanged;
			tmpJuncSupNumInEachSampleStr += "\t";
			string tmpSJsupNumStr = int_to_str(tmpSJsupNum);
			tmpJuncSupNumInEachSampleStr += tmpSJsupNumStr;
		}
		junc_merged_ofs << "\t" << tmpJunc_existNumInSJvec << "\t" << tmpJunc_totalSupNum << tmpJuncSupNumInEachSampleStr << endl;
	}
	cout << "Total junc #: " << totalJuncNum << endl;
	log_ofs << "Total junc #: " << totalJuncNum << endl;	
	cout << "All jobs done !" << endl;
	log_ofs << "All jobs done !" << endl;
	junc_merged_ofs.close();
	log_ofs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}