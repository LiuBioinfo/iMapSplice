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
#include "../../../general/splice_info.h"
#include "../../../general/transcript_set.h"
#include "../general/SNPhash_info.h"

using namespace std;

//typedef set<int> 

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputSNPfile_1 inputSNPfile_2 outputFolderStr" << endl;
		exit(1);
	}
	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputFile_log = outputFolderStr + "/log";
	ofstream log_ofs(outputFile_log.c_str());

	log_ofs << "Command: \n" << argv[0] << endl << argv[1] << endl << argv[2] << endl << argv[3] << endl << argv[4] << endl;

	cout << "initiate indexInfo ..." << endl;
	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	indexInfo->initiate();
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	///////////////////////  load SNPs //////////////////////////////////
	string inputSNPfile_1 = argv[2];
	string inputSNPfile_2 = argv[3];
	SNPhash_Info tmpSNPhashInfo_1;
	SNPhash_Info tmpSNPhashInfo_2;
	cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec ..." << endl;
	tmpSNPhashInfo_1.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	tmpSNPhashInfo_2.initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
	cout << "start to do generateSNPhash_formattedSNPfile ..." << endl;
	tmpSNPhashInfo_1.generateSNPhash_formattedSNPfile(inputSNPfile_1, indexInfo);
	tmpSNPhashInfo_2.generateSNPhash_formattedSNPfile(inputSNPfile_2, indexInfo);

	string unique_1_SNP_file = outputFolderStr + "unique_1_SNP.txt";
	string unique_2_SNP_file = outputFolderStr + "unique_2_SNP.txt";
	string shared_SNP_file = outputFolderStr + "shared_SNP.txt";
	string sharedPos_difBase_1_SNP_file = outputFolderStr + "sharedPos_difBase_1_SNP.txt";
	string sharedPos_difBase_2_SNP_file = outputFolderStr + "sharedPos_difBase_2_SNP.txt";
	string cmp_file = outputFolderStr + "cmp.txt";
	ofstream unique_1_ofs(unique_1_SNP_file.c_str());
	ofstream unique_2_ofs(unique_2_SNP_file.c_str());
	ofstream shared_SNP_ofs(shared_SNP_file.c_str());
	ofstream sharedPos_difBase_1_SNP_ofs(sharedPos_difBase_1_SNP_file.c_str());
	ofstream sharedPos_difBase_2_SNP_ofs(sharedPos_difBase_2_SNP_file.c_str());
	ofstream cmp_ofs(cmp_file.c_str());

	int total_1_SNP_num = 0;
	int total_2_SNP_num = 0;
	int unique_1_SNP_num = 0;
	int unique_2_SNP_num = 0;
	int shared_SNP_num = 0;
	int sharedPos_difBase_SNP_num = 0;

	cout << "start to cmp SNP_1 to SNP_2 " << endl;
	int SNP_num_1 = tmpSNPhashInfo_1.returnSNPnum();
	total_1_SNP_num = SNP_num_1;
	for(int tmp = 0; tmp < SNP_num_1; tmp++)
	{
		int tmpSNP_chrNameInt = tmpSNPhashInfo_1.returnSNP_chrNameInt(tmp);
		string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
		int tmpSNP_chrPos = tmpSNPhashInfo_1.returnSNP_chrPos(tmp);
		string tmpSNP_alterBase = tmpSNPhashInfo_1.returnSNP_alterBase(tmp);
		string tmpSNP_referBase = tmpSNPhashInfo_1.returnSNP_referBase(tmp);
	
		string tmpSNP_alterBase_in_SNPhash_2;
		bool tmpSNP_found_in_SNPhash_2_bool = tmpSNPhashInfo_2.searchAndReturnSNPbase(
			tmpSNP_chrNameInt, tmpSNP_chrPos, tmpSNP_alterBase_in_SNPhash_2);
		if(tmpSNP_found_in_SNPhash_2_bool)
		{
			if(tmpSNP_alterBase == tmpSNP_alterBase_in_SNPhash_2)
			{
				shared_SNP_num ++;
				shared_SNP_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t" << tmpSNP_referBase << "\t" << tmpSNP_alterBase << endl;
			}
			else
			{
				sharedPos_difBase_SNP_num ++;
				sharedPos_difBase_1_SNP_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t" << tmpSNP_referBase << "\t" << tmpSNP_alterBase << endl;
				sharedPos_difBase_2_SNP_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t" << tmpSNP_referBase << "\t" << tmpSNP_alterBase_in_SNPhash_2 << endl;
			}
		}
		else
		{
			unique_1_SNP_num ++;
			unique_1_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t" << tmpSNP_referBase << "\t" << tmpSNP_alterBase << endl;
		}
	}

	cout << "start to cmp SNP_2 to SNP_1 " << endl;
	int SNP_num_2 = tmpSNPhashInfo_2.returnSNPnum();
	total_2_SNP_num = SNP_num_2;
	for(int tmp = 0; tmp < SNP_num_2; tmp++)
	{
		int tmpSNP_chrNameInt = tmpSNPhashInfo_2.returnSNP_chrNameInt(tmp);
		string tmpSNP_chrNameStr = indexInfo->returnChrNameStr(tmpSNP_chrNameInt);
		int tmpSNP_chrPos = tmpSNPhashInfo_2.returnSNP_chrPos(tmp);
		string tmpSNP_alterBase = tmpSNPhashInfo_2.returnSNP_alterBase(tmp);
		string tmpSNP_referBase = tmpSNPhashInfo_2.returnSNP_referBase(tmp);

		string tmpSNP_alterBase_in_SNPhash_1;
		bool tmpSNP_found_in_SNPhash_1_bool = tmpSNPhashInfo_1.searchAndReturnSNPbase(
			tmpSNP_chrNameInt, tmpSNP_chrPos, tmpSNP_alterBase_in_SNPhash_1);
		if(!tmpSNP_found_in_SNPhash_1_bool)
		{
			unique_2_SNP_num ++;
			unique_2_ofs << tmpSNP_chrNameStr << "\t" << tmpSNP_chrPos << "\t" << tmpSNP_referBase << "\t" << tmpSNP_alterBase << endl;
		}
	}

	cmp_ofs << "total_1_SNP_num: " << total_1_SNP_num << endl;
	cmp_ofs << "total_2_SNP_num: " << total_2_SNP_num << endl << endl;

	cmp_ofs << "shared_SNP_num: " << shared_SNP_num << endl;
	cmp_ofs << "sharedPos_difBase_SNP_num: " << sharedPos_difBase_SNP_num << endl << endl;

	cmp_ofs << "unique_1_SNP_num: " << unique_1_SNP_num << endl;
	cmp_ofs << "unique_2_SNP_num: " << unique_2_SNP_num << endl << endl;

	cout << "all jobs done ..." << endl;
	unique_1_ofs.close();
	unique_2_ofs.close();
	shared_SNP_ofs.close();
	sharedPos_difBase_1_SNP_ofs.close();
	sharedPos_difBase_2_SNP_ofs.close();
	cmp_ofs.close();
	log_ofs.close();
	chrom_bit_file_ifs.close();
	parameter_ifs.close();
	delete indexInfo;
	return 0;
}