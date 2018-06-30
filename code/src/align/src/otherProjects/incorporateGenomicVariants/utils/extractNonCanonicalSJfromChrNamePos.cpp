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

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputIndexFolderPath inputSJlistWithChrNamePos outputFolder SJ_size_min SJ_size_max" << endl;
		exit(1);
	}
	string SJ_size_min_str = argv[4];
	string SJ_size_max_str = argv[5];
	int SJ_size_max = atoi(SJ_size_max_str.c_str());
	int SJ_size_min = atoi(SJ_size_min_str.c_str());

	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "inputSJlistWithChrNamePos: " << argv[2] << endl;
	log_ofs << "outputFolder: " << argv[3] << endl;
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

	string total_file = outputFolderStr + "totalValid_junc.txt";
	string canonical_file = outputFolderStr + "canonical_junc.txt";
	string semiNonCanonical_file = outputFolderStr + "semi_non_junc.txt";
	string semiCanonical_file = outputFolderStr + "semi_junc.txt";
	string nonCanonical_file = outputFolderStr + "non_junc.txt";
	ofstream total_ofs(total_file.c_str());
	ofstream canonical_ofs(canonical_file.c_str());
	ofstream semiNonCanonical_ofs(semiNonCanonical_file.c_str());
	ofstream semiCanonical_ofs(semiCanonical_file.c_str());
	ofstream nonCanonical_ofs(nonCanonical_file.c_str());
	string inputSJ_file = argv[2];
	ifstream SJ_ifs(inputSJ_file.c_str());

	int invalid_junc_num = 0;
	int totalValid_junc_num = 0;
	int canonical_junc_num = 0;
	int semi_junc_num = 0;
	int non_junc_num = 0;
	while(!SJ_ifs.eof())
	{
		string tmpStr;
		getline(SJ_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int tabLoc_1 = tmpStr.find("\t");
		int tabLoc_2 = tmpStr.find("\t", tabLoc_1 + 1);
		int tabLoc_3 = tmpStr.find("\t", tabLoc_2 + 1);
		string tmpSJ_chrNameStr = tmpStr.substr(0, tabLoc_1);
		string tmpSJ_startPosStr = tmpStr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
		string tmpSJ_endPosStr;
		if(tabLoc_3 != string::npos)
 			tmpSJ_endPosStr = tmpStr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);
 		else
 			tmpSJ_endPosStr = tmpStr.substr(tabLoc_2 + 1);
 		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrNameStr);
 		int tmpSJ_startPos = atoi(tmpSJ_startPosStr.c_str());
 		int tmpSJ_endPos = atoi(tmpSJ_endPosStr.c_str());
 		int tmpSJ_size = tmpSJ_endPos - tmpSJ_startPos - 1;
 		if((tmpSJ_chrNameInt < 0)||(tmpSJ_size < SJ_size_min)||(tmpSJ_size > SJ_size_max))
 		{
 			invalid_junc_num ++;
 			continue;
 		}
 		totalValid_junc_num ++; 		
 		string tmpSJ_flankString = indexInfo->returnFlankString(tmpSJ_chrNameInt,
 			tmpSJ_startPos, tmpSJ_endPos);
 		total_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t" << tmpSJ_flankString << endl;
 		if((tmpSJ_flankString == "GTAG")||(tmpSJ_flankString == "CTAC"))
 		{
 			canonical_junc_num ++;
 			canonical_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t" << tmpSJ_flankString << endl;
 		}
 		else if((tmpSJ_flankString == "GCAG")||(tmpSJ_flankString == "CTGC")||(tmpSJ_flankString == "ATAC")||(tmpSJ_flankString == "GTAT"))
 		{	
 			semi_junc_num ++;
 			semiCanonical_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t" << tmpSJ_flankString << endl;
 			semiNonCanonical_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t" << tmpSJ_flankString << endl;
 		}
 		else
 		{
 			non_junc_num ++;
 			nonCanonical_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t" << tmpSJ_flankString << endl;
 			semiNonCanonical_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos << "\t" << tmpSJ_flankString << endl;
 		}
 	}
	SJ_ifs.close();
	log_ofs << endl << "invalid_junc_num: " << invalid_junc_num << endl << "totalValid_junc_num: " << totalValid_junc_num << endl
		<< "canonical_junc_num: " << canonical_junc_num << endl << "semi_junc_num: " << semi_junc_num << endl
		<< "non_junc_num: " << non_junc_num << endl;

	log_ofs.close();
	total_ofs.close();
	canonical_ofs.close();
	semiNonCanonical_ofs.close();
	semiCanonical_ofs.close();
	nonCanonical_ofs.close();
	return 0;
}