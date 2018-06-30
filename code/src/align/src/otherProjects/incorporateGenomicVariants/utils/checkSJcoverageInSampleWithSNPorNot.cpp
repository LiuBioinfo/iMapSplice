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
//#include "../../../general/splice_info.h"
//#include "../../../general/transcript_set.h"
#include "../../../general/alignInferJunctionHash_info.h"
#include "../general/SNPhash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc <= 8)
	{
		cout << "Executable inputIndexFolder inputMergedSJ outputFolder populationSize SNP_1 SNP_2 (SNP_3 ... SNP_n)";
		cout << " SJ_1 SJ_2 (SJ_3 ... SJ_n)" << endl;
		exit(1);
	}
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "InputIndexFolderPath: " << argv[1] << endl;
	log_ofs << "inputMergedSJ: " << argv[2] << endl;
	log_ofs << "outputFolder: " << argv[3] << endl;
	log_ofs << "populationSize: " << argv[4] << endl;
	int SNPfileNum = (argc - 5)/2;
	int juncFileNum = (argc - 5)/2;
	log_ofs << "SNPfileNum == juncFileNum: " << SNPfileNum << endl;
	if(SNPfileNum + juncFileNum != argc - 5)
	{
		cout << "SNPfileNum + juncFileNum != argc - 5" << endl;
		exit(1);
	}
	string populationSizeStr = argv[4];
	int populationSize = atoi(populationSizeStr.c_str());
	if(populationSize != SNPfileNum)
	{
		cout << "populationSize != SNPfileNum " << endl;
		cout << "populationSize: " << populationSize << endl;
		cout << "SNPfileNum: " << SNPfileNum << endl;
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

	cout << "loading SNP file paths......" << endl;
	log_ofs << "loading SNP file paths......" << endl;
	vector<string> SNPfilePathVec;
	for(int tmp = 0; tmp < populationSize; tmp++)
	{
		string tmpSNPfile = argv[5 + tmp];
		log_ofs << "SNPfile[" << tmp + 1 << "]:" << tmpSNPfile << endl;
		SNPfilePathVec.push_back(tmpSNPfile);
	}

	vector<string> juncFilePathVec;
	for(int tmp = 0; tmp < populationSize; tmp++)
	{
		string tmpJuncFile = argv[5 + tmp + populationSize];
		log_ofs << "juncFile[" << tmp+1 << "]:" << tmpJuncFile << endl;
		juncFilePathVec.push_back(tmpJuncFile);
	}	

	vector<SNPhash_Info*> snpHashVec;
	for(int tmp = 0; tmp < populationSize; tmp++)
	{
		string tmpSNPfile = SNPfilePathVec[tmp];
		cout << "start to do initiate_SNPinfoIndexMapVec_SNPareaMapVec: " << tmp + 1 << endl;
		SNPhash_Info* tmpSNPhashInfo = new SNPhash_Info();
		tmpSNPhashInfo->initiate_SNPinfoIndexMapVec_SNPareaMapVec(chromNum);
		tmpSNPhashInfo->generateSNPhash_formattedSNPfile(tmpSNPfile, indexInfo);
		snpHashVec.push_back(tmpSNPhashInfo);
	}

	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 0; tmp < populationSize; tmp++)
	{
		string tmpJuncFile = juncFilePathVec[tmp];
		cout << "start to do initiate_juncHash: " << tmp + 1 << endl;
		AlignInferJunctionHash_Info* tmpJuncHash = new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum(tmpJuncFile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
	}	

	string outputFile = outputFolderStr + "SJcoverageInSampleWithSNPorNot.txt";
	ofstream SJcov_ofs(outputFile.c_str());
	string inputMergedSJ = argv[2];
	ifstream mergedSJ_ifs(inputMergedSJ.c_str());
	while(!mergedSJ_ifs.eof())
	{
		string tmpStr;
		getline(mergedSJ_ifs, tmpStr);
		if(tmpStr == "")
			break;
		int fieldNum = 8 + 2*populationSize;
		int startLoc = 0;
		vector<string> fieldStrVec;
		for(int tmp = 0; tmp < fieldNum - 1; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpField = tmpStr.substr(startLoc, tabLoc - startLoc);
			fieldStrVec.push_back(tmpField);
			startLoc = tabLoc + 1;
		}
		string tmpLastFieldStr = tmpStr.substr(startLoc);
		fieldStrVec.push_back(tmpLastFieldStr);

		string tmpSJ_chrNameStr = fieldStrVec[0];
		string tmpSJ_startPosStr = fieldStrVec[1];
		string tmpSJ_endPosStr = fieldStrVec[2];
		//string tmpSJ_supNumStr = fieldStrVec[7]
		int tmpSJ_chrNameInt = indexInfo->convertStringToInt(tmpSJ_chrNameStr);
		int tmpSJ_startPos = atoi(tmpSJ_startPosStr.c_str());
		int tmpSJ_endPos = atoi(tmpSJ_endPosStr.c_str());
		//cout << "tmpSJ_chrNameInt: " << tmpSJ_chrNameInt << endl;
		//cout << "tmpSJ_startPos: " << tmpSJ_startPos << endl;
		//int tmpSJ_supNum = atoi(tmpSJ_supNumStr.c_str());
		/*vector<int> SJcovVec;
		for(int tmp = 0; tmp < populationSize; tmp++)
		{
			string tmpSJcovStr = fieldStrVec[9 + tmp*2];
			int tmpSJcov = atoi(tmpSJcovStr.c_str());
			//cout << "tmpSJcov: " << tmpSJcov << endl;
			SJcovVec.push_back(tmpSJcov);
		}*/
		int SNPsampleNum = 0;
		int nonSNPsampleNum = 0;
		int SNPsample_SJcov = 0;
		int nonSNPsample_SJcov = 0;
		string tmpSJcovInEachSampleStr = "";
		for(int tmp = 0; tmp < populationSize; tmp++)
		{
			int candiSNP_1 = tmpSJ_startPos + 1;
			int candiSNP_2 = tmpSJ_startPos + 2;
			int candiSNP_3 = tmpSJ_endPos - 2;
			int candiSNP_4 = tmpSJ_endPos - 1;
			int candiSNP_1_index = snpHashVec[tmp]->searchAndReturnSNPinfoVecIndex(tmpSJ_chrNameInt, candiSNP_1);
			int candiSNP_2_index = snpHashVec[tmp]->searchAndReturnSNPinfoVecIndex(tmpSJ_chrNameInt, candiSNP_2);
			int candiSNP_3_index = snpHashVec[tmp]->searchAndReturnSNPinfoVecIndex(tmpSJ_chrNameInt, candiSNP_3);
			int candiSNP_4_index = snpHashVec[tmp]->searchAndReturnSNPinfoVecIndex(tmpSJ_chrNameInt, candiSNP_4);			
			//int tmpSJ_supNum = SJcovVec[tmp];
			int tmpSJ_supNum = juncHashVec[tmp]->searchAndReturnAlignInferJuncHashSupNum(tmpSJ_chrNameInt,
				tmpSJ_startPos, tmpSJ_endPos);
			//cout << "tmpSJ_supNum: " << tmpSJ_supNum << endl;
			if((candiSNP_1_index != -1)||(candiSNP_2_index != -1)||(candiSNP_3_index != -1)||(candiSNP_4_index != -1))
			{
				//cout << "SNP exists " << endl;
				SNPsampleNum ++;
				SNPsample_SJcov += tmpSJ_supNum;
				tmpSJcovInEachSampleStr += "\tY\t";
				string tmpSJ_supNum_str = int_to_str(tmpSJ_supNum);
				tmpSJcovInEachSampleStr += tmpSJ_supNum_str;
			}	
			else
			{
				//cout << "nonSNP exists " << endl;
				nonSNPsampleNum ++;
				nonSNPsample_SJcov += tmpSJ_supNum;
				tmpSJcovInEachSampleStr += "\tN\t";
				string tmpSJ_supNum_str = int_to_str(tmpSJ_supNum);
				tmpSJcovInEachSampleStr += tmpSJ_supNum_str;
			}
		}
		SJcov_ofs << tmpSJ_chrNameStr << "\t" << tmpSJ_startPosStr << "\t" << tmpSJ_endPosStr << "\t"
			<< fieldStrVec[3] << "\t" << fieldStrVec[4] << "\t" << fieldStrVec[5] << "\t" << fieldStrVec[6];
		SJcov_ofs << tmpSJcovInEachSampleStr;
		SJcov_ofs << "\t" << SNPsampleNum << "\t" << SNPsample_SJcov 
			<< "\t" << nonSNPsampleNum << "\t" << nonSNPsample_SJcov << endl;
	}
	SJcov_ofs.close();
	mergedSJ_ifs.close();
	return 0;
}