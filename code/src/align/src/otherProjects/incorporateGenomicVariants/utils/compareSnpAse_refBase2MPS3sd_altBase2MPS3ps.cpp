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

using namespace std;

int returnMax(int a, int b, int c)
{
	int tmpMax = a;
	if(b > tmpMax)
		tmpMax = b;
	if(c > tmpMax)
		tmpMax = c;
	return tmpMax;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputMergedSnpAseResults outputFolderPath supNumThres" << endl;
		exit(1);
	}
	string supNumThresStr = argv[3];
	int supNumThres = atoi(supNumThresStr.c_str());

	string inputMergedSnpAseResults = argv[1];
	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string output_prefix = outputFolderStr;
	string outputFile_cmp2sd_refBase = output_prefix + "/cmp2sd_refBase.txt";
	ofstream cmp2sd_refBase_ofs(outputFile_cmp2sd_refBase.c_str());
	string outputFile_cmp2ps_altBase = output_prefix + "/cmp2ps_altBase.txt";
	ofstream cmp2ps_altBase_ofs(outputFile_cmp2ps_altBase.c_str());

	string outputFile_cmp2sd_refBase_2plotInR = output_prefix + "/cmp2sd_refBase_2plotInR.txt";
	ofstream cmp2sd_refBase_2plotInR_ofs(outputFile_cmp2sd_refBase_2plotInR.c_str());
	string outputFile_cmp2ps_altBase_2plotInR = output_prefix + "/cmp2ps_altBase_2plotInR.txt";
	ofstream cmp2ps_altBase_2plotInR_ofs(outputFile_cmp2ps_altBase_2plotInR.c_str());


	string outputFile_refBaseReadCountOnSdAboveThres = output_prefix + "/aboveThres_refBaseReadCountOnSd.txt";
	string outputFile_altBaseReadCountOnPsAboveThres = output_prefix + "/aboveThres_altBaseReadCountOnPs.txt";
	ofstream refBaseReadCountOnSd_ofs(outputFile_refBaseReadCountOnSdAboveThres.c_str());
	ofstream altBaseReadCountOnPs_ofs(outputFile_altBaseReadCountOnPsAboveThres.c_str());
	string outputFile_log = output_prefix + "/log";
	ofstream log_ofs(outputFile_log.c_str());

	int totalNumOfResults2compare = 3;
	int targetResultsIDtoCompare = 3;

	int snpNum_total = 0;
	int snpNum_refBaseReadCountOnSdAboveThres = 0;
	int snpNum_altBaseReadCountOnPsAboveThres = 0;
	ifstream snpASE_ifs(inputMergedSnpAseResults.c_str());
	while(!snpASE_ifs.eof())
	{
		string tmpStr;
		getline(snpASE_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> tmpSnpAseFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < totalNumOfResults2compare * 8 - 1; tmp ++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpSnpAseFieldStr = tmpStr.substr(startLoc, tabLoc-startLoc);
			tmpSnpAseFieldVec.push_back(tmpSnpAseFieldStr);
			startLoc = tabLoc + 1;			
		}
		tmpSnpAseFieldVec.push_back(tmpStr.substr(startLoc));
		string tmpChrName = tmpSnpAseFieldVec[0];
		string tmpChrPos = tmpSnpAseFieldVec[1];
		string tmpRefBase = tmpSnpAseFieldVec[2];
		string tmpAltBase = tmpSnpAseFieldVec[3];
		string tmpReadCountStr_rf = tmpSnpAseFieldVec[4];
		string tmpReadCountStr_ps = tmpSnpAseFieldVec[12];
		string tmpReadCountStr_i = tmpSnpAseFieldVec[20];
		string tmpRefBaseReadCountStr_rf = tmpSnpAseFieldVec[5];
		string tmpRefBaseReadCountStr_ps = tmpSnpAseFieldVec[13];
		string tmpRefBaseReadCountStr_i = tmpSnpAseFieldVec[21];
		string tmpAltBaseReadCountStr_rf = tmpSnpAseFieldVec[6];
		string tmpAltBaseReadCountStr_ps = tmpSnpAseFieldVec[14];
		string tmpAltBaseReadCountStr_i = tmpSnpAseFieldVec[22];
		int tmpReadCount_sd = atoi(tmpReadCountStr_rf.c_str());
		int tmpReadCount_ps = atoi(tmpReadCountStr_ps.c_str());
		int tmpReadCount_i = atoi(tmpReadCountStr_i.c_str());
		int tmpRefBaseReadCount_sd = atoi(tmpRefBaseReadCountStr_rf.c_str());
		int tmpRefBaseReadCount_ps = atoi(tmpRefBaseReadCountStr_ps.c_str());
		int tmpRefBaseReadCount_i = atoi(tmpRefBaseReadCountStr_i.c_str());
		int tmpAltBaseReadCount_sd = atoi(tmpAltBaseReadCountStr_rf.c_str());
		int tmpAltBaseReadCount_ps = atoi(tmpAltBaseReadCountStr_ps.c_str());
		int tmpAltBaseReadCount_i = atoi(tmpAltBaseReadCountStr_i.c_str());
		if(tmpRefBaseReadCount_sd >= supNumThres)
		{
			snpNum_refBaseReadCountOnSdAboveThres ++;
			refBaseReadCountOnSd_ofs << tmpStr << endl;
			double refBaseReadCountCmpRatio_MPS3i_MPS3sd = (double)tmpRefBaseReadCount_i / (double)tmpRefBaseReadCount_sd;
			double refBaseReadCountCmpRatio_MPS3ps_MPS3sd = (double)tmpRefBaseReadCount_ps / (double)tmpRefBaseReadCount_sd;
			cmp2sd_refBase_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
				<< refBaseReadCountCmpRatio_MPS3i_MPS3sd << "\t" << refBaseReadCountCmpRatio_MPS3ps_MPS3sd << endl;
			cmp2sd_refBase_2plotInR_ofs << refBaseReadCountCmpRatio_MPS3i_MPS3sd << "\tMPS3-i vs. MPS3-sd" << endl;
			cmp2sd_refBase_2plotInR_ofs << refBaseReadCountCmpRatio_MPS3ps_MPS3sd << "\tMPS3-ps vs. MPS3-sd" << endl;
		}
		if(tmpAltBaseReadCount_ps >= supNumThres)
		{
			snpNum_altBaseReadCountOnPsAboveThres ++;
			altBaseReadCountOnPs_ofs << tmpStr << endl;
			double altBaseReadCountCmpRatio_MPS3i_MPS3ps = (double)tmpAltBaseReadCount_i / (double)tmpAltBaseReadCount_ps;
			double altBaseReadCountCmpRatio_MPS3sd_MPS3ps = (double)tmpAltBaseReadCount_sd / (double)tmpAltBaseReadCount_ps;
			cmp2ps_altBase_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
				<< altBaseReadCountCmpRatio_MPS3i_MPS3ps << "\t" << altBaseReadCountCmpRatio_MPS3sd_MPS3ps << endl;
			cmp2ps_altBase_2plotInR_ofs << altBaseReadCountCmpRatio_MPS3i_MPS3ps << "\tMPS3-i vs. MPS3-ps" << endl;
			cmp2ps_altBase_2plotInR_ofs << altBaseReadCountCmpRatio_MPS3sd_MPS3ps << "\tMPS3-sd vs. MPS3-ps" << endl;			
		}
		snpNum_total ++;
	}
	log_ofs << "snpNum_total: " << snpNum_total << endl;
	log_ofs << "snpNum_refBaseReadCountOnSdAboveThres: " << snpNum_refBaseReadCountOnSdAboveThres << endl;
	log_ofs << "snpNum_altBaseReadCountOnPsAboveThres: " << snpNum_altBaseReadCountOnPsAboveThres << endl;
	cmp2sd_refBase_ofs.close();
	cmp2ps_altBase_ofs.close();
	refBaseReadCountOnSd_ofs.close();
	altBaseReadCountOnPs_ofs.close();
	snpASE_ifs.close();
	log_ofs.close();
	return 0;
}