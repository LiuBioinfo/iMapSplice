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
	string outputFile_addMax = output_prefix + "/snpAse_addMax.txt";
	ofstream addMax_ofs(outputFile_addMax.c_str());
	string outputFile_cmp2max = output_prefix + "/cmp2max.txt";
	ofstream cmp2max_ofs(outputFile_cmp2max.c_str());
	// string outputFile_cmp2max_aboveThres = output_prefix + "/snpAse_cmp2max_aboveThres.txt";
	// ofstream cmp2max_aboveThres_ofs(outputFile_cmp2max_aboveThres.c_str());
	string outputFile_cmp2i = output_prefix + "/cmp2i.txt";
	ofstream cmp2i_ofs(outputFile_cmp2i.c_str());
	// string outputFile_cmp2i_aboveThres = output_prefix + "/snpAse_cmp2i_aboveThres.txt";
	// ofstream cmp2i_aboveThres_ofs(outputFile_cmp2i_aboveThres.c_str());
	string outputFile_cmp2max_refMaxAboveThres = output_prefix + "/cmp2max_refMaxAboveThres.txt";
	ofstream cmp2max_refMaxAboveThres_ofs(outputFile_cmp2max_refMaxAboveThres.c_str());
	string outputFile_cmp2max_altMaxAboveThres = output_prefix + "/cmp2max_altMaxAboveThres.txt";
	ofstream cmp2max_altMaxAboveThres_ofs(outputFile_cmp2max_altMaxAboveThres.c_str());

	string outputFile_log = output_prefix + "/log";
	ofstream log_ofs(outputFile_log.c_str());

	int totalNumOfResults2compare = 3;
	int targetResultsIDtoCompare = 3;

	int snpNum_total = 0;
	int snpNum_aboveThres_cmp2max = 0;
	int snpNum_aboveThres_cmp2i = 0;
	log_ofs << "chrName chrPos refBase altBase readcountTotal_max refBaseCountTotal_max ";
	log_ofs << "altBaseCountTotal_max readcountTotal_sd refBaseCountTotal_sd altBaseCountTotal_sd";
	log_ofs << "readcountTotal_ps refBaseCountTotal_ps altBaseCountTotal_ps";
	log_ofs << "readcountTotal_i refBaseCountTotal_i altBaseCountTotal_i" << endl;
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
		int tmpReadCount_rf = atoi(tmpReadCountStr_rf.c_str());
		int tmpReadCount_ps = atoi(tmpReadCountStr_ps.c_str());
		int tmpReadCount_i = atoi(tmpReadCountStr_i.c_str());
		int tmpRefBaseReadCount_rf = atoi(tmpRefBaseReadCountStr_rf.c_str());
		int tmpRefBaseReadCount_ps = atoi(tmpRefBaseReadCountStr_ps.c_str());
		int tmpRefBaseReadCount_i = atoi(tmpRefBaseReadCountStr_i.c_str());
		int tmpAltBaseReadCount_rf = atoi(tmpAltBaseReadCountStr_rf.c_str());
		int tmpAltBaseReadCount_ps = atoi(tmpAltBaseReadCountStr_ps.c_str());
		int tmpAltBaseReadCount_i = atoi(tmpAltBaseReadCountStr_i.c_str());
		int tmpReadCount_max = returnMax(tmpReadCount_rf, tmpReadCount_ps, tmpReadCount_i);
		int tmpRefBaseReadCount_max = returnMax(tmpRefBaseReadCount_rf, tmpRefBaseReadCount_ps, tmpRefBaseReadCount_i);
		int tmpAltBaseReadCount_max = returnMax(tmpAltBaseReadCount_rf, tmpAltBaseReadCount_ps, tmpAltBaseReadCount_i);
		addMax_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
			<< tmpReadCount_max << "\t" << tmpRefBaseReadCount_max << "\t" << tmpAltBaseReadCount_max << "\t"
			<< tmpReadCount_rf << "\t" << tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
			<< tmpReadCount_ps << "\t" << tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t"
			<< tmpReadCount_i << "\t" << tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;
		double tmpReadCount_rf_cmp2max = (double)tmpReadCount_rf / (double)tmpReadCount_max;
		double tmpReadCount_rf_cmp2i = (double)tmpReadCount_rf / (double)tmpReadCount_i;
		double tmpReadCount_ps_cmp2max = (double)tmpReadCount_ps / (double)tmpReadCount_max;
		double tmpReadCount_ps_cmp2i = (double)tmpReadCount_ps / (double)tmpReadCount_i;
		double tmpReadCount_i_cmp2max = (double)tmpReadCount_i / (double)tmpReadCount_max;
		double tmpReadCount_i_cmp2i = 1;
		double tmpRefBaseReadCount_rf_cmp2max = (double)tmpRefBaseReadCount_rf / (double)tmpRefBaseReadCount_max;
		double tmpRefBaseReadCount_rf_cmp2i = (double)tmpRefBaseReadCount_rf / (double)tmpRefBaseReadCount_i;
		double tmpRefBaseReadCount_ps_cmp2max = (double)tmpRefBaseReadCount_ps / (double)tmpRefBaseReadCount_max;
		double tmpRefBaseReadCount_ps_cmp2i = (double)tmpRefBaseReadCount_ps / (double)tmpRefBaseReadCount_i;
		double tmpRefBaseReadCount_i_cmp2max = (double)tmpRefBaseReadCount_i / (double)tmpRefBaseReadCount_max;
		double tmpRefBaseReadCount_i_cmp2i = 1;
		double tmpAltBaseReadCount_rf_cmp2max = (double)tmpAltBaseReadCount_rf / (double)tmpAltBaseReadCount_max;
		double tmpAltBaseReadCount_rf_cmp2i = (double)tmpAltBaseReadCount_rf / (double)tmpAltBaseReadCount_i;
		double tmpAltBaseReadCount_ps_cmp2max = (double)tmpAltBaseReadCount_ps / (double)tmpAltBaseReadCount_max;
		double tmpAltBaseReadCount_ps_cmp2i = (double)tmpAltBaseReadCount_ps / (double)tmpAltBaseReadCount_i;
		double tmpAltBaseReadCount_i_cmp2max = (double)tmpAltBaseReadCount_i / (double)tmpAltBaseReadCount_max;
		double tmpAltBaseReadCount_i_cmp2i = 1;
		cmp2max_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
			<< tmpRefBaseReadCount_rf_cmp2max << "\t" << tmpAltBaseReadCount_rf_cmp2max << "\t"
			<< tmpRefBaseReadCount_ps_cmp2max << "\t" << tmpAltBaseReadCount_ps_cmp2max << "\t" 
			<< tmpRefBaseReadCount_i_cmp2max << "\t" << tmpAltBaseReadCount_i_cmp2max << "\t" 
			<< tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
			<< tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t"
			<< tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;
		cmp2i_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
			<< tmpRefBaseReadCount_rf_cmp2i << "\t" << tmpAltBaseReadCount_rf_cmp2i << "\t"
			<< tmpRefBaseReadCount_ps_cmp2i << "\t" << tmpAltBaseReadCount_ps_cmp2i << "\t" 
			<< tmpRefBaseReadCount_i_cmp2i << "\t" << tmpAltBaseReadCount_i_cmp2i 
			<< tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
			<< tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t"
			<< tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;
		if(tmpRefBaseReadCount_max >= supNumThres)
			cmp2max_refMaxAboveThres_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
				<< tmpRefBaseReadCount_rf_cmp2max << "\t" << tmpAltBaseReadCount_rf_cmp2max << "\t" << tmpRefBaseReadCount_ps_cmp2max << "\t" << tmpAltBaseReadCount_ps_cmp2max << "\t" 
				<< tmpRefBaseReadCount_i_cmp2max << "\t" << tmpAltBaseReadCount_i_cmp2max << "\t" << tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
				<< tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t" << tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;
		if(tmpAltBaseReadCount_max >= supNumThres)
			cmp2max_altMaxAboveThres_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
				<< tmpRefBaseReadCount_rf_cmp2max << "\t" << tmpAltBaseReadCount_rf_cmp2max << "\t" << tmpRefBaseReadCount_ps_cmp2max << "\t" << tmpAltBaseReadCount_ps_cmp2max << "\t" 
				<< tmpRefBaseReadCount_i_cmp2max << "\t" << tmpAltBaseReadCount_i_cmp2max << "\t" << tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
				<< tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t" << tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;			
		// if((tmpRefBaseReadCount_rf >= supNumThres)&&(tmpRefBaseReadCount_ps >= supNumThres)&&(tmpRefBaseReadCount_i >= supNumThres)
		// 	&&(tmpAltBaseReadCount_rf >= supNumThres)&&(tmpAltBaseReadCount_ps >= supNumThres)&&(tmpAltBaseReadCount_i >= supNumThres))
		// {
		// 	snpNum_aboveThres_cmp2max ++;
		// 	cmp2max_aboveThres_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
		// 		<< tmpRefBaseReadCount_rf_cmp2max << "\t" << tmpAltBaseReadCount_rf_cmp2max << "\t"
		// 		<< tmpRefBaseReadCount_ps_cmp2max << "\t" << tmpAltBaseReadCount_ps_cmp2max << "\t" 
		// 		<< tmpRefBaseReadCount_i_cmp2max << "\t" << tmpAltBaseReadCount_i_cmp2max << "\t" 
		// 		<< tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
		// 		<< tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t"
		// 		<< tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;			
		// 	snpNum_aboveThres_cmp2i ++;
		// 	cmp2i_aboveThres_ofs << tmpChrName << "\t" << tmpChrPos << "\t" << tmpRefBase << "\t" << tmpAltBase << "\t"
		// 		<< tmpRefBaseReadCount_rf_cmp2i << "\t" << tmpAltBaseReadCount_rf_cmp2i << "\t"
		// 		<< tmpRefBaseReadCount_ps_cmp2i << "\t" << tmpAltBaseReadCount_ps_cmp2i << "\t" 
		// 		<< tmpRefBaseReadCount_i_cmp2i << "\t" << tmpAltBaseReadCount_i_cmp2i 
		// 		<< tmpRefBaseReadCount_rf << "\t" << tmpAltBaseReadCount_rf << "\t"
		// 		<< tmpRefBaseReadCount_ps << "\t" << tmpAltBaseReadCount_ps << "\t"
		// 		<< tmpRefBaseReadCount_i << "\t" << tmpAltBaseReadCount_i << endl;	
		// }




		snpNum_total ++;
	}
	log_ofs << "snpNum_total: " << snpNum_total << endl;
	// log_ofs << "snpNum_aboveThres_cmp2max: " << snpNum_aboveThres_cmp2max << endl;
	// log_ofs << "snpNum_aboveThres_cmp2i: " << snpNum_aboveThres_cmp2i << endl;
	addMax_ofs.close();
	snpASE_ifs.close();
	log_ofs.close();
	return 0;
}