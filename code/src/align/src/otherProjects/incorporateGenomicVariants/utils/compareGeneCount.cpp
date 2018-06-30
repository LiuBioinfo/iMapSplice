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

using namespace std;

bool geneAboveThresholdBool(int gc_1, int gc_2, bool supNumIn1stData_or_averageSupNum_bool, int thres_supNum)
{
	if(supNumIn1stData_or_averageSupNum_bool)
	{
		if(gc_1 >= thres_supNum)
			return true;
		else
			return false;
	}
	else
	{
		double gc_average = ((double)gc_1 + (double)gc_2)/2;
		if(gc_average >= thres_supNum)
			return true;
		else
			return false;
	}
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable inputGeneCount_1 inputGeneCount_2 outputFolder thres_supNum supNumIn1stData_or_averageSupNum_bool" << endl;
		exit(1);
	}
	int diff_ratio_max = 100;

	string supNumIn1stData_or_averageSupNum_bool_str = argv[5];
	bool supNumIn1stData_or_averageSupNum_bool;
	if((supNumIn1stData_or_averageSupNum_bool_str == "true")||(supNumIn1stData_or_averageSupNum_bool_str == "True")
		||(supNumIn1stData_or_averageSupNum_bool_str == "TRUE"))
		supNumIn1stData_or_averageSupNum_bool = true;
	else if((supNumIn1stData_or_averageSupNum_bool_str == "false")||(supNumIn1stData_or_averageSupNum_bool_str == "False")
		||(supNumIn1stData_or_averageSupNum_bool_str == "FALSE"))
		supNumIn1stData_or_averageSupNum_bool = false;
	else
	{
		cout << "supNumIn1stData_or_averageSupNum_bool_str should be set as true or false" << endl;
		cout << "exiting ......" << endl;
		exit(1);
	}
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	log_ofs << "inputGeneCount_1: " << argv[1] << endl;
	log_ofs << "inputGeneCount_2: " << argv[2] << endl;
	log_ofs << "outputFolder: " << argv[3] << endl;
	log_ofs << "thres_supNum: " << argv[4] << endl;
	log_ofs << "supNumIn1stData_or_averageSupNum_bool: " << argv[5] << endl; 

	string compare_file = outputFolderStr + "compare.txt";
	string geneCountMerged_file = outputFolderStr + "geneCountMerged.txt";
	ofstream compare_ofs(compare_file.c_str());
	ofstream geneCountMerged_ofs(geneCountMerged_file.c_str());

	string thres_supNum_str = argv[4];
	int thres_supNum = atoi(thres_supNum_str.c_str());

	int totalGeneNum = 0;
	int aboveThresholdGeneNum = 0;
	vector<int> geneCountDiffRatioGeneNumVec;
	for(int tmp = 0; tmp <= diff_ratio_max; tmp++)
		geneCountDiffRatioGeneNumVec.push_back(0);

	string inputGeneCount_1 = argv[1];
	string inputGeneCount_2 = argv[2];
	ifstream gc_ifs_1(inputGeneCount_1.c_str());
	ifstream gc_ifs_2(inputGeneCount_2.c_str());
	while(!gc_ifs_1.eof())
	{
		string tmpGeneCountStr_1;
		string tmpGeneCountStr_2;
		getline(gc_ifs_1, tmpGeneCountStr_1);
		getline(gc_ifs_2, tmpGeneCountStr_2);
		if((tmpGeneCountStr_1 == "")||(tmpGeneCountStr_2 == ""))
			break;
		int tabLoc_inFile1 = tmpGeneCountStr_1.find("\t");
		int tabLoc_inFile2 = tmpGeneCountStr_2.find("\t");
		string tmpGeneName_1 = tmpGeneCountStr_1.substr(0, tabLoc_inFile1);
		string tmpGeneName_2 = tmpGeneCountStr_2.substr(0, tabLoc_inFile2);
		if(tmpGeneName_1 != tmpGeneName_2)
		{
			cout << "tmpGeneName_1 != tmpGeneName_2" << endl;
			cout << "exiting ......" << endl;
			break;
		}
		totalGeneNum ++;
		string tmpGcStr_1 = tmpGeneCountStr_1.substr(tabLoc_inFile1 + 1);
		string tmpGcStr_2 = tmpGeneCountStr_2.substr(tabLoc_inFile2 + 1);
		int tmpGc_1 = atoi(tmpGcStr_1.c_str());
		int tmpGc_2 = atoi(tmpGcStr_2.c_str());	
		bool aboveThresholdGeneNum_bool = geneAboveThresholdBool(tmpGc_1, tmpGc_2, 
			supNumIn1stData_or_averageSupNum_bool, thres_supNum);
		if(aboveThresholdGeneNum_bool)
			aboveThresholdGeneNum ++;
		else
			continue;
		int tmpDiff_geneCount = tmpGc_2 - tmpGc_1;
		if(tmpDiff_geneCount < 0)
			tmpDiff_geneCount = 0 - tmpDiff_geneCount;
		int tmpGeneCountDiffRatio;
		if(tmpGc_1 == 0)
		{
			if(tmpGc_2 == 0)
				tmpGeneCountDiffRatio = 0;
			else
				tmpGeneCountDiffRatio = 100;
		}
		else
		{
			double tmpDiffRatio = ((double)tmpDiff_geneCount/(double)tmpGc_1)*100;
			int tmpDiffRatio_int = (int)tmpDiffRatio;
			if(tmpDiffRatio_int >= diff_ratio_max)
				tmpGeneCountDiffRatio = diff_ratio_max;
			else
				tmpGeneCountDiffRatio = tmpDiffRatio;
		}
		geneCountDiffRatioGeneNumVec[tmpGeneCountDiffRatio] ++;
		geneCountMerged_ofs << tmpGeneName_1 << "\t" << tmpGc_1 << "\t" << tmpGc_2 << "\t" 
			<< tmpGeneCountDiffRatio << "%" << endl;
	}
	int tmpGeneCountSum = 0;
	for(int tmp = diff_ratio_max; tmp >= 0; tmp--)
	{	
		tmpGeneCountSum += geneCountDiffRatioGeneNumVec[tmp];
		double perc_in_totalGene = ((double)tmpGeneCountSum/(double)totalGeneNum)*100;
		double perc_in_validGene = ((double)tmpGeneCountSum/(double)aboveThresholdGeneNum)*100;
		compare_ofs << "geneCount_diff_ratio[>=" << tmp << "]\t" << tmpGeneCountSum 
			<< "\t" << perc_in_validGene << "%\t" << perc_in_totalGene << "%" << endl;
	}

	log_ofs << "totalGeneNum: " << totalGeneNum << endl;
	log_ofs << "aboveThresholdGeneNum: " << aboveThresholdGeneNum << endl;
	compare_ofs.close();
	geneCountMerged_ofs.close();
	gc_ifs_1.close();
	gc_ifs_2.close();
	log_ofs.close();
	return 0;
}