// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

using namespace std;

bool compare2geneNameVec_trueWhenAtLeastOneEndMatch(vector<string>& geneNameVec_groundTruth, vector<string>& geneNameVec_toCompare)
{
	int geneNameVec_groundTruth_size = geneNameVec_groundTruth.size();
	int geneNameVec_toCompare_size = geneNameVec_toCompare.size();
	for(int tmp = 0; tmp < geneNameVec_groundTruth_size; tmp++)
	{
		//bool tmpGeneName_groundTruth_in_geneNameVec_toCompare = false;
		string tmpGeneName_groundTruth = geneNameVec_groundTruth[tmp];
		for(int tmp2 = 0; tmp2 < geneNameVec_toCompare_size; tmp2++)
		{
			string tmpGeneName_toCompare = geneNameVec_toCompare[tmp2];
			if(tmpGeneName_groundTruth == tmpGeneName_toCompare)
			{
				return true;
				//tmpGeneName_groundTruth_in_geneNameVec_toCompare = true;
				//break;
			}	
		}
		//if(tmpGeneName_groundTruth_in_geneNameVec_toCompare == true)
		//	return true;
	}
	return false;
}

bool compare2geneNameVec_trueWhenBothEndsMatch(vector<string>& geneNameVec_groundTruth, vector<string>& geneNameVec_toCompare)
{
	int geneNameVec_groundTruth_size = geneNameVec_groundTruth.size();
	int geneNameVec_toCompare_size = geneNameVec_toCompare.size();
	for(int tmp = 0; tmp < geneNameVec_groundTruth_size; tmp++)
	{
		bool tmpGeneName_groundTruth_in_geneNameVec_toCompare = false;
		string tmpGeneName_groundTruth = geneNameVec_groundTruth[tmp];
		for(int tmp2 = 0; tmp2 < geneNameVec_toCompare_size; tmp2++)
		{
			string tmpGeneName_toCompare = geneNameVec_toCompare[tmp2];
			if(tmpGeneName_groundTruth == tmpGeneName_toCompare)
			{
				//return true;
				tmpGeneName_groundTruth_in_geneNameVec_toCompare = true;
				break;
			}	
		}
		if(tmpGeneName_groundTruth_in_geneNameVec_toCompare == false)
			return false;
	}
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputGroundTruthGenePairList inputToCompareGenePairList outputFolder" << endl;
		exit(1);
	}
	string inputGroundTruthGenePairList = argv[1];
	string inputToCompareGenePairList = argv[2];
	string outputFolder = argv[3];
	string mkdirCmd = "mkdir " + outputFolder;
	system(mkdirCmd.c_str());

	vector<string> groundTruthGeneNameStrVec;
	vector< vector<string> > groundTruthGeneNameVecVec;
	vector<bool> detectedOrNotBoolVec_trueWhenBothEndsMatch;
	vector<bool> detectedOrNotBoolVec_trueWhenAtLeastOneEndMatch;
	vector<string> toCompareGeneNameStrVec;
	vector< vector<string> > toCompareGeneNameVecVec;
	vector<bool> correctOrNotBoolVec_trueWhenBothEndsMatch;
	vector<bool> correctOrNotBoolVec_trueWhenAtLeastOneEndMatch;

	ifstream fus_gt_ifs(inputGroundTruthGenePairList.c_str());
	ifstream fus_2compare_ifs(inputToCompareGenePairList.c_str());

	string detectedTrueFus_path_bothEndsMatch = outputFolder + "/detected_true_fusion_bothEndsMatch.txt";
	string missedTrueFus_path_bothEndsMatch = outputFolder + "/missed_true_fusion_bothEndsMatch.txt";
	string correctlyDetectedFus_path_bothEndsMatch = outputFolder + "/correctlyDetected_fusion_bothEndsMatch.txt";
	string incorrectlyDetectedFus_path_bothEndsMatch = outputFolder + "/incorrectlyDetected_fusion_bothEndsMatch.txt";

	ofstream detectedTrueFus_trueWhenBothEndsMatch_ofs(detectedTrueFus_path_bothEndsMatch.c_str());
	ofstream missedTrueFus_trueWhenBothEndsMatch_ofs(missedTrueFus_path_bothEndsMatch.c_str());
	ofstream correctlyDetectedFus_trueWhenBothEndsMatch_ofs(correctlyDetectedFus_path_bothEndsMatch.c_str());	
	ofstream incorrectlyDetectedFus_trueWhenBothEndsMatch_ofs(incorrectlyDetectedFus_path_bothEndsMatch.c_str());

	string detectedTrueFus_path_atLeastOneEndMatch = outputFolder + "/detected_true_fusion_atLeastOneEndMatch.txt";
	string missedTrueFus_path_atLeastOneEndMatch = outputFolder + "/missed_true_fusion_atLeastOneEndMatch.txt";
	string correctlyDetectedFus_path_atLeastOneEndMatch = outputFolder + "/correctlyDetected_fusion_atLeastOneEndMatch.txt";
	string incorrectlyDetectedFus_path_atLeastOneEndMatch = outputFolder + "/incorrectlyDetected_fusion_atLeastOneEndMatch.txt";

	ofstream detectedTrueFus_trueWhenAtLeastOneEndMatch_ofs(detectedTrueFus_path_atLeastOneEndMatch.c_str());
	ofstream missedTrueFus_trueWhenAtLeastOneEndMatch_ofs(missedTrueFus_path_atLeastOneEndMatch.c_str());
	ofstream correctlyDetectedFus_trueWhenAtLeastOneEndMatch_ofs(correctlyDetectedFus_path_atLeastOneEndMatch.c_str());
	ofstream incorrectlyDetectedFus_trueWhenAtLeastOneEndMatch_ofs(incorrectlyDetectedFus_path_atLeastOneEndMatch.c_str());

	int gt_num = 0;
	int detected_true_num = 0;
	int missed_true_num = 0;
	int detected_false_num = 0;

	cout << "start to load groundtruth gene list" << endl;
	while(!fus_gt_ifs.eof())
	{
		string tmpFusStr;
		getline(fus_gt_ifs, tmpFusStr);
		if(tmpFusStr == "")
			break;
		gt_num ++;
		//cout << "tmpFusStr" << tmpFusStr << endl;
		groundTruthGeneNameStrVec.push_back(tmpFusStr);
		vector<string> tmpGeneNameVec;
		int startLoc = 0;
		//cout << gt_num - 1 << " tmpGene: ";
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpFusStr.find(",", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpGeneNameField = tmpFusStr.substr(startLoc, tabLoc-startLoc);
			tmpGeneNameVec.push_back(tmpGeneNameField);
			startLoc = tabLoc + 1;
			cout << tmpGeneNameField << ",";
		}
		//cout << endl;
		groundTruthGeneNameVecVec.push_back(tmpGeneNameVec);
		detectedOrNotBoolVec_trueWhenBothEndsMatch.push_back(false);
		detectedOrNotBoolVec_trueWhenAtLeastOneEndMatch.push_back(false);
	}
	cout << "start to load 2compare gene list" << endl;
	int tmp_gt_num = 0;
	while(!fus_2compare_ifs.eof())
	{
		string tmpFusStr_ori;
		getline(fus_2compare_ifs, tmpFusStr_ori);
		if(tmpFusStr_ori == "")
			break;
		tmp_gt_num ++;		
		toCompareGeneNameStrVec.push_back(tmpFusStr_ori);
		vector<string> tmpFusStrOriFieldVec;
		int startLoc = 0;
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpFusStr_ori.find("\t", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpFusField = tmpFusStr_ori.substr(startLoc, tabLoc-startLoc);
			tmpFusStrOriFieldVec.push_back(tmpFusField);
			startLoc = tabLoc + 1;
		}
		tmpFusStrOriFieldVec.push_back(tmpFusStr_ori.substr(startLoc));
		int tmpFusStrOriFieldVecSize = tmpFusStrOriFieldVec.size();

		string tmpFusStr = tmpFusStrOriFieldVec[tmpFusStrOriFieldVecSize-2] 
			+ tmpFusStrOriFieldVec[tmpFusStrOriFieldVecSize-1];
		//cout << "tmpFusStr:" << tmpFusStr << endl;
		vector<string> tmpGeneNameVec;
		startLoc = 0;
		//cout << tmp_gt_num - 1 << " gene: ";
		for(int tmp = 0; ; tmp++)
		{
			int tabLoc = tmpFusStr.find(",", startLoc);
			if(tabLoc == string::npos)
				break;
			string tmpGeneNameField = tmpFusStr.substr(startLoc, tabLoc-startLoc);
			tmpGeneNameVec.push_back(tmpGeneNameField);
			startLoc = tabLoc + 1;
			//cout << tmpGeneNameField << ",";
		}
		//cout << endl;
		toCompareGeneNameVecVec.push_back(tmpGeneNameVec);
		correctOrNotBoolVec_trueWhenBothEndsMatch.push_back(false);
		correctOrNotBoolVec_trueWhenAtLeastOneEndMatch.push_back(false);
	}
	cout << "start to compare" << endl;
	for(int tmp = 0; tmp < toCompareGeneNameVecVec.size(); tmp++)
	{
		//cout << "toCompareGeneNameVecVec: " << tmp << endl;
		for(int tmp2 = 0; tmp2 < groundTruthGeneNameVecVec.size(); tmp2 ++)
		{
			//cout << "groundTruthGeneNameVecVec: " << tmp2 << endl;
			bool matchBool_trueWhenBothEndsMatch = compare2geneNameVec_trueWhenBothEndsMatch(
				groundTruthGeneNameVecVec[tmp2], toCompareGeneNameVecVec[tmp]);
			bool matchBool_trueWhenAtLeastOneEndMatch = compare2geneNameVec_trueWhenAtLeastOneEndMatch(
				groundTruthGeneNameVecVec[tmp2], toCompareGeneNameVecVec[tmp]);			
			//cout << "matchBool_trueWhenBothEndsMatch: " << matchBool_trueWhenBothEndsMatch << endl;
			//cout << "matchBool_trueWhenAtLeastOneEndMatch: " << matchBool_trueWhenAtLeastOneEndMatch << endl;
			if(matchBool_trueWhenBothEndsMatch)
			{
				detectedOrNotBoolVec_trueWhenBothEndsMatch[tmp2] = true;
				correctOrNotBoolVec_trueWhenBothEndsMatch[tmp] = true;
			}
			if(matchBool_trueWhenAtLeastOneEndMatch)
			{
				detectedOrNotBoolVec_trueWhenAtLeastOneEndMatch[tmp2] = true;
				correctOrNotBoolVec_trueWhenAtLeastOneEndMatch[tmp] = true;				
			}
		}
	}
	cout << "start to output resutls ..." << endl;
	for(int tmp = 0; tmp < groundTruthGeneNameVecVec.size(); tmp++)
	{
		bool detectedOrNotBool_trueWhenBothEndsMatch = detectedOrNotBoolVec_trueWhenBothEndsMatch[tmp];
		if(detectedOrNotBool_trueWhenBothEndsMatch)
			detectedTrueFus_trueWhenBothEndsMatch_ofs << groundTruthGeneNameStrVec[tmp] << endl;
		else
			missedTrueFus_trueWhenBothEndsMatch_ofs << groundTruthGeneNameStrVec[tmp] << endl;
		bool detectedOrNotBool_trueWhenAtLeastOneEndMatch = detectedOrNotBoolVec_trueWhenAtLeastOneEndMatch[tmp];
		if(detectedOrNotBool_trueWhenAtLeastOneEndMatch)
			detectedTrueFus_trueWhenAtLeastOneEndMatch_ofs << groundTruthGeneNameStrVec[tmp] << endl;
		else
			missedTrueFus_trueWhenAtLeastOneEndMatch_ofs << groundTruthGeneNameStrVec[tmp] << endl;
	}
	for(int tmp = 0; tmp < toCompareGeneNameVecVec.size(); tmp++)
	{
		bool correctOrNotBool_trueWhenBothEndsMatch = correctOrNotBoolVec_trueWhenBothEndsMatch[tmp];
		if(correctOrNotBool_trueWhenBothEndsMatch)
			correctlyDetectedFus_trueWhenBothEndsMatch_ofs << toCompareGeneNameStrVec[tmp] << endl;
		else
			incorrectlyDetectedFus_trueWhenBothEndsMatch_ofs << toCompareGeneNameStrVec[tmp] << endl;
		bool correctOrNotBool_trueWhenAtLeastOneEndMatch = correctOrNotBoolVec_trueWhenAtLeastOneEndMatch[tmp];
		if(correctOrNotBool_trueWhenAtLeastOneEndMatch)
			correctlyDetectedFus_trueWhenAtLeastOneEndMatch_ofs << toCompareGeneNameStrVec[tmp] << endl;
		else
			incorrectlyDetectedFus_trueWhenAtLeastOneEndMatch_ofs << toCompareGeneNameStrVec[tmp] << endl;
	}


	fus_gt_ifs.close();
	fus_2compare_ifs.close();

	detectedTrueFus_trueWhenAtLeastOneEndMatch_ofs.close();
	missedTrueFus_trueWhenAtLeastOneEndMatch_ofs.close();
	correctlyDetectedFus_trueWhenBothEndsMatch_ofs.close();
	incorrectlyDetectedFus_trueWhenBothEndsMatch_ofs.close();

	detectedTrueFus_trueWhenAtLeastOneEndMatch_ofs.close();
	missedTrueFus_trueWhenAtLeastOneEndMatch_ofs.close();
	correctlyDetectedFus_trueWhenAtLeastOneEndMatch_ofs.close();
	incorrectlyDetectedFus_trueWhenAtLeastOneEndMatch_ofs.close();
	return 0;
}