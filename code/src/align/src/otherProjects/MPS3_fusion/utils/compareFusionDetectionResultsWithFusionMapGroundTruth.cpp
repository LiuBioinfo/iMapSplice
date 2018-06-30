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

#include "../../../general/read_block_test.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputGroundTruth inputFusionDetectionResults outputFolderPath" << endl;
		exit(1);
	}
	string inputGroundTruth = argv[1];
	string inputFusionDetectionResults = argv[2];
	string outputFolderPath = argv[3];
	outputFolderPath += "/";
	string mkdir= "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());

	string detected_correctFusion = outputFolderPath + "detectedFusion_correct.txt";
	string detected_incorrectFusion = outputFolderPath + "detectedFusion_incorrect.txt";
	string groundTruthFusion_detected = outputFolderPath + "groundTruthFusion_detected.txt";
	string groundTruthFusion_missed = outputFolderPath + "groundTruthFusion_missed.txt";
	ofstream detected_correctFusion_ofs(detected_correctFusion.c_str());
	ofstream detected_incorrectFusion_ofs(detected_incorrectFusion.c_str());
	ofstream groundTruthFusion_detected_ofs(groundTruthFusion_detected.c_str());
	ofstream groundTruthFusion_missed_ofs(groundTruthFusion_missed.c_str());
	cout << "start to read ground truth fusionJunc" << endl;
	// read ground truth fusions
	vector<string> chrNameVec_groundTruth_left;
	vector<string> chrNameVec_groundTruth_right;
	vector<int> chrPosVec_groundTruth_left;
	vector<int> chrPosVec_groundTruth_right;
	vector<string> strandVec_groundTruth_left;
	vector<string> strandVec_groundTruth_right;
	vector<string> otherStrVec_groundTruth;

	ifstream groundTruth_ifs(inputGroundTruth.c_str());
	while(!groundTruth_ifs.eof())
	{
		string tmpGroundTruthStr;
		getline(groundTruth_ifs, tmpGroundTruthStr);
		if(tmpGroundTruthStr == "")
			break;
		vector<string> tmpFieldStrVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = tmpGroundTruthStr.find("\t", startLoc);
			string tmpFieldStr = tmpGroundTruthStr.substr(startLoc, tabLoc - startLoc);
			tmpFieldStrVec.push_back(tmpFieldStr);
			startLoc = tabLoc + 1;
		}
		string tmpOtherStr = tmpGroundTruthStr.substr(startLoc);
		otherStrVec_groundTruth.push_back(tmpOtherStr);
		string tmpChrName_left = tmpFieldStrVec[0];
		string tmpChrName_right = tmpFieldStrVec[1];
		string tmpChrPosStr_left = tmpFieldStrVec[2];
		int tmpChrPosInt_left = atoi(tmpChrPosStr_left.c_str());
		string tmpChrPosStr_right = tmpFieldStrVec[3];
		int tmpChrPosInt_right = atoi(tmpChrPosStr_right.c_str());
		string tmpStrand_left = tmpFieldStrVec[4];
		string tmpStrand_right = tmpFieldStrVec[5];
		chrNameVec_groundTruth_left.push_back(tmpChrName_left);
		chrNameVec_groundTruth_right.push_back(tmpChrName_right);
		chrPosVec_groundTruth_left.push_back(tmpChrPosInt_left);
		chrPosVec_groundTruth_right.push_back(tmpChrPosInt_right);
		strandVec_groundTruth_left.push_back(tmpStrand_left);
		strandVec_groundTruth_right.push_back(tmpStrand_right);
	}
	groundTruth_ifs.close();

	cout << "start to read detected fusionJunc" << endl;
	// read detected fusions
	vector<string> chrNameVec_detected_left;
	vector<string> chrNameVec_detected_right;
	vector<int> chrPosVec_detected_left;
	vector<int> chrPosVec_detected_right;
	vector<string> strandVec_detected_left;
	vector<string> strandVec_detected_right;
	vector<string> otherStrVec_detected;

	ifstream detected_ifs(inputFusionDetectionResults.c_str());
	while(!detected_ifs.eof())
	{
		string tmpDetectedStr;
		getline(detected_ifs, tmpDetectedStr);
		if(tmpDetectedStr == "")
			break;
		vector<string> tmpFieldStrVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 6; tmp++)
		{
			int tabLoc = tmpDetectedStr.find("\t", startLoc);
			string tmpFieldStr = tmpDetectedStr.substr(startLoc, tabLoc - startLoc);
			tmpFieldStrVec.push_back(tmpFieldStr);
			startLoc = tabLoc + 1;
		}
		string tmpOtherStr = tmpDetectedStr.substr(startLoc);
		otherStrVec_detected.push_back(tmpOtherStr);
		string tmpChrName_left = tmpFieldStrVec[0];
		string tmpChrName_right = tmpFieldStrVec[1];
		string tmpChrPosStr_left = tmpFieldStrVec[2];
		int tmpChrPosInt_left = atoi(tmpChrPosStr_left.c_str());
		string tmpChrPosStr_right = tmpFieldStrVec[3];
		int tmpChrPosInt_right = atoi(tmpChrPosStr_right.c_str());
		string tmpStrand_left = tmpFieldStrVec[4];
		string tmpStrand_right = tmpFieldStrVec[5];
		chrNameVec_detected_left.push_back(tmpChrName_left);
		chrNameVec_detected_right.push_back(tmpChrName_right);
		chrPosVec_detected_left.push_back(tmpChrPosInt_left);
		chrPosVec_detected_right.push_back(tmpChrPosInt_right);
		strandVec_detected_left.push_back(tmpStrand_left);
		strandVec_detected_right.push_back(tmpStrand_right);
	}
	detected_ifs.close();

	cout << "start to compare detected results to ground truth" << endl;
	// compare detected results to ground truth
	cout << "chrNameVec_detected_left.size(): " << chrNameVec_detected_left.size() << endl;
	for(int tmp = 0; tmp < chrNameVec_detected_left.size(); tmp++)
	{
		string tmpChrName_left = chrNameVec_detected_left[tmp];
		cout << "tmpChrName_left: " << tmpChrName_left << endl;
		string tmpChrName_right = chrNameVec_detected_right[tmp];
		cout << "tmpChrName_right: " << tmpChrName_right << endl;
		int tmpChrPos_left = chrPosVec_detected_left[tmp];
		int tmpChrPos_right = chrPosVec_detected_right[tmp];
		cout << "tmpChrPos_left: " <<  tmpChrPos_left << endl;
		string tmpStrand_left = strandVec_detected_left[tmp];
		string tmpStrand_right = strandVec_detected_right[tmp];
		cout << "tmpStrand_left: " << tmpStrand_left << endl;
		string tmpOtherStr = otherStrVec_detected[tmp];
		cout << "tmpChrPos_left: " << tmpChrPos_left << endl;
		bool tmpStrandedBool = (!((tmpStrand_left == "N")&&(tmpStrand_right == "N")));

		bool detectedInGroundTruth_bool = false;
		cout << "chrNameVec_groundTruth_left.size(): " << chrNameVec_groundTruth_left.size() << endl;
		for(int tmp2 = 0; tmp2 < chrNameVec_groundTruth_left.size(); tmp2++)
		{
			string tmpChrName_groundTruth_left = chrNameVec_groundTruth_left[tmp2];
			string tmpChrName_groundTruth_right = chrNameVec_groundTruth_right[tmp2];
			int tmpChrPos_groundTruth_left = chrPosVec_groundTruth_left[tmp2];
			int tmpChrPos_groundTruth_right = chrPosVec_groundTruth_right[tmp2];
			string tmpStrand_groundTruth_left = strandVec_groundTruth_left[tmp2];
			string tmpStrand_groundTruth_right = strandVec_groundTruth_right[tmp2];
			cout << "tmpChrPos_groundTruth_left: " << tmpChrPos_groundTruth_left << endl; 
			if(tmpStrandedBool)
			{
				if((tmpChrName_left == tmpChrName_groundTruth_left)
					&&(tmpChrName_right == tmpChrName_groundTruth_right)
					&&(tmpStrand_left == tmpStrand_groundTruth_left)
					&&(tmpStrand_right == tmpStrand_groundTruth_right)
					&&((tmpChrPos_left >= tmpChrPos_groundTruth_left - 50)&&(tmpChrPos_left <= tmpChrPos_groundTruth_left + 50))
					&&((tmpChrPos_right >= tmpChrPos_groundTruth_right - 50)&&(tmpChrPos_right <= tmpChrPos_groundTruth_right + 50)))
				{
					detectedInGroundTruth_bool = true;
					break;
				}
			}
			else
			{
				if(((tmpChrName_left == tmpChrName_groundTruth_left)
					&&(tmpChrName_right == tmpChrName_groundTruth_right)
					&&(tmpStrand_left == tmpStrand_groundTruth_left)
					&&(tmpStrand_right == tmpStrand_groundTruth_right)
					&&((tmpChrPos_left >= tmpChrPos_groundTruth_left - 50)&&(tmpChrPos_left <= tmpChrPos_groundTruth_left + 50))
					&&((tmpChrPos_right >= tmpChrPos_groundTruth_right - 50)&&(tmpChrPos_right <= tmpChrPos_groundTruth_right + 50)))
					||
					((tmpChrName_right == tmpChrName_groundTruth_left)
					&&(tmpChrName_left == tmpChrName_groundTruth_right)
					&&(tmpStrand_right == tmpStrand_groundTruth_left)
					&&(tmpStrand_left == tmpStrand_groundTruth_right)
					&&((tmpChrPos_right >= tmpChrPos_groundTruth_left - 50)&&(tmpChrPos_right <= tmpChrPos_groundTruth_left + 50))
					&&((tmpChrPos_left >= tmpChrPos_groundTruth_right - 50)&&(tmpChrPos_left <= tmpChrPos_groundTruth_right + 50))))
				{
					detectedInGroundTruth_bool = true;
					break;
				}
			}
		}
		string tmpDetectedStr = tmpChrName_left + "\t" + tmpChrName_right + "\t" 
			+ int_to_str(tmpChrPos_left) + "\t" + int_to_str(tmpChrPos_right) + "\t"
			+ tmpStrand_left + "\t" + tmpStrand_right + "\t" + tmpOtherStr;
		if(detectedInGroundTruth_bool)
			detected_correctFusion_ofs << tmpDetectedStr << endl;
		else
			detected_incorrectFusion_ofs << tmpDetectedStr << endl;
	}

	cout << "start to compare ground truth to detected fusions" << endl;
	// compare ground truth to detected results
	for(int tmp = 0; tmp < chrNameVec_groundTruth_left.size(); tmp++)
	{
		string tmpChrName_groundTruth_left = chrNameVec_groundTruth_left[tmp];
		string tmpChrName_groundTruth_right = chrNameVec_groundTruth_right[tmp];
		int tmpChrPos_groundTruth_left = chrPosVec_groundTruth_left[tmp];
		int tmpChrPos_groundTruth_right = chrPosVec_groundTruth_right[tmp];
		string tmpStrand_groundTruth_left = strandVec_groundTruth_left[tmp];
		string tmpStrand_groundTruth_right = strandVec_groundTruth_right[tmp];
		string tmpOtherStr_groundTruth = otherStrVec_groundTruth[tmp];
		bool groundTruthDetected_bool = false;
		for(int tmp2 = 0; tmp2 < chrNameVec_detected_left.size(); tmp2++)
		{
			string tmpChrName_left = chrNameVec_detected_left[tmp2];
			string tmpChrName_right = chrNameVec_detected_right[tmp2];
			int tmpChrPos_left = chrPosVec_detected_left[tmp2];
			int tmpChrPos_right = chrPosVec_detected_right[tmp2];
			string tmpStrand_left = strandVec_detected_left[tmp2];
			string tmpStrand_right = strandVec_detected_right[tmp2];

			bool tmpStrandedBool = (!((tmpStrand_left == "N")&&(tmpStrand_right == "N")));
			
			if(tmpStrandedBool)
			{
				if((tmpChrName_left == tmpChrName_groundTruth_left)
					&&(tmpChrName_right == tmpChrName_groundTruth_right)
					&&(tmpStrand_left == tmpStrand_groundTruth_left)
					&&(tmpStrand_right == tmpStrand_groundTruth_right)
					&&((tmpChrPos_left >= tmpChrPos_groundTruth_left - 50)&&(tmpChrPos_left <= tmpChrPos_groundTruth_left + 50))
					&&((tmpChrPos_right >= tmpChrPos_groundTruth_right - 50)&&(tmpChrPos_right <= tmpChrPos_groundTruth_right + 50)))
				{
					groundTruthDetected_bool = true;
					break;
				}
			}
			else
			{
				if(((tmpChrName_left == tmpChrName_groundTruth_left)
					&&(tmpChrName_right == tmpChrName_groundTruth_right)
					&&(tmpStrand_left == tmpStrand_groundTruth_left)
					&&(tmpStrand_right == tmpStrand_groundTruth_right)
					&&((tmpChrPos_left >= tmpChrPos_groundTruth_left - 50)&&(tmpChrPos_left <= tmpChrPos_groundTruth_left + 50))
					&&((tmpChrPos_right >= tmpChrPos_groundTruth_right - 50)&&(tmpChrPos_right <= tmpChrPos_groundTruth_right + 50)))
					||
					((tmpChrName_right == tmpChrName_groundTruth_left)
					&&(tmpChrName_left == tmpChrName_groundTruth_right)
					&&(tmpStrand_right == tmpStrand_groundTruth_left)
					&&(tmpStrand_left == tmpStrand_groundTruth_right)
					&&((tmpChrPos_right >= tmpChrPos_groundTruth_left - 50)&&(tmpChrPos_right <= tmpChrPos_groundTruth_left + 50))
					&&((tmpChrPos_left >= tmpChrPos_groundTruth_right - 50)&&(tmpChrPos_left <= tmpChrPos_groundTruth_right + 50))))
				{
					groundTruthDetected_bool = true;
					break;
				}
			}
		}
		string tmpGroundTruthStr = tmpChrName_groundTruth_left + "\t" + tmpChrName_groundTruth_right + "\t" 
			+ int_to_str(tmpChrPos_groundTruth_left) + "\t" + int_to_str(tmpChrPos_groundTruth_right) + "\t"
			+ tmpStrand_groundTruth_left + "\t" + tmpStrand_groundTruth_right + "\t" + tmpOtherStr_groundTruth;
		if(groundTruthDetected_bool)
			groundTruthFusion_detected_ofs << tmpGroundTruthStr << endl;
		else
			groundTruthFusion_missed_ofs << tmpGroundTruthStr << endl;
	}

	detected_correctFusion_ofs.close();
	detected_incorrectFusion_ofs.close();
	groundTruthFusion_detected_ofs.close();
	groundTruthFusion_missed_ofs.close();
	return 0;
}