// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to compare two alignment results, especailly the aligner's results and grount truth
// 0 -- exactly the same, perfect
// 1 -- totally different
// 2 -- overlapped, but contradictory in some bases
// 3 -- without contradictory, not totally the same, covered bases are a parent set with another
// 4 -- without contradictory, not totally the same, covered bases are a child set with another
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
//#include <omp.h>
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/splice_info.h"
#include "alignment_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable inputIndexFolderPath inputGroundTruthAlignmentFile inputAlignerResults outputFolderPath readPairNum" << endl;
		exit(1);
	}
	string inputIndexFolderPath = argv[1];
	string inputGroundTruthSamPath = argv[2];
	string inputAlignerResultsPath = argv[3];
	string outputFolderPath = argv[4];
	string readPairNumStr = argv[5];
	int readPairNum = atoi(readPairNumStr.c_str());

	outputFolderPath += "/";
	string mkdir = "mkdir -p " + outputFolderPath;
	system(mkdir.c_str());	
	string outputLogStr = outputFolderPath + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string indexParameterFileStr = inputIndexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);

	string outputPath_perfect = outputFolderPath + "perfect.sam";
	string outputPath_completelyIncorrect = outputFolderPath + "completelyIncorrect.sam";
	string outputPath_overlapButDifferent = outputFolderPath + "overlapButDifferent.sam";
	string outputPath_parent = outputFolderPath + "parent.sam";
	string outputPath_child = outputFolderPath + "child.sam";
	string outputPath_missedSamInFile1 = outputFolderPath + "missed_1.sam";
	string outputPath_missedSamInFile2 = outputFolderPath + "missed_2.sam";
	string outputPath_missedSamInBothFiles = outputFolderPath + "missed_both.sam";

	ofstream perfect_ofs(outputPath_perfect.c_str());
	ofstream completelyIncorrect_ofs(outputPath_completelyIncorrect.c_str());
	ofstream overlapButDifferent_ofs(outputPath_overlapButDifferent.c_str());
	ofstream parent_ofs(outputPath_parent.c_str());
	ofstream child_ofs(outputPath_child.c_str());
	ofstream missedSam_1_ofs(outputPath_missedSamInFile1.c_str());
	ofstream missedSam_2_ofs(outputPath_missedSamInFile2.c_str());
	ofstream missedSam_both_ofs(outputPath_missedSamInBothFiles.c_str());

	vector< Alignment_Info* > groundTruthSamInfoVec_end1;
	vector< bool > groundTruthExistBoolVec_end1;
	vector< Alignment_Info* > groundTruthSamInfoVec_end2;
	vector< bool > groundTruthExistBoolVec_end2;

	vector< Alignment_Info* > alignerSamInfoVec_end1;
	vector< bool > alignerSamExistBoolVec_end1;
	vector< Alignment_Info* > alignerSamInfoVec_end2;
	vector< bool > alignerSamExistBoolVec_end2;

	for(int tmp = 0; tmp < readPairNum; tmp++)
	{
		Alignment_Info* tmpAlignInfo_gt_end1 = new Alignment_Info();
		Alignment_Info* tmpAlignInfo_gt_end2 = new Alignment_Info();
		Alignment_Info* tmpAlignInfo_rs_end1 = new Alignment_Info();
		Alignment_Info* tmpAlignInfo_rs_end2 = new Alignment_Info();
		groundTruthSamInfoVec_end1.push_back(tmpAlignInfo_gt_end1);
		groundTruthSamInfoVec_end2.push_back(tmpAlignInfo_gt_end2);
		alignerSamInfoVec_end1.push_back(tmpAlignInfo_rs_end1);
		alignerSamInfoVec_end2.push_back(tmpAlignInfo_rs_end2);
		groundTruthExistBoolVec_end1.push_back(false);
		groundTruthExistBoolVec_end2.push_back(false);
		alignerSamExistBoolVec_end1.push_back(false);
		alignerSamExistBoolVec_end2.push_back(false);
	}

	// generate groundTruthSamInfoVec from groundTruthAlignmentFile


	// generate alignerSamInfoVec from alignerResultFile



	// start to compare groundTruhtSamInfo with alignerSamInfo
	for(int tmp = 0; tmp < readPairNum; tmp++)
	{
		// compare two samInfo for end_1

		// compare two samInfo for end_2
	
	}


	for(int tmp = 0; tmp < readPairNum; tmp++)
	{
		delete groundTruthSamInfoVec_end1[tmp];
		delete groundTruthSamInfoVec_end2[tmp];
		delete alignerSamInfoVec_end1[tmp];
		delete alignerSamInfoVec_end2[tmp];
	}

	missedSam_both_ofs.close();
	missedSam_1_ofs.close();
	missedSam_2_ofs.close();
	perfect_ofs.close();
	completelyIncorrect_ofs.close();
	overlapButDifferent_ofs.close();
	parent_ofs.close();
	child_ofs.close();
	delete indexInfo;
	parameter_ifs.close();
	log_ofs.close();
	return 0;
}