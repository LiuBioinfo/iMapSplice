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

#include "../../general/index_info.h"

using namespace std;

bool searchBackSpliceLeftSiteInExonStartPosSetVec(
	vector< set<int> >& exonStartPosSetVec,
	int tmpBackSpliceChrNameInt,
	int tmpBackSpliceSite_left)
{
 	set<int>::iterator tmpSetIter_1 
 		= exonStartPosSetVec[tmpBackSpliceChrNameInt].find(tmpBackSpliceSite_left);
 	if(tmpSetIter_1 != exonStartPosSetVec[tmpBackSpliceChrNameInt].end())	
 		return true;
 	else 
 		return false;
}

bool searchBackSpliceRightSiteInExonEndPosSetVec(
	vector< set<int> >& exonEndPosSetVec,
	int tmpBackSpliceChrNameInt,
	int tmpBackSpliceSite_right)
{
 	set<int>::iterator tmpSetIter_2
 		= exonEndPosSetVec[tmpBackSpliceChrNameInt].find(tmpBackSpliceSite_right);
 	if(tmpSetIter_2 != exonEndPosSetVec[tmpBackSpliceChrNameInt].end())
 		return true;
 	else
 		return false;
}

bool searchBackSpliceSiteInExonEndSiteSetVec(
	vector< set<int> >& exonStartPosSetVec,
	vector< set<int> >& exonEndPosSetVec,
	int tmpBackSpliceChrNameInt,
	int tmpBackSpliceSite_left,
	int tmpBackSpliceSite_right)
{
 	set<int>::iterator tmpSetIter_1 
 		= exonStartPosSetVec[tmpBackSpliceChrNameInt].find(tmpBackSpliceSite_left); 
 	set<int>::iterator tmpSetIter_2
 		= exonEndPosSetVec[tmpBackSpliceChrNameInt].find(tmpBackSpliceSite_right);
 	if((tmpSetIter_1 != exonStartPosSetVec[tmpBackSpliceChrNameInt].end())
 		&&(tmpSetIter_2 != exonEndPosSetVec[tmpBackSpliceChrNameInt].end()))
 		return true;
 	else
 		return false;
}


int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <IndexInput> <inputExonInferFile(gtf)> <toCheckBackSpliceSJfile> <outputFilePrefix>" << endl;
		exit(1);
	}
	cout << "start to initiate indexInfo" << endl;
	string indexFolderPath = argv[1];
	string indexStr = indexFolderPath;
	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputBackSpliceSJ_prefix = outputFolderStr;
	string logFilePath = outputBackSpliceSJ_prefix + "log.txt";
	string outputBackSpliceSJ_annotated = outputBackSpliceSJ_prefix + "annotated.SJ";
	string outputBackSpliceSJ_annotated_oneEnd = outputBackSpliceSJ_prefix + "annotated_oneEnd.SJ";
	string outputBackSpliceSJ_annotated_bothEnds = outputBackSpliceSJ_prefix + "annotated_bothEnds.SJ";
	string outputBackSpliceSJ_unannotated = outputBackSpliceSJ_prefix + "unannotated.SJ";
	ofstream log_ofs(logFilePath.c_str());
	ofstream backSpliceSJ_annotated_ofs(outputBackSpliceSJ_annotated.c_str());
	ofstream backSpliceSJ_annotated_oneEnd_ofs(outputBackSpliceSJ_annotated_oneEnd.c_str());
	ofstream backSpliceSJ_annotated_bothEnds_ofs(outputBackSpliceSJ_annotated_bothEnds.c_str());
	ofstream backSpliceSJ_unannotated_ofs(outputBackSpliceSJ_unannotated.c_str());

	cout << "start to initiate exonStartPosSetVec and exonEndPosSetVec " << endl; 
	log_ofs << "start to initiate exonStartPosSetVec and exonEndPosSetVec " << endl;
	vector< set<int> > exonStartPosSetVec;
	vector< set<int> > exonEndPosSetVec;
	for(int tmpChr = 0; tmpChr < chromNum; tmpChr++)
	{
		set<int> tmpExonStartPosSet;
		set<int> tmpExonEndPosSet;
		exonStartPosSetVec.push_back(tmpExonStartPosSet);
		exonEndPosSetVec.push_back(tmpExonEndPosSet);
	}
	cout << "start to extract exonEnd pos from gtf file" << endl;
	log_ofs << "start to extract exonEnd pos from gtf file" << endl;
	string inputGTFfile = argv[2];
	ifstream gtf_ifs(inputGTFfile.c_str());
	int exonNum = 0;
	while(!gtf_ifs.eof())
	{
		string gtfStr;
		getline(gtf_ifs, gtfStr);
		if((gtf_ifs.eof())||(gtfStr == ""))
			break;
		//cout << "gtfStr: " << endl << gtfStr << endl;
		vector<string> gtfFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 5; tmp++)
		{
			int tabLoc = gtfStr.find("\t", startLoc);
			string tmpGtfField = gtfStr.substr(startLoc, tabLoc-startLoc);
			gtfFieldVec.push_back(tmpGtfField);
			startLoc = tabLoc + 1;
		}
		string gtfChrNameStr = gtfFieldVec[0];
		if(gtfFieldVec[2] == "exon")
		{
			//cout << "exon exists ...." << endl;
			//cout << "gtfChrNameStr: " << gtfChrNameStr << endl;
			int gtfChrNameInt = indexInfo->convertStringToInt(gtfChrNameStr);
			if(gtfChrNameInt < 0)
				continue;
			string exonStartPosStr = gtfFieldVec[3];
			string exonEndPosStr = gtfFieldVec[4];

			int exonStartPosInt = atoi(exonStartPosStr.c_str());
			int exonEndPosInt = atoi(exonEndPosStr.c_str());
			exonStartPosSetVec[gtfChrNameInt].insert(exonStartPosInt);
			exonEndPosSetVec[gtfChrNameInt].insert(exonEndPosInt);
			exonNum ++;
		}
	}
	gtf_ifs.close();
	cout << "exon num: " << exonNum << endl;
	log_ofs << "exon num: " << exonNum << endl;
	cout << "start to check BackSplice SJ annotatedOrNot ..." << endl;
	log_ofs << "start to check BackSplice SJ annotatedOrNot ..." << endl;
	int annotated_oneEnd_backSpliceSJ = 0;
	int annotated_bothEnds_backSpliceSJ = 0;
	int unannotated_backSpliceSJ = 0;
	int total_backSpliceSJ = 0;
	string inputJuncFile = argv[3];
	ifstream junc_ifs(inputJuncFile.c_str());
	while(!(junc_ifs.eof()))
	{	
		string juncStr;
		getline(junc_ifs, juncStr);
		 if(junc_ifs.eof()||(juncStr == ""))
		 	break;

		vector<string> juncFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 3; tmp++)
		{
			int tabLoc = juncStr.find("\t", startLoc);
			string tmpJuncField = juncStr.substr(startLoc, tabLoc-startLoc);
			juncFieldVec.push_back(tmpJuncField);
			startLoc = tabLoc + 1;
		}
		juncFieldVec.push_back(juncStr.substr(startLoc));
		string tmpSJchrName = juncFieldVec[0];
		int tmpSJchrNameInt = indexInfo->convertStringToInt(tmpSJchrName);
		string tmpSJstartPosStr = juncFieldVec[1];
		int tmpSJstartPosInt = atoi(tmpSJstartPosStr.c_str());
		string tmpSJendPosStr = juncFieldVec[2];
		int tmpSJendPosInt = atoi(tmpSJendPosStr.c_str());

		if(tmpSJstartPosInt > tmpSJendPosInt)
		{
			int tmpSJsite_right = tmpSJstartPosInt;
			int tmpSJsite_left = tmpSJendPosInt;
			total_backSpliceSJ ++;
			//bool annotatedOrNot = searchBackSpliceSiteInExonEndSiteSetVec(
			//	tmpSJchrNameInt, tmpSJsite_left, tmpSJsite_right);
			bool annotatedOrNot_leftSite = searchBackSpliceLeftSiteInExonStartPosSetVec(
				exonStartPosSetVec, tmpSJchrNameInt, tmpSJsite_left);
			bool annotatedOrNot_rightSite = searchBackSpliceRightSiteInExonEndPosSetVec(
				exonEndPosSetVec, tmpSJchrNameInt, tmpSJsite_right);

			if(annotatedOrNot_leftSite && annotatedOrNot_rightSite)
			{
				annotated_bothEnds_backSpliceSJ ++;
				backSpliceSJ_annotated_ofs << juncStr << endl;
				backSpliceSJ_annotated_bothEnds_ofs << juncStr << endl;
			}
			else if(annotatedOrNot_leftSite || annotatedOrNot_rightSite)
			{
				annotated_oneEnd_backSpliceSJ ++;
				backSpliceSJ_annotated_ofs << juncStr << endl;
				backSpliceSJ_annotated_oneEnd_ofs << juncStr << endl;
			}
			else
			{
				unannotated_backSpliceSJ ++;
				backSpliceSJ_unannotated_ofs << juncStr << endl;
			}
		}
		else
		{
			cout << "Not backSplice junctions... error !" << endl;
			exit(1);
		}
	}	
	junc_ifs.close();
	cout << "annotated_oneEnd_backSpliceSJ #: " << annotated_oneEnd_backSpliceSJ << endl;
	cout << "annotated_bothEnds_backSpliceSJ #: " << annotated_bothEnds_backSpliceSJ << endl;
	cout << "unannotated_backSpliceSJ #: " << unannotated_backSpliceSJ << endl;
	cout << "total_backSpliceSJ #: " << total_backSpliceSJ << endl;
	log_ofs << "annotated_oneEnd_backSpliceSJ #: " << annotated_oneEnd_backSpliceSJ << endl;
	log_ofs << "annotated_bothEnds_backSpliceSJ #: " << annotated_bothEnds_backSpliceSJ << endl;
	log_ofs << "unannotated_backSpliceSJ #: " << unannotated_backSpliceSJ << endl;
	log_ofs << "total_backSpliceSJ #: " << total_backSpliceSJ << endl;
	delete indexInfo;
	junc_ifs.close();
	backSpliceSJ_annotated_oneEnd_ofs.close();
	backSpliceSJ_annotated_bothEnds_ofs.close();
	backSpliceSJ_unannotated_ofs.close();
	return 0;
}