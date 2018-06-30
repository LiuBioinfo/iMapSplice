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
		cout << "Executable <IndexInput> <inputExonInferFile(gtf)> <toCheckCirRNAreadCountfile> <outputFilePrefix>" << endl;
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
	string outputBackSpliceSJ_annotated = outputBackSpliceSJ_prefix + "CirRNA_withAnnotatedSJ";
	string outputBackSpliceSJ_annotated_oneEnd = outputBackSpliceSJ_prefix + "CirRNA_withAnnotatedSJ_oneEnd";
	string outputBackSpliceSJ_annotated_bothEnds = outputBackSpliceSJ_prefix + "CirRNA_withAnnotatedSJ_bothEnds";
	string outputBackSpliceSJ_unannotated = outputBackSpliceSJ_prefix + "CirRNA_withUnannotatedSJ";
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
	string inputCirRNAreadCountFile = argv[3];
	ifstream cirRNAreadCount_ifs(inputCirRNAreadCountFile.c_str());
	while(!(cirRNAreadCount_ifs.eof()))
	{	
		string cirRNAreadCountStr;
		getline(cirRNAreadCount_ifs, cirRNAreadCountStr);
		 if(cirRNAreadCount_ifs.eof()||(cirRNAreadCountStr == ""))
		 	break;

		int startLoc = 0;
		int tabLoc = cirRNAreadCountStr.find("\t", startLoc);
		string tmpCirRNAreadCount_chrNameAndPosStr = cirRNAreadCountStr.substr(startLoc, tabLoc-startLoc);
		startLoc = 0;
		int lineLoc = tmpCirRNAreadCount_chrNameAndPosStr.find("_", startLoc);
		string tmpCirRNAreadCount_chrNameStr = tmpCirRNAreadCount_chrNameAndPosStr.substr(0, lineLoc - startLoc);
		startLoc = lineLoc + 1;
		lineLoc = tmpCirRNAreadCount_chrNameAndPosStr.find("_", startLoc);
		string tmpCirRNAreadCount_leftPosStr = tmpCirRNAreadCount_chrNameAndPosStr.substr(startLoc, lineLoc - startLoc);
		string tmpCirRNAreadCount_rightPosStr = tmpCirRNAreadCount_chrNameAndPosStr.substr(lineLoc+1);

		int tmpBackSpliceSJchrNameInt = indexInfo->convertStringToInt(tmpCirRNAreadCount_chrNameStr);
		int tmpBackSpliceSJsite_left = atoi(tmpCirRNAreadCount_leftPosStr.c_str());
		int tmpBackSpliceSJsite_right = atoi(tmpCirRNAreadCount_rightPosStr.c_str());

		total_backSpliceSJ ++;
		//bool annotatedOrNot = searchBackSpliceSiteInExonEndSiteSetVec(
		//	tmpSJchrNameInt, tmpSJsite_left, tmpSJsite_right);
		bool annotatedOrNot_leftSite = searchBackSpliceLeftSiteInExonStartPosSetVec(
			exonStartPosSetVec, tmpBackSpliceSJchrNameInt, tmpBackSpliceSJsite_left);
		bool annotatedOrNot_rightSite = searchBackSpliceRightSiteInExonEndPosSetVec(
			exonEndPosSetVec, tmpBackSpliceSJchrNameInt, tmpBackSpliceSJsite_right);
		cout << "tmpBackSpliceSJsite_left: " << tmpBackSpliceSJsite_left << endl;
		cout << "tmpBackSpliceSJsite_right: " << tmpBackSpliceSJsite_right << endl; 
		if(annotatedOrNot_leftSite && annotatedOrNot_rightSite)
		{
			annotated_bothEnds_backSpliceSJ ++;
			backSpliceSJ_annotated_ofs << cirRNAreadCountStr << endl;
			backSpliceSJ_annotated_bothEnds_ofs << cirRNAreadCountStr << endl;
		}
		else if(annotatedOrNot_leftSite || annotatedOrNot_rightSite)
		{
			annotated_oneEnd_backSpliceSJ ++;
			backSpliceSJ_annotated_ofs << cirRNAreadCountStr << endl;
			backSpliceSJ_annotated_oneEnd_ofs << cirRNAreadCountStr << endl;
		}
		else
		{
			unannotated_backSpliceSJ ++;
			backSpliceSJ_unannotated_ofs << cirRNAreadCountStr << endl;
		}
	}
	cirRNAreadCount_ifs.close();
	cout << "CirRNA with annotated_oneEnd_backSpliceSJ #: " << annotated_oneEnd_backSpliceSJ << endl;
	cout << "CirRNA with annotated_bothEnds_backSpliceSJ #: " << annotated_bothEnds_backSpliceSJ << endl;
	cout << "CirRNA with unannotated_backSpliceSJ #: " << unannotated_backSpliceSJ << endl;
	cout << "CirRNA with total_backSpliceSJ #: " << total_backSpliceSJ << endl;
	log_ofs << "CirRNA with annotated_oneEnd_backSpliceSJ #: " << annotated_oneEnd_backSpliceSJ << endl;
	log_ofs << "CirRNA with annotated_bothEnds_backSpliceSJ #: " << annotated_bothEnds_backSpliceSJ << endl;
	log_ofs << "CirRNA with unannotated_backSpliceSJ #: " << unannotated_backSpliceSJ << endl;
	log_ofs << "CirRNA with total_backSpliceSJ #: " << total_backSpliceSJ << endl;
	delete indexInfo;
	backSpliceSJ_annotated_oneEnd_ofs.close();
	backSpliceSJ_annotated_bothEnds_ofs.close();
	backSpliceSJ_unannotated_ofs.close();
	return 0;
}