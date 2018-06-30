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

#include "../../../general/index_info.h"

using namespace std;

bool searchSpliceSiteInExonBoundaryPosSetVec(
	vector< set<int> >& exonBoundaryPosSetVec, int tmpSpliceChrNameInt, int tmpSpliceSite)
{
 	set<int>::iterator tmpSetIter = exonBoundaryPosSetVec[tmpSpliceChrNameInt].find(tmpSpliceSite);
 	if(tmpSetIter != exonBoundaryPosSetVec[tmpSpliceChrNameInt].end())	
 		return true;
 	else 
 		return false;
}


int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <IndexInput> <inputExonInferFile(gtf)> <toCheckSpliceSJfile> <outputFileFolder>" << endl;
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
	string outputSJ_prefix = outputFolderStr;
	string logFilePath = outputSJ_prefix + "log.txt";
	string outputSJ_annotated_oneEnd = outputSJ_prefix + "annotated_oneEnd.SJ";
	string outputSJ_annotated_bothEnds = outputSJ_prefix + "annotated_bothEnds.SJ";
	string outputSJ_unannotated = outputSJ_prefix + "unannotated.SJ";

	ofstream log_ofs(logFilePath.c_str());
	ofstream annotated_oneEnd_ofs(outputSJ_annotated_oneEnd.c_str());
	ofstream annotated_bothEnds_ofs(outputSJ_annotated_bothEnds.c_str());
	ofstream unannotated_ofs(outputSJ_unannotated.c_str());

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
	cout << "start to extract exon boudary pos from gtf file" << endl;
	log_ofs << "start to extract exon boudary pos from gtf file" << endl;
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
	cout << "start to check  SJ annotatedOrNot ..." << endl;
	log_ofs << "start to check  SJ annotatedOrNot ..." << endl;
	int annotated_oneEnd_SJ = 0;
	int annotated_bothEnds_SJ = 0;
	int unannotated_SJ = 0;
	int total_SJ = 0;
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

		total_SJ ++;
		bool annotatedOrNot_spliceDonerSite = searchSpliceSiteInExonBoundaryPosSetVec(
			exonEndPosSetVec, tmpSJchrNameInt, tmpSJstartPosInt);
		bool annotatedOrNot_spliceAcceptorSite = searchSpliceSiteInExonBoundaryPosSetVec(
			exonStartPosSetVec, tmpSJchrNameInt, tmpSJendPosInt);
		if(annotatedOrNot_spliceDonerSite && annotatedOrNot_spliceAcceptorSite)
		{
			annotated_bothEnds_SJ ++;
			annotated_bothEnds_ofs << juncStr << endl;
		}
		else if(annotatedOrNot_spliceDonerSite || annotatedOrNot_spliceAcceptorSite)
		{
			annotated_oneEnd_SJ ++;
			annotated_oneEnd_ofs << juncStr << endl;
		}
		else
		{
			unannotated_SJ ++;
			unannotated_ofs << juncStr << endl;
		}
	}	
	junc_ifs.close();
	cout << "annotated_oneEnd_SJ #: " << annotated_oneEnd_SJ << endl;
	cout << "annotated_bothEnds_SJ #: " << annotated_bothEnds_SJ << endl;
	cout << "unannotated_SJ #: " << unannotated_SJ << endl;
	cout << "total_SJ #: " << total_SJ << endl;
	log_ofs << "annotated_oneEnd_SJ #: " << annotated_oneEnd_SJ << endl;
	log_ofs << "annotated_bothEnds_SJ #: " << annotated_bothEnds_SJ << endl;
	log_ofs << "unannotated_SJ #: " << unannotated_SJ << endl;
	log_ofs << "total_SJ #: " << total_SJ << endl;
	delete indexInfo;
	junc_ifs.close();
	annotated_oneEnd_ofs.close();
	annotated_bothEnds_ofs.close();
	unannotated_ofs.close();
	return 0;
}