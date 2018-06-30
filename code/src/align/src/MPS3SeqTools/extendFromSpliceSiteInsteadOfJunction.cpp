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

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"
#include "../general/alignmentToJunc_supportNum.h"
#include "../general/fixSingleAnchorNWDP_info.h"

using namespace std;

string jumpCodeVec2cigarString(vector<Jump_Code>& jumpCodeVec)
{
	string tmpCigarString;
	//cout << "********" << "jumpCodeVecSize: " << jumpCodeVec.size() << endl;
	for(int tmp = 0; tmp < jumpCodeVec.size(); tmp++)
	{
		tmpCigarString += jumpCodeVec[tmp].toString();
	}
	return tmpCigarString;
}

pair<int,int> returnAnchorLengthPair(string& anchorLengthStr)
{
	//cout << "anchorLengthStr: " << anchorLengthStr << endl;
	int firstCommaLoc = anchorLengthStr.find(",");
	//cout << "firstCommaLoc: " << firstCommaLoc << endl;
	string donerAnchorLengthStr = anchorLengthStr.substr(0, firstCommaLoc);
	//cout << "donerAnchorLengthStr: " << donerAnchorLengthStr << endl;
	int donerAnchorLength = atoi(donerAnchorLengthStr.c_str());
	//cout << "donerAnchorLength: " << donerAnchorLength << endl;
	int secondCommaLoc = anchorLengthStr.find(",", firstCommaLoc+1);
	//cout << "secondCommaLoc: " << secondCommaLoc << endl;
	string acceptorAnchorLengthStr = anchorLengthStr.substr(firstCommaLoc+1, secondCommaLoc-firstCommaLoc-1);
	//cout << "acceptorAnchorLengthStr: " << acceptorAnchorLengthStr << endl;
	int acceptorAnchorLength = atoi(acceptorAnchorLengthStr.c_str());
	//cout << "acceptorAnchorLength: " << acceptorAnchorLength << endl;
	return pair<int,int>(donerAnchorLength, acceptorAnchorLength);
}

void extendSJ_AcceptorSite(int chrNameInt, int donerEndPos, int acceptorStartPos,
	Index_Info* indexInfo, int donerAnchorLength, 
	vector<Jump_Code>& extendedSJjumpCodeVec_backward, int& extendedSJpenalty_backward)
{
	string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(chrNameInt, donerEndPos-donerAnchorLength+1, donerAnchorLength);
	string SJextendAtAcceptorSiteStr = indexInfo->returnChromStrSubstr(chrNameInt, acceptorStartPos-donerAnchorLength, donerAnchorLength);
	//cout << "SJdonerAnchorStr: " << SJdonerAnchorStr << endl;
	//cout << "SJextendAtAcceptorSiteStr: " << SJextendAtAcceptorSiteStr << endl;

	FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
	fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(SJdonerAnchorStr, SJextendAtAcceptorSiteStr);
	fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(extendedSJjumpCodeVec_backward);
	extendedSJpenalty_backward = fixSingleAnchorNWDPinfo.getPenalty();
}


void extendSJ_DonerSite(int chrNameInt, int donerEndPos, int acceptorStartPos,
	Index_Info* indexInfo, int acceptorAnchorLength, 
	vector<Jump_Code>& extendedSJjumpCodeVec_forward, int& extendedSJpenalty_forward)
{
	string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(chrNameInt, acceptorStartPos, acceptorAnchorLength);
	string SJextendAtDonerSiteStr = indexInfo->returnChromStrSubstr(chrNameInt, donerEndPos+1, acceptorAnchorLength);
	//cout << "SJacceptorAnchorStr: " << SJacceptorAnchorStr << endl;
	//cout << "SJextendAtDonerSiteStr: " << SJextendAtDonerSiteStr << endl;

	FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
	fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(SJacceptorAnchorStr, SJextendAtDonerSiteStr);
	fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(extendedSJjumpCodeVec_forward);
	extendedSJpenalty_forward = fixSingleAnchorNWDPinfo.getPenalty();
}


void SJextendAndOutput(string& tmpJuncStr, 
	Index_Info* indexInfo, ofstream& SJinfo_ofs, ofstream& SJextension_ofs)
{
	vector<string> juncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 14; tmp++)
	{
		int tabLoc = tmpJuncStr.find("\t", startLoc);
		string tmpJuncField = tmpJuncStr.substr(startLoc, tabLoc-startLoc);
		juncFieldVec.push_back(tmpJuncField);
		startLoc = tabLoc + 1;
	}
	string chrNameStr = juncFieldVec[0];
	int chrNameInt = indexInfo->convertStringToInt(chrNameStr);
	string donerEndPosStr = juncFieldVec[1];
	int donerEndPos = atoi(donerEndPosStr.c_str());
	string acceptorStartPosStr = juncFieldVec[2];
	int acceptorStartPos = atoi(acceptorStartPosStr.c_str());
	string juncSupNumStr = juncFieldVec[4];
	int juncSupNum = atoi(juncSupNumStr.c_str());
	string anchorLengthStr = juncFieldVec[10];
	pair<int,int> anchorLengthPair = returnAnchorLengthPair(anchorLengthStr);
	int donerAnchorLength = anchorLengthPair.first;
	int acceptorAnchorLength = anchorLengthPair.second;
	//cout << "donerAnchorLength: " << donerAnchorLength << endl;
	//cout << "acceptorAnchorLength: " << acceptorAnchorLength << endl;
	string flankStringCaseStr = juncFieldVec[13];
	int flankStringCase = atoi(flankStringCaseStr.c_str());

	int extendedSJpenalty_backward;
	int extendedSJpenalty_forward;
	vector<Jump_Code> extendedSJjumpCodeVec_backward;
	vector<Jump_Code> extendedSJjumpCodeVec_forward;

	// extend backward at acceptorSite
	extendSJ_AcceptorSite(chrNameInt, donerEndPos, acceptorStartPos, indexInfo, donerAnchorLength,
		extendedSJjumpCodeVec_backward, extendedSJpenalty_backward);
	// extend forward at donerSite
	extendSJ_DonerSite(chrNameInt, donerEndPos, acceptorStartPos, indexInfo, acceptorAnchorLength,
		extendedSJjumpCodeVec_forward, extendedSJpenalty_forward);
	string extendedSJjumpCodeVec_backward_cigarString = jumpCodeVec2cigarString(extendedSJjumpCodeVec_backward);
	string extendedSJjumpCodeVec_forward_cigarString = jumpCodeVec2cigarString(extendedSJjumpCodeVec_forward);
	//cout << "*******************************" << endl << extendedSJjumpCodeVec_forward_cigarString << endl
	//	<< "***********************************" << endl;
	SJinfo_ofs << chrNameStr << "\t" << donerEndPos << "\t" << acceptorStartPos << "\t"
		<< juncSupNumStr << "\t" << donerAnchorLength << "\t" << acceptorAnchorLength << "\t"
		<< flankStringCaseStr << "\t" 
		<< extendedSJjumpCodeVec_backward_cigarString << "\t" 
		<< extendedSJjumpCodeVec_forward_cigarString << "\t" 
		<< extendedSJpenalty_backward << "\t"
		<< extendedSJpenalty_forward << "\t" << endl;
}



int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable <InputIndexInfo> <inputSJfile> <outputFolderPath>" << endl;
		exit(0);
	}
	string indexFolderPath = argv[1];
	string inputSJpath = argv[2];
	string outputFolder = argv[3];

	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";
	string outputSJinfoFile = outputFolder + "SJinfo.txt";
	string outputSJextensionFile = outputFolder + "SJextension.txt";

	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;

	ofstream SJinfo_ofs(outputSJinfoFile.c_str());
	ofstream SJextension_ofs(outputSJextensionFile.c_str());
	ifstream inputSJ_ifs(inputSJpath.c_str());
	while(1)
	{
		if(inputSJ_ifs.eof())
			break;
		string tmpAlignStr;
		getline(inputSJ_ifs, tmpAlignStr);
		if(inputSJ_ifs.eof())
			break;
		if(tmpAlignStr.substr(0,3) != "chr")
			continue;
		SJextendAndOutput(tmpAlignStr, indexInfo, SJinfo_ofs, SJextension_ofs);
	}

	return 0;
}