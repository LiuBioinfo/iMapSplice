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

void returnSpliceSite(string& tmpJuncStr, int& chrNameInt,
	 int& donerEndPos, int& acceptorStartPos, Index_Info* indexInfo)
{
	vector<string> juncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = tmpJuncStr.find("\t", startLoc);
		string tmpJuncField = tmpJuncStr.substr(startLoc, tabLoc-startLoc);
		juncFieldVec.push_back(tmpJuncField);
		startLoc = tabLoc + 1;
	}
	string chrNameStr = juncFieldVec[0];
	chrNameInt = indexInfo->convertStringToInt(chrNameStr);
	string donerEndPosStr = juncFieldVec[1];
	donerEndPos = atoi(donerEndPosStr.c_str());
	string acceptorStartPosStr = juncFieldVec[2];
	acceptorStartPos = atoi(acceptorStartPosStr.c_str());
}

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

void checkAnchorSimilarity(Index_Info* indexInfo,
	int targetSJchrNameInt,
	int targetSJdonerEndPos, int targetSJacceptorStartPos,
	int tmpDonerAnchorLength, int tmpAcceptorAnchorLength,
	vector<int>& tmpMultiDonerSiteVec, vector<int>& tmpMultiAcceptorSiteVec,
	vector< vector<Jump_Code> >& tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec,
	vector< vector<Jump_Code> >& tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec,
	vector<int>& tmpMultiDonerSiteAnchorNWDPpenaltyVec,
	vector<int>& tmpMultiAcceptorSiteAnchorNWDPpenaltyVec
	)
{
	//cout << "checkAnchorSimilarity starts ...." << endl;
	string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(
		targetSJchrNameInt, targetSJdonerEndPos - tmpDonerAnchorLength + 1, tmpDonerAnchorLength);
	//cout << "SJdonerAnchorStr: " << SJdonerAnchorStr << endl;
	for(int tmp = 0; tmp < tmpMultiDonerSiteVec.size(); tmp++)
	{
		//cout << "tmpDonerSiteVecIndex: " << tmp << endl;
		int tmpNewDonerSite = tmpMultiDonerSiteVec[tmp];
		//cout << "tmpNewDonerSite: " << tmpNewDonerSite << endl;
		string tmpNewDonerAnchorStr = indexInfo->returnChromStrSubstr(
			targetSJchrNameInt, tmpNewDonerSite - tmpDonerAnchorLength + 1, tmpDonerAnchorLength);
		//cout << "targetSJchrNameInt: " << targetSJchrNameInt << endl;

		//cout << "tmpNewDonerAnchorStr: " << tmpNewDonerAnchorStr << endl;
		vector<Jump_Code> tmpNewDonerNWDPjumpCodeVec;
		FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
		fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
			SJdonerAnchorStr, tmpNewDonerAnchorStr);
		fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewDonerNWDPjumpCodeVec);
		int tmpNewDonerNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
		tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewDonerNWDPjumpCodeVec);
		tmpMultiDonerSiteAnchorNWDPpenaltyVec.push_back(tmpNewDonerNWDPpenalty);
	}

	string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(
		targetSJchrNameInt, targetSJacceptorStartPos, tmpAcceptorAnchorLength);
	//cout << "SJacceptorAnchorStr: " << SJacceptorAnchorStr << endl;
	for(int tmp = 0; tmp < tmpMultiAcceptorSiteVec.size(); tmp++)
	{
		int tmpNewAcceptorSite = tmpMultiAcceptorSiteVec[tmp];
		string tmpNewAcceptorAnchorStr = indexInfo->returnChromStrSubstr(
			targetSJchrNameInt, tmpNewAcceptorSite, tmpAcceptorAnchorLength);
		vector<Jump_Code> tmpNewAcceptorNWDPjumpCodeVec;
		FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
		fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
			SJacceptorAnchorStr, tmpNewAcceptorAnchorStr);
		fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewAcceptorNWDPjumpCodeVec);
		int tmpNewAcceptorNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
		tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewAcceptorNWDPjumpCodeVec);
		tmpMultiAcceptorSiteAnchorNWDPpenaltyVec.push_back(tmpNewAcceptorNWDPpenalty);
	}
}


bool classifySJwithSpliceSiteAnchorSimilarity(
	Index_Info* indexInfo,
	int tmpDonerAnchorLength, int tmpAcceptorAnchorLength,
	vector<int>& tmpMultiDonerSiteAnchorNWDPpenaltyVec,
	vector<int>& tmpMultiAcceptorSiteAnchorNWDPpenaltyVec
	)
{
	int donerAnchorNWDPpenalty_max = 1 + tmpDonerAnchorLength/8;
	int acceptorAnchorNWDPpenalty_max = 1 + tmpAcceptorAnchorLength/8;
	for(int tmp = 0; tmp < tmpMultiDonerSiteAnchorNWDPpenaltyVec.size(); tmp++)
	{
		int tmpPenalty = tmpMultiDonerSiteAnchorNWDPpenaltyVec[tmp];
		if(tmpPenalty <= donerAnchorNWDPpenalty_max)
			return true;
	}
	for(int tmp = 0; tmp < tmpMultiAcceptorSiteAnchorNWDPpenaltyVec.size(); tmp++)
	{
		int tmpPenalty = tmpMultiAcceptorSiteAnchorNWDPpenaltyVec[tmp];
		if(tmpPenalty <= acceptorAnchorNWDPpenalty_max)
			return true;
	}
	return false;
}

void returnMultiDonerEndPosVec(
	int tmpChrNameInt, int tmpAcceptorStartPos, int tmpDonerEndPos,
 	vector< map<int, vector<int> > >& tmpAcceptorStartPosMapVec,
	vector<int>& tmpMultiDonerEndPosVec, int offset, int tmpDonerAnchorLength, Index_Info* indexInfo)
{
	tmpMultiDonerEndPosVec.push_back(tmpAcceptorStartPos - 1);
	map<int, vector<int> >::iterator tmpIntMapIter;
	for(int tmpNewAcceptorStartPos = tmpAcceptorStartPos - offset; 
		tmpNewAcceptorStartPos <= tmpAcceptorStartPos + offset; tmpNewAcceptorStartPos ++)
	{
		tmpIntMapIter = tmpAcceptorStartPosMapVec[tmpChrNameInt].find(tmpNewAcceptorStartPos);
		if(tmpIntMapIter != tmpAcceptorStartPosMapVec[tmpChrNameInt].end())
		{
			for(int tmpValVecIndex = 0; tmpValVecIndex <= (tmpIntMapIter->second).size(); tmpValVecIndex++)
			{
				int tmpNewDonerEndPos = (tmpIntMapIter->second)[tmpValVecIndex];
				if((tmpNewDonerEndPos != tmpDonerEndPos)
					&&((tmpNewDonerEndPos-tmpDonerAnchorLength+1) > 0)
					&&(tmpNewDonerEndPos < indexInfo->returnChromLength(tmpChrNameInt) ))
				{	
					tmpMultiDonerEndPosVec.push_back(tmpNewDonerEndPos);
				}
			}
		}
	}
}

void returnMultiAcceptorStartPosVec(
	int tmpChrNameInt, int tmpDonerEndPos, int tmpAcceptorStartPos,
	vector< map<int, vector<int> > >& tmpDonerEndPosMapVec,
	vector<int>& tmpMultiAcceptorStartPosVec, int offset, int tmpAcceptorAnchorLength, Index_Info* indexInfo)
{
	tmpMultiAcceptorStartPosVec.push_back(tmpDonerEndPos + 1);
	map<int, vector<int> >::iterator tmpIntMapIter;
	for(int tmpNewDonerEndPos = tmpDonerEndPos - offset; 
		tmpNewDonerEndPos <= tmpDonerEndPos + offset; tmpNewDonerEndPos ++)
	{
		tmpIntMapIter = tmpDonerEndPosMapVec[tmpChrNameInt].find(tmpNewDonerEndPos);
		if(tmpIntMapIter != tmpDonerEndPosMapVec[tmpChrNameInt].end())
		{
			for(int tmpValVecIndex = 0; tmpValVecIndex <= (tmpIntMapIter->second).size(); tmpValVecIndex++)
			{
				int tmpNewAcceptorStartPos = (tmpIntMapIter->second)[tmpValVecIndex];
				if((tmpNewAcceptorStartPos != tmpAcceptorStartPos)
					&&(tmpNewAcceptorStartPos > 0)
					&&(tmpNewAcceptorStartPos + tmpAcceptorAnchorLength - 1 
						<= indexInfo->returnChromLength(tmpChrNameInt)))
				{	
					tmpMultiAcceptorStartPosVec.push_back(tmpNewAcceptorStartPos);
				}
			}
		}
	}
}

void returnMultiSpliceSiteVecFromSJstr(
	string& tmpJuncStr, Index_Info* indexInfo,
	int& tmpDonerAnchorLength,
	int& tmpAcceptorAnchorLength,
	vector<int>& tmpMultiDonerEndPosVec,
	vector<int>& tmpMultiAcceptorStartPosVec,
	vector< map<int, vector<int> > >& donerEndPosMapVec,
	vector< map<int, vector<int> > >& acceptorStartPosMapVec,
	int offset)
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
	string anchorLengthStr = juncFieldVec[10];
	pair<int,int> anchorLengthPair = returnAnchorLengthPair(anchorLengthStr);

	tmpDonerAnchorLength = anchorLengthPair.first;
	tmpAcceptorAnchorLength = anchorLengthPair.second;
	// cout << "chrNameInt: " << chrNameInt << endl;
	// cout << "donerEndPos: " << donerEndPos << endl;
	// cout << "acceptorStartPos: " << acceptorStartPos << endl;
	// cout << "tmpDonerAnchorLength: " << tmpDonerAnchorLength << endl;
	// cout << "tmpAcceptorAnchorLength: " << tmpAcceptorAnchorLength << endl;
	// cout << "start to returnMultiAcceptorStartPosVec ..." << endl;
	returnMultiAcceptorStartPosVec(
		chrNameInt, donerEndPos, acceptorStartPos,
		donerEndPosMapVec, tmpMultiAcceptorStartPosVec, offset, tmpAcceptorAnchorLength, indexInfo);
	//cout << "start to returnMultiDonerEndPosVec ... " << endl;
	returnMultiDonerEndPosVec(
		chrNameInt, acceptorStartPos, donerEndPos,
		acceptorStartPosMapVec, tmpMultiDonerEndPosVec, offset, tmpDonerAnchorLength, indexInfo);
	//cout << "end of returnMultiSpliceSiteVec ..." << endl;
}


void insert2mapVec(int tmpMapIndex, int tmpKey, int tmpVal, 
	vector< map<int, vector<int> > >& targetMapVec)
{
	map<int, vector<int> >::iterator tmpIntMapIter;
	tmpIntMapIter = targetMapVec[tmpMapIndex].find(tmpKey);
	if(tmpIntMapIter == targetMapVec[tmpMapIndex].end())
	{
		vector<int> tmpValVec;
		tmpValVec.push_back(tmpVal);
		targetMapVec[tmpMapIndex].insert(pair<int, vector<int> >(tmpKey, tmpValVec));
	}
	else // tmpKey found, then check valVec
	{
		int tmpValVecSize = (tmpIntMapIter->second).size();
		for(int tmp = 0; tmp < tmpValVecSize; tmp++)
		{
			if(tmpVal == (tmpIntMapIter->second)[tmp])
				return; 	
		}
		(tmpIntMapIter->second).push_back(tmpVal);
	}
}


void outputMultiSpliceSiteAnchorSimilarityStr(
	int tmpComparedSJchrNameInt, 
	int tmpComparedSJdonerEndPos, int tmpComparedSJacceptorStartPos,
	int tmpDonerAnchorLength, int tmpAcceptorAnchorLength,
	vector<int>& tmpMultiDonerSiteVec, vector<int>& tmpMultiAcceptorSiteVec,
	vector< vector<Jump_Code> >& tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec,
	vector< vector<Jump_Code> >& tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec,
	vector<int>& tmpMultiDonerSiteAnchorNWDPpenaltyVec,
	vector<int>& tmpMultiAcceptorSiteAnchorNWDPpenaltyVec,
	ofstream& output_ofs, Index_Info* indexInfo)
{
	//cout << "outputMultiSpliceSiteAnchorSimilarityStr starts ..." << endl;
	string tmpComparedSJchrNameStr = indexInfo->returnChrNameStr(tmpComparedSJchrNameInt);
	output_ofs << endl << tmpComparedSJchrNameStr 
		<< "\tdonerSite: " << tmpComparedSJdonerEndPos 
		<< "\tacceptorSite: " << tmpComparedSJacceptorStartPos
		<< "\tdonerAnchorLength: " << tmpDonerAnchorLength 
		<< "\tacceptorAnchorLength: " << tmpAcceptorAnchorLength << endl;
	//cout << "Doner splice site anchor similarity " << endl;
	output_ofs << "Doner splice site anchor similarity: " << endl;
	for(int tmp = 0; tmp < tmpMultiDonerSiteVec.size(); tmp++)	
	{
		int tmpNewDonerAnchorSpliceSite = tmpMultiDonerSiteVec[tmp];
		string tmpNewDonerAnchorNWDPcigarString 
			= jumpCodeVec2cigarString(tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec[tmp]);
		int tmpNewDonerAnchorNWDPpenalty = tmpMultiDonerSiteAnchorNWDPpenaltyVec[tmp];
		output_ofs << tmpNewDonerAnchorSpliceSite << "\t" 
			<< tmpNewDonerAnchorNWDPcigarString << "\t"
			<< tmpNewDonerAnchorNWDPpenalty << endl; 
	}
	//cout << "Acceptor splice site anchor similarity " << endl;
	output_ofs << "Acceptor splice site anchor similarity: " << endl;
	for(int tmp = 0; tmp < tmpMultiAcceptorSiteVec.size(); tmp++)	
	{
		//cout << "tmpIndex in tmpMultiDonerSiteVec: " << tmp << endl;
		int tmpNewAcceptorAnchorSpliceSite = tmpMultiAcceptorSiteVec[tmp];
		//cout << "tmpNewAcceptorAnchorSpliceSite: " << tmpNewAcceptorAnchorSpliceSite << endl;
		string tmpNewAcceptorAnchorNWDPcigarString 
			= jumpCodeVec2cigarString(tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec[tmp]);
		//cout << "tmpNewAcceptorAnchorNWDPcigarString: " << tmpNewAcceptorAnchorNWDPcigarString << endl;
		int tmpNewAcceptorAnchorNWDPpenalty = tmpMultiAcceptorSiteAnchorNWDPpenaltyVec[tmp];
		//cout << "tmpNewAcceptorAnchorNWDPpenalty: " << tmpNewAcceptorAnchorNWDPpenalty << endl;
		output_ofs << tmpNewAcceptorAnchorSpliceSite << "\t" 
			<< tmpNewAcceptorAnchorNWDPcigarString << "\t"
			<< tmpNewAcceptorAnchorNWDPpenalty << endl; 
	}
	output_ofs << "*****************************************************" << endl;
}	


int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable indexFolderPath groundTruthSJfilePath toCompareSJfilePath outputFolderPath offset " << endl;
		exit(1);
	}
	string indexFolderPath = argv[1];
	string inputGroundTruthSJfilePath = argv[2];
	string toCompareSJfilePath = argv[3];
	string outputFolder = argv[4];
	string offsetStr = argv[5];
	int offset = atoi(offsetStr.c_str());

	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	outputFolder += "/";	
	string anchorComparison_file = outputFolder + "anchorComparison.txt";
	ofstream anchorComparison_ofs(anchorComparison_file.c_str());
	string keptSJ_anchorComparison_file = outputFolder + "keptSJ_anchorComparison.txt";
	ofstream keptSJ_anchorComparison_ofs(keptSJ_anchorComparison_file.c_str());	
	string filterOutSJ_anchorComparison_file = outputFolder + "filterOutSJ_anchorComparison.txt";
	ofstream filterOutSJ_anchorComparison_ofs(filterOutSJ_anchorComparison_file.c_str());	

	string log_file = outputFolder + "log";
	ofstream log_ofs(log_file.c_str());
	string keptSJ_file = outputFolder + "keptSJ.junc";
	string filterOutSJ_file = outputFolder + "filterOutSJ.junc";
	ofstream keptSJ_ofs(keptSJ_file.c_str());
	ofstream filterOutSJ_ofs(filterOutSJ_file.c_str());


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

	vector< map<int, vector<int> > > donerEndPosMapVec_groundTruth; // <doner, vector<acceptor> >
	vector< map<int, vector<int> > > acceptorStartPosMapVec_groundTruth; // <acceptor, vector<doner> >
	int chromNum = indexInfo->returnChromNum();
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		map<int, vector<int> > newDonerEndPosMap;
		map<int, vector<int> > newAcceptorEndPosMap;
		donerEndPosMapVec_groundTruth.push_back(newDonerEndPosMap);
		acceptorStartPosMapVec_groundTruth.push_back(newAcceptorEndPosMap);
	}

	cout << "start to insert ground truth SJ into posMapVec_groundTruth" << endl;
	int groundTruthSJ_num = 0;
	ifstream GroundTruthSJ_ifs(inputGroundTruthSJfilePath.c_str());
 	while(1)
	{
		if(GroundTruthSJ_ifs.eof())
			break;
		string tmpGroundTruthSJstr;
		getline(GroundTruthSJ_ifs, tmpGroundTruthSJstr);
		if(GroundTruthSJ_ifs.eof())
			break;
		if(tmpGroundTruthSJstr.substr(0, 3) != "chr")
			continue;
		int tmpGroundTruthSJchrNameInt, 
			tmpGroundTruthSJdonerEndPos, tmpGroundTruthSJacceptorStartPos;
		returnSpliceSite(tmpGroundTruthSJstr, tmpGroundTruthSJchrNameInt, 
			tmpGroundTruthSJdonerEndPos, tmpGroundTruthSJacceptorStartPos, indexInfo);
		int tmpJuncSize = tmpGroundTruthSJacceptorStartPos - tmpGroundTruthSJdonerEndPos - 1;
		if((tmpJuncSize < 50)||(tmpJuncSize > 300000))
			continue;
		groundTruthSJ_num ++;
		//cout << "tmpGroundTruthSJchrNameInt: " << tmpGroundTruthSJchrNameInt << endl;
		//cout << "tmpGroundTruthSJacceptorStartPos: " << tmpGroundTruthSJacceptorStartPos << endl;
		//cout << "tmpGroundTruthSJdonerEndPos: " << tmpGroundTruthSJdonerEndPos << endl;
		insert2mapVec(tmpGroundTruthSJchrNameInt, tmpGroundTruthSJdonerEndPos, 
			tmpGroundTruthSJacceptorStartPos, donerEndPosMapVec_groundTruth);
		insert2mapVec(tmpGroundTruthSJchrNameInt, tmpGroundTruthSJacceptorStartPos, 
			tmpGroundTruthSJdonerEndPos, acceptorStartPosMapVec_groundTruth);
	}

	cout << "start to compare splice site anchor string similarity ...." << endl;
	//int toCompareSJ_num = 0;
	int keptSJ_num = 0;
	int filterOutSJ_num = 0;
	ifstream toCompareSJ_ifs(toCompareSJfilePath.c_str());
	while(1)
	{
		if(toCompareSJ_ifs.eof())
			break;
		string tmpComparedSJstr;
		getline(toCompareSJ_ifs, tmpComparedSJstr);
		if(toCompareSJ_ifs.eof())
			break;
		if(tmpComparedSJstr.substr(0, 3) != "chr")
			continue;
		int tmpComparedSJchrNameInt, 
			tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos;
		returnSpliceSite(tmpComparedSJstr, tmpComparedSJchrNameInt, 
			tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos, indexInfo);
		int tmpJunctionSize = tmpComparedSJacceptorStartPos - tmpComparedSJdonerEndPos - 1;
		if((tmpJunctionSize < 50)||(tmpJunctionSize > 300000))
			continue;
		//cout << "tmpComparedSJchrNameInt: " << tmpComparedSJchrNameInt << endl;
		//cout << "tmpComparedSJdonerEndPos: " << tmpComparedSJdonerEndPos << endl;
		//cout << "tmpComparedSJacceptorStartPos: " << tmpComparedSJacceptorStartPos << endl;
		int tmpDonerAnchorLength, tmpAcceptorAnchorLength;
		vector<int> tmpMultiDonerSiteVec, tmpMultiAcceptorSiteVec;
		returnMultiSpliceSiteVecFromSJstr(
			tmpComparedSJstr, indexInfo,
			tmpDonerAnchorLength, tmpAcceptorAnchorLength,
			tmpMultiDonerSiteVec, tmpMultiAcceptorSiteVec,
			donerEndPosMapVec_groundTruth,
			acceptorStartPosMapVec_groundTruth, offset);

		//cout << "start to check anchor similarity ..." << endl;
		vector< vector<Jump_Code> > tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec;
		vector< vector<Jump_Code> > tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec;
		vector<int> tmpMultiDonerSiteAnchorNWDPpenaltyVec;
		vector<int> tmpMultiAcceptorSiteAnchorNWDPpenaltyVec;
		checkAnchorSimilarity(indexInfo,
			tmpComparedSJchrNameInt,
			tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos,
			tmpDonerAnchorLength, tmpAcceptorAnchorLength,
			tmpMultiDonerSiteVec, tmpMultiAcceptorSiteVec,
			tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec,
			tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec,
			tmpMultiDonerSiteAnchorNWDPpenaltyVec,
			tmpMultiAcceptorSiteAnchorNWDPpenaltyVec
			);

		//cout << "start to classify SJ with spliceSiteAnchorSimilarity " << endl;
		bool multiSpliceSiteAnchorSimilar_bool  
			= classifySJwithSpliceSiteAnchorSimilarity(indexInfo,  
				tmpDonerAnchorLength, tmpAcceptorAnchorLength,
				tmpMultiDonerSiteAnchorNWDPpenaltyVec,
				tmpMultiAcceptorSiteAnchorNWDPpenaltyVec
				);
		//cout << "multiSpliceSiteAnchorSimilar_bool: " << multiSpliceSiteAnchorSimilar_bool << endl;
		if(multiSpliceSiteAnchorSimilar_bool)
		{
			filterOutSJ_num ++;
			filterOutSJ_ofs << tmpComparedSJstr << endl;
			//cout << "start to outputMultiSpliceSiteAnchorSimilarityStr ..." << endl;
			outputMultiSpliceSiteAnchorSimilarityStr(
				tmpComparedSJchrNameInt, 
				tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos,
				tmpDonerAnchorLength, tmpAcceptorAnchorLength,
				tmpMultiDonerSiteVec, tmpMultiAcceptorSiteVec,
				tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec,
				tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec,
				tmpMultiDonerSiteAnchorNWDPpenaltyVec,
				tmpMultiAcceptorSiteAnchorNWDPpenaltyVec,
				filterOutSJ_anchorComparison_ofs, indexInfo
				);
			//cout << "end of outputMultiSpliceSiteAnchorSimilarityStr ..." << endl;
		}
		else
		{
			keptSJ_num ++;
			keptSJ_ofs << tmpComparedSJstr << endl;
			outputMultiSpliceSiteAnchorSimilarityStr(
				tmpComparedSJchrNameInt, 
				tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos,
				tmpDonerAnchorLength, tmpAcceptorAnchorLength,
				tmpMultiDonerSiteVec, tmpMultiAcceptorSiteVec,
				tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec,
				tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec,
				tmpMultiDonerSiteAnchorNWDPpenaltyVec,
				tmpMultiAcceptorSiteAnchorNWDPpenaltyVec,
				keptSJ_anchorComparison_ofs, indexInfo
				);
		}
		//cout << "start to output to anchorSimilarityFile ..." << endl;
		outputMultiSpliceSiteAnchorSimilarityStr(
			tmpComparedSJchrNameInt, 
			tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos,
			tmpDonerAnchorLength, tmpAcceptorAnchorLength,
			tmpMultiDonerSiteVec, tmpMultiAcceptorSiteVec,
			tmpMultiDonerSiteAnchorNWDPjumpCodeVecVec,
			tmpMultiAcceptorSiteAnchorNWDPjumpCodeVecVec,
			tmpMultiDonerSiteAnchorNWDPpenaltyVec,
			tmpMultiAcceptorSiteAnchorNWDPpenaltyVec,
			anchorComparison_ofs, indexInfo
			);		
	}	

	int toCompareSJ_num = keptSJ_num + filterOutSJ_num;
	double keptSJ_num_perc = ((double)keptSJ_num / (double)toCompareSJ_num) * 100;
	double filterOutSJ_num_perc = ((double)filterOutSJ_num / (double)toCompareSJ_num) * 100;
	log_ofs << "groundTruthSJ_num: " << groundTruthSJ_num << endl;
	log_ofs << "toCompareSJ_num: " << toCompareSJ_num << endl;
	log_ofs << endl;
	log_ofs << "keptSJ_num: " << keptSJ_num << " -- " << keptSJ_num_perc << "%" << endl;
	log_ofs << "filterOutSJ_num: " << filterOutSJ_num << " -- " << filterOutSJ_num_perc << "%" << endl;
	free(chrom);
	delete indexInfo;
	anchorComparison_ofs.close();
	keptSJ_anchorComparison_ofs.close();
	filterOutSJ_anchorComparison_ofs.close();
	log_ofs.close();
	keptSJ_ofs.close();
	filterOutSJ_ofs.close();
	parameter_ifs.close();
	GroundTruthSJ_ifs.close();
	toCompareSJ_ifs.close();
	return 0;
}
