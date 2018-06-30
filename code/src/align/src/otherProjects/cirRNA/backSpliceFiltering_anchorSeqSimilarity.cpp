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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"

using namespace std;

int minAmongThreeNum(int a, int b, int c)
{
	if(a <= b)
	{
		if(a <= c)
			return a;
		else
			return c;
	}
	else // b < a
	{
		if(b <= c)
			return b;
		else
			return c;
	}
}

void checkAlterSpliceAnchorSeqSimilarity_backSpliceJunc(
	SJhash_Info* SJ, AlignInferJunctionHash_Info* alignInferJuncHashInfo,
	int tmpChrNameInt, int tmpBackSpliceStartPos, int tmpBackSpliceEndPos,
	int tmpAnchorSize_doner, int tmpAnchorSize_acceptor, 
	Index_Info* indexInfo, int toCheckAnchorLengthMax,	
	vector<int>& validAnchorLengthVec_alterAcceptor, vector<int>& validAnchorLengthVec_alterDoner,
	vector<int>& penaltyVec_alterAcceptor, vector<int>& penaltyVec_alterDoner,
	vector<int>& alterSiteVec_alterAcceptor, vector<int>& alterSiteVec_alterDoner)
{
	// check alternative acceptor start pos vec
	vector<int> tmpAlterAcceptorSpliceSitePosVec;
	SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChrNameInt,
		tmpBackSpliceStartPos, tmpAlterAcceptorSpliceSitePosVec);
	for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
		if(tmpAlterAcceptorSpliceSitePos != tmpBackSpliceEndPos)
		{
			int tmpIndexInAlignInferJuncVec
				= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
					tmpChrNameInt, tmpBackSpliceStartPos, tmpAlterAcceptorSpliceSitePos);
			if(tmpIndexInAlignInferJuncVec < 0)
			{
				cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
				exit(1);
			}
			else
			{
				int tmpAlterSplice_anchorLength 
					= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(tmpIndexInAlignInferJuncVec);
				int tmpValidAnchorLength_acceptor = minAmongThreeNum(
					tmpAnchorSize_acceptor, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);
				
				string tmpAcceptorSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt,
					tmpBackSpliceEndPos, tmpValidAnchorLength_acceptor);
				string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt,
					tmpAlterAcceptorSpliceSitePos, tmpValidAnchorLength_acceptor);
				
				FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
				tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(
					tmpAcceptorSeq, tmpAlterAcceptorSeq);
				int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
				
				validAnchorLengthVec_alterAcceptor.push_back(tmpValidAnchorLength_acceptor);
				penaltyVec_alterAcceptor.push_back(tmpPenalty_alterAcceptor);
				alterSiteVec_alterAcceptor.push_back(tmpAlterAcceptorSpliceSitePos);
			}
		}
	}

	// check alternative doner end pos vec
	vector<int> tmpAlterDonerSpliceSitePosVec;
	SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChrNameInt,
		tmpBackSpliceEndPos, tmpAlterDonerSpliceSitePosVec);
	for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
		if(tmpAlterDonerSpliceSitePos != tmpBackSpliceStartPos)
		{
			int tmpIndexInAlignInferJuncVec
				= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
					tmpChrNameInt, tmpAlterDonerSpliceSitePos, tmpBackSpliceEndPos);
			if(tmpIndexInAlignInferJuncVec < 0)
			{
				cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
				exit(1);
			}
			else
			{
				int tmpAlterSplice_anchorLength 
					= alignInferJuncHashInfo->returnAnchorSizeMax_doner(tmpIndexInAlignInferJuncVec);
				int tmpValidAnchorLength_doner = minAmongThreeNum(
					tmpAnchorSize_doner, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

				string tmpDonerSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt,
					tmpBackSpliceStartPos - tmpValidAnchorLength_doner + 1, 
					tmpValidAnchorLength_doner);
				string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt,
					tmpAlterDonerSpliceSitePos - tmpValidAnchorLength_doner + 1,
					tmpValidAnchorLength_doner);

				FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
				tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
					tmpDonerSeq, tmpAlterDonerSeq);
				int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();

				validAnchorLengthVec_alterDoner.push_back(tmpValidAnchorLength_doner);
				penaltyVec_alterDoner.push_back(tmpPenalty_alterDoner);
				alterSiteVec_alterDoner.push_back(tmpAlterDonerSpliceSitePos);
			}
		}
	}
}

void checkMatchThroughSeqSimilarity_backSpliceJunc(
	int tmpChrNameInt, int tmpBackSpliceStartPos, int tmpBackSpliceEndPos,
	int tmpAnchorSize_doner, int tmpAnchorSize_acceptor, 
	Index_Info* indexInfo, int toCheckAnchorLengthMax,
	int& validAnchorLength_matchThroughAtDonerSite, 
	int& validAnchorLength_matchThroughAtAcceptorSite,
	int& penalty_matchThroughAtDonerSite, 
	int& penalty_matchThroughAtAcceptorSite)
{
	int tmpChromSeqLength = indexInfo->returnChromLength(tmpChrNameInt);
	validAnchorLength_matchThroughAtDonerSite = tmpAnchorSize_acceptor;
	validAnchorLength_matchThroughAtAcceptorSite = tmpAnchorSize_doner;
	if(validAnchorLength_matchThroughAtDonerSite > toCheckAnchorLengthMax)
		validAnchorLength_matchThroughAtDonerSite = toCheckAnchorLengthMax;
	if(validAnchorLength_matchThroughAtAcceptorSite > toCheckAnchorLengthMax)
		validAnchorLength_matchThroughAtAcceptorSite = toCheckAnchorLengthMax;

	int matchThroughSeq_atDonerSite_startPos = tmpBackSpliceStartPos + 1;
	//int matchThroughSeq_atDonerSite_endPos 
	//	= matchThroughSeq_atDonerSite_startPos + validAnchorLength_acceptor - 1;
	string matchThroughSeq_atDonerSite
		= indexInfo->returnChromStrSubstr(tmpChrNameInt, 
			matchThroughSeq_atDonerSite_startPos, 
			validAnchorLength_matchThroughAtDonerSite);
	string anchorSeq_acceptorStart
		= indexInfo->returnChromStrSubstr(tmpChrNameInt,
			tmpBackSpliceEndPos, validAnchorLength_matchThroughAtDonerSite);

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_atDonerSite;
	tmpFixNWDPinfo_atDonerSite.doNWDP_withMismatchJumpCode(
		matchThroughSeq_atDonerSite, anchorSeq_acceptorStart);
	penalty_matchThroughAtDonerSite = tmpFixNWDPinfo_atDonerSite.getPenalty();

	int matchThroughSeq_atAcceptorSite_endPos = tmpBackSpliceEndPos - 1;
	int matchThroughSeq_atAcceptorSite_startPos 
		= matchThroughSeq_atAcceptorSite_endPos 
			- validAnchorLength_matchThroughAtAcceptorSite + 1;
	string matchThroughSeq_atAcceptorSite
		= indexInfo->returnChromStrSubstr(tmpChrNameInt,
			matchThroughSeq_atAcceptorSite_startPos, 
			validAnchorLength_matchThroughAtAcceptorSite);
	string anchorSeq_donerEnd
		= indexInfo->returnChromStrSubstr(tmpChrNameInt,
			tmpBackSpliceStartPos - validAnchorLength_matchThroughAtAcceptorSite + 1, 
			validAnchorLength_matchThroughAtAcceptorSite);	

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_atAcceptorSite;
	tmpFixNWDPinfo_atAcceptorSite.doNWDP_withMismatchJumpCode(
		matchThroughSeq_atAcceptorSite, anchorSeq_donerEnd);
	penalty_matchThroughAtAcceptorSite = tmpFixNWDPinfo_atAcceptorSite.getPenalty();
}

void extractBackSpliceInfoFromStr(int& tmpChrNameInt, 
	int& tmpBackSpliceStartPos, int& tmpBackSpliceEndPos,
	int& tmpAnchorSize_doner, int& tmpAnchorSize_acceptor,
	string& tmpBackSpliceStr, Index_Info* indexInfo)
{
	vector<string> tmpBackSpliceFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 6; tmp++)
	{
		int tabLoc = tmpBackSpliceStr.find("\t", startLoc);
		string tmpBackSpliceField = tmpBackSpliceStr.substr(startLoc, tabLoc - startLoc);
		tmpBackSpliceFieldVec.push_back(tmpBackSpliceField);
		startLoc = tabLoc + 1;
	}
	tmpBackSpliceFieldVec.push_back(tmpBackSpliceStr.substr(startLoc));
	string tmpChrNameStr = tmpBackSpliceFieldVec[0];
	string tmpBackSpliceStartPosStr = tmpBackSpliceFieldVec[1];
	string tmpBackSpliceEndPosStr = tmpBackSpliceFieldVec[2];
	//string tmpBackSpliceSupportNum = tmpBackSpliceFieldVec[4];
	string tmpBackSpliceDonerAnchorSizeStr = tmpBackSpliceFieldVec[5];
	string tmpBackSpliceAcceptorAnchorSizeStr = tmpBackSpliceFieldVec[6];
	tmpChrNameInt = indexInfo->convertStringToInt(tmpChrNameStr);
	tmpBackSpliceStartPos = atoi(tmpBackSpliceStartPosStr.c_str());
	tmpBackSpliceEndPos = atoi(tmpBackSpliceEndPosStr.c_str());
	tmpAnchorSize_doner = atoi(tmpBackSpliceDonerAnchorSizeStr.c_str());
	tmpAnchorSize_acceptor = atoi(tmpBackSpliceAcceptorAnchorSizeStr.c_str());
}

int main(int argc, char** argv)
{
	if((argc != 4)&&(argc != 5))
	{
		cout << "Executable inputIndexInfoPath inputBackSpliceJuncPath outputFolderPath";
		cout << " (inputNormalAlignInferJuncHashWithAnchorSizePath)" << endl;
		exit(1);
	}
	int toCheckAnchorLengthMax = 30;
	bool checkAlterSpliceBool;
	if(argc == 4)
		checkAlterSpliceBool = false;
	else
		checkAlterSpliceBool = true;

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string chrom_bit_file = indexFolderPath; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);	
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo ..." << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char));
	indexInfo->readGenome(chrom);
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	cout << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);\
	cout << "end of initiating indexInfo" << endl;


	string inputBackSpliceJuncPath = argv[2];
	ifstream backSplice_ifs(inputBackSpliceJuncPath.c_str());

	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	string outputPathPrefix = outputFolderStr + "/backSpliceJunc";
	string outputClassifiedBackSpliceJuncPath = outputPathPrefix + "_classified.txt";
	ofstream backSpliceJunc_classified_ofs(outputClassifiedBackSpliceJuncPath.c_str());
	string outputClassifiedBackSpliceJuncPath_pass = outputPathPrefix + "_classified_pass.txt";
	ofstream backSpliceJunc_pass_ofs(outputClassifiedBackSpliceJuncPath_pass.c_str());
	string outputClassifiedBackSpliceJuncPath_filterOut = outputPathPrefix + "_classified_filterOut.txt";
	ofstream backSpliceJunc_filterOut_ofs(outputClassifiedBackSpliceJuncPath_filterOut.c_str());
	
	cout << "start to initiateAlignInferJuncHashInfo" << endl;
	AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo = new AlignInferJunctionHash_Info();
	if(checkAlterSpliceBool)
	{
		string inputNormalAlignInferJuncHashWithAnchorSizePath = argv[4];
		tmpAlignInferJuncHashInfo->initiateAlignInferJunctionInfo(chromNum);
		tmpAlignInferJuncHashInfo->insertJuncFromJuncFile_chrNamePos_supportNum_anchorSize(
			inputNormalAlignInferJuncHashWithAnchorSizePath, indexInfo);
		tmpAlignInferJuncHashInfo->insertJuncFromJuncFile_chrNamePos_supportNum_anchorSize(
			inputBackSpliceJuncPath, indexInfo);
	}
	cout << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	tmpAlignInferJuncHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);

	while(!backSplice_ifs.eof())
	{
		string tmpBackSpliceStr;
		getline(backSplice_ifs, tmpBackSpliceStr);
		if(//(backSplice_ifs.eof())||
			(tmpBackSpliceStr == ""))
			break;
		int tmpChrNameInt;
		int tmpBackSpliceStartPos;
		int tmpBackSpliceEndPos;
		int tmpAnchorSize_doner;
		int tmpAnchorSize_acceptor;
		extractBackSpliceInfoFromStr(tmpChrNameInt, 
			tmpBackSpliceStartPos, tmpBackSpliceEndPos,
			tmpAnchorSize_doner, tmpAnchorSize_acceptor,
			tmpBackSpliceStr, indexInfo);

		string tmpBackSpliceJuncStr_afterClassify = "";
		int validAnchorLength_matchThroughAtDonerSite; 
		int validAnchorLength_matchThroughAtAcceptorSite;
		int penalty_matchThroughAtDonerSite; 
		int penalty_matchThroughAtAcceptorSite;
		checkMatchThroughSeqSimilarity_backSpliceJunc(
			tmpChrNameInt, tmpBackSpliceStartPos, tmpBackSpliceEndPos,
			tmpAnchorSize_doner, tmpAnchorSize_acceptor, indexInfo, 
			toCheckAnchorLengthMax,
			validAnchorLength_matchThroughAtDonerSite, 
			validAnchorLength_matchThroughAtAcceptorSite,
			penalty_matchThroughAtDonerSite, 
			penalty_matchThroughAtAcceptorSite);

		vector<int> validAnchorLengthVec_alterAcceptor;
		vector<int> validAnchorLengthVec_alterDoner;
		vector<int> penaltyVec_alterAcceptor;
		vector<int> penaltyVec_alterDoner;
		vector<int> alterSiteVec_alterAcceptor;
		vector<int> alterSiteVec_alterDoner;
		if(checkAlterSpliceBool)
		{
			checkAlterSpliceAnchorSeqSimilarity_backSpliceJunc(
				SJhashInfo, tmpAlignInferJuncHashInfo, tmpChrNameInt, 
				tmpBackSpliceStartPos, tmpBackSpliceEndPos,
				tmpAnchorSize_doner, tmpAnchorSize_acceptor, 
				indexInfo, toCheckAnchorLengthMax,	
				validAnchorLengthVec_alterAcceptor, validAnchorLengthVec_alterDoner,
				penaltyVec_alterAcceptor, penaltyVec_alterDoner,
				alterSiteVec_alterAcceptor, alterSiteVec_alterDoner);
		}
		// backSpliceClassified_ofs << tmpBackSpliceStr << "\t" 
		// 	<< validAnchorLength_acceptor << ":" << penalty_matchThroughAtDonerSite << "\t"
		// 	<< validAnchorLength_doner << ":" << penalty_matchThroughAtAcceptorSite << "\t";
		tmpBackSpliceJuncStr_afterClassify = tmpBackSpliceJuncStr_afterClassify
			+ tmpBackSpliceStr + "\t" 
			+ int_to_str(validAnchorLength_matchThroughAtDonerSite) + ":"
			+ int_to_str(penalty_matchThroughAtDonerSite) + "\t"
			+ int_to_str(validAnchorLength_matchThroughAtAcceptorSite) + ":"
			+ int_to_str(penalty_matchThroughAtAcceptorSite) + "\t";

		bool matchThroughSeqSimilar_bool = false;
		if((penalty_matchThroughAtDonerSite <= (validAnchorLength_matchThroughAtDonerSite/4))
			||(penalty_matchThroughAtAcceptorSite <= (validAnchorLength_matchThroughAtAcceptorSite/4)))
		{
			tmpBackSpliceJuncStr_afterClassify += "filterOut_matchThroughSeqSimi";
			matchThroughSeqSimilar_bool = true;
		}
		else
		{	
			tmpBackSpliceJuncStr_afterClassify += "kept_noMatchThroughSeqSimi";
			matchThroughSeqSimilar_bool = false;
		}

		bool anchorSeqSimilarAlterSpliceExists_alterAcceptor_bool = false;
		bool anchorSeqSimilarAlterSpliceExists_alterDoner_bool = false;
		if(checkAlterSpliceBool)
		{
			if(validAnchorLengthVec_alterAcceptor.size() == 0)
				tmpBackSpliceJuncStr_afterClassify += "\tnoAS\t";
			else
			{
				tmpBackSpliceJuncStr_afterClassify += "\t";
				for(int tmp = 0; tmp < validAnchorLengthVec_alterAcceptor.size(); tmp++)
				{
					int tmpValidAnchorLength = validAnchorLengthVec_alterAcceptor[tmp];
					int tmpPenalty = penaltyVec_alterAcceptor[tmp];
					tmpBackSpliceJuncStr_afterClassify = tmpBackSpliceJuncStr_afterClassify
						+ int_to_str(alterSiteVec_alterAcceptor[tmp]) + ":"
						+ int_to_str(tmpValidAnchorLength) + ":" + int_to_str(tmpPenalty) + ",";
					if(tmpPenalty <= (tmpValidAnchorLength/4))
						anchorSeqSimilarAlterSpliceExists_alterAcceptor_bool = true;
				}
				tmpBackSpliceJuncStr_afterClassify += "\t";
			}

			if(validAnchorLengthVec_alterDoner.size() == 0)
				tmpBackSpliceJuncStr_afterClassify += "\tnoAS\t";			
			else
			{
				tmpBackSpliceJuncStr_afterClassify +="\t";
				for(int tmp = 0; tmp < validAnchorLengthVec_alterDoner.size(); tmp++)
				{
					int tmpValidAnchorLength = validAnchorLengthVec_alterDoner[tmp];
					int tmpPenalty = penaltyVec_alterDoner[tmp];
					tmpBackSpliceJuncStr_afterClassify = tmpBackSpliceJuncStr_afterClassify
						+ int_to_str(alterSiteVec_alterDoner[tmp]) + ":"
						+ int_to_str(tmpValidAnchorLength) + ":" + int_to_str(tmpPenalty) + ",";
					if(tmpPenalty <= (tmpValidAnchorLength/4))
						anchorSeqSimilarAlterSpliceExists_alterDoner_bool = true;
				}
				tmpBackSpliceJuncStr_afterClassify += "\t";
			}

			if(anchorSeqSimilarAlterSpliceExists_alterAcceptor_bool
				|| anchorSeqSimilarAlterSpliceExists_alterDoner_bool)
				tmpBackSpliceJuncStr_afterClassify += "filterOut_alterSpliceAnchorSimi";
			else
				tmpBackSpliceJuncStr_afterClassify += "kept_noAlterSpliceAnchorSimi";
		}

		if(matchThroughSeqSimilar_bool
			||anchorSeqSimilarAlterSpliceExists_alterAcceptor_bool
			|| anchorSeqSimilarAlterSpliceExists_alterDoner_bool)
		{	
			tmpBackSpliceJuncStr_afterClassify += "\tFILTER_OUT";
			backSpliceJunc_filterOut_ofs << tmpBackSpliceJuncStr_afterClassify << endl;
		}
		else
		{	
			tmpBackSpliceJuncStr_afterClassify += "\tPASS";
			backSpliceJunc_pass_ofs << tmpBackSpliceJuncStr_afterClassify << endl;
		}
		backSpliceJunc_classified_ofs << tmpBackSpliceJuncStr_afterClassify << endl;
	}
	backSpliceJunc_classified_ofs.close();
	backSpliceJunc_pass_ofs.close();
	backSpliceJunc_filterOut_ofs.close();
	delete tmpAlignInferJuncHashInfo;
	delete indexInfo;
	return 0;
}