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

time_t nowtime;
struct tm *local;

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/alignInferJunctionHash_info.h"
#include "../general/alignInferJunctionHash_info_vec.h"

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

void checkAlterSpliceAnchorSeqSimilarity(
	AlignInferJunctionHash_Info* alignInferJuncHashInfo, 
	SJhash_Info* SJ, int tmpSJ_chrNameInt,
	int tmpSJ_donerEndPos, int tmpSJ_acceptorStartPos,
	int tmpSJ_anchorSize_doner, int tmpSJ_anchorSize_acceptor,
	vector<int>& validAnchorLengthVec_alterAcceptor, 
	vector<int>& validAnchorLengthVec_alterDoner,
	vector<int>& penaltyVec_alterAcceptor, 
	vector<int>& penaltyVec_alterDoner,
	vector<int>& alterSiteVec_alterAcceptor, 
	vector<int>& alterSiteVec_alterDoner,
	int toCheckAnchorLengthMax, Index_Info* indexInfo)
{
	vector<int> tmpAlterAcceptorSpliceSitePosVec;
	SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpSJ_chrNameInt,
		tmpSJ_donerEndPos, tmpAlterAcceptorSpliceSitePosVec);
	for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
		if(tmpAlterAcceptorSpliceSitePos == tmpSJ_acceptorStartPos)
			continue;
		int tmpIndexInAlignInferJuncVec
			= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
				tmpSJ_chrNameInt, tmpSJ_donerEndPos, tmpAlterAcceptorSpliceSitePos);
		if(tmpIndexInAlignInferJuncVec < 0)
		{
			cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
			exit(1);
		}
		else
		{
			int tmpAlterAcceptorSite_anchorLength 
				= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(tmpIndexInAlignInferJuncVec); 
			int tmpValidAnchorLength_acceptor = minAmongThreeNum(tmpSJ_anchorSize_acceptor, 
					tmpAlterAcceptorSite_anchorLength, toCheckAnchorLengthMax);			
			string tmpAcceptorSeq = indexInfo->returnChromStrSubstr(tmpSJ_chrNameInt,
					tmpSJ_acceptorStartPos, tmpValidAnchorLength_acceptor);
			string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(tmpSJ_chrNameInt,
					tmpAlterAcceptorSpliceSitePos, tmpValidAnchorLength_acceptor);
			FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
			tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(tmpAcceptorSeq, tmpAlterAcceptorSeq);
			int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
			validAnchorLengthVec_alterAcceptor.push_back(tmpValidAnchorLength_acceptor);
			penaltyVec_alterAcceptor.push_back(tmpPenalty_alterAcceptor);
			alterSiteVec_alterAcceptor.push_back(tmpAlterAcceptorSpliceSitePos);
		}
	}

	vector<int> tmpAlterDonerSpliceSitePosVec;
	SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpSJ_chrNameInt,
		tmpSJ_acceptorStartPos, tmpAlterDonerSpliceSitePosVec);	
	for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
	{
		int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
		if(tmpAlterDonerSpliceSitePos == tmpSJ_donerEndPos)
			continue;
		int tmpIndexInAlignInferJuncVec
			= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
				tmpSJ_chrNameInt, tmpAlterDonerSpliceSitePos, tmpSJ_acceptorStartPos);
		if(tmpIndexInAlignInferJuncVec < 0)
		{
			cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
			exit(1);
		}
		else
		{
			int tmpAlterSplice_anchorLength 
				= alignInferJuncHashInfo->returnAnchorSizeMax_doner(tmpIndexInAlignInferJuncVec); 
			int tmpValidAnchorLength_doner = minAmongThreeNum(tmpSJ_anchorSize_doner, 
					tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);			
			string tmpDonerSeq = indexInfo->returnChromStrSubstr(tmpSJ_chrNameInt,
					tmpSJ_donerEndPos - tmpValidAnchorLength_doner + 1, tmpValidAnchorLength_doner);
			string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(tmpSJ_chrNameInt,
					tmpAlterDonerSpliceSitePos - tmpValidAnchorLength_doner + 1, tmpValidAnchorLength_doner);
			FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
			tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(tmpDonerSeq, tmpAlterDonerSeq);
			int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
			validAnchorLengthVec_alterDoner.push_back(tmpValidAnchorLength_doner);
			penaltyVec_alterDoner.push_back(tmpPenalty_alterDoner);
			alterSiteVec_alterDoner.push_back(tmpAlterDonerSpliceSitePos);
		}
	}
}		

void checkMatchTroughSeqSimilarity(int tmpSJ_chrNameInt, 
	int tmpSJ_donerEndPos, int tmpSJ_acceptorStartPos,
	int tmpSJ_anchorSize_doner, int tmpSJ_anchorSize_acceptor,
	int& validAnchorLength_matchThrough_doner, int& validAnchorLength_matchThrough_acceptor,
	int& penalty_matchThrough_doner, int& penalty_matchThrough_acceptor, 
	int toCheckAnchorLength, Index_Info* indexInfo)
{
	validAnchorLength_matchThrough_doner = tmpSJ_anchorSize_acceptor;
	validAnchorLength_matchThrough_acceptor = tmpSJ_anchorSize_doner;
	if(validAnchorLength_matchThrough_doner > toCheckAnchorLength)
		validAnchorLength_matchThrough_doner = toCheckAnchorLength;
	if(validAnchorLength_matchThrough_acceptor > toCheckAnchorLength)
		validAnchorLength_matchThrough_acceptor = toCheckAnchorLength;	
	
	string matchThroughSeq_atDoner;
	string anchorSeq_atAcceptor;
	string matchThroughSeq_atAcceptor;
	string anchorSeq_atDoner;

	matchThroughSeq_atDoner = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_donerEndPos+1, validAnchorLength_matchThrough_doner);
	anchorSeq_atAcceptor = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_acceptorStartPos, validAnchorLength_matchThrough_doner);
	matchThroughSeq_atAcceptor = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_acceptorStartPos - 1 - validAnchorLength_matchThrough_acceptor + 1,
		validAnchorLength_matchThrough_acceptor);
	anchorSeq_atDoner = indexInfo->returnChromStrSubstr(
		tmpSJ_chrNameInt, tmpSJ_donerEndPos - validAnchorLength_matchThrough_acceptor + 1,
		validAnchorLength_matchThrough_acceptor);

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtDoner;
	tmpFixNWDPinfo_matchThroughAtDoner.doNWDP_withMismatchJumpCode(
		matchThroughSeq_atDoner, anchorSeq_atAcceptor);
	penalty_matchThrough_doner = tmpFixNWDPinfo_matchThroughAtDoner.getPenalty();

	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_matchThroughAtAcceptor;
	tmpFixNWDPinfo_matchThroughAtAcceptor.doNWDP_withMismatchJumpCode(
		matchThroughSeq_atAcceptor, anchorSeq_atDoner);
	penalty_matchThrough_acceptor = tmpFixNWDPinfo_matchThroughAtAcceptor.getPenalty();
}

bool inconfidentSJ_anchorSeqSimilarity(int penalty, int anchorSize)
{
	if(anchorSize <= 3)
	{
		if(penalty == 0)
			return true;
		else
			return false;
	}
	else
	{
		if(penalty <= 1 + anchorSize/8)
			return true;
		else
			return false;
	}
}

bool inconfidentSJ_lowSupportNumNonCanonical_bool(string& flankString, int supportNum)
{
	if((flankString != "GTAG")&&(flankString != "CTAC"))
	{
		if(supportNum <= 2)
			return true;
		else
			return false;
	}
	else
		return false;
}

int main(int argc, char**argv)
{
	if(argc < 4)
	{
		cout << "Executable <InputIndexInfo> <threadNum> <OutputSJ> <inputSAM_1> ... (other input SAM files)" << endl;
		exit(1);
	}
	//int maxReadBaseNumInPathStructure = 30;
	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string outputSJstr = outputFolderStr + "output.alignInferJunc.txt";
	string outputLogStr = outputFolderStr + "log";
	ofstream log_ofs(outputLogStr.c_str());

	string outputSJstr_withAnchorSeqSimilarity 
		= outputFolderStr + "output.alignInferJunc_classified.txt";
	ofstream SJ_ofs(outputSJstr_withAnchorSeqSimilarity.c_str());
	string outputSJstr_withAnchorSeqSimilarity_filterOut 
		= outputFolderStr + "output.alignInferJunc_classified_filterOut.txt";
	ofstream filterOut_ofs(outputSJstr_withAnchorSeqSimilarity_filterOut.c_str());
	string outputSJstr_withAnchorSeqSimilarity_kept 
		= outputFolderStr + "output.alignInferJunc_classified_kept.txt";
	ofstream kept_ofs(outputSJstr_withAnchorSeqSimilarity_kept.c_str());

	string threadNumStr = argv[2];
	int threadNum_int = atoi(threadNumStr.c_str());

	int alignInferJuncHashInfoVecSize = threadNum_int;

	vector<string> inputSAMfileVec;
	for(int tmp = 4; tmp < argc; tmp++)
	{
		string tmpSAMstr = argv[tmp];
		inputSAMfileVec.push_back(tmpSAMstr);
	}

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); 
	ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	cout << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initiate chrNameIndexArray" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initiate chrNameIndexArray" << endl;	
	indexInfo->initiateChrNameIndexArray(1000);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite merged alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite merged alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo_merged 
		= new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo_merged->initiateAlignInferJunctionHashInfo(chromNum);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite alignInferJunctionHashInfo vec " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite alignInferJunctionHashInfo vec " << endl;	
	AlignInferJunctionHash_Info_Vec* alignInferJunctionHashInfoVec 
		= new AlignInferJunctionHash_Info_Vec();
	alignInferJunctionHashInfoVec->initiateAlignInferJunctionHashInfoVec(alignInferJuncHashInfoVecSize, chromNum);
	
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to insert SJ into SJmap from SAM file" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to insert SJ into SJmap from SAM file" << endl;
	// insert SJ into SJmap
	alignInferJunctionHashInfoVec->insertJuncFromAlignmentFileVec_chrNamePos_supportNum_anchorSize_XM_parallel(
		inputSAMfileVec, indexInfo, alignInferJuncHashInfoVecSize, log_ofs);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to merge alignInferJuncHashInfo in vec" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to merge alignInferJuncHashInfo in vec" << endl;
	alignInferJunctionHashInfoVec->mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum_anchorSize_XM(
		alignInferJunctionHashInfo_merged, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initiate SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initiate SJhashInfo" << endl;	
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to convert alignment2juncInfoHash 2 SJhashInfo" << endl;
	alignInferJunctionHashInfo_merged->convert2SJhashInfo(SJ, indexInfo);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to output SJ map" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to output SJ map" << endl;	
	// output SJmap
	alignInferJunctionHashInfo_merged->outputAlignInferInfoHashInfo_chrNamePos_supportNum_anchorSize_XM(
		indexInfo, outputSJstr);

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "end of running getAlignInferJuncFromSAMfile_chrNamePos_supportNum_anchorSize_XM" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "end of running getAlignInferJuncFromSAMfile_chrNamePos_supportNum_anchorSize_XM" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to do classification" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to do classification" << endl;

	int toCheckAnchorLengthMax = 30;
	int alignInferInfoVecSize = alignInferJunctionHashInfo_merged->returnAlignInferInfoVecSize();
	for(int tmp = 0; tmp < alignInferInfoVecSize; tmp++)
	{
		int tmpSJ_chrNameInt = alignInferJunctionHashInfo_merged->returnAlignInferInfo_chrNameInt(tmp);
		int tmpSJ_donerEndPos = alignInferJunctionHashInfo_merged->returnAlignInferInfo_donerEndPos(tmp);
		int tmpSJ_acceptorStartPos = alignInferJunctionHashInfo_merged->returnAlignInferInfo_acceptorStartPos(tmp);
		int tmpSJ_anchorSize_doner = alignInferJunctionHashInfo_merged->returnAnchorSizeMax_doner(tmp);
		int tmpSJ_anchorSize_acceptor = alignInferJunctionHashInfo_merged->returnAnchorSizeMax_acceptor(tmp);
		int tmpSJ_supportNum = alignInferJunctionHashInfo_merged->returnAlignInferInfo_supportNum(tmp);
		int tmpSJ_XMmin = alignInferJunctionHashInfo_merged->returnAlignInferInfo_XMmin(tmp);
		int tmpSJ_XMmax = alignInferJunctionHashInfo_merged->returnAlignInferInfo_XMmax(tmp);
		string tmpSJ_flankString = alignInferJunctionHashInfo_merged->returnAlignInferInfo_flankString(tmp, indexInfo);

		int validAnchorLength_matchThrough_doner;
		int validAnchorLength_matchThrough_acceptor;
		int penalty_matchThrough_doner;
		int penalty_matchThrough_acceptor;
		checkMatchTroughSeqSimilarity(tmpSJ_chrNameInt, 
			tmpSJ_donerEndPos, tmpSJ_acceptorStartPos,
			tmpSJ_anchorSize_doner, tmpSJ_anchorSize_acceptor,
			validAnchorLength_matchThrough_doner, validAnchorLength_matchThrough_acceptor,
			penalty_matchThrough_doner, penalty_matchThrough_acceptor, 
			toCheckAnchorLengthMax, indexInfo);

		vector<int> validAnchorLengthVec_alterAcceptor;
		vector<int> validAnchorLengthVec_alterDoner;
		vector<int> penaltyVec_alterAcceptor;
		vector<int> penaltyVec_alterDoner;
		vector<int> alterSiteVec_alterAcceptor;
		vector<int> alterSiteVec_alterDoner;
		checkAlterSpliceAnchorSeqSimilarity(
			alignInferJunctionHashInfo_merged, SJ, tmpSJ_chrNameInt, 
			tmpSJ_donerEndPos, tmpSJ_acceptorStartPos,
			tmpSJ_anchorSize_doner, tmpSJ_anchorSize_acceptor,
			validAnchorLengthVec_alterAcceptor, validAnchorLengthVec_alterDoner,
			penaltyVec_alterAcceptor, penaltyVec_alterDoner,
			alterSiteVec_alterAcceptor, alterSiteVec_alterDoner,
			toCheckAnchorLengthMax, indexInfo);

		string tmpJuncStr = indexInfo->returnChrNameStr(tmpSJ_chrNameInt) + "\t"
			+ int_to_str(tmpSJ_donerEndPos) + "\t" + int_to_str(tmpSJ_acceptorStartPos) + "\t"
			+ "JUNC_" + int_to_str(tmp+1) + "\t" + int_to_str(tmpSJ_supportNum) + "\t"
			+ int_to_str(tmpSJ_anchorSize_doner) + "\t" + int_to_str(tmpSJ_anchorSize_acceptor) + "\t"
			+ int_to_str(tmpSJ_XMmin) + "\t" + int_to_str(tmpSJ_XMmax) + "\t" + tmpSJ_flankString + "\t";
		bool lowSupportNumNonCanonical_bool = inconfidentSJ_lowSupportNumNonCanonical_bool(tmpSJ_flankString, tmpSJ_supportNum);
		if(lowSupportNumNonCanonical_bool)
			tmpJuncStr += "filterOut_signal\t";
		else
			tmpJuncStr += "kept_signal\t";
		
		tmpJuncStr = tmpJuncStr	
			+ int_to_str(validAnchorLength_matchThrough_doner) + ":" + int_to_str(penalty_matchThrough_doner) + "\t"
			+ int_to_str(validAnchorLength_matchThrough_acceptor) + ":" + int_to_str(penalty_matchThrough_acceptor) + "\t";
		bool matchThroughSeqSimilar_bool = ((inconfidentSJ_anchorSeqSimilarity(penalty_matchThrough_doner, validAnchorLength_matchThrough_doner))
			||(inconfidentSJ_anchorSeqSimilarity(penalty_matchThrough_acceptor, validAnchorLength_matchThrough_acceptor)));
		if(matchThroughSeqSimilar_bool)
			tmpJuncStr += "filterOut_MT\t";
		else
			tmpJuncStr += "kept_MT\t";	

		bool alterSpliceSiteSeqSimilar_bool = false;
		if(validAnchorLengthVec_alterAcceptor.size() == 0)
			tmpJuncStr += "noAS\t";
		else
		{
			//tmpJuncStr += "\t";
			for(int tmp = 0; tmp < validAnchorLengthVec_alterAcceptor.size(); tmp++)
			{
				int tmpValidAnchorLength = validAnchorLengthVec_alterAcceptor[tmp];
				int tmpPenalty = penaltyVec_alterAcceptor[tmp];
				int tmpAlterSite = alterSiteVec_alterAcceptor[tmp];
				if(inconfidentSJ_anchorSeqSimilarity(tmpPenalty, tmpValidAnchorLength))
					alterSpliceSiteSeqSimilar_bool = true;
				tmpJuncStr = tmpJuncStr + int_to_str(tmpAlterSite) + ":"
					+ int_to_str(tmpValidAnchorLength) + ":"
					+ int_to_str(tmpPenalty) + ",";
			}
			tmpJuncStr += "\t";
		}
		if(validAnchorLengthVec_alterDoner.size() == 0)
			tmpJuncStr += "noAS\t";
		else
		{
			//tmpJuncStr += "\t";
			for(int tmp = 0; tmp < validAnchorLengthVec_alterDoner.size(); tmp++)
			{
				int tmpValidAnchorLength = validAnchorLengthVec_alterDoner[tmp];
				int tmpPenalty = penaltyVec_alterDoner[tmp];
				int tmpAlterSite = alterSiteVec_alterDoner[tmp];
				if(inconfidentSJ_anchorSeqSimilarity(tmpPenalty, tmpValidAnchorLength))
					alterSpliceSiteSeqSimilar_bool = true;
				tmpJuncStr = tmpJuncStr + int_to_str(tmpAlterSite) + ":"
					+ int_to_str(tmpValidAnchorLength) + ":"
					+ int_to_str(tmpPenalty) + ",";
			}
			tmpJuncStr += "\t";
		}	
		if(alterSpliceSiteSeqSimilar_bool)
			tmpJuncStr += "filterOut_AS\t";
		else
			tmpJuncStr += "kept_AS\t";	


		if(lowSupportNumNonCanonical_bool || matchThroughSeqSimilar_bool || alterSpliceSiteSeqSimilar_bool)
		{
			tmpJuncStr += "FILTER_OUT";
			SJ_ofs << tmpJuncStr << endl;
			filterOut_ofs << tmpJuncStr << endl;
		}	
		else
		{
			tmpJuncStr += "PASS";
			SJ_ofs << tmpJuncStr << endl;
			kept_ofs << tmpJuncStr << endl;
		}
	}

	SJ_ofs.close();
	filterOut_ofs.close();
	kept_ofs.close();
	alignInferJunctionHashInfoVec->freeMemory();
	delete alignInferJunctionHashInfoVec;
	delete alignInferJunctionHashInfo_merged;
	delete indexInfo;
	//sj_ofs.close();
	parameter_ifs.close();
	return 0;
}