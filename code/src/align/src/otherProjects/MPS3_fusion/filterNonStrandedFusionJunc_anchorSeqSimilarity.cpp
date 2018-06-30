// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
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

void checkAlterSpliceAnchorSeqSimilarity_fusionJunc(
	AlignInferJunctionHash_Info* alignInferJuncHashInfo, SJhash_Info* SJ,
	string& strand_1, string& strand_2, 
	int chrNameInt_1, int chrNameInt_2, int breakPointPos_1, int breakPointPos_2,
	int anchorLength_1, int anchorLength_2, int toCheckAnchorLengthMax, Index_Info* indexInfo, 
	vector<int>& validAnchorLengthVec_alterSiteAtGene1toReplaceGene2, 
	vector<int>& validAnchorLengthVec_alterSiteAtGene2toReplaceGene1,
	vector<int>& penaltyVec_alterSiteAtGene1toReplaceGene2, 
	vector<int>& penaltyVec_alterSiteAtGene2toReplaceGene1,
	vector<int>& alterSiteVec_atGene1toReplaceGene2, 
	vector<int>& alterSiteVec_atGene2toReplaceGene1)
{
	if(((strand_1 == "+")&&(strand_2 == "+"))
		||((strand_1 == "-")&&(strand_2 == "-")))
	{
		int tmpChrNameInt_doner;
		int tmpChrNameInt_acceptor;
		int tmpChrPos_doner; 
		int tmpChrPos_acceptor;
		int tmpAnchorSize_doner;
		int tmpAnchorSize_acceptor;

		if((strand_1 == "+")&&(strand_2 == "+"))
		{
			tmpChrNameInt_doner = chrNameInt_1;
			tmpChrNameInt_acceptor = chrNameInt_2;
			tmpChrPos_doner = breakPointPos_1;
			tmpChrPos_acceptor = breakPointPos_2;
			tmpAnchorSize_doner = anchorLength_1;
			tmpAnchorSize_acceptor = anchorLength_2;
		}
		else
		{
			tmpChrNameInt_doner = chrNameInt_2;
			tmpChrNameInt_acceptor = chrNameInt_1;
			tmpChrPos_doner = breakPointPos_2;
			tmpChrPos_acceptor = breakPointPos_1;
			tmpAnchorSize_doner = anchorLength_2;
			tmpAnchorSize_acceptor = anchorLength_1;
		}
		// cout << "tmpChrNameInt_doner: " << tmpChrNameInt_doner << endl;
		// cout << "tmpChrNameInt_acceptor: " << tmpChrNameInt_acceptor << endl;
		// cout << "tmpChrPos_doner: " << tmpChrPos_doner << endl;
		// cout << "tmpChrPos_acceptor: " << tmpChrPos_acceptor << endl;
		// cout << "tmpAnchorSize_doner: " << tmpAnchorSize_doner << endl;
		// cout << "tmpAnchorSize_acceptor: " << tmpAnchorSize_acceptor << endl; 
		// check alternative acceptor start pos vec
		vector<int> tmpAlterAcceptorSpliceSitePosVec;
		SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(tmpChrNameInt_doner,
			tmpChrPos_doner, tmpAlterAcceptorSpliceSitePosVec);
		for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
		{
			int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
			if((tmpAlterAcceptorSpliceSitePos != tmpChrPos_acceptor)
				||(tmpChrNameInt_doner != tmpChrNameInt_acceptor))
			{
				int tmpIndexInAlignInferJuncVec
					= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
						tmpChrNameInt_doner, tmpChrPos_doner, tmpAlterAcceptorSpliceSitePos);
				if(tmpIndexInAlignInferJuncVec < 0)
				{
					cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
					exit(1);
				}
				else
				{
					int tmpAlterAcceptorSite_anchorLength 
						= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(
							tmpIndexInAlignInferJuncVec); 
					int tmpValidAnchorLength_acceptor = minAmongThreeNum(
						tmpAnchorSize_acceptor, tmpAlterAcceptorSite_anchorLength, toCheckAnchorLengthMax);
					string tmpAcceptorSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_acceptor,
						tmpChrPos_acceptor, tmpValidAnchorLength_acceptor);
					string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_doner,
						tmpAlterAcceptorSpliceSitePos, tmpValidAnchorLength_acceptor);
					FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
					tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(tmpAcceptorSeq, tmpAlterAcceptorSeq);
					int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
					if((strand_1 == "+")&&(strand_2 == "+"))
					{
						validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_acceptor);
						penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterAcceptor);
						alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterAcceptorSpliceSitePos);
					}
					else
					{
						validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_acceptor);
						penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterAcceptor);
						alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterAcceptorSpliceSitePos);						
					}
				}
			}
		}
		// check alternative doner end pos vec
		vector<int> tmpAlterDonerSpliceSitePosVec;
		SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(tmpChrNameInt_acceptor,
			tmpChrPos_acceptor, tmpAlterDonerSpliceSitePosVec);
		for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
		{
			int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
			if((tmpAlterDonerSpliceSitePos != tmpChrPos_doner)
				||(tmpChrNameInt_doner != tmpChrNameInt_acceptor))
			{
				int tmpIndexInAlignInferJuncVec
					= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(
						tmpChrNameInt_acceptor, tmpAlterDonerSpliceSitePos, 
						tmpChrPos_acceptor);
				if(tmpIndexInAlignInferJuncVec < 0)
				{
					cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
					exit(1);
				}
				else
				{
					int tmpAlterSplice_anchorLength 
						= alignInferJuncHashInfo->returnAnchorSizeMax_doner(
							tmpIndexInAlignInferJuncVec); 
					int tmpValidAnchorLength_doner = minAmongThreeNum(
						tmpAnchorSize_doner, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);
					string tmpDonerSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_doner,
						tmpChrPos_doner - tmpValidAnchorLength_doner + 1, tmpValidAnchorLength_doner);
					string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(tmpChrNameInt_acceptor,
						tmpAlterDonerSpliceSitePos - tmpValidAnchorLength_doner + 1, 
						tmpValidAnchorLength_doner);
					FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
					tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
						tmpDonerSeq, tmpAlterDonerSeq);
					int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
					if((strand_1 == "+")&&(strand_2 == "+"))
					{
						validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_doner);
						penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterDoner);
						alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterDonerSpliceSitePos);
					}
					else
					{
						validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_doner);
						penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterDoner);
						alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterDonerSpliceSitePos);					
					}
				}
			}
		}
	}
	else if((strand_1 == "+")&&(strand_2 == "-"))
	{
		// donerEndPos = breakPointPos_1
		vector<int> tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2;
		SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt_1,
			breakPointPos_1, tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2);		
		for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2.size(); tmp++)
		{
			int tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2 
				= tmpAlterAcceptorSpliceSitePosVec_atGene1toReplaceGene2[tmp];
			int tmpIndexInAlignInferJuncVec 
				= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_1,
					breakPointPos_1, tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2);
			if(tmpIndexInAlignInferJuncVec < 0)
			{
				cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
				exit(1);
			}
			else
			{
				int tmpAlterSplice_anchorLength
					= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(
						tmpIndexInAlignInferJuncVec);
				int tmpValidAnchorLength_acceptor = minAmongThreeNum(
					anchorLength_2, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

				string tmpAcceptorSeq = convertStringToReverseComplement(
					indexInfo->returnChromStrSubstr(chrNameInt_2,
						breakPointPos_2 - tmpValidAnchorLength_acceptor + 1, 
						tmpValidAnchorLength_acceptor));
				string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(
					chrNameInt_1, tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2, 
					tmpValidAnchorLength_acceptor);
				FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
				tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(
					tmpAcceptorSeq, tmpAlterAcceptorSeq);
				int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
				validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_acceptor);
				penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterAcceptor);
				alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterAcceptorSpliceSitePos_atGene1toReplaceGene2);
			}	
		}		
		// donerEndPos = breakPointPos_2		
		vector<int> tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1;
		SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt_2,
			breakPointPos_2, tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1);		
		for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1.size(); tmp++)
		{	
			int tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1
				= tmpAlterAcceptorSpliceSitePosVec_atGene2toReplaceGene1[tmp];
			int tmpIndexInAlignInferJuncVec 
				= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_2, 
					breakPointPos_2, tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1);
			if(tmpIndexInAlignInferJuncVec < 0)
			{
				cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
				exit(1);
			}
			else
			{
				int tmpAlterSplice_anchorLength
					= alignInferJuncHashInfo->returnAnchorSizeMax_acceptor(
						tmpIndexInAlignInferJuncVec);
				int tmpValidAnchorLength_acceptor = minAmongThreeNum(
					anchorLength_1, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);
				string tmpAcceptorSeq = convertStringToReverseComplement(
					indexInfo->returnChromStrSubstr(chrNameInt_1,
						breakPointPos_1 - tmpValidAnchorLength_acceptor + 1, 
						tmpValidAnchorLength_acceptor));
				string tmpAlterAcceptorSeq = indexInfo->returnChromStrSubstr(
					chrNameInt_2, tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1, 
					tmpValidAnchorLength_acceptor);
				FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterAcceptor;
				tmpFixNWDPinfo_alterAcceptor.doNWDP_withMismatchJumpCode(
					tmpAcceptorSeq, tmpAlterAcceptorSeq);
				int tmpPenalty_alterAcceptor = tmpFixNWDPinfo_alterAcceptor.getPenalty();
				validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_acceptor);
				penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterAcceptor);
				alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterAcceptorSpliceSitePos_atGene2toReplaceGene1);
			}
		}	
	}
	else // strand_1 == "-" && strand_2 == "+"
	{
		// acceptorStartPos = breakPointPos_1
		vector<int> tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2;
		SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt_1,
			breakPointPos_1, tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2);
		for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2.size(); tmp++)
		{
			int tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2
				= tmpAlterDonerSpliceSitePosVec_atGene1toReplaceGene2[tmp];
			int tmpIndexInAlignInferJuncVec
				= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_1,
					tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2, breakPointPos_1);
			if(tmpIndexInAlignInferJuncVec < 0)
			{
				cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
				exit(1);
			}
			else
			{
				int tmpAlterSplice_anchorLength
					= alignInferJuncHashInfo->returnAnchorSizeMax_doner(
						tmpIndexInAlignInferJuncVec);
				int tmpValidAnchorLength_doner = minAmongThreeNum(
					anchorLength_2, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

				string tmpDonerSeq = convertStringToReverseComplement(
					indexInfo->returnChromStrSubstr(chrNameInt_2,
						breakPointPos_2, tmpValidAnchorLength_doner));
				string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(chrNameInt_1, 
					tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2 - tmpValidAnchorLength_doner + 1, 
					tmpValidAnchorLength_doner);
			
				FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
				tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
					tmpDonerSeq, tmpAlterDonerSeq);
				int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
				validAnchorLengthVec_alterSiteAtGene1toReplaceGene2.push_back(tmpValidAnchorLength_doner);
				penaltyVec_alterSiteAtGene1toReplaceGene2.push_back(tmpPenalty_alterDoner);
				alterSiteVec_atGene1toReplaceGene2.push_back(tmpAlterDonerSpliceSitePos_atGene1toReplaceGene2);			
			}
		}
		// acceptorStartPos = breakPointPos_2
		vector<int> tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1;
		SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt_2,
			breakPointPos_2, tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1);
		for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1.size(); tmp++)
		{
			int tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1
				= tmpAlterDonerSpliceSitePosVec_atGene2toReplaceGene1[tmp];
			int tmpIndexInAlignInferJuncVec
				= alignInferJuncHashInfo->searchAndReturnAlignInferInfoVecIndex(chrNameInt_2,
					tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1, breakPointPos_2);
			if(tmpIndexInAlignInferJuncVec < 0)
			{
				cout << "SJ does not exist in alignInferJuncHashInfo" << endl;
				exit(1);
			}
			else
			{
				int tmpAlterSplice_anchorLength
					= alignInferJuncHashInfo->returnAnchorSizeMax_doner(
						tmpIndexInAlignInferJuncVec);
				int tmpValidAnchorLength_doner = minAmongThreeNum(
					anchorLength_1, tmpAlterSplice_anchorLength, toCheckAnchorLengthMax);

				string tmpDonerSeq = convertStringToReverseComplement(
					indexInfo->returnChromStrSubstr(chrNameInt_1,
						breakPointPos_1, tmpValidAnchorLength_doner));
				string tmpAlterDonerSeq = indexInfo->returnChromStrSubstr(chrNameInt_2,
					tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1 - tmpValidAnchorLength_doner + 1,
					tmpValidAnchorLength_doner);

				FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_alterDoner;
				tmpFixNWDPinfo_alterDoner.doNWDP_withMismatchJumpCode(
					tmpDonerSeq, tmpAlterDonerSeq);
				int tmpPenalty_alterDoner = tmpFixNWDPinfo_alterDoner.getPenalty();
				validAnchorLengthVec_alterSiteAtGene2toReplaceGene1.push_back(tmpValidAnchorLength_doner);
				penaltyVec_alterSiteAtGene2toReplaceGene1.push_back(tmpPenalty_alterDoner);
				alterSiteVec_atGene2toReplaceGene1.push_back(tmpAlterDonerSpliceSitePos_atGene2toReplaceGene1);
			}
		}		
	}
}

void checkMatchTroughSeqSimilarity_fusionJunc(string& strand_1, string& strand_2, 
	int chrNameInt_1, int chrNameInt_2, int breakPointPos_1, int breakPointPos_2,
	int anchorLength_1, int anchorLength_2, int toCheckAnchorLengthMax, Index_Info* indexInfo, 
	int& validAnchorLength_matchThroughAtGene1, int& validAnchorLength_matchThroughAtGene2,
	int& penalty_matchThroughAtGene1, int& penalty_matchThroughAtGene2)
{
	// if((strand_1 == "N") || (strand_2 == "N"))
	// 	return false;

	validAnchorLength_matchThroughAtGene1 = anchorLength_2;
	validAnchorLength_matchThroughAtGene2 = anchorLength_1;
	if(validAnchorLength_matchThroughAtGene1 > toCheckAnchorLengthMax)
		validAnchorLength_matchThroughAtGene1 = toCheckAnchorLengthMax;
	if(validAnchorLength_matchThroughAtGene2 > toCheckAnchorLengthMax)
		validAnchorLength_matchThroughAtGene2 = toCheckAnchorLengthMax;

	string matchThroughSeq_gene1;
	string anchorSeq_gene2;

	string matchThroughSeq_gene2;
	string anchorSeq_gene1;
	// cout << "validAnchorLength_1: " << validAnchorLength_1 << endl;
	// cout << "validAnchorLength_2: " << validAnchorLength_2 << endl;
	// cout << "strand_1: " << strand_1 << endl;
	// cout << "strand_2: " << strand_2 << endl;
	if((strand_1 == "+")&&(strand_2 == "+"))
	{
		matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1+1, validAnchorLength_matchThroughAtGene1);
		anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2, validAnchorLength_matchThroughAtGene1);
		matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_matchThroughAtGene2 + 1, 
			validAnchorLength_matchThroughAtGene2);
		anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - validAnchorLength_matchThroughAtGene2 + 1, 
			validAnchorLength_matchThroughAtGene2);
	}
	else if((strand_1 == "-")&&(strand_2 == "-"))
	{
		matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_matchThroughAtGene1 + 1, 
			validAnchorLength_matchThroughAtGene1));
		anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - validAnchorLength_matchThroughAtGene1 + 1, 
			validAnchorLength_matchThroughAtGene1));
		matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_matchThroughAtGene2));
		anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1, validAnchorLength_matchThroughAtGene2));
	}
	else if((strand_1 == "+")&&(strand_2 == "-"))
	{
		matchThroughSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1+1, validAnchorLength_matchThroughAtGene1);
		anchorSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - validAnchorLength_matchThroughAtGene1 + 1, 
			validAnchorLength_matchThroughAtGene1));
		matchThroughSeq_gene2 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 + 1, validAnchorLength_matchThroughAtGene2));
		anchorSeq_gene1 = indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - validAnchorLength_matchThroughAtGene2 + 1, 
			validAnchorLength_matchThroughAtGene2);
	}
	else
	{
		matchThroughSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1 - 1 - validAnchorLength_matchThroughAtGene1 + 1, 
			validAnchorLength_matchThroughAtGene1));
		anchorSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2, validAnchorLength_matchThroughAtGene1);
		matchThroughSeq_gene2 = indexInfo->returnChromStrSubstr(
			chrNameInt_2, breakPointPos_2 - 1 - validAnchorLength_matchThroughAtGene2 + 1, 
			validAnchorLength_matchThroughAtGene2);
		anchorSeq_gene1 = convertStringToReverseComplement(indexInfo->returnChromStrSubstr(
			chrNameInt_1, breakPointPos_1, validAnchorLength_matchThroughAtGene2));
	}
	//cout << "matchThroughSeq_gene1: " << matchThroughSeq_gene1 << endl;
	//cout << "matchThroughSeq_gene2: " << matchThroughSeq_gene2 << endl;
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_1;
	tmpFixNWDPinfo_1.doNWDP_withMismatchJumpCode(
		matchThroughSeq_gene1, anchorSeq_gene2);
	penalty_matchThroughAtGene1 = tmpFixNWDPinfo_1.getPenalty();
	//delete tmpFixNWDPinfo_1;
	FixSingleAnchor_NWDP_Info tmpFixNWDPinfo_2;
	tmpFixNWDPinfo_2.doNWDP_withMismatchJumpCode(
		matchThroughSeq_gene2, anchorSeq_gene1);
	penalty_matchThroughAtGene2 = tmpFixNWDPinfo_2.getPenalty();
	//delete tmpFixNWDPinfo_2;
	//return true;
}

void extractFusionJuncInfoFromStr(
	int& tmpChrNameInt_1, int& tmpChrNameInt_2,
	int& tmpBreakPointPos_1, int& tmpBreakPointPos_2,
	string& tmpStrand_1, string& tmpStrand_2, 
	int& tmpAnchorSize_1, int& tmpAnchorSize_2,
	string& tmpCandiFusionJuncCaseStr, string& tmpFusionJuncStr,
	Index_Info* indexInfo)
{
	vector<string> tmpFusionJuncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 13; tmp++)
	{
		int tabLoc = tmpFusionJuncStr.find("\t", startLoc);
		string tmpFusionJuncField = tmpFusionJuncStr.substr(startLoc, tabLoc-startLoc);
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncField);
		startLoc = tabLoc + 1;
	}
	int nextTabLoc = tmpFusionJuncStr.find("\t", startLoc);
	if(nextTabLoc == string::npos)
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncStr.substr(startLoc));
	else
		tmpFusionJuncFieldVec.push_back(tmpFusionJuncStr.substr(startLoc, nextTabLoc - startLoc));
	//cout << "tmpCandiFusionJuncCaseStr: " << tmpFusionJuncStr.substr(startLoc) << endl;
	tmpChrNameInt_1 = indexInfo->convertStringToInt(tmpFusionJuncFieldVec[0]);
	tmpChrNameInt_2 = indexInfo->convertStringToInt(tmpFusionJuncFieldVec[1]);	
	string tmpChrPosStr_1 = tmpFusionJuncFieldVec[2];
	string tmpChrPosStr_2 = tmpFusionJuncFieldVec[3];
	tmpBreakPointPos_1 = atoi(tmpChrPosStr_1.c_str());
	tmpBreakPointPos_2 = atoi(tmpChrPosStr_2.c_str());
	tmpStrand_1 = tmpFusionJuncFieldVec[4];
	tmpStrand_2 = tmpFusionJuncFieldVec[5];
	tmpAnchorSize_1 = atoi(tmpFusionJuncFieldVec[7].c_str());
	tmpAnchorSize_2 = atoi(tmpFusionJuncFieldVec[8].c_str());
	tmpCandiFusionJuncCaseStr = tmpFusionJuncFieldVec[13];
}

int main(int argc, char** argv)
{
	if((argc != 4)&&(argc != 5))
	{
		cout << "Executable inputIndexInfoPath inputNonStrandedFusionJuncPath outputFolderPath";
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
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "end of initiating indexInfo" << endl;

	string inputNonStrandedFusionJuncPath = argv[2];
	ifstream nonStrandedFusionJunc_ifs(inputNonStrandedFusionJuncPath);

	string outputFolderStr = argv[3];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());

	string outputPathPrefix = outputFolderStr + "/fusionJunc";
	string outputNonStrandedFusionJuncPath_classified = outputPathPrefix + "_classified.txt";
	ofstream nonStrandedFusionJunc_classified_ofs(outputNonStrandedFusionJuncPath_classified.c_str());
	string outputNonStrandedFusionJuncPath_pass = outputPathPrefix + "_classified_pass.txt";
	ofstream nonStrandedFusionJunc_pass_ofs(outputNonStrandedFusionJuncPath_pass.c_str());
	string outputNonStrandedFusionJuncPath_filterOut = outputPathPrefix + "_classified_filterOut.txt";
	ofstream nonStrandedFusionJunc_filterOut_ofs(outputNonStrandedFusionJuncPath_filterOut.c_str());	

	cout << "start to initiateAlignInferJuncHashInfo" << endl;
	AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo = new AlignInferJunctionHash_Info();
	if(checkAlterSpliceBool)
	{
		string inputNormalAlignInferJuncHashWithAnchorSizePath = argv[4];
		tmpAlignInferJuncHashInfo->initiateAlignInferJunctionInfo(chromNum);
		tmpAlignInferJuncHashInfo->insertJuncFromJuncFile_chrNamePos_supportNum_anchorSize(
			inputNormalAlignInferJuncHashWithAnchorSizePath, indexInfo);
	}
	cout << "start to initaite SJhashInfo" << endl;
	SJhash_Info* SJhashInfo = new SJhash_Info();
	SJhashInfo->initiateAreaAndStringHash(chromNum);
	cout << "start to convert 2 SJhashInfo" << endl;
	tmpAlignInferJuncHashInfo->convert2SJhashInfo(SJhashInfo, indexInfo);
	cout << "end of converting 2 SJhashInfo" << endl;
	cout << "start to do filtering based on anchorSeq similarity" << endl;
	while(!nonStrandedFusionJunc_ifs.eof())
	{
		string tmpFusionJuncStr;
		getline(nonStrandedFusionJunc_ifs, tmpFusionJuncStr);
		//cout << "tmpFusionJuncStr: " << tmpFusionJuncStr << endl;
		if(tmpFusionJuncStr == "")
			break;
		int detectedChrNameInt_left, detectedChrNameInt_right;
		int detectedBreakPoint_left, detectedBreakPoint_right;
		string detectedFusionJuncStrand_left, detectedFusionJuncStrand_right;
		int detectedFusionAnchorLength_left, detectedFusionAnchorLength_right;
		string candiFusionCaseStr;
		extractFusionJuncInfoFromStr(
			detectedChrNameInt_left, detectedChrNameInt_right,
			detectedBreakPoint_left, detectedBreakPoint_right,
			detectedFusionJuncStrand_left, detectedFusionJuncStrand_right,
			detectedFusionAnchorLength_left, detectedFusionAnchorLength_right,
			candiFusionCaseStr, tmpFusionJuncStr, indexInfo);
		if(!((detectedFusionJuncStrand_left == "N")&&(detectedFusionJuncStrand_right == "N")))
		{
			cout << "incorrect strand, should be nonStranded ones" << endl;
			exit(1);
		}

		////////////   set candi fusion strands with candiFusion case ////////////////
		string detectedFusionJuncStrand_1_leftGeneAsGene1, detectedFusionJuncStrand_2_leftGeneAsGene1;
		string detectedFusionJuncStrand_1_leftGeneAsGene2, detectedFusionJuncStrand_2_leftGeneAsGene2;
		if(candiFusionCaseStr == "1,4,")
		{
			detectedFusionJuncStrand_1_leftGeneAsGene1 = "+";
			detectedFusionJuncStrand_2_leftGeneAsGene1 = "+";
			detectedFusionJuncStrand_1_leftGeneAsGene2 = "-";
			detectedFusionJuncStrand_2_leftGeneAsGene2 = "-";			
		}
		else if(candiFusionCaseStr == "2,5,")
		{
			detectedFusionJuncStrand_1_leftGeneAsGene1 = "+";
			detectedFusionJuncStrand_2_leftGeneAsGene1 = "+";
			detectedFusionJuncStrand_1_leftGeneAsGene2 = "-";
			detectedFusionJuncStrand_2_leftGeneAsGene2 = "-";
		}
		else if(candiFusionCaseStr == "7,8,")
		{
			detectedFusionJuncStrand_1_leftGeneAsGene1 = "+";
			detectedFusionJuncStrand_2_leftGeneAsGene1 = "-";
			detectedFusionJuncStrand_1_leftGeneAsGene2 = "+";
			detectedFusionJuncStrand_2_leftGeneAsGene2 = "-";
		}
		else if(candiFusionCaseStr == "10,11,")
		{
			detectedFusionJuncStrand_1_leftGeneAsGene1 = "-";
			detectedFusionJuncStrand_2_leftGeneAsGene1 = "+";
			detectedFusionJuncStrand_1_leftGeneAsGene2 = "-";
			detectedFusionJuncStrand_2_leftGeneAsGene2 = "+";
		}
		else
		{
			cout << "incorrect candi fusion cases....: " << candiFusionCaseStr << endl;
			exit(1);
		}
		//cout << "start to check anchor seq similarity " << endl;
		////////////////////////   start to check anchor seq similarity  //////////////////////////////
		string tmpFusionJuncStr_afterClassify = ""; 
		//cout << "left gene as gene1 " << endl;
		///////////////////    left gene as gene1 /////////////////////////////////////////////////////
		int detectedChrNameInt_1_leftGeneAsGene1 = detectedChrNameInt_left, 
			detectedChrNameInt_2_leftGeneAsGene1 = detectedChrNameInt_right;
		int detectedBreakPoint_1_leftGeneAsGene1 = detectedBreakPoint_left, 
			detectedBreakPoint_2_leftGeneAsGene1 = detectedBreakPoint_right;
		int detectedFusionAnchorLength_1_leftGeneAsGene1 = detectedFusionAnchorLength_left, 
			detectedFusionAnchorLength_2_leftGeneAsGene1 = detectedFusionAnchorLength_right;

		int validAnchorLength_matchThroughAtGene1_leftGeneAsGene1, 
			validAnchorLength_matchThroughAtGene2_leftGeneAsGene1, 
			penalty_matchThroughAtGene1_leftGeneAsGene1, 
			penalty_matchThroughAtGene2_leftGeneAsGene1;
 
		checkMatchTroughSeqSimilarity_fusionJunc(
			detectedFusionJuncStrand_1_leftGeneAsGene1, detectedFusionJuncStrand_2_leftGeneAsGene1, 
			detectedChrNameInt_1_leftGeneAsGene1, detectedChrNameInt_2_leftGeneAsGene1, 
			detectedBreakPoint_1_leftGeneAsGene1, detectedBreakPoint_2_leftGeneAsGene1,
			detectedFusionAnchorLength_1_leftGeneAsGene1, detectedFusionAnchorLength_2_leftGeneAsGene1, 
			toCheckAnchorLengthMax, indexInfo, 
			validAnchorLength_matchThroughAtGene1_leftGeneAsGene1, 
			validAnchorLength_matchThroughAtGene2_leftGeneAsGene1,
			penalty_matchThroughAtGene1_leftGeneAsGene1, penalty_matchThroughAtGene2_leftGeneAsGene1);			

		vector<int> validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1;
		vector<int> validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1;
		vector<int> penaltyVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1;
		vector<int> penaltyVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1;
		vector<int> alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene1;
		vector<int> alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene1;
		if(checkAlterSpliceBool)
		{	
			checkAlterSpliceAnchorSeqSimilarity_fusionJunc(
				tmpAlignInferJuncHashInfo, SJhashInfo,
				detectedFusionJuncStrand_1_leftGeneAsGene1, detectedFusionJuncStrand_2_leftGeneAsGene1,  
				detectedChrNameInt_1_leftGeneAsGene1, detectedChrNameInt_2_leftGeneAsGene1, 
				detectedBreakPoint_1_leftGeneAsGene1, detectedBreakPoint_2_leftGeneAsGene1,
				detectedFusionAnchorLength_1_leftGeneAsGene1, detectedFusionAnchorLength_2_leftGeneAsGene1, 
				toCheckAnchorLengthMax, indexInfo, 
				validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1,
				validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1,
				penaltyVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1,
				penaltyVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1,
				alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene1,
				alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene1);
		}

		//cout << "left gene as gene2 " << endl;
		///////////////////    left gene as gene2 /////////////////////////////////////////////////////
		int detectedChrNameInt_1_leftGeneAsGene2 = detectedChrNameInt_right, 
			detectedChrNameInt_2_leftGeneAsGene2 = detectedChrNameInt_left;
		int detectedBreakPoint_1_leftGeneAsGene2 = detectedBreakPoint_right, 
			detectedBreakPoint_2_leftGeneAsGene2 = detectedBreakPoint_left;
		int detectedFusionAnchorLength_1_leftGeneAsGene2 = detectedFusionAnchorLength_right, 
			detectedFusionAnchorLength_2_leftGeneAsGene2 = detectedFusionAnchorLength_left;

		int validAnchorLength_matchThroughAtGene1_leftGeneAsGene2,
			validAnchorLength_matchThroughAtGene2_leftGeneAsGene2, 
			penalty_matchThroughAtGene1_leftGeneAsGene2, 
			penalty_matchThroughAtGene2_leftGeneAsGene2;
		checkMatchTroughSeqSimilarity_fusionJunc(
			detectedFusionJuncStrand_1_leftGeneAsGene2, detectedFusionJuncStrand_2_leftGeneAsGene2, 
			detectedChrNameInt_1_leftGeneAsGene2, detectedChrNameInt_2_leftGeneAsGene2, 
			detectedBreakPoint_1_leftGeneAsGene2, detectedBreakPoint_2_leftGeneAsGene2,
			detectedFusionAnchorLength_1_leftGeneAsGene2, detectedFusionAnchorLength_2_leftGeneAsGene2, 
			toCheckAnchorLengthMax, indexInfo, 
			validAnchorLength_matchThroughAtGene1_leftGeneAsGene2, 
			validAnchorLength_matchThroughAtGene2_leftGeneAsGene2,
			penalty_matchThroughAtGene1_leftGeneAsGene2, penalty_matchThroughAtGene2_leftGeneAsGene2);	

		vector<int> validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2;
		vector<int> validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2;
		vector<int> penaltyVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2;
		vector<int> penaltyVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2;
		vector<int> alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene2;
		vector<int> alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene2;
		if(checkAlterSpliceBool)
		{	
			checkAlterSpliceAnchorSeqSimilarity_fusionJunc(
				tmpAlignInferJuncHashInfo, SJhashInfo,
				detectedFusionJuncStrand_1_leftGeneAsGene2, detectedFusionJuncStrand_2_leftGeneAsGene2,  
				detectedChrNameInt_1_leftGeneAsGene2, detectedChrNameInt_2_leftGeneAsGene2, 
				detectedBreakPoint_1_leftGeneAsGene2, detectedBreakPoint_2_leftGeneAsGene2,
				detectedFusionAnchorLength_1_leftGeneAsGene2, detectedFusionAnchorLength_2_leftGeneAsGene2, 
				toCheckAnchorLengthMax, indexInfo, 
				validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2,
				validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2,
				penaltyVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2,
				penaltyVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2,
				alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene2,
				alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene2);
		}

		tmpFusionJuncStr_afterClassify = tmpFusionJuncStr_afterClassify + tmpFusionJuncStr + "\t" 
			+ int_to_str(validAnchorLength_matchThroughAtGene1_leftGeneAsGene1) + ":"
			+ int_to_str(penalty_matchThroughAtGene1_leftGeneAsGene1) + "\t" 
			+ int_to_str(validAnchorLength_matchThroughAtGene2_leftGeneAsGene1) + ":"
			+ int_to_str(penalty_matchThroughAtGene2_leftGeneAsGene1) + "\t"
			+ int_to_str(validAnchorLength_matchThroughAtGene1_leftGeneAsGene2) + ":"
			+ int_to_str(penalty_matchThroughAtGene1_leftGeneAsGene2) + "\t" 
			+ int_to_str(validAnchorLength_matchThroughAtGene2_leftGeneAsGene2) + ":"
			+ int_to_str(penalty_matchThroughAtGene2_leftGeneAsGene2) + "\t";
		bool matchThroughSeqSimilar_bool = false;
		if((penalty_matchThroughAtGene1_leftGeneAsGene1 <= (validAnchorLength_matchThroughAtGene1_leftGeneAsGene1/4))
			||(penalty_matchThroughAtGene2_leftGeneAsGene1 <= (validAnchorLength_matchThroughAtGene2_leftGeneAsGene1/4))
			||(penalty_matchThroughAtGene1_leftGeneAsGene2 <= (validAnchorLength_matchThroughAtGene1_leftGeneAsGene2/4))
			||(penalty_matchThroughAtGene2_leftGeneAsGene2 <= (validAnchorLength_matchThroughAtGene2_leftGeneAsGene2/4)))
		{
			tmpFusionJuncStr_afterClassify += "filterOut_MT";
			matchThroughSeqSimilar_bool = true;
		}
		else
		{	
			tmpFusionJuncStr_afterClassify += "kept_MT";
			matchThroughSeqSimilar_bool = false;
		}

		bool anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene1 = false;
		bool anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene1 = false;
		bool anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene2 = false;
		bool anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene2 = false;		
		if(checkAlterSpliceBool)
		{
			if(alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene1.size() == 0)
				tmpFusionJuncStr_afterClassify += "\tnoAS\t";
			else
			{
				tmpFusionJuncStr_afterClassify += "\t";
				for(int tmp = 0; tmp < validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1.size(); tmp++)
				{	
					int tmpValidAnchorLength 
						= validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1[tmp];
					int tmpPenalty = penaltyVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene1[tmp];	
					tmpFusionJuncStr_afterClassify = tmpFusionJuncStr_afterClassify 
						+ int_to_str(alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene1[tmp]) + ":"
						+ int_to_str(tmpValidAnchorLength) + ":" + int_to_str(tmpPenalty) + ",";
					if(tmpPenalty <= (tmpValidAnchorLength/4))
						anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene1 = true;
				}
				tmpFusionJuncStr_afterClassify += "\t";
			}	

			if(alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene1.size() == 0)
				tmpFusionJuncStr_afterClassify += "\tnoAS\t";
			else
			{
				tmpFusionJuncStr_afterClassify += "\t";
				for(int tmp = 0; tmp < validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1.size(); tmp++)
				{	
					int tmpValidAnchorLength 
						= validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1[tmp];
					int tmpPenalty = penaltyVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene1[tmp];	
					tmpFusionJuncStr_afterClassify = tmpFusionJuncStr_afterClassify
						+ int_to_str(alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene1[tmp]) + ":"
						+ int_to_str(tmpValidAnchorLength) + ":" + int_to_str(tmpPenalty) + ",";
					if(tmpPenalty <= (tmpValidAnchorLength/4))
						anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene1 = true;
				}
				tmpFusionJuncStr_afterClassify += "\t";
			}	

			if(alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene2.size() == 0)
				tmpFusionJuncStr_afterClassify += "\tnoAS\t";
			else
			{
				tmpFusionJuncStr_afterClassify += "\t";
				for(int tmp = 0; tmp < validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2.size(); tmp++)
				{	
					int tmpValidAnchorLength 
						= validAnchorLengthVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2[tmp];
					int tmpPenalty = penaltyVec_alterSiteAtGene1toReplaceGene2_leftGeneAsGene2[tmp];	
					tmpFusionJuncStr_afterClassify = tmpFusionJuncStr_afterClassify 
						+ int_to_str(alterSiteVec_atGene1toReplaceGene2_leftGeneAsGene2[tmp]) + ":"
						+ int_to_str(tmpValidAnchorLength) + ":" + int_to_str(tmpPenalty) + ",";
					if(tmpPenalty <= (tmpValidAnchorLength/4))
						anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene2 = true;
				}
				tmpFusionJuncStr_afterClassify += "\t";
			}	

			if(alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene2.size() == 0)
				tmpFusionJuncStr_afterClassify += "\tnoAS\t";
			else
			{
				tmpFusionJuncStr_afterClassify += "\t";
				for(int tmp = 0; tmp < validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2.size(); tmp++)
				{	
					int tmpValidAnchorLength 
						= validAnchorLengthVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2[tmp];
					int tmpPenalty = penaltyVec_alterSiteAtGene2toReplaceGene1_leftGeneAsGene2[tmp];	
					tmpFusionJuncStr_afterClassify = tmpFusionJuncStr_afterClassify
						+ int_to_str(alterSiteVec_atGene2toReplaceGene1_leftGeneAsGene2[tmp]) + ":"
						+ int_to_str(tmpValidAnchorLength) + ":" + int_to_str(tmpPenalty) + ",";
					if(tmpPenalty <= (tmpValidAnchorLength/4))
						anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene2 = true;
				}
				tmpFusionJuncStr_afterClassify += "\t";
			}

			if(anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene1
				||anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene1
				||anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene2
				||anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene2)
				tmpFusionJuncStr_afterClassify += "filterOut_AS\t";
			else
				tmpFusionJuncStr_afterClassify += "kept_AS\t";
		}

		if(matchThroughSeqSimilar_bool 
			|| anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene1
			|| anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene1
			|| anchorSeqSimilarAlterSpliceExists_atGene1toReplaceGene2_bool_leftGeneAsGene2
			|| anchorSeqSimilarAlterSpliceExists_atGene2toReplaceGene1_bool_leftGeneAsGene2)
		{
			tmpFusionJuncStr_afterClassify += "FILTER_OUT";
			nonStrandedFusionJunc_filterOut_ofs << tmpFusionJuncStr_afterClassify << endl;
		}
		else
		{
			tmpFusionJuncStr_afterClassify += "PASS";
			nonStrandedFusionJunc_pass_ofs << tmpFusionJuncStr_afterClassify << endl;
		}
		nonStrandedFusionJunc_classified_ofs << tmpFusionJuncStr_afterClassify << endl;
	}

	nonStrandedFusionJunc_ifs.close();
	nonStrandedFusionJunc_classified_ofs.close();
	nonStrandedFusionJunc_pass_ofs.close();
	nonStrandedFusionJunc_filterOut_ofs.close();
	delete SJhashInfo;
	delete tmpAlignInferJuncHashInfo;
	delete indexInfo;
	return 0;
	return 0;
}	