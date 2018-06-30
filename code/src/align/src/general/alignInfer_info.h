// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGNINFER_INFO_H
#define ALIGNINFER_INFO_H

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
#include <sstream>

#include "../constantDefinitions.h"
#include "splice_info.h"
#include "../general/fixSingleAnchorNWDP_info.h"
#include "../phase2/spliceJunctionHash_info.h"
#include "refineSAM_info.h"

#define SHORT_ANCHOR_LEN 20
#define LONG_ANCHOR_LEN 50
#define SHORT_SJ_DISTANCE_LEN 10000
#define MID_SJ_DISTANCE_LEN 50000
#define LONG_SJ_DISTANCE_LEN 100000
#define STORED_SUPPORT_READ_NUM_MAX 5
#define DEFAULT_SEQ_SIMILARITY_PENALTY_HIGH_SUPPORT_SJ 1000;

using namespace std;

class AlignInfer_Info
{
private:
	int chrNameInt;
	int donerEndPos;
	int acceptorStartPos;

	string strand;

	int flankStringCase;
	string flankStringStr;

	string flankStringChanged;

	int pathNum_backward;
	int pathNum_forward;

	vector< vector<Jump_Code> > pathVec_backward;
	vector< vector<Jump_Code> > pathVec_forward;

	// Note: support numbers for paths are not that accurate as some alignment with 
	// small(er) number bases beside SJ may comply with muliple paths, but only the first
	// complied one got extended.
	vector< int > pathSupportNumVec_backward; 
	vector< int > pathSupportNumVec_forward;

	vector< int > pathAnchorSizeMax_backward;
	vector< int > pathAnchorSizeMax_forward;

	vector< pair<int,int> > alterDonerSpliceSitePairVec;
	vector< pair<int,int> > alterAcceptorSpliceSitePairVec;

	vector< vector< Jump_Code > > alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec;
	vector< vector< Jump_Code > > alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec;	
	vector< int > alterDonerSpliceSiteAnchorNWDPpenaltyVec;
	vector< int > alterAcceptorSpliceSiteAnchorNWDPpenaltyVec;

	int anchorSizeMax_doner;
	int anchorSizeMax_acceptor;

	vector< Jump_Code > extensionJumpCodeVec_doner;
	int extension_penalty_doner;
	vector< Jump_Code > extensionJumpCodeVec_acceptor;
	int extension_penalty_acceptor;

	int supportNum;
	int encompassingNum;

	int donerAnchorNWDPpenalty_max;
	int acceptorAnchorNWDPpenalty_max;

	int XM_min;
	int XM_max;

	// store RefineSAMinfo with that SJ, supportNumType: L, F, H
	//char SJsupportNumType;
	vector<int> indexInRefineSAMinfoVec;
public:
	AlignInfer_Info()
	{
		pathNum_backward = 0;
		pathNum_forward = 0;
		anchorSizeMax_doner = 0;
		anchorSizeMax_acceptor = 0;
		extension_penalty_doner = 0;
		extension_penalty_acceptor = 0;
		supportNum = 0;
		encompassingNum = 0;
		donerAnchorNWDPpenalty_max = 0;
		acceptorAnchorNWDPpenalty_max = 0;
		//SJsupportNumType = 'L';
	}

	string returnFlankStringChanged()
	{
		return flankStringChanged;
	}

	int returnXMmin()
	{
		return XM_min;
	}

	int returnXMmax()
	{
		return XM_max;
	}

	int returnExtension_penalty_doner()
	{
		return extension_penalty_doner;
	}

	int returnExtension_penalty_acceptor()
	{
		return extension_penalty_acceptor;
	}

	int returnAnchorSizeMax_doner()
	{
		return anchorSizeMax_doner;
	}

	int returnAnchorSizeMax_acceptor()
	{
		return anchorSizeMax_acceptor;
	}

	int returnStoredRefineSAMinfoVecSize()
	{
		return indexInRefineSAMinfoVec.size();
	}

	bool supportNumNoHigherThanStoredReadNumMaxBool()
	{
		// cout << "start to do supportNumNoHigherThanStoredReadNumMaxBool" << endl;
		// cout << "donerEndPos: " << donerEndPos << endl;
		// cout << "acceptorStartPos: " << acceptorStartPos << endl;
		// cout << "supportNum: " << supportNum << endl;
		if(supportNum > STORED_SUPPORT_READ_NUM_MAX)
			return false;
		else
			return true;
	}

	int returnRefineSAMinfoIndex(int index_refineSAMinfoVec)
	{
		return indexInRefineSAMinfoVec[index_refineSAMinfoVec];
	}

	// char returnSupportNum()
	// {
	// 	return SJsupportNumType;
	// }

	// void setNewSupportNumType(char newSupportNumType)
	// {
	// 	SJsupportNumType = newSupportNumType;
	// }

	int returnSJdistanceLevel()
	{
		int shortSJdistanceLen = SHORT_SJ_DISTANCE_LEN;
		int midSJdistanceLen = MID_SJ_DISTANCE_LEN;
		int longSJdistanceLen = LONG_SJ_DISTANCE_LEN;
		int tmpSJdistance = acceptorStartPos - donerEndPos - 1;
		if(tmpSJdistance <= shortSJdistanceLen)
			return 3;
		else if(tmpSJdistance <= midSJdistanceLen)
			return 2;
		else if(tmpSJdistance <= longSJdistanceLen)
			return 1;
		else
			return 0;
	}

	int returnAnchorSizeLevel()
	{
		int shortAnchorLen = SHORT_ANCHOR_LEN;
		int longAnchorLen= LONG_ANCHOR_LEN;
		if((anchorSizeMax_doner <= shortAnchorLen)||(anchorSizeMax_acceptor <= shortAnchorLen))
			return 0;
		else if((anchorSizeMax_doner < longAnchorLen)||(anchorSizeMax_acceptor < longAnchorLen))
			return 1;
		else
			return 2;
	}

	bool nonCanonical_bool()
	{
		return (flankStringCase == 0);
	}

	bool canonical_bool()
	{
		return (flankStringCase > 4);
	}

	bool semiCanonical_bool()
	{
		return ((flankStringCase > 0)&&(flankStringCase <= 4));
	}

	bool lowSupportNum_bool(int lowSupNum_max)
	{
		if(supportNum <= lowSupNum_max)
			return true;
		else
			return false;
	} 

	bool shortAnchor_bool(int shortAnchor_max)
	{
		if((anchorSizeMax_doner <= shortAnchor_max) || (anchorSizeMax_acceptor <= shortAnchor_max))
			return true;
		else 
			return false;
	}

	// bool validSJ()
	// {
	// 	bool canBeExtendedBool = this->canBeExtended_bool();
	// 	//bool exactTheSame2AlterSpliceSiteBool = this->exactTheSame2someAlterSpliceSite();
	// 	bool similar2someAlterSpliceSiteBool = this->similar2someAlterSpliceSite();
	// 	if(canBeExtendedBool || similar2someAlterSpliceSiteBool)
	// 		return false;
	// 	else
	// 		return true;
	// }

	bool lowSupportNumNonCanonical_bool(int lowSupNum_max)
	{
		if((supportNum <= lowSupNum_max)&&(this->nonCanonical_bool()))
			return true;
		else
			return false;
	}

	bool validSJ()
	{
		bool canBeExtendedBool = this->canBeExtended_bool();
		int lowSupNum_max = 2;
		bool lowSupportNumNonCanonicalBool = this->lowSupportNumNonCanonical_bool(2);
		bool exactTheSame2AlterSpliceSiteBool = this->exactTheSame2someAlterSpliceSite();
		//bool similar2someAlterSpliceSiteBool = this->similar2someAlterSpliceSite();
		if(canBeExtendedBool || 
			exactTheSame2AlterSpliceSiteBool
			//similar2someAlterSpliceSiteBool
			)
			return false;
		else
		{	
			if(lowSupportNumNonCanonicalBool)
				return false;
			else
				return true;
		}
	}

	bool similar2someAlterSpliceSite_doner_bool()
	{
		if(anchorSizeMax_doner <= 5)
		{
			for(int tmp = 0; tmp < alterDonerSpliceSiteAnchorNWDPpenaltyVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSiteAnchorNWDPpenalty = alterDonerSpliceSiteAnchorNWDPpenaltyVec[tmp];
				if(tmpAlterDonerSpliceSiteAnchorNWDPpenalty == 0)
					return true;
			}
		}
		else
		{		
			for(int tmp = 0; tmp < alterDonerSpliceSiteAnchorNWDPpenaltyVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSiteAnchorNWDPpenalty = alterDonerSpliceSiteAnchorNWDPpenaltyVec[tmp];
				if(tmpAlterDonerSpliceSiteAnchorNWDPpenalty <= donerAnchorNWDPpenalty_max)
					return true;
			}
		}
		return false;
	}

	bool similar2someAlterSpliceSite_acceptor_bool()
	{
		if(anchorSizeMax_acceptor <= 5)
		{
			for(int tmp = 0; tmp < alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSiteAnchorNWDPpenalty = alterAcceptorSpliceSiteAnchorNWDPpenaltyVec[tmp];
				if(tmpAlterAcceptorSpliceSiteAnchorNWDPpenalty == 0)
					return true;
			}
		}
		else
		{	
			for(int tmp = 0; tmp < alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSiteAnchorNWDPpenalty = alterAcceptorSpliceSiteAnchorNWDPpenaltyVec[tmp];
				if(tmpAlterAcceptorSpliceSiteAnchorNWDPpenalty <= acceptorAnchorNWDPpenalty_max)
					return true;
			}
		}
		return false;
	}	

	bool similar2someAlterSpliceSite()
	{
		bool similar2someAlterDonerSpliceSite_bool = this->similar2someAlterSpliceSite_doner_bool();
		bool similar2someAlterAcceptorSpliceSite_bool = this->similar2someAlterSpliceSite_acceptor_bool();
		if(similar2someAlterDonerSpliceSite_bool || similar2someAlterAcceptorSpliceSite_bool)
			return true;
		else
			return false;
	}

	bool exactTheSame2someAlterSpliceSite_doner()
	{
		for(int tmp = 0; tmp < alterDonerSpliceSiteAnchorNWDPpenaltyVec.size(); tmp++)
		{
			int tmpAlterDonerSpliceSiteAnchorNWDPpenalty = alterDonerSpliceSiteAnchorNWDPpenaltyVec[tmp];
			if(tmpAlterDonerSpliceSiteAnchorNWDPpenalty == 0)
				return true;
		}
		return false;		
	}

	bool exactTheSame2someAlterSpliceSite_acceptor()
	{
		for(int tmp = 0; tmp < alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.size(); tmp++)
		{
			int tmpAlterAcceptorSpliceSiteAnchorNWDPpenalty = alterAcceptorSpliceSiteAnchorNWDPpenaltyVec[tmp];
			if(tmpAlterAcceptorSpliceSiteAnchorNWDPpenalty == 0)
				return true;
		}
		return false;		
	}

	bool exactTheSame2someAlterSpliceSite()
	{
		bool exactTheSame2someAlterSpliceSite_doner_bool = this->exactTheSame2someAlterSpliceSite_doner();
		bool exactTheSame2someAlterSpliceSite_acceptor_bool = this->exactTheSame2someAlterSpliceSite_acceptor();
		if(exactTheSame2someAlterSpliceSite_doner_bool || exactTheSame2someAlterSpliceSite_acceptor_bool)
			return true;
		else
			return false;
	}	

	bool canBeExtended_doner_bool()
	{
		if(anchorSizeMax_doner <= 3)
		{
			if(extension_penalty_doner == 0)
				return true;
			else
				return false;
		}
		else
		{	
			if(extension_penalty_doner <= donerAnchorNWDPpenalty_max)
				return true;
			else
				return false;
		}
	}

	bool canBeExtended_acceptor_bool()
	{
		if(anchorSizeMax_acceptor <= 3)
		{
			if(extension_penalty_acceptor == 0)
				return true;
			else
				return false;
		}	
		else
		{ 
			if(extension_penalty_acceptor <= acceptorAnchorNWDPpenalty_max)
				return true;
			else
				return false;
		}
	}

	bool canBeExtended_bool()
	{
		bool extend_doner_bool = this->canBeExtended_doner_bool();
		bool extend_acceptor_bool = this->canBeExtended_acceptor_bool();
		if(extend_doner_bool || extend_acceptor_bool)
			return true;
		else
			return false;
	}

	int returnDonerEndPos()
	{
		return donerEndPos;
	} 

	int returnAcceptorStartPos()
	{
		return acceptorStartPos;
	}

	int returnSupportNum()
	{
		return supportNum;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	void checkExtension_compareAnchorSimilarity(Index_Info* indexInfo)
	{
		// generateAlterDonerSpliceSitePairVec
		if(donerEndPos-anchorSizeMax_doner < 0)
			anchorSizeMax_doner = donerEndPos;
		string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		// try direct extension
		string extensionFromAcceptorSite = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos - anchorSizeMax_doner, anchorSizeMax_doner);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDPinfo_extension_doner;
		tmpFixSingleAnchorNWDPinfo_extension_doner.doNWDP_withMismatchJumpCode(
			SJdonerAnchorStr, extensionFromAcceptorSite);
		extension_penalty_doner = tmpFixSingleAnchorNWDPinfo_extension_doner.getPenalty();
		if(extension_penalty_doner <= donerAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDPinfo_extension_doner.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_doner);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
			extensionJumpCodeVec_doner.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}		

		// generateAlterAcceptorSpliceSitePairVec
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		if(acceptorStartPos + anchorSizeMax_acceptor - 1 > tmpChromLength)
			anchorSizeMax_acceptor = tmpChromLength + 1 - acceptorStartPos;

		string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos, anchorSizeMax_acceptor);
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		// try direct extension
		string extensionFromDonerSite = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos + 1, anchorSizeMax_acceptor);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDP_extension_acceptor;
		tmpFixSingleAnchorNWDP_extension_acceptor.doNWDP_withMismatchJumpCode(
			SJacceptorAnchorStr, extensionFromDonerSite);
		extension_penalty_acceptor = tmpFixSingleAnchorNWDP_extension_acceptor.getPenalty();
		if(extension_penalty_acceptor <= acceptorAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDP_extension_acceptor.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_acceptor);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
			extensionJumpCodeVec_acceptor.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

	}

	void checkAlterSpliceSites_compareAnchorSimilarity(SJhash_Info* SJ, int offset, Index_Info* indexInfo)
	{
		// generateAlterDonerSpliceSitePairVec
		if(donerEndPos-anchorSizeMax_doner < 0)
			anchorSizeMax_doner = donerEndPos;
		string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos - anchorSizeMax_doner + 1, anchorSizeMax_doner);		
		// check alternative acceptor splice sites ......
		for(int tmpAcceptorSpliceSite = acceptorStartPos-offset; 
			tmpAcceptorSpliceSite <= acceptorStartPos+offset;
			tmpAcceptorSpliceSite++)
		{
			vector< int > tmpAlterDonerSpliceSitePosVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt,
				tmpAcceptorSpliceSite, tmpAlterDonerSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
				if(tmpAlterDonerSpliceSitePos != donerEndPos)
				{
					alterDonerSpliceSitePairVec.push_back(pair<int,int>(
						tmpAlterDonerSpliceSitePos, tmpAcceptorSpliceSite));
					int tmpNewDonerNWDPpenalty;
					vector<Jump_Code> tmpNewDonerNWDPjumpCodeVec;
					if(tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1 >= 1)
					{	
						string tmpAlterDonerSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJdonerAnchorStr, tmpAlterDonerSpliceSiteAnchorStr);
						tmpNewDonerNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewDonerNWDPpenalty <= donerAnchorNWDPpenalty_max)
						{
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewDonerNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
							tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewDonerNWDPpenalty = anchorSizeMax_doner;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
						tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterDonerSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewDonerNWDPpenalty);
					alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewDonerNWDPjumpCodeVec);
				}
			}
 		}
		
		// generateAlterAcceptorSpliceSitePairVec
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		if(acceptorStartPos + anchorSizeMax_acceptor - 1 > tmpChromLength)
			anchorSizeMax_acceptor = tmpChromLength + 1 - acceptorStartPos;

		string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos, anchorSizeMax_acceptor);
		int acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		// check alternative doner splice sites
		for(int tmpDonerSpliceSite = donerEndPos-offset;
			tmpDonerSpliceSite <= donerEndPos+offset;
			tmpDonerSpliceSite++)
		{
			vector< int > tmpAlterAcceptorSpliceSitePosVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt,
				tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
				if(tmpAlterAcceptorSpliceSitePos != acceptorStartPos)
				{
					alterAcceptorSpliceSitePairVec.push_back(pair<int,int>(
						tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePos));
					int tmpNewAcceptorNWDPpenalty;
					vector<Jump_Code> tmpNewAcceptorNWDPjumpCodeVec;
					
					if(tmpAlterAcceptorSpliceSitePos + anchorSizeMax_acceptor - 1 <= tmpChromLength)
					{	
						string tmpAlterAcceptorSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterAcceptorSpliceSitePos, anchorSizeMax_acceptor);					
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJacceptorAnchorStr, tmpAlterAcceptorSpliceSiteAnchorStr);
						tmpNewAcceptorNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewAcceptorNWDPpenalty <= acceptorAnchorNWDPpenalty_max)
						{	
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewAcceptorNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
							tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewAcceptorNWDPpenalty = anchorSizeMax_acceptor;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
						tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewAcceptorNWDPpenalty);
					alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewAcceptorNWDPjumpCodeVec);		
				}
			}
		}
	}

	void getAlterSpliceSites_compareAnchorSimilarity(SJhash_Info* SJ, int offset, Index_Info* indexInfo)
	{
		//cout << "chrNameInt: " << chrNameInt << " donerEnd: " << donerEndPos << " acceptorStart: " << acceptorStartPos << endl;
		// generateAlterDonerSpliceSitePairVec
		if(donerEndPos-anchorSizeMax_doner < 0)
			anchorSizeMax_doner = donerEndPos;
		string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
		
		// try direct extension
		string extensionFromAcceptorSite = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos - anchorSizeMax_doner, anchorSizeMax_doner);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDPinfo_extension_doner;
		tmpFixSingleAnchorNWDPinfo_extension_doner.doNWDP_withMismatchJumpCode(
			SJdonerAnchorStr, extensionFromAcceptorSite);
		extension_penalty_doner = tmpFixSingleAnchorNWDPinfo_extension_doner.getPenalty();
		if(extension_penalty_doner <= donerAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDPinfo_extension_doner.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_doner);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
			extensionJumpCodeVec_doner.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

		// check alternative acceptor splice sites ......
		for(int tmpAcceptorSpliceSite = acceptorStartPos-offset; 
			tmpAcceptorSpliceSite <= acceptorStartPos+offset;
			tmpAcceptorSpliceSite++)
		{
			vector< int > tmpAlterDonerSpliceSitePosVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt,
				tmpAcceptorSpliceSite, tmpAlterDonerSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
				if(tmpAlterDonerSpliceSitePos != donerEndPos)
				{
					alterDonerSpliceSitePairVec.push_back(pair<int,int>(
						tmpAlterDonerSpliceSitePos, tmpAcceptorSpliceSite));
					int tmpNewDonerNWDPpenalty;
					vector<Jump_Code> tmpNewDonerNWDPjumpCodeVec;
					if(tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1 >= 1)
					{	
						string tmpAlterDonerSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJdonerAnchorStr, tmpAlterDonerSpliceSiteAnchorStr);
						tmpNewDonerNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewDonerNWDPpenalty <= donerAnchorNWDPpenalty_max)
						{
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewDonerNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
							tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewDonerNWDPpenalty = anchorSizeMax_doner;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
						tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterDonerSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewDonerNWDPpenalty);
					alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewDonerNWDPjumpCodeVec);
				}
			}
 		}
	
		// generateAlterAcceptorSpliceSitePairVec
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		if(acceptorStartPos + anchorSizeMax_acceptor - 1 > tmpChromLength)
			anchorSizeMax_acceptor = tmpChromLength + 1 - acceptorStartPos;

		string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos, anchorSizeMax_acceptor);
		int acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		// try direct extension
		string extensionFromDonerSite = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos + 1, anchorSizeMax_acceptor);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDP_extension_acceptor;
		tmpFixSingleAnchorNWDP_extension_acceptor.doNWDP_withMismatchJumpCode(
			SJacceptorAnchorStr, extensionFromDonerSite);
		extension_penalty_acceptor = tmpFixSingleAnchorNWDP_extension_acceptor.getPenalty();
		if(extension_penalty_acceptor <= acceptorAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDP_extension_acceptor.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_acceptor);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
			extensionJumpCodeVec_acceptor.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

		// check alternative doner splice sites
		for(int tmpDonerSpliceSite = donerEndPos-offset;
			tmpDonerSpliceSite <= donerEndPos+offset;
			tmpDonerSpliceSite++)
		{
			vector< int > tmpAlterAcceptorSpliceSitePosVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt,
				tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
				if(tmpAlterAcceptorSpliceSitePos != acceptorStartPos)
				{
					alterAcceptorSpliceSitePairVec.push_back(pair<int,int>(
						tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePos));
					int tmpNewAcceptorNWDPpenalty;
					vector<Jump_Code> tmpNewAcceptorNWDPjumpCodeVec;
					
					if(tmpAlterAcceptorSpliceSitePos + anchorSizeMax_acceptor - 1 <= tmpChromLength)
					{	
						string tmpAlterAcceptorSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterAcceptorSpliceSitePos, anchorSizeMax_acceptor);					
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJacceptorAnchorStr, tmpAlterAcceptorSpliceSiteAnchorStr);
						tmpNewAcceptorNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewAcceptorNWDPpenalty <= acceptorAnchorNWDPpenalty_max)
						{	
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewAcceptorNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
							tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewAcceptorNWDPpenalty = anchorSizeMax_acceptor;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
						tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewAcceptorNWDPpenalty);
					alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewAcceptorNWDPjumpCodeVec);		
				}
			}
		}
	}

	void getAlterSpliceSites_compareAnchorSimilarity_allSJ_withDefaultAnchorSize(SJhash_Info* SJ, int offset, 
		int defaultAnchorSize, Index_Info* indexInfo)
	{
		anchorSizeMax_doner = defaultAnchorSize;
		anchorSizeMax_acceptor = defaultAnchorSize;
		
		//cout << "chrNameInt: " << chrNameInt << " donerEnd: " << donerEndPos << " acceptorStart: " << acceptorStartPos << endl;
		// generateAlterDonerSpliceSitePairVec
		if(donerEndPos-anchorSizeMax_doner < 0)
			anchorSizeMax_doner = donerEndPos;
		string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
		
		// try direct extension
		string extensionFromAcceptorSite = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos - anchorSizeMax_doner, anchorSizeMax_doner);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDPinfo_extension_doner;
		tmpFixSingleAnchorNWDPinfo_extension_doner.doNWDP_withMismatchJumpCode(
			SJdonerAnchorStr, extensionFromAcceptorSite);
		extension_penalty_doner = tmpFixSingleAnchorNWDPinfo_extension_doner.getPenalty();
		if(extension_penalty_doner <= donerAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDPinfo_extension_doner.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_doner);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
			extensionJumpCodeVec_doner.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

		// check alternative acceptor splice sites ......
		for(int tmpAcceptorSpliceSite = acceptorStartPos-offset; 
			tmpAcceptorSpliceSite <= acceptorStartPos+offset;
			tmpAcceptorSpliceSite++)
		{
			vector< int > tmpAlterDonerSpliceSitePosVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt,
				tmpAcceptorSpliceSite, tmpAlterDonerSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
				if(tmpAlterDonerSpliceSitePos != donerEndPos)
				{
					alterDonerSpliceSitePairVec.push_back(pair<int,int>(
						tmpAlterDonerSpliceSitePos, tmpAcceptorSpliceSite));
					int tmpNewDonerNWDPpenalty;
					vector<Jump_Code> tmpNewDonerNWDPjumpCodeVec;
					if(tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1 >= 1)
					{	
						string tmpAlterDonerSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJdonerAnchorStr, tmpAlterDonerSpliceSiteAnchorStr);
						tmpNewDonerNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewDonerNWDPpenalty <= donerAnchorNWDPpenalty_max)
						{
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewDonerNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
							tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewDonerNWDPpenalty = anchorSizeMax_doner;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
						tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterDonerSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewDonerNWDPpenalty);
					alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewDonerNWDPjumpCodeVec);
				}
			}
 		}
	
		// generateAlterAcceptorSpliceSitePairVec
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		if(acceptorStartPos + anchorSizeMax_acceptor - 1 > tmpChromLength)
			anchorSizeMax_acceptor = tmpChromLength + 1 - acceptorStartPos;

		string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos, anchorSizeMax_acceptor);
		int acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		// try direct extension
		string extensionFromDonerSite = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos + 1, anchorSizeMax_acceptor);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDP_extension_acceptor;
		tmpFixSingleAnchorNWDP_extension_acceptor.doNWDP_withMismatchJumpCode(
			SJacceptorAnchorStr, extensionFromDonerSite);
		extension_penalty_acceptor = tmpFixSingleAnchorNWDP_extension_acceptor.getPenalty();
		if(extension_penalty_acceptor <= acceptorAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDP_extension_acceptor.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_acceptor);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
			extensionJumpCodeVec_acceptor.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

		// check alternative doner splice sites
		for(int tmpDonerSpliceSite = donerEndPos-offset;
			tmpDonerSpliceSite <= donerEndPos+offset;
			tmpDonerSpliceSite++)
		{
			vector< int > tmpAlterAcceptorSpliceSitePosVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt,
				tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
				if(tmpAlterAcceptorSpliceSitePos != acceptorStartPos)
				{
					alterAcceptorSpliceSitePairVec.push_back(pair<int,int>(
						tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePos));
					int tmpNewAcceptorNWDPpenalty;
					vector<Jump_Code> tmpNewAcceptorNWDPjumpCodeVec;
					
					if(tmpAlterAcceptorSpliceSitePos + anchorSizeMax_acceptor - 1 <= tmpChromLength)
					{	
						string tmpAlterAcceptorSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterAcceptorSpliceSitePos, anchorSizeMax_acceptor);					
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJacceptorAnchorStr, tmpAlterAcceptorSpliceSiteAnchorStr);
						tmpNewAcceptorNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewAcceptorNWDPpenalty <= acceptorAnchorNWDPpenalty_max)
						{	
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewAcceptorNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
							tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewAcceptorNWDPpenalty = anchorSizeMax_acceptor;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
						tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewAcceptorNWDPpenalty);
					alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewAcceptorNWDPjumpCodeVec);		
				}
			}
		}

	}

	void getAlterSpliceSites_compareAnchorSimilarity_onlyLowSupportSJ(SJhash_Info* SJ, int offset, Index_Info* indexInfo)
	{
		// return if SJ high support
		if(supportNum > STORED_SUPPORT_READ_NUM_MAX)
		{
			if(donerEndPos-anchorSizeMax_doner < 0)
				anchorSizeMax_doner = donerEndPos;
			extension_penalty_doner = DEFAULT_SEQ_SIMILARITY_PENALTY_HIGH_SUPPORT_SJ;
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode_doner(anchorSizeMax_doner, "X");
			extensionJumpCodeVec_doner.push_back(tmpExtensionDefaultNoSimilarityJumpCode_doner);

			int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
			if(acceptorStartPos + anchorSizeMax_acceptor - 1 > tmpChromLength)
				anchorSizeMax_acceptor = tmpChromLength + 1 - acceptorStartPos;
			extension_penalty_acceptor = DEFAULT_SEQ_SIMILARITY_PENALTY_HIGH_SUPPORT_SJ;							
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode_acceptor(anchorSizeMax_acceptor, "X");
			extensionJumpCodeVec_acceptor.push_back(tmpExtensionDefaultNoSimilarityJumpCode_acceptor);
			return;
		}

		//cout << "chrNameInt: " << chrNameInt << " donerEnd: " << donerEndPos << " acceptorStart: " << acceptorStartPos << endl;
		// generateAlterDonerSpliceSitePairVec
		if(donerEndPos-anchorSizeMax_doner < 0)
			anchorSizeMax_doner = donerEndPos;
		string SJdonerAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
		
		// try direct extension
		string extensionFromAcceptorSite = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos - anchorSizeMax_doner, anchorSizeMax_doner);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDPinfo_extension_doner;
		tmpFixSingleAnchorNWDPinfo_extension_doner.doNWDP_withMismatchJumpCode(
			SJdonerAnchorStr, extensionFromAcceptorSite);
		extension_penalty_doner = tmpFixSingleAnchorNWDPinfo_extension_doner.getPenalty();
		if(extension_penalty_doner <= donerAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDPinfo_extension_doner.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_doner);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
			extensionJumpCodeVec_doner.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

		// check alternative acceptor splice sites ......
		for(int tmpAcceptorSpliceSite = acceptorStartPos-offset; 
			tmpAcceptorSpliceSite <= acceptorStartPos+offset;
			tmpAcceptorSpliceSite++)
		{
			vector< int > tmpAlterDonerSpliceSitePosVec;
			SJ->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(chrNameInt,
				tmpAcceptorSpliceSite, tmpAlterDonerSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterDonerSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterDonerSpliceSitePos = tmpAlterDonerSpliceSitePosVec[tmp];
				if(tmpAlterDonerSpliceSitePos != donerEndPos)
				{
					alterDonerSpliceSitePairVec.push_back(pair<int,int>(
						tmpAlterDonerSpliceSitePos, tmpAcceptorSpliceSite));
					int tmpNewDonerNWDPpenalty;
					vector<Jump_Code> tmpNewDonerNWDPjumpCodeVec;
					if(tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1 >= 1)
					{	
						string tmpAlterDonerSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterDonerSpliceSitePos - anchorSizeMax_doner + 1, anchorSizeMax_doner);
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJdonerAnchorStr, tmpAlterDonerSpliceSiteAnchorStr);
						tmpNewDonerNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewDonerNWDPpenalty <= donerAnchorNWDPpenalty_max)
						{
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewDonerNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
							tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewDonerNWDPpenalty = anchorSizeMax_doner;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_doner, "X");
						tmpNewDonerNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterDonerSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewDonerNWDPpenalty);
					alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewDonerNWDPjumpCodeVec);
				}
			}
 		}
	
		// generateAlterAcceptorSpliceSitePairVec
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		if(acceptorStartPos + anchorSizeMax_acceptor - 1 > tmpChromLength)
			anchorSizeMax_acceptor = tmpChromLength + 1 - acceptorStartPos;

		string SJacceptorAnchorStr = indexInfo->returnChromStrSubstr(
			chrNameInt, acceptorStartPos, anchorSizeMax_acceptor);
		int acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		// try direct extension
		string extensionFromDonerSite = indexInfo->returnChromStrSubstr(
			chrNameInt, donerEndPos + 1, anchorSizeMax_acceptor);
		FixSingleAnchor_NWDP_Info tmpFixSingleAnchorNWDP_extension_acceptor;
		tmpFixSingleAnchorNWDP_extension_acceptor.doNWDP_withMismatchJumpCode(
			SJacceptorAnchorStr, extensionFromDonerSite);
		extension_penalty_acceptor = tmpFixSingleAnchorNWDP_extension_acceptor.getPenalty();
		if(extension_penalty_acceptor <= acceptorAnchorNWDPpenalty_max)
		{
			tmpFixSingleAnchorNWDP_extension_acceptor.copyJumpCodeVec2TargetVec(extensionJumpCodeVec_acceptor);
		}
		else
		{
			Jump_Code tmpExtensionDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
			extensionJumpCodeVec_acceptor.push_back(tmpExtensionDefaultNoSimilarityJumpCode);
		}

		// check alternative doner splice sites
		for(int tmpDonerSpliceSite = donerEndPos-offset;
			tmpDonerSpliceSite <= donerEndPos+offset;
			tmpDonerSpliceSite++)
		{
			vector< int > tmpAlterAcceptorSpliceSitePosVec;
			SJ->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(chrNameInt,
				tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePosVec);
			for(int tmp = 0; tmp < tmpAlterAcceptorSpliceSitePosVec.size(); tmp++)
			{
				int tmpAlterAcceptorSpliceSitePos = tmpAlterAcceptorSpliceSitePosVec[tmp];
				if(tmpAlterAcceptorSpliceSitePos != acceptorStartPos)
				{
					alterAcceptorSpliceSitePairVec.push_back(pair<int,int>(
						tmpDonerSpliceSite, tmpAlterAcceptorSpliceSitePos));
					int tmpNewAcceptorNWDPpenalty;
					vector<Jump_Code> tmpNewAcceptorNWDPjumpCodeVec;
					
					if(tmpAlterAcceptorSpliceSitePos + anchorSizeMax_acceptor - 1 <= tmpChromLength)
					{	
						string tmpAlterAcceptorSpliceSiteAnchorStr = indexInfo->returnChromStrSubstr(
							chrNameInt, tmpAlterAcceptorSpliceSitePos, anchorSizeMax_acceptor);					
						FixSingleAnchor_NWDP_Info fixSingleAnchorNWDPinfo;
						fixSingleAnchorNWDPinfo.doNWDP_withMismatchJumpCode(
							SJacceptorAnchorStr, tmpAlterAcceptorSpliceSiteAnchorStr);
						tmpNewAcceptorNWDPpenalty = fixSingleAnchorNWDPinfo.getPenalty();
						if(tmpNewAcceptorNWDPpenalty <= acceptorAnchorNWDPpenalty_max)
						{	
							fixSingleAnchorNWDPinfo.copyJumpCodeVec2TargetVec(tmpNewAcceptorNWDPjumpCodeVec);
						}
						else
						{
							Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
							tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);
						}
					}
					else
					{
						tmpNewAcceptorNWDPpenalty = anchorSizeMax_acceptor;
						Jump_Code tmpDefaultNoSimilarityJumpCode(anchorSizeMax_acceptor, "X");
						tmpNewAcceptorNWDPjumpCodeVec.push_back(tmpDefaultNoSimilarityJumpCode);						
					}
					alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.push_back(tmpNewAcceptorNWDPpenalty);
					alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(tmpNewAcceptorNWDPjumpCodeVec);		
				}
			}
		}
	}

	int selectBestCandiInferedPath_minMismatch(
		vector< vector<Jump_Code> >& tmpCandiInferedPathJumpCodeVecVec, 
		vector< vector<int> >& tmpCandiInferedPathMismatchPosVecVec, 
		vector< vector<char> >& tmpCandiInferedPathMismatchCharVecVec)
	{
		int mismatch_num_min = 10000;
		int index = 0;
		for(int tmp = 0; tmp < tmpCandiInferedPathMismatchPosVecVec.size(); tmp++)
		{
			int tmpCandiInferedPathMismatchNum = (tmpCandiInferedPathMismatchPosVecVec[tmp]).size();
			if(tmpCandiInferedPathMismatchNum < mismatch_num_min)
			{
				index = tmp;
			}
		}
		return index;
	}

	int return_readBaseNum_path_backward(int index_pathVec_backward)
	{
		int path_backward_jumpCodeVecSize = pathVec_backward[index_pathVec_backward].size();
		int tmpReadBaseNum = 0;
		for(int tmp = 0; tmp < path_backward_jumpCodeVecSize; tmp++)
		{
			int path_backward_jumpCodeLength = (pathVec_backward[index_pathVec_backward])[tmp].len;
			string path_backward_jumpCodeType = (pathVec_backward[index_pathVec_backward])[tmp].type;
			if((path_backward_jumpCodeType == "M")||(path_backward_jumpCodeType == "I"))
			{
				tmpReadBaseNum += path_backward_jumpCodeLength;
			}
		}
		return tmpReadBaseNum;
	}

	int return_readBaseNum_path_forward(int index_pathVec_forward)
	{
		int path_forward_jumpCodeVecSize = pathVec_forward[index_pathVec_forward].size();
		int tmpReadBaseNum = 0;
		for(int tmp = 0; tmp < path_forward_jumpCodeVecSize; tmp++)
		{
			int path_forward_jumpCodeLength = (pathVec_forward[index_pathVec_forward])[tmp].len;
			string path_forward_jumpCodeType = (pathVec_forward[index_pathVec_forward])[tmp].type;
			if((path_forward_jumpCodeType == "M")||(path_forward_jumpCodeType == "I"))
			{
				tmpReadBaseNum += path_forward_jumpCodeLength;
			}
		}
		return tmpReadBaseNum;
	}


	bool generateCandiInferedUnfixedHeadPath(
		vector<Jump_Code>& tmpCandiInferedPathJumpCodeVec, 
		vector<int>& tmpCandiInferedPathMismatchPosVec, 
		vector<char>& tmpCandiInferedPathMismatchCharVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedHeadLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr,
		int index_pathVec_backward)
	{
		//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
		//	<< "generateCandiInferedUnfixedHeadPath starts ......" << endl;

		int max_mismatch = unfixedHeadLength / LengthOfSeqPerMismatchAllowed_REMAPPING + 1;
		bool generateCandiInferedUnfixedHeadPath_bool = true;
		int tmpMismatchNum = 0;

		vector<Jump_Code> tmpReverseCandiJumpCodeVec;
		vector<int> tmpCandiMismatchPosVec;
		vector<char> tmpCandiMismatchCharVec;
		//cout << "start to check acceptor side seq consistent or not ..." << endl;
		// check acceptor side seq conosistent or not
		if(tmpSpliceAcceptorStartPosInRead < unfixedHeadLength)
		{	
			for(int tmpLocInRead = tmpSpliceAcceptorStartPosInRead; 
				tmpLocInRead <= unfixedHeadLength;
				tmpLocInRead ++)
			{
				char tmpCharInRead = readSeqWithDirection.at(tmpLocInRead-1);
				char tmpCharInChrom = indexInfo->returnOneBaseCharInGenome(
					tmpChromNameInt, tmpSpliceAcceptorStartPosInChr + tmpLocInRead - tmpSpliceAcceptorStartPosInRead);
				if(tmpCharInRead != tmpCharInChrom)
				{
					tmpMismatchNum ++;
					tmpCandiInferedPathMismatchPosVec.push_back(tmpLocInRead);
					tmpCandiInferedPathMismatchCharVec.push_back(tmpCharInRead);
					if(tmpMismatchNum > max_mismatch)
						return false;
				}			
			}
		}
		Jump_Code acceptorMatchJumpCode(unfixedHeadLength - tmpSpliceAcceptorStartPosInRead + 1, "M");
		Jump_Code spliceJunctionJumpCode(tmpSpliceAcceptorStartPosInChr - tmpSpliceDonerEndPosInChr - 1, "N");
		tmpReverseCandiJumpCodeVec.push_back(acceptorMatchJumpCode);	
		tmpReverseCandiJumpCodeVec.push_back(spliceJunctionJumpCode);	
		//cout << "start to check doner side path ......" << endl;
		// check doner side path
		int path_backward_jumpCodeVecSize = pathVec_backward[index_pathVec_backward].size();
		int tmp_JumpCode_startLoc_inRead = tmpSpliceAcceptorStartPosInRead;// = 0;
		int tmp_JumpCode_endLoc_inRead;
		int tmp_JumpCode_startLoc_inChr = tmpSpliceDonerEndPosInChr + 1;
		int tmp_JumpCode_endLoc_inChr;
		//cout << "path_backward_jumpCodeVecSize: " << path_backward_jumpCodeVecSize << endl;
		for(int tmp = 0; tmp < path_backward_jumpCodeVecSize; tmp++)
		{
			int path_backward_jumpCodeLength = (pathVec_backward[index_pathVec_backward])[tmp].len;
			string path_backward_jumpCodeType = (pathVec_backward[index_pathVec_backward])[tmp].type;
			//cout << "path_backward_jumpCodeLength: " << path_backward_jumpCodeLength << endl;
			//cout << "path_backward_jumpCodeType: " << path_backward_jumpCodeType << endl;

			if(path_backward_jumpCodeType == "M")
			{
				tmp_JumpCode_endLoc_inRead = tmp_JumpCode_startLoc_inRead - 1;
				tmp_JumpCode_startLoc_inRead = tmp_JumpCode_endLoc_inRead - path_backward_jumpCodeLength + 1;
				tmp_JumpCode_endLoc_inChr = tmp_JumpCode_startLoc_inChr - 1;
				tmp_JumpCode_startLoc_inChr = tmp_JumpCode_endLoc_inChr - path_backward_jumpCodeLength + 1;
				//cout << "tmp_JumpCode_endLoc_inRead: " << tmp_JumpCode_endLoc_inRead << endl;
				//cout << "tmp_JumpCode_endLoc_inChr: " << tmp_JumpCode_endLoc_inChr << endl;
				//cout << "tmp_JumpCode_startLoc_inRead: " << tmp_JumpCode_startLoc_inRead << endl;
				//cout << "tmp_JumpCode_startLoc_inChr: " << tmp_JumpCode_startLoc_inChr << endl;

				if(tmp_JumpCode_startLoc_inRead > 1)
				{
					for(int tmpLocInRead = tmp_JumpCode_startLoc_inRead; 
						tmpLocInRead <= tmp_JumpCode_endLoc_inRead; tmpLocInRead ++)
					{
						char tmpCharInRead = readSeqWithDirection.at(tmpLocInRead-1);
						char tmpCharInChrom = indexInfo->returnOneBaseCharInGenome(
							tmpChromNameInt, tmp_JumpCode_startLoc_inChr + tmpLocInRead - tmp_JumpCode_startLoc_inRead);
						if(tmpCharInRead != tmpCharInChrom)
						{
							tmpMismatchNum ++;
							tmpCandiMismatchPosVec.push_back(tmpLocInRead);
							tmpCandiMismatchCharVec.push_back(tmpCharInRead);
							if(tmpMismatchNum > max_mismatch)
							{
								return false;
							}
						}
					}
					tmpReverseCandiJumpCodeVec.push_back((pathVec_backward[index_pathVec_backward])[tmp]);
				}
				else
				{
					for(int tmpLocInRead = 1; 
						tmpLocInRead <= tmp_JumpCode_endLoc_inRead; tmpLocInRead ++)
					{
						char tmpCharInRead = readSeqWithDirection.at(tmpLocInRead-1);
						char tmpCharInChrom = indexInfo->returnOneBaseCharInGenome(
							tmpChromNameInt, tmp_JumpCode_startLoc_inChr + tmpLocInRead - tmp_JumpCode_startLoc_inRead);
						if(tmpCharInRead != tmpCharInChrom)
						{
							tmpMismatchNum ++;
							tmpCandiMismatchPosVec.push_back(tmpLocInRead);
							tmpCandiMismatchCharVec.push_back(tmpCharInRead);
							if(tmpMismatchNum > max_mismatch)
							{
								return false;
							}
						}
					}
					int tmpLastMatchLength = tmp_JumpCode_endLoc_inRead;
					Jump_Code tmpLastMatchJumpCode(tmpLastMatchLength, "M");
					tmpReverseCandiJumpCodeVec.push_back(tmpLastMatchJumpCode);
					break;
					//return true;				
				}
			}
			else if(path_backward_jumpCodeType == "N")
			{
				tmp_JumpCode_endLoc_inChr = tmp_JumpCode_startLoc_inChr - 1;
				tmp_JumpCode_startLoc_inChr = tmp_JumpCode_endLoc_inChr - path_backward_jumpCodeLength + 1;
				tmpReverseCandiJumpCodeVec.push_back((pathVec_backward[index_pathVec_backward])[tmp]);			
			}
			else if(path_backward_jumpCodeType == "I")
			{
				tmp_JumpCode_endLoc_inRead = tmp_JumpCode_startLoc_inRead - 1;
				tmp_JumpCode_startLoc_inRead = tmp_JumpCode_endLoc_inRead - path_backward_jumpCodeLength + 1;
				if(tmp_JumpCode_startLoc_inRead <= 1)
				{
					int finalInsertionLength = tmp_JumpCode_endLoc_inRead;
					Jump_Code finalInsertionJumpCode(finalInsertionLength, "I");
					tmpReverseCandiJumpCodeVec.push_back(finalInsertionJumpCode);
					break;
				}
				else
				{
					tmpReverseCandiJumpCodeVec.push_back((pathVec_backward[index_pathVec_backward])[tmp]);
				}		
			}
			else if(path_backward_jumpCodeType == "D")
			{
				tmp_JumpCode_endLoc_inChr = tmp_JumpCode_startLoc_inChr - 1;
				tmp_JumpCode_startLoc_inChr = tmp_JumpCode_endLoc_inChr - path_backward_jumpCodeLength + 1;				
				tmpReverseCandiJumpCodeVec.push_back((pathVec_backward[index_pathVec_backward])[tmp]);	
			}
			else
			{	
				cout << "invalid jumpcode type in path_backward : " << path_backward_jumpCodeType << endl;
				exit(1);
			}				
		}

		//cout << "tmp_JumpCode_startLoc_inRead: " << tmp_JumpCode_startLoc_inRead << endl;
		int softclippingLength = tmp_JumpCode_startLoc_inRead - 1;
		if(softclippingLength > 0)
		{
			Jump_Code softClipJumpCode(softclippingLength, "S");
			tmpCandiInferedPathJumpCodeVec.push_back(softClipJumpCode);
		}
		for(int tmp = tmpReverseCandiJumpCodeVec.size()-1; tmp >= 0; tmp--)
		{
			//cout << "tmpReverseCandiJumpCode.type: " << tmpReverseCandiJumpCodeVec[tmp].type << endl;
			//cout << "tmpReverseCandiJumpCode.len: " << tmpReverseCandiJumpCodeVec[tmp].len << endl; 
			tmpCandiInferedPathJumpCodeVec.push_back(tmpReverseCandiJumpCodeVec[tmp]);
		}
		//cout << "start to push back mismatch pos and char" << endl;
		for(int tmp = 0; tmp < tmpCandiMismatchPosVec.size(); tmp++)
		{
			//cout << "tmp_index in mismatchVec: " << tmp << endl;
			//cout << "mismatchPso: " << tmpCandiMismatchPosVec[tmp] << endl;
			//cout << "mismatchChar: " << tmpCandiMismatchCharVec[tmp] << endl;
			tmpCandiInferedPathMismatchPosVec.push_back(tmpCandiMismatchPosVec[tmp]);
			tmpCandiInferedPathMismatchCharVec.push_back(tmpCandiMismatchCharVec[tmp]);
		}
		//cout << "generateCandiInferedUnfixedHeadPath_bool: " << generateCandiInferedUnfixedHeadPath_bool << endl;
		//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		return generateCandiInferedUnfixedHeadPath_bool;
	}


	bool generateCandiInferedUnfixedTailPath(
		vector<Jump_Code>& tmpCandiInferedPathJumpCodeVec, 
		vector<int>& tmpCandiInferedPathMismatchPosVec, 
		vector<char>& tmpCandiInferedPathMismatchCharVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedTailLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr,
		int index_pathVec_forward)
	{
		int readSeqWithDirectionLength = readSeqWithDirection.length();
		int midPartEndLocInRead = readSeqWithDirectionLength - unfixedTailLength;

		int max_mismatch = unfixedTailLength / LengthOfSeqPerMismatchAllowed_REMAPPING + 1;
		bool generateCandiInferedUnfixedTailPath_bool = true;
		int tmpMismatchNum = 0;

		//vector<Jump_Code> tmpForwardCandiJumpCodeVec;
		vector<int> tmpCandiMismatchPosVec;
		vector<char> tmpCandiMismatchCharVec;
		//cout << "tmpSpliceDonerEndPosInRead: " << tmpSpliceDonerEndPosInRead << endl;
		//cout << "midPartEndLocInRead: " << midPartEndLocInRead << endl;
		// check doner side seq conosistent or not
		if(tmpSpliceDonerEndPosInRead > midPartEndLocInRead)
		{	
			for(int tmpLocInRead = midPartEndLocInRead + 1; 
				tmpLocInRead <= tmpSpliceDonerEndPosInRead;
				tmpLocInRead ++)
			{
				char tmpCharInRead = readSeqWithDirection.at(tmpLocInRead-1);
				char tmpCharInChrom = indexInfo->returnOneBaseCharInGenome(
					tmpChromNameInt, tmpSpliceDonerEndPosInChr - (tmpSpliceDonerEndPosInRead - tmpLocInRead));
				if(tmpCharInRead != tmpCharInChrom)
				{
					tmpMismatchNum ++;
					tmpCandiInferedPathMismatchPosVec.push_back(tmpLocInRead);
					tmpCandiInferedPathMismatchCharVec.push_back(tmpCharInRead);
					if(tmpMismatchNum > max_mismatch)
						return false;
				}			
			}
		}
		//cout << "finishs doner side seq checking ..." << endl;
		Jump_Code donerMatchJumpCode(tmpSpliceDonerEndPosInRead - midPartEndLocInRead, "M");
		Jump_Code spliceJunctionJumpCode(tmpSpliceAcceptorStartPosInChr - tmpSpliceDonerEndPosInChr - 1, "N");
		tmpCandiInferedPathJumpCodeVec.push_back(donerMatchJumpCode);	
		tmpCandiInferedPathJumpCodeVec.push_back(spliceJunctionJumpCode);	


		// check acceptor side path
		int path_forward_jumpCodeVecSize = pathVec_forward[index_pathVec_forward].size();
		int tmp_JumpCode_startLoc_inRead;
		int tmp_JumpCode_endLoc_inRead = tmpSpliceDonerEndPosInRead;// = 0;
		int tmp_JumpCode_startLoc_inChr;
		int tmp_JumpCode_endLoc_inChr = tmpSpliceAcceptorStartPosInChr - 1;

		for(int tmp = 0; tmp < path_forward_jumpCodeVecSize; tmp++)
		{
			int path_forward_jumpCodeLength = (pathVec_forward[index_pathVec_forward])[tmp].len;
			string path_forward_jumpCodeType = (pathVec_forward[index_pathVec_forward])[tmp].type;
			if(path_forward_jumpCodeType == "M")
			{
				tmp_JumpCode_startLoc_inRead = tmp_JumpCode_endLoc_inRead + 1;
				tmp_JumpCode_endLoc_inRead = tmp_JumpCode_startLoc_inRead + path_forward_jumpCodeLength - 1;
				tmp_JumpCode_startLoc_inChr = tmp_JumpCode_endLoc_inChr + 1;
				tmp_JumpCode_endLoc_inChr = tmp_JumpCode_startLoc_inChr + path_forward_jumpCodeLength - 1;

				if(tmp_JumpCode_endLoc_inRead < readSeqWithDirectionLength)
				{
					for(int tmpLocInRead = tmp_JumpCode_startLoc_inRead; 
						tmpLocInRead <= tmp_JumpCode_endLoc_inRead; tmpLocInRead ++)
					{
						char tmpCharInRead = readSeqWithDirection.at(tmpLocInRead-1);
						char tmpCharInChrom = indexInfo->returnOneBaseCharInGenome(
							tmpChromNameInt, tmp_JumpCode_startLoc_inChr + tmpLocInRead - tmp_JumpCode_startLoc_inRead);
						if(tmpCharInRead != tmpCharInChrom)
						{
							tmpMismatchNum ++;
							tmpCandiMismatchPosVec.push_back(tmpLocInRead);
							tmpCandiMismatchCharVec.push_back(tmpCharInRead);
							if(tmpMismatchNum > max_mismatch)
							{
								return false;
							}
						}
					}
					tmpCandiInferedPathJumpCodeVec.push_back((pathVec_forward[index_pathVec_forward])[tmp]);
				}
				else
				{
					for(int tmpLocInRead = tmp_JumpCode_startLoc_inRead; 
						tmpLocInRead <= readSeqWithDirectionLength; tmpLocInRead ++)
					{
						char tmpCharInRead = readSeqWithDirection.at(tmpLocInRead-1);
						char tmpCharInChrom = indexInfo->returnOneBaseCharInGenome(
							tmpChromNameInt, tmp_JumpCode_startLoc_inChr + tmpLocInRead - tmp_JumpCode_startLoc_inRead);
						if(tmpCharInRead != tmpCharInChrom)
						{
							tmpMismatchNum ++;
							tmpCandiMismatchPosVec.push_back(tmpLocInRead);
							tmpCandiMismatchCharVec.push_back(tmpCharInRead);
							if(tmpMismatchNum > max_mismatch)
							{
								return false;
							}
						}
					}
					int tmpLastMatchLength = readSeqWithDirectionLength - tmp_JumpCode_startLoc_inRead + 1;
					Jump_Code tmpLastMatchJumpCode(tmpLastMatchLength, "M");
					tmpCandiInferedPathJumpCodeVec.push_back(tmpLastMatchJumpCode);
					break;				
				}
			}
			else if(path_forward_jumpCodeType == "N")
			{
				tmp_JumpCode_startLoc_inChr = tmp_JumpCode_endLoc_inChr + 1;
				tmp_JumpCode_endLoc_inChr = tmp_JumpCode_startLoc_inChr + path_forward_jumpCodeLength - 1;
				tmpCandiInferedPathJumpCodeVec.push_back((pathVec_forward[index_pathVec_forward])[tmp]);		
			}
			else if(path_forward_jumpCodeType == "I")
			{
				tmp_JumpCode_startLoc_inRead = tmp_JumpCode_endLoc_inRead + 1;
				tmp_JumpCode_endLoc_inRead = tmp_JumpCode_startLoc_inRead + path_forward_jumpCodeLength - 1;
				if(tmp_JumpCode_endLoc_inRead >= readSeqWithDirectionLength)
				{
					int finalInsertionLength = readSeqWithDirectionLength - tmp_JumpCode_startLoc_inRead + 1;
					Jump_Code finalInsertionJumpCode(finalInsertionLength, "I");
					tmpCandiInferedPathJumpCodeVec.push_back(finalInsertionJumpCode);
					break;
				}
				else
				{
					tmpCandiInferedPathJumpCodeVec.push_back((pathVec_forward[index_pathVec_forward])[tmp]);
				}
			}
			else if(path_forward_jumpCodeType == "D")
			{
				tmp_JumpCode_startLoc_inChr = tmp_JumpCode_endLoc_inChr + 1;
				tmp_JumpCode_endLoc_inChr = tmp_JumpCode_startLoc_inChr + path_forward_jumpCodeLength - 1;
				tmpCandiInferedPathJumpCodeVec.push_back((pathVec_forward[index_pathVec_forward])[tmp]);
			}
			else
			{	
				cout << "invalid jumpcode type in path_forward : " << path_forward_jumpCodeType << endl;
				exit(1);
			}				
		}

		int softclippingLength = readSeqWithDirectionLength - tmp_JumpCode_endLoc_inRead;
		if(softclippingLength > 0)
		{
			Jump_Code softClipJumpCode(softclippingLength, "S");
			tmpCandiInferedPathJumpCodeVec.push_back(softClipJumpCode);
		}
		// for(int tmp = tmpReverseCandiJumpCodeVec.size()-1; tmp >= 0; tmp--)
		// {
		// 	tmpCandiInferedPathJumpCodeVec.push_back(tmpReverseCandiJumpCodeVec[tmp]);
		// }

		for(int tmp = 0; tmp < tmpCandiMismatchPosVec.size(); tmp++)
		{
			tmpCandiInferedPathMismatchPosVec.push_back(tmpCandiMismatchPosVec[tmp]);
			tmpCandiInferedPathMismatchCharVec.push_back(tmpCandiMismatchCharVec[tmp]);
		}
		return generateCandiInferedUnfixedTailPath_bool;
	}


	void generateCandiInferedUnfixedHeadPathVec(
		vector< vector<Jump_Code> >& tmpCandiInferedPathJumpCodeVecVec, 
		vector< vector<int> >& tmpCandiInferedPathMismatchPosVecVec, 
		vector< vector<char> >& tmpCandiInferedPathMismatchCharVecVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedHeadLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr)
	{
		//cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl
		//	<< "generateCandiInferedUnfixedHeadPathVec starts ......" << endl;
		for(int tmp = 0; tmp < pathVec_backward.size(); tmp++)
		{
			//cout << "index_pathVec_backward: " << tmp << endl;
			vector<Jump_Code> tmpCandiInferedPathJumpCodeVec;
			vector<int> tmpCandiInferedPathMismatchPosVec;
			vector<char> tmpCandiInferedPathMismatchCharVec;
			bool candiPathGenerated_bool = this->generateCandiInferedUnfixedHeadPath(
				tmpCandiInferedPathJumpCodeVec,
				tmpCandiInferedPathMismatchPosVec,
				tmpCandiInferedPathMismatchCharVec,
				indexInfo, tmpChromNameInt, unfixedHeadLength, readSeqWithDirection,
				tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
				tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr, tmp);
			//cout << "candiPathGenerated_bool: " << candiPathGenerated_bool << endl;
			if(candiPathGenerated_bool)
			{
				// for(int tmp = 0; tmp < tmpCandiInferedPathJumpCodeVec.size(); tmp++)
				// {
				// 	cout << "tmpCandiPathJumpCodeType: " << tmpCandiInferedPathJumpCodeVec[tmp].type << endl;
				// 	cout << "tmpCandiPathJumpCodeLength: " << tmpCandiInferedPathJumpCodeVec[tmp].len << endl;
				// }
				tmpCandiInferedPathJumpCodeVecVec.push_back(tmpCandiInferedPathJumpCodeVec);
				tmpCandiInferedPathMismatchPosVecVec.push_back(tmpCandiInferedPathMismatchPosVec);
				tmpCandiInferedPathMismatchCharVecVec.push_back(tmpCandiInferedPathMismatchCharVec);
			}
		}
		//cout << "generateCandiInferedUnfixedHeadPathVec ends ......" << endl
		//	<< "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;		
	}

	
	void generateCandiInferedUnfixedTailPathVec(
		vector< vector<Jump_Code> >& tmpCandiInferedPathJumpCodeVecVec, 
		vector< vector<int> >& tmpCandiInferedPathMismatchPosVecVec, 
		vector< vector<char> >& tmpCandiInferedPathMismatchCharVecVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedTailLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr)
	{
		//cout << "pathVec_forward.size(): " << pathVec_forward.size() << endl;
		for(int tmp = 0; tmp < pathVec_forward.size(); tmp++)
		{
			vector<Jump_Code> tmpCandiInferedPathJumpCodeVec;
			vector<int> tmpCandiInferedPathMismatchPosVec;
			vector<char> tmpCandiInferedPathMismatchCharVec;
			bool candiPathGenerated_bool = this->generateCandiInferedUnfixedTailPath(
				tmpCandiInferedPathJumpCodeVec,
				tmpCandiInferedPathMismatchPosVec,
				tmpCandiInferedPathMismatchCharVec,
				indexInfo, tmpChromNameInt, unfixedTailLength, readSeqWithDirection,
				tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
				tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr, tmp);
			if(candiPathGenerated_bool)
			{
				tmpCandiInferedPathJumpCodeVecVec.push_back(tmpCandiInferedPathJumpCodeVec);
				tmpCandiInferedPathMismatchPosVecVec.push_back(tmpCandiInferedPathMismatchPosVec);
				tmpCandiInferedPathMismatchCharVecVec.push_back(tmpCandiInferedPathMismatchCharVec);
			}
		}
	}



	bool generateInferedUnfixedHeadPath_alignInfer(
		vector<Jump_Code>& tmpInferedPathJumpCodeVec, 
		vector<int>& tmpInferedPathMismatchPosVec, 
		vector<char>& tmpInferedPathMismatchCharVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedHeadLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr)
	{
		//cout << " :::::::::::::::::::::::::::::::::::::::::" << endl 
		//	<< "generateInferedUnfixedHeadPath_alignInfer starts" << endl;

		vector< vector<Jump_Code> > tmpCandiInferedPathJumpCodeVecVec;
		vector< vector<int> > tmpCandiInferedPathMismatchPosVecVec;
		vector< vector<char> > tmpCandiInferedPathMismatchCharVecVec; 		
		this->generateCandiInferedUnfixedHeadPathVec(
			tmpCandiInferedPathJumpCodeVecVec, 
			tmpCandiInferedPathMismatchPosVecVec, 
			tmpCandiInferedPathMismatchCharVecVec, 
			indexInfo, tmpChromNameInt, unfixedHeadLength, readSeqWithDirection,
			tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
			tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
		
		//cout << "tmpCandiInferedPathMismatchCharVecVec.size(): " 
		//	<< tmpCandiInferedPathMismatchCharVecVec.size() << endl
		//	 << ":::::::::::::::::::::::::::::::::::::" << endl;
		if(tmpCandiInferedPathMismatchCharVecVec.size() == 0)
			return false;
		int index_bestSelectedInferedPath = this->selectBestCandiInferedPath_minMismatch(
			tmpCandiInferedPathJumpCodeVecVec, tmpCandiInferedPathMismatchPosVecVec, 
			tmpCandiInferedPathMismatchCharVecVec);
		//cout << "index_bestSelectedInferedPath: " << index_bestSelectedInferedPath << endl;
		for(int tmp = 0; 
			tmp < tmpCandiInferedPathJumpCodeVecVec[index_bestSelectedInferedPath].size();
			tmp++)
		{
			tmpInferedPathJumpCodeVec.push_back(
				(tmpCandiInferedPathJumpCodeVecVec[index_bestSelectedInferedPath])[tmp]);
		}
		for(int tmp = 0; 
			tmp < tmpCandiInferedPathMismatchPosVecVec[index_bestSelectedInferedPath].size();
			tmp++)
		{
			tmpInferedPathMismatchPosVec.push_back(
				(tmpCandiInferedPathMismatchPosVecVec[index_bestSelectedInferedPath])[tmp]);
			tmpInferedPathMismatchCharVec.push_back(
				(tmpCandiInferedPathMismatchCharVecVec[index_bestSelectedInferedPath])[tmp]);
		}
		return true;
	}


	bool generateInferedUnfixedTailPath_alignInfer(
		vector<Jump_Code>& tmpInferedPathJumpCodeVec, 
		vector<int>& tmpInferedPathMismatchPosVec, 
		vector<char>& tmpInferedPathMismatchCharVec, 
		Index_Info* indexInfo, int tmpChromNameInt, 
		int unfixedTailLength, const string& readSeqWithDirection,
		int tmpSpliceDonerEndPosInRead, int tmpSpliceAcceptorStartPosInRead,
		int	tmpSpliceDonerEndPosInChr, int tmpSpliceAcceptorStartPosInChr)
	{
		vector< vector<Jump_Code> > tmpCandiInferedPathJumpCodeVecVec;
		vector< vector<int> > tmpCandiInferedPathMismatchPosVecVec;
		vector< vector<char> > tmpCandiInferedPathMismatchCharVecVec; 		
		this->generateCandiInferedUnfixedTailPathVec(
			tmpCandiInferedPathJumpCodeVecVec, 
			tmpCandiInferedPathMismatchPosVecVec, 
			tmpCandiInferedPathMismatchCharVecVec, 
			indexInfo, tmpChromNameInt, unfixedTailLength, readSeqWithDirection,
			tmpSpliceDonerEndPosInRead, tmpSpliceAcceptorStartPosInRead,
			tmpSpliceDonerEndPosInChr, tmpSpliceAcceptorStartPosInChr);
		//cout << "tmpCandiInferedPathMismatchCharVecVec.size(): " << tmpCandiInferedPathMismatchCharVecVec.size() << endl;
		if(tmpCandiInferedPathMismatchCharVecVec.size() == 0)
			return false;
		int index_bestSelectedInferedPath = this->selectBestCandiInferedPath_minMismatch(
			tmpCandiInferedPathJumpCodeVecVec, tmpCandiInferedPathMismatchPosVecVec, 
			tmpCandiInferedPathMismatchCharVecVec);
		//cout << "index_bestSelectedInferedPath: " << index_bestSelectedInferedPath << endl;
		for(int tmp = 0; 
			tmp < tmpCandiInferedPathJumpCodeVecVec[index_bestSelectedInferedPath].size();
			tmp++)
		{
			tmpInferedPathJumpCodeVec.push_back(
				(tmpCandiInferedPathJumpCodeVecVec[index_bestSelectedInferedPath])[tmp]);
		}
		for(int tmp = 0; 
			tmp < tmpCandiInferedPathMismatchPosVecVec[index_bestSelectedInferedPath].size();
			tmp++)
		{
			tmpInferedPathMismatchPosVec.push_back(
				(tmpCandiInferedPathMismatchPosVecVec[index_bestSelectedInferedPath])[tmp]);
			tmpInferedPathMismatchCharVec.push_back(
				(tmpCandiInferedPathMismatchCharVecVec[index_bestSelectedInferedPath])[tmp]);
		}
		return true;
	}

	/*bool extendHeadWithAlignInferInfo_backward(const string& readSeq_inProcess,
		int unfixedHeadSeqLen, Index_Info* indexInfo, vector<Jump_Code>& extendedHeadJumpCodeVec,
		vector<int>& extendedHeadMismatchPosVecInRead, vector<char>& extendedHeadMismatchCharVec) // FIX ME: 04 / 15 / 2015
	{
		int alreadyExtendedHeadLen_inRead = 0;
		int alreadyExtendedHeadLen_inChr = 0;
		for(int tmp = 0; tmp < jumpCodeVec_backward.size(); tmp++)
		{
			string tmpJumpCodeType = jumpCodeVec_backward[tmp].type;
			int tmpJumpCodeLen = jumpCodeVec_backward[tmp].len;
			if(tmpJumpCodeType == "M")
			{
				if(tmpJumpCodeLen >= unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead)
				{
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					string tmpChromSeq = indexInfo->returnChromStrSubstr(chrNameInt, 
						donerEndPos - alreadyExtendedHeadLen_inChr - (unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead) + 1,
					    (unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead));
					string tmpReadSeq = readSeq_inProcess.substr(0, (unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead));
					int max_mismatch = (unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead)/MATCH_BASE_PER_MISMATCH_BASE;
					bool matchBool = fixMatchInfo->fixMatch(tmpReadSeq, tmpChromSeq, max_mismatch, 1);
					if(matchBool)
					{
						Jump_Code tmpMatchJumpCode((unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead), "M");
						extendedHeadJumpCodeVec.push_back(tmpMatchJumpCode);
						alreadyExtendedHeadLen_inRead = unfixedHeadSeqLen;
						alreadyExtendedHeadLen_inChr = unfixedHeadSeqLen;
						delete fixMatchInfo;
						return true;
					}
					else
					{
						Jump_Code tmpSoftClipJumpCode((unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead), "S");
						extendedHeadJumpCodeVec.push_back(tmpSoftClipJumpCode);
						delete fixMatchInfo;
						return true;
					}
				}
				else //(tmpJumpCodeLen < unfixedHeadSeqLen)
				{
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					string tmpChromSeq = indexInfo->returnChromStrSubstr(chrNameInt, donerEndPos - alreadyExtendedHeadLen_inChr - tmpJumpCodeLen+ 1,
						  tmpJumpCodeLen);
					string tmpReadSeq = readSeq_inProcess.substr(unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead - tmpJumpCodeLen, tmpJumpCodeLen);
					int max_mismatch = (tmpJumpCodeLen)/MATCH_BASE_PER_MISMATCH_BASE;
					bool matchBool = fixMatchInfo->fixMatch(tmpReadSeq, tmpChromSeq, max_mismatch, 1);
					if(matchBool)
					{
						Jump_Code tmpMatchJumpCode(tmpJumpCodeLen, "M");
						extendedHeadJumpCodeVec.push_back(tmpMatchJumpCode);
						alreadyExtendedHeadLen_inRead += tmpJumpCodeLen;
						alreadyExtendedHeadLen_inChr += tmpJumpCodeLen;
						delete fixMatchInfo;
						//return true;
					}
					else
					{
						Jump_Code tmpSoftClipJumpCode((unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead), "S");
						extendedHeadJumpCodeVec.push_back(tmpSoftClipJumpCode);
						delete fixMatchInfo;
						return true;
					}					
				}
			}
			else if(tmpJumpCodeType == "I")
			{
				if(tmpJumpCodeLen >= unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead)
				{
					Jump_Code tmpInsJumpCode(unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead, "I");
					extendedHeadJumpCodeVec.push_back(tmpInsJumpCode);
					alreadyExtendedHeadLen_inRead = unfixedHeadSeqLen;
					return true;
				}
				else
				{
					Jump_Code tmpInsJumpCode(tmpJumpCodeLen, "I");
					extendedHeadJumpCodeVec.push_back(tmpInsJumpCode);
					alreadyExtendedHeadLen_inRead += tmpJumpCodeLen;
				}
			}
			else if(tmpJumpCodeType == "D")
			{
				Jump_Code tmpDelJumpCode(tmpJumpCodeLen, "D");
				extendedHeadJumpCodeVec.push_back(tmpDelJumpCode);
				alreadyExtendedHeadLen_inChr += tmpJumpCodeLen;
			}
			else if(tmpJumpCodeType == "N")
			{
				Jump_Code tmpSpliceJumpCode(tmpJumpCodeLen, "N");
				extendedHeadJumpCodeVec.push_back(tmpSpliceJumpCode);
				alreadyExtendedHeadLen_inChr += tmpJumpCodeLen;				
			}
			else if(tmpJumpCodeType == "S")
			{
				Jump_Code tmpSoftClipJumpCode(unfixedHeadSeqLen - alreadyExtendedHeadLen_inRead, "S");
				extendedHeadJumpCodeVec.push_back(tmpSoftClipJumpCode);
				alreadyExtendedHeadLen_inRead = unfixedHeadSeqLen;
				return true;
			}
			else
			{
				cout << "error in extendHeadWithAlignInferInfo_backward ..." << endl;
				cout << "invaid jumpCode type: " << tmpJumpCodeType << endl;
				exit(1);
			}
		}
		return true;
	}*/

	string returnAlignInferInfoStr_onlyExtensionAnchorSimilarity(Index_Info* indexInfo, int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) 
			+ "\t" + int_to_str(anchorSizeMax_doner) + "\t" + int_to_str(anchorSizeMax_acceptor);

		// extension results
		tmpAlignInferInfoStr += "\t";
		string extensionJumpCodeVecStr_doner = jumpCodeVec2cigarString(extensionJumpCodeVec_doner);
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + extensionJumpCodeVecStr_doner + ":" + int_to_str(extension_penalty_doner); 
		tmpAlignInferInfoStr += "\t";
		string extensionJumpCodeVecStr_acceptor = jumpCodeVec2cigarString(extensionJumpCodeVec_acceptor);
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + extensionJumpCodeVecStr_acceptor + ":" + int_to_str(extension_penalty_acceptor);

		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_chrNamePosOnly(Index_Info* indexInfo, int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum);
		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_chrNamePos_strand(Index_Info* indexInfo, int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) + "\t" + strand;
		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_chrNamePos_supportNum(Index_Info* indexInfo,
		int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) + "\t" + int_to_str(supportNum);
		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_chrNamePos_supportNum_encompassingNum(Index_Info* indexInfo,
		int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) + "\t" + int_to_str(supportNum+encompassingNum)
			+ "\t" + int_to_str(supportNum) + "\t" + int_to_str(encompassingNum);
		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_chrNamePos_supportNum_anchorSize(Index_Info* indexInfo,
		int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) + "\t" + int_to_str(supportNum) + "\t"
			+ int_to_str(anchorSizeMax_doner) + "\t" + int_to_str(anchorSizeMax_acceptor);
		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_chrNamePos_supportNum_anchorSize_XM(Index_Info* indexInfo,
		int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) + "\t" + int_to_str(supportNum) + "\t"
			+ int_to_str(anchorSizeMax_doner) + "\t" + int_to_str(anchorSizeMax_acceptor) + "\t"
			+ int_to_str(XM_min) + "\t" + int_to_str(XM_max);
		return tmpAlignInferInfoStr;
	}

	string returnAlignInferInfoStr_maxAnchorSize(Index_Info* indexInfo, int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\t" + int_to_str(supportNum) + "\t" + flankStringStr + "\t" + int_to_str(flankStringCase)
			+ "\tJUNC_" + int_to_str(tmpJuncNum) 
			+ "\t" + int_to_str(anchorSizeMax_doner) + "\t" + int_to_str(anchorSizeMax_acceptor) + "\t";
		//cout << "tmpAlignInferInfoStr: " << tmpAlignInferInfoStr << endl;
		for(int tmpPath = 0; tmpPath < pathVec_backward.size(); tmpPath++)
		{
			for(int tmp = 0; tmp < pathVec_backward[tmpPath].size(); tmp++)
			{
				tmpAlignInferInfoStr = tmpAlignInferInfoStr + (pathVec_backward[tmpPath])[tmp].toString();
			}
			tmpAlignInferInfoStr += ",";
		}
		//cout << "tmpAlignInferInfoStr: " << tmpAlignInferInfoStr << endl;
		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathVec_forward.size(); tmpPath++)
		{
			for(int tmp = 0; tmp < pathVec_forward[tmpPath].size(); tmp++)
			{
				tmpAlignInferInfoStr = tmpAlignInferInfoStr + (pathVec_forward[tmpPath])[tmp].toString();
			}
			tmpAlignInferInfoStr += ",";
		}

		tmpAlignInferInfoStr += "\t";
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + int_to_str(pathNum_backward);
		tmpAlignInferInfoStr += "\t";
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + int_to_str(pathNum_forward);

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathSupportNumVec_backward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathSupportNumVec_backward[tmpPath]) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathSupportNumVec_forward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathSupportNumVec_forward[tmpPath]) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathAnchorSizeMax_backward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(pathAnchorSizeMax_backward[tmpPath]) + ","; 
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathAnchorSizeMax_forward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathAnchorSizeMax_forward[tmpPath]) + ",";
		}

		return tmpAlignInferInfoStr;
	}

	/*
	string returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites(Index_Info* indexInfo, int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\t" + int_to_str(supportNum) + "\t" + flankStringStr + "\t" + int_to_str(flankStringCase) 
			+ "\tJUNC_" + int_to_str(tmpJuncNum) 
			+ "\t" + int_to_str(anchorSizeMax_doner) + "\t" + int_to_str(anchorSizeMax_acceptor) + "\t";
		//cout << "tmpAlignInferInfoStr: " << tmpAlignInferInfoStr << endl;
		for(int tmpPath = 0; tmpPath < pathVec_backward.size(); tmpPath++)
		{
			for(int tmp = 0; tmp < pathVec_backward[tmpPath].size(); tmp++)
			{
				tmpAlignInferInfoStr = tmpAlignInferInfoStr + (pathVec_backward[tmpPath])[tmp].toString();
			}
			tmpAlignInferInfoStr += ",";
		}
		//cout << "tmpAlignInferInfoStr: " << tmpAlignInferInfoStr << endl;
		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathVec_forward.size(); tmpPath++)
		{
			for(int tmp = 0; tmp < pathVec_forward[tmpPath].size(); tmp++)
			{
				tmpAlignInferInfoStr = tmpAlignInferInfoStr + (pathVec_forward[tmpPath])[tmp].toString();
			}
			tmpAlignInferInfoStr += ",";
		}

		tmpAlignInferInfoStr += "\t";
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + int_to_str(pathNum_backward);
		tmpAlignInferInfoStr += "\t";
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + int_to_str(pathNum_forward);

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathSupportNumVec_backward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathSupportNumVec_backward[tmpPath]) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathSupportNumVec_forward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathSupportNumVec_forward[tmpPath]) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathAnchorSizeMax_backward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(pathAnchorSizeMax_backward[tmpPath]) + ","; 
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathAnchorSizeMax_forward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathAnchorSizeMax_forward[tmpPath]) + ","; 
		}

		tmpAlignInferInfoStr += "\t";
		if(alterDonerSpliceSitePairVec.size() == 0)
			tmpAlignInferInfoStr += "*";
		for(int tmpAlterSpliceSitePair = 0; 
			tmpAlterSpliceSitePair < alterDonerSpliceSitePairVec.size();
			tmpAlterSpliceSitePair++)
		{
			int tmpAlterDonerSpliceSite = alterDonerSpliceSitePairVec[tmpAlterSpliceSitePair].first;
			int tmpAlterAcceptorSpliceSite = alterDonerSpliceSitePairVec[tmpAlterSpliceSitePair].second;
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(tmpAlterDonerSpliceSite) + ":" + int_to_str(tmpAlterAcceptorSpliceSite) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		if(alterAcceptorSpliceSitePairVec.size() == 0)
			tmpAlignInferInfoStr += "*";
		for(int tmpAlterSpliceSitePair = 0; 
			tmpAlterSpliceSitePair < alterAcceptorSpliceSitePairVec.size();
			tmpAlterSpliceSitePair++)
		{
			int tmpAlterDonerSpliceSite = alterAcceptorSpliceSitePairVec[tmpAlterSpliceSitePair].first;
			int tmpAlterAcceptorSpliceSite = alterAcceptorSpliceSitePairVec[tmpAlterSpliceSitePair].second;
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(tmpAlterDonerSpliceSite) + ":" + int_to_str(tmpAlterAcceptorSpliceSite) + ",";
		}

		return tmpAlignInferInfoStr;
	}*/

	string returnAlignInferInfoStr_maxAnchorSize_alterSpliceSites_compareAnchorSimilarity(Index_Info* indexInfo, int tmpJuncNum)
	{
		string tmpAlignInferInfoStr;
		tmpAlignInferInfoStr = indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(donerEndPos) + "\t" + int_to_str(acceptorStartPos)
			+ "\t" + int_to_str(supportNum) + "\t" + flankStringStr + "\t" + int_to_str(flankStringCase) 
			+ "\tJUNC_" + int_to_str(tmpJuncNum) 
			+ "\t" + int_to_str(anchorSizeMax_doner) + "\t" + int_to_str(anchorSizeMax_acceptor) + "\t";
		//cout << "tmpAlignInferInfoStr: " << tmpAlignInferInfoStr << endl;
		for(int tmpPath = 0; tmpPath < pathVec_backward.size(); tmpPath++)
		{
			for(int tmp = 0; tmp < pathVec_backward[tmpPath].size(); tmp++)
			{
				tmpAlignInferInfoStr = tmpAlignInferInfoStr + (pathVec_backward[tmpPath])[tmp].toString();
			}
			tmpAlignInferInfoStr += ",";
		}
		//cout << "tmpAlignInferInfoStr: " << tmpAlignInferInfoStr << endl;
		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathVec_forward.size(); tmpPath++)
		{
			for(int tmp = 0; tmp < pathVec_forward[tmpPath].size(); tmp++)
			{
				tmpAlignInferInfoStr = tmpAlignInferInfoStr + (pathVec_forward[tmpPath])[tmp].toString();
			}
			tmpAlignInferInfoStr += ",";
		}

		tmpAlignInferInfoStr += "\t";
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + int_to_str(pathNum_backward);
		tmpAlignInferInfoStr += "\t";
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + int_to_str(pathNum_forward);

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathSupportNumVec_backward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathSupportNumVec_backward[tmpPath]) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathSupportNumVec_forward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathSupportNumVec_forward[tmpPath]) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathAnchorSizeMax_backward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(pathAnchorSizeMax_backward[tmpPath]) + ","; 
		}

		tmpAlignInferInfoStr += "\t";
		for(int tmpPath = 0; tmpPath < pathAnchorSizeMax_forward.size(); tmpPath++)
		{
			tmpAlignInferInfoStr = tmpAlignInferInfoStr 
				+ int_to_str(pathAnchorSizeMax_forward[tmpPath]) + ","; 
		}

		// extension results
		tmpAlignInferInfoStr += "\t";
		string extensionJumpCodeVecStr_doner = jumpCodeVec2cigarString(extensionJumpCodeVec_doner);
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + extensionJumpCodeVecStr_doner + ":" + int_to_str(extension_penalty_doner); 
		tmpAlignInferInfoStr += "\t";
		string extensionJumpCodeVecStr_acceptor = jumpCodeVec2cigarString(extensionJumpCodeVec_acceptor);
		tmpAlignInferInfoStr = tmpAlignInferInfoStr + extensionJumpCodeVecStr_acceptor + ":" + int_to_str(extension_penalty_acceptor);

		tmpAlignInferInfoStr += "\t";
		if(alterDonerSpliceSitePairVec.size() == 0)
			tmpAlignInferInfoStr += "*";
		for(int tmpAlterSpliceSitePair = 0; 
			tmpAlterSpliceSitePair < alterDonerSpliceSitePairVec.size();
			tmpAlterSpliceSitePair++)
		{
			int tmpAlterDonerSpliceSite = alterDonerSpliceSitePairVec[tmpAlterSpliceSitePair].first;
			int tmpAlterAcceptorSpliceSite = alterDonerSpliceSitePairVec[tmpAlterSpliceSitePair].second;
			string tmpNWDPjumpCodeVecStr = jumpCodeVec2cigarString(alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec[tmpAlterSpliceSitePair]);
			int tmpNWDPpenalty = alterDonerSpliceSiteAnchorNWDPpenaltyVec[tmpAlterSpliceSitePair];
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(tmpAlterDonerSpliceSite) + ":" + int_to_str(tmpAlterAcceptorSpliceSite) + ":"
				+ tmpNWDPjumpCodeVecStr + ":" + int_to_str(tmpNWDPpenalty) + ",";
		}

		tmpAlignInferInfoStr += "\t";
		if(alterAcceptorSpliceSitePairVec.size() == 0)
			tmpAlignInferInfoStr += "*";
		for(int tmpAlterSpliceSitePair = 0; 
			tmpAlterSpliceSitePair < alterAcceptorSpliceSitePairVec.size();
			tmpAlterSpliceSitePair++)
		{
			int tmpAlterDonerSpliceSite = alterAcceptorSpliceSitePairVec[tmpAlterSpliceSitePair].first;
			int tmpAlterAcceptorSpliceSite = alterAcceptorSpliceSitePairVec[tmpAlterSpliceSitePair].second;
			string tmpNWDPjumpCodeVecStr = jumpCodeVec2cigarString(alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec[tmpAlterSpliceSitePair]);
			int tmpNWDPpenalty = alterAcceptorSpliceSiteAnchorNWDPpenaltyVec[tmpAlterSpliceSitePair];
			tmpAlignInferInfoStr = tmpAlignInferInfoStr
				+ int_to_str(tmpAlterDonerSpliceSite) + ":" + int_to_str(tmpAlterAcceptorSpliceSite) + ":"
				+ tmpNWDPjumpCodeVecStr + ":" + int_to_str(tmpNWDPpenalty) + ",";
		}

		return tmpAlignInferInfoStr;
	}

	string returnFlankString(Index_Info* indexInfo)
	{
		string tmpFlankString = indexInfo->returnFlankString(chrNameInt, donerEndPos, acceptorStartPos);
		return tmpFlankString;
	}

	int getFlankStringCase(Index_Info* indexInfo)
	{
		string tmpFlankString = indexInfo->returnFlankString(chrNameInt, donerEndPos, acceptorStartPos);
		flankStringStr = tmpFlankString;
		if(tmpFlankString == "GTAG")
			return 5;
		else if(tmpFlankString == "CTAC")
			return 6;
		else if(tmpFlankString == "ATAC")
			return 1;
		else if(tmpFlankString == "GTAT")
			return 2;
		else if(tmpFlankString == "CTGC")
			return 3;
		else if(tmpFlankString == "GCAG")
			return 4;
		else
			return 0;
	}

	void initiateAlignInferInfo_maxAnchorSize(int mapChrNameInt, int tmpSJposDonerEnd, 
		int tmpSJposAcceptorStart, 
		vector<Jump_Code>& tmpSJjumpCodeVec_backward, 
		vector<Jump_Code>& tmpSJjumpCodeVec_forward,
		int tmpDonerAnchorSize, int tmpAcceptorAnchorSize, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;

		//flankStringStr = 
		flankStringCase = getFlankStringCase(indexInfo);
		supportNum = 1;
		// copy backward jumpCodeVec
		/*for(int tmp = 0; tmp < tmpSJjumpCodeVec_backward.size(); tmp++)
		{
			jumpCodeVec_backward.push_back(tmpSJjumpCodeVec_backward[tmp]);
			//cout << "jumpCodeStr: " << tmpSJjumpCodeVec_backward[tmp].toString() << endl;
		}*/
		pathVec_backward.push_back(tmpSJjumpCodeVec_backward);
		pathSupportNumVec_backward.push_back(1);
		pathAnchorSizeMax_backward.push_back(tmpDonerAnchorSize);
		pathNum_backward = 1;	
		// copy forward jumpCodeVec
		/*for(int tmp = 0; tmp < tmpSJjumpCodeVec_forward.size(); tmp++)
		{
			jumpCodeVec_forward.push_back(tmpSJjumpCodeVec_forward[tmp]);
			//cout << "jumpCodeStr: " << tmpSJjumpCodeVec_forward[tmp].toString() << endl;
		}*/
		pathVec_forward.push_back(tmpSJjumpCodeVec_forward);
		pathSupportNumVec_forward.push_back(1);
		pathAnchorSizeMax_forward.push_back(tmpAcceptorAnchorSize);
		pathNum_forward = 1;

		anchorSizeMax_doner = tmpDonerAnchorSize;
		anchorSizeMax_acceptor = tmpAcceptorAnchorSize;

		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;
		//cout << "end of initiateAlignInferInfo_maxAnchorSize" << endl;
	}

	void initiateAlignInferInfo_onlyChrNamePos(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, Index_Info* indexInfo,
		int tmpAnchorSizeMax_doner, int tmpAnchorSizeMax_acceptor) 
		// used to genearte the histogram graph of annotated SJs
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;

		anchorSizeMax_doner = tmpAnchorSizeMax_doner;
		anchorSizeMax_acceptor = tmpAnchorSizeMax_acceptor;

		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;
	}

	void initiateAlignInferInfo_chrNamePosOnly(int mapChrNameInt,
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
	}

	void initiateAlignInferInfo_chrNamePos_supportNum(int mapChrNameInt,
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = 1;
	}	

	void initiateAlignInferInfo_chrNamePos_supportNum(int mapChrNameInt,
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, string& tmpStrand, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = 1;
		strand = tmpStrand;
	}	

	void initiateAlignInferInfo_chrNamePos_supportNum(int mapChrNameInt,
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJsupportNum, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = tmpSJsupportNum;
	}	

	void initiateAlignInferInfo_chrNamePos_supportNum_flankStringChange(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJsupportNum, Index_Info* indexInfo, string& tmpSJflankStringChange)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = tmpSJsupportNum;
		flankStringChanged = tmpSJflankStringChange.substr(6, 4);
	}	

	void initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJsupportNum,
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = tmpSJsupportNum;
		anchorSizeMax_doner = tmpSJanchorSizeDonerEnd;
		anchorSizeMax_acceptor = tmpSJanchorSizeAcceptorStart;
	}

	void initiateAlignInferInfo_chrNamePos_supportNum_anchorSize(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = 1;
		anchorSizeMax_doner = tmpSJanchorSizeDonerEnd;
		anchorSizeMax_acceptor = tmpSJanchorSizeAcceptorStart;
	}

	void initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		int tmpXM, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = 1;
		anchorSizeMax_doner = tmpSJanchorSizeDonerEnd;
		anchorSizeMax_acceptor = tmpSJanchorSizeAcceptorStart;
		XM_min = tmpXM;
		XM_max = tmpXM;
	}

	void initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		int tmpXMmin, int tmpXMmax, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = 1;
		anchorSizeMax_doner = tmpSJanchorSizeDonerEnd;
		anchorSizeMax_acceptor = tmpSJanchorSizeAcceptorStart;
		XM_min = tmpXMmin;
		XM_max = tmpXMmax;
	}

	void initiateAlignInferInfo_chrNamePos_supportNum_anchorSize_XM(
		int mapChrNameInt, int tmpSJposDonerEnd, int tmpSJposAcceptorStart, 
		int tmpSJsupportNum,
		int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart,
		int tmpXMmin, int tmpXMmax, Index_Info* indexInfo)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;		
		supportNum = tmpSJsupportNum;
		anchorSizeMax_doner = tmpSJanchorSizeDonerEnd;
		anchorSizeMax_acceptor = tmpSJanchorSizeAcceptorStart;
		XM_min = tmpXMmin;
		XM_max = tmpXMmax;
	}	

	void addSupportNum(int supportNum2add)
	{
		supportNum = supportNum + supportNum2add;
	}

	void initiateAlignInferInfoWithMultiPathVec(int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, int tmpSJsupportNum,
		vector< vector<Jump_Code> >& tmpSJdonerPathVec, 
		vector< vector<Jump_Code> >& tmpSJacceptorPathVec, 
		int tmpSJdonerPathVecNum, int tmpSJacceptorPathVecNum, 
		vector<int>& tmpSJdonerPathSupportNumVec,
		vector<int>& tmpSJacceptorPathSupportNumVec,
		vector<int>& tmpDonerPathMaxAnchorSizeVec,
		vector<int>& tmpAcceptorPathMaxAnchorSizeVec,
		int tmpDonerAnchorSizeMax,
		int tmpAcceptorAnchorSizeMax, Index_Info* indexInfo,
		string& tmpFlankString, int tmpFlankStringCase)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;
		flankStringStr = tmpFlankString;
		flankStringCase = tmpFlankStringCase;
		supportNum = tmpSJsupportNum;
		// copy backward jumpCodeVec
		/*for(int tmp = 0; tmp < tmpSJjumpCodeVec_backward.size(); tmp++)
		{
			jumpCodeVec_backward.push_back(tmpSJjumpCodeVec_backward[tmp]);
			//cout << "jumpCodeStr: " << tmpSJjumpCodeVec_backward[tmp].toString() << endl;
		}*/
		for(int tmp = 0; tmp < tmpSJdonerPathVecNum; tmp++)
		{	
			pathVec_backward.push_back(tmpSJdonerPathVec[tmp]);
			pathSupportNumVec_backward.push_back(tmpSJdonerPathSupportNumVec[tmp]);
		}
		pathNum_backward = tmpSJdonerPathVecNum;	
		// copy forward jumpCodeVec
		/*for(int tmp = 0; tmp < tmpSJjumpCodeVec_forward.size(); tmp++)
		{
			jumpCodeVec_forward.push_back(tmpSJjumpCodeVec_forward[tmp]);
			//cout << "jumpCodeStr: " << tmpSJjumpCodeVec_forward[tmp].toString() << endl;
		}*/
		for(int tmp = 0; tmp < tmpSJacceptorPathVecNum; tmp++)
		{	
			pathVec_forward.push_back(tmpSJacceptorPathVec[tmp]);
			pathSupportNumVec_forward.push_back(tmpSJacceptorPathSupportNumVec[tmp]);
		}
		pathNum_forward = tmpSJacceptorPathVecNum;

		anchorSizeMax_doner = tmpDonerAnchorSizeMax;
		anchorSizeMax_acceptor = tmpAcceptorAnchorSizeMax;

		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		for(int tmp = 0; tmp < tmpDonerPathMaxAnchorSizeVec.size(); tmp++)
		{
			pathAnchorSizeMax_backward.push_back(tmpDonerPathMaxAnchorSizeVec[tmp]);
		}
		for(int tmp = 0; tmp < tmpAcceptorPathMaxAnchorSizeVec.size(); tmp++)
		{
			pathAnchorSizeMax_forward.push_back(tmpAcceptorPathMaxAnchorSizeVec[tmp]);
		}
	}

	void initiateAlignInferInfoWithMultiPathVec_withAlterSpliceSiteAnchorSimilarity(
		int mapChrNameInt, 
		int tmpSJposDonerEnd, int tmpSJposAcceptorStart, int tmpSJsupportNum,
		vector< vector<Jump_Code> >& tmpSJdonerPathVec, 
		vector< vector<Jump_Code> >& tmpSJacceptorPathVec, 
		int tmpSJdonerPathVecNum, int tmpSJacceptorPathVecNum, 
		vector<int>& tmpSJdonerPathSupportNumVec,
		vector<int>& tmpSJacceptorPathSupportNumVec,
		vector<int>& tmpDonerPathMaxAnchorSizeVec,
		vector<int>& tmpAcceptorPathMaxAnchorSizeVec,
		int tmpDonerAnchorSizeMax,
		int tmpAcceptorAnchorSizeMax, Index_Info* indexInfo,
		string& tmpFlankString, int tmpFlankStringCase,
		vector<Jump_Code>& tmpExtensionJumpCodeVec_doner,
		int& tmpExtensionPenalty_doner,
		vector<Jump_Code>& tmpExtensionJumpCodeVec_acceptor,
		int& tmpExtensionPenalty_acceptor,
		vector< pair<int,int> >& tmpAlterSpliceSitePairVec_doner,
		vector< vector<Jump_Code> >& tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner,
		vector< int >& tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner,
		vector< pair<int,int> >& tmpAlterSpliceSitePairVec_acceptor,
		vector< vector<Jump_Code> >& tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor,
		vector< int >& tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor)
	{
		chrNameInt = mapChrNameInt;
		//cout << "chrNameInt: " << mapChrNameInt << endl;
		donerEndPos = tmpSJposDonerEnd;
		//cout << "donerEndPos: " << donerEndPos << endl;
		acceptorStartPos = tmpSJposAcceptorStart;
		flankStringStr = tmpFlankString;
		flankStringCase = tmpFlankStringCase;
		supportNum = tmpSJsupportNum;
		// copy backward jumpCodeVec
		/*for(int tmp = 0; tmp < tmpSJjumpCodeVec_backward.size(); tmp++)
		{
			jumpCodeVec_backward.push_back(tmpSJjumpCodeVec_backward[tmp]);
			//cout << "jumpCodeStr: " << tmpSJjumpCodeVec_backward[tmp].toString() << endl;
		}*/
		for(int tmp = 0; tmp < tmpSJdonerPathVecNum; tmp++)
		{	
			pathVec_backward.push_back(tmpSJdonerPathVec[tmp]);
			pathSupportNumVec_backward.push_back(tmpSJdonerPathSupportNumVec[tmp]);
		}
		pathNum_backward = tmpSJdonerPathVecNum;	
		// copy forward jumpCodeVec
		/*for(int tmp = 0; tmp < tmpSJjumpCodeVec_forward.size(); tmp++)
		{
			jumpCodeVec_forward.push_back(tmpSJjumpCodeVec_forward[tmp]);
			//cout << "jumpCodeStr: " << tmpSJjumpCodeVec_forward[tmp].toString() << endl;
		}*/
		for(int tmp = 0; tmp < tmpSJacceptorPathVecNum; tmp++)
		{	
			pathVec_forward.push_back(tmpSJacceptorPathVec[tmp]);
			pathSupportNumVec_forward.push_back(tmpSJacceptorPathSupportNumVec[tmp]);
		}
		pathNum_forward = tmpSJacceptorPathVecNum;

		anchorSizeMax_doner = tmpDonerAnchorSizeMax;
		anchorSizeMax_acceptor = tmpAcceptorAnchorSizeMax;

		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;

		for(int tmp = 0; tmp < tmpDonerPathMaxAnchorSizeVec.size(); tmp++)
		{
			pathAnchorSizeMax_backward.push_back(tmpDonerPathMaxAnchorSizeVec[tmp]);
		}
		for(int tmp = 0; tmp < tmpAcceptorPathMaxAnchorSizeVec.size(); tmp++)
		{
			pathAnchorSizeMax_forward.push_back(tmpAcceptorPathMaxAnchorSizeVec[tmp]);
		}

		// extension -- doner 
		for(int tmp = 0; tmp < tmpExtensionJumpCodeVec_doner.size(); tmp++)
		{
			extensionJumpCodeVec_doner.push_back(tmpExtensionJumpCodeVec_doner[tmp]);
			//cout << "tmpDonerExtensionJumpCode: " << tmpExtensionJumpCodeVec_doner[tmp].toString() << endl;
		}
		extension_penalty_doner = tmpExtensionPenalty_doner;
		//cout << "extension_penalty_doner: " << extension_penalty_doner << endl;
		// extension -- acceptor
		for(int tmp = 0; tmp < tmpExtensionJumpCodeVec_acceptor.size(); tmp++)
		{
			extensionJumpCodeVec_acceptor.push_back(tmpExtensionJumpCodeVec_acceptor[tmp]);
			//cout << "tmpAcceptorExtensionJumpCode: " << tmpExtensionJumpCodeVec_acceptor[tmp].toString() << endl;
		}
		extension_penalty_acceptor = tmpExtensionPenalty_acceptor;
		//cout << "extension_penalty_acceptor: " << extension_penalty_acceptor << endl;
		// alterSpliceSite -- doner
		for(int tmp = 0; tmp < tmpAlterSpliceSitePairVec_doner.size(); tmp++)
		{
			alterDonerSpliceSitePairVec.push_back(tmpAlterSpliceSitePairVec_doner[tmp]);
			alterDonerSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(
				tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_doner[tmp]);
			alterDonerSpliceSiteAnchorNWDPpenaltyVec.push_back(
				tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_doner[tmp]);
		}
		// alterSpliceSite -- acceptor
		for(int tmp = 0; tmp < tmpAlterSpliceSitePairVec_acceptor.size(); tmp++)
		{
			alterAcceptorSpliceSitePairVec.push_back(tmpAlterSpliceSitePairVec_acceptor[tmp]);
			alterAcceptorSpliceSiteAnchorNWDPjumpCodeVecVec.push_back(
				tmpAlterSpliceSiteAnchorSimilarityJumpCodeVecVec_acceptor[tmp]);
			alterAcceptorSpliceSiteAnchorNWDPpenaltyVec.push_back(
				tmpAlterSpliceSiteAnchorSimilarityPenaltyVec_acceptor[tmp]);
		}	
	}


	void updateMaxAnchorSize_NWDPpenaltyMax(int tmpDonerAnchorSize, int tmpAcceptorAnchorSize)
	{
		if(tmpDonerAnchorSize > anchorSizeMax_doner)
		{	
			anchorSizeMax_doner = tmpDonerAnchorSize;
			donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;				
		}
		if(tmpAcceptorAnchorSize > anchorSizeMax_acceptor)
		{	
			anchorSizeMax_acceptor = tmpAcceptorAnchorSize;
			acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;
		}
	}

	void updateOldAlignInfo_maxAnchorSize(
		vector<Jump_Code>& tmpSJjumpCodeVec_backward,
		vector<Jump_Code>& tmpSJjumpCodeVec_forward, int maxReadBaseNum,
		int tmpDonerAnchorSize, int tmpAcceptorAnchorSize)
	{
		pair<int,int> checkWithOldPath_backward_relation_pair 
			= this->checkWithOldPathVec_backward(tmpSJjumpCodeVec_backward);
		int compliedOldPathIndex_backward = checkWithOldPath_backward_relation_pair.second;
		int checkWithOldPath_backward_relation = checkWithOldPath_backward_relation_pair.first;			
		if(checkWithOldPath_backward_relation == -1) // no need to update path
		{
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
			int tmpPathMaxAnchorSize_backward = pathAnchorSizeMax_backward[compliedOldPathIndex_backward];
			if(tmpDonerAnchorSize > tmpPathMaxAnchorSize_backward)
				pathAnchorSizeMax_backward[compliedOldPathIndex_backward] = tmpDonerAnchorSize;
		}
		else if(checkWithOldPath_backward_relation == -2) // new path
		{
			pathVec_backward.push_back(tmpSJjumpCodeVec_backward);
			pathSupportNumVec_backward.push_back(1);
			pathAnchorSizeMax_backward.push_back(tmpDonerAnchorSize);
			pathNum_backward ++;
		}
		else if(checkWithOldPath_backward_relation == 0)
		{
			pathVec_backward[compliedOldPathIndex_backward] = tmpSJjumpCodeVec_backward;
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
			int tmpPathMaxAnchorSize_backward = pathAnchorSizeMax_backward[compliedOldPathIndex_backward];
			if(tmpDonerAnchorSize > tmpPathMaxAnchorSize_backward)
				pathAnchorSizeMax_backward[compliedOldPathIndex_backward] = tmpDonerAnchorSize;
		}
		else
		{
			cout << "checkWithOldPath_backward_relation: " << checkWithOldPath_backward_relation << endl;
			cout << "error in checkWithOldPath_backward_relation ..." << endl;
			exit(1);
		}

		pair<int,int> checkWithOldPath_forward_relation_pair
			= this->checkWithOldPathVec_forward(tmpSJjumpCodeVec_forward);
		int compliedOldPathIndex_forward = checkWithOldPath_forward_relation_pair.second;
		int checkWithOldPath_forward_relation = checkWithOldPath_forward_relation_pair.first; 
		if(checkWithOldPath_forward_relation == -1) // no need to update path
		{
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
			int tmpPathMaxAnchorSize_forward = pathAnchorSizeMax_forward[compliedOldPathIndex_forward];
			if(tmpAcceptorAnchorSize > tmpPathMaxAnchorSize_forward)
				pathAnchorSizeMax_forward[compliedOldPathIndex_forward] = tmpAcceptorAnchorSize;
		}
		else if(checkWithOldPath_forward_relation == -2) // new path
		{
			pathVec_forward.push_back(tmpSJjumpCodeVec_forward);
			pathSupportNumVec_forward.push_back(1);
			pathAnchorSizeMax_forward.push_back(tmpAcceptorAnchorSize);
			pathNum_forward ++;
		}
		else if(checkWithOldPath_forward_relation == 0)
		{
			pathVec_forward[compliedOldPathIndex_forward] = tmpSJjumpCodeVec_forward;
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
			int tmpPathMaxAnchorSize_forward = pathAnchorSizeMax_forward[compliedOldPathIndex_forward];
			if(tmpAcceptorAnchorSize > tmpPathMaxAnchorSize_forward)
				pathAnchorSizeMax_forward[compliedOldPathIndex_forward] = tmpAcceptorAnchorSize;
		}
		else
		{
			cout << "checkWithOldPath_forward_relation: " << checkWithOldPath_forward_relation << endl;
			cout << "error in checkWithOldPath_forward_relation ..." << endl;
			exit(1);
		}
		this->updateMaxAnchorSize_NWDPpenaltyMax(tmpDonerAnchorSize, tmpAcceptorAnchorSize);
		this->updateSupportNum();
	}

	pair<int,int> checkWithOldPathVec_backward(
		vector<Jump_Code>& newSJjumpCodeVec_backward)
	{
		// check with the old paths and return the relation checking results:
		// Output:
		//	<0, indexInOldPath(>=0)>, -- comply with an old path and need to extend
		//  or
		//	<-1, indexInOldPath(>=0)> -- comply with an old path but no need to extend
		//	or
		//	<-2, -2> -- new path	
		// Methods:				 
		// 1. comply with one of the old paths
		//		1.1 longer than the old path
		//			1.1.1 bases num in old path < ALIGNINFER_BASE_NUM_MAX -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX), 
		//				 return <0,indexInOldPath>
		//			1.1.2 bases num in old path == ALIGNINFER_BASE_NUM_MAX -- skip it 
		//				 return <-1,indexInOldPath>;
		//		1.2 shorter than the old path or the same with it -- skip it
		//				return <-1,indexInOldPath>;
		// 2. got a new path -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX)
		//			return <-2,-2>
		for(int tmp = 0; tmp < pathVec_backward.size(); tmp++)
		{
			int checkOldPath_relation 
				= this->checkWithOldPath_backward(newSJjumpCodeVec_backward, tmp);
			if(checkOldPath_relation == 0)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(-1,tmp);
			}
			else if(checkOldPath_relation == 1)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(0,tmp);
			}
			else
			{}
		}
		return pair<int,int>(-2,-2);
	}

	int checkWithOldPath_backward(
		vector<Jump_Code>& newSJjumpCodeVec_backward, int index_oldPathVec_backward)
		// return: 1 -- comply with old path, need to update
		// 		   0 -- comply with old path, no need to update
		//        -1 -- new path
	{
		int newSJjumpCodeVecSize = newSJjumpCodeVec_backward.size();
		int oldSJjumpCodeVecSize = (pathVec_backward[index_oldPathVec_backward]).size();
		int minSJjumpCodeVecSize = newSJjumpCodeVecSize;
		if(oldSJjumpCodeVecSize < minSJjumpCodeVecSize)
			minSJjumpCodeVecSize = oldSJjumpCodeVecSize;

		for(int tmp = 0; tmp < minSJjumpCodeVecSize-1; tmp++) // check first N-1 JumpCode
		{
			int newSJjumpCodeLength = newSJjumpCodeVec_backward[tmp].len;
			string newSJjumpCodeType = newSJjumpCodeVec_backward[tmp].type;
			int oldSJjumpCodeLength = (pathVec_backward[index_oldPathVec_backward])[tmp].len;
			string oldSJjumpCodeType = (pathVec_backward[index_oldPathVec_backward])[tmp].type;			

			if( (newSJjumpCodeLength != oldSJjumpCodeLength) || (newSJjumpCodeType != oldSJjumpCodeType))
				return -1;
		}

		// check last shared JumpCode
		int newSJjumpCodeLength_lastShared = newSJjumpCodeVec_backward[minSJjumpCodeVecSize-1].len;
		string newSJjumpCodeType_lastShared = newSJjumpCodeVec_backward[minSJjumpCodeVecSize-1].type;
		int oldSJjumpCodeLength_lastShared = (pathVec_backward[index_oldPathVec_backward])[minSJjumpCodeVecSize-1].len;
		string oldSJjumpCodeType_lastShared = (pathVec_backward[index_oldPathVec_backward])[minSJjumpCodeVecSize-1].type;	
		if(newSJjumpCodeType_lastShared == oldSJjumpCodeType_lastShared)
		{
			if(newSJjumpCodeVecSize < oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return -1;
			}
			else if(newSJjumpCodeVecSize == oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return 1;				
			}
			else // newSJjumpCodeVecSize > oldSJjumpCodeVecSize
			{
				if(newSJjumpCodeLength_lastShared < oldSJjumpCodeLength_lastShared)
					return -1;
				else
					return 1;				
			}
		}
		else // last shared JumpCode inconsistent
		{
			return -1;
		}
	}

	pair<int,int> checkWithOldPathVec_forward(
		vector<Jump_Code>& newSJjumpCodeVec_forward)
	{
		// check with the old paths and return the relation checking results:
		// Output:
		//	<0, indexInOldPath(>=0)>, -- comply with an old path and need to extend
		//  or
		//	<-1, indexInOldPath(>=0)> -- comply with an old path but no need to extend
		//	or
		//	<-2, -2> -- new path	
		// Methods:				 
		// 1. comply with one of the old paths
		//		1.1 longer than the old path
		//			1.1.1 bases num in old path < ALIGNINFER_BASE_NUM_MAX -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX), 
		//				 return <0,indexInOldPath>
		//			1.1.2 bases num in old path == ALIGNINFER_BASE_NUM_MAX -- skip it 
		//				 return <-1,indexInOldPath>;
		//		1.2 shorter than the old path or the same with it -- skip it
		//				return <-1,indexInOldPath>;
		// 2. got a new path -- do extension (base num <= ALIGNINFER_BASE_NUM_MAX)
		//			return <-2,-1>
		for(int tmp = 0; tmp < pathVec_forward.size(); tmp++)
		{
			int checkOldPath_relation 
				= this->checkWithOldPath_forward(newSJjumpCodeVec_forward, tmp);
			if(checkOldPath_relation == 0)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(-1,tmp);
			}
			else if(checkOldPath_relation == 1)
			{
				//compliedOldPathIndex = tmp;
				return pair<int,int>(0,tmp);
			}
			else
			{}
		}
		return pair<int,int>(-2,-2);
	}

	int checkWithOldPath_forward(
		vector<Jump_Code>& newSJjumpCodeVec_forward, int index_oldPathVec_forward)
		// return: 1 -- comply with old path, need to update
		// 		   0 -- comply with old path, no need to update
		//        -1 -- new path
	{
		int newSJjumpCodeVecSize = newSJjumpCodeVec_forward.size();
		int oldSJjumpCodeVecSize = (pathVec_forward[index_oldPathVec_forward]).size();
		int minSJjumpCodeVecSize = newSJjumpCodeVecSize;
		if(oldSJjumpCodeVecSize < minSJjumpCodeVecSize)
			minSJjumpCodeVecSize = oldSJjumpCodeVecSize;

		for(int tmp = 0; tmp < minSJjumpCodeVecSize-1; tmp++) // check first N-1 JumpCode
		{
			int newSJjumpCodeLength = newSJjumpCodeVec_forward[tmp].len;
			string newSJjumpCodeType = newSJjumpCodeVec_forward[tmp].type;
			int oldSJjumpCodeLength = (pathVec_forward[index_oldPathVec_forward])[tmp].len;
			string oldSJjumpCodeType = (pathVec_forward[index_oldPathVec_forward])[tmp].type;			

			if( (newSJjumpCodeLength != oldSJjumpCodeLength) || (newSJjumpCodeType != oldSJjumpCodeType))
				return -1;
		}

		// check last shared JumpCode
		int newSJjumpCodeLength_lastShared = newSJjumpCodeVec_forward[minSJjumpCodeVecSize-1].len;
		string newSJjumpCodeType_lastShared = newSJjumpCodeVec_forward[minSJjumpCodeVecSize-1].type;
		int oldSJjumpCodeLength_lastShared = (pathVec_forward[index_oldPathVec_forward])[minSJjumpCodeVecSize-1].len;
		string oldSJjumpCodeType_lastShared = (pathVec_forward[index_oldPathVec_forward])[minSJjumpCodeVecSize-1].type;	
		if(newSJjumpCodeType_lastShared == oldSJjumpCodeType_lastShared)
		{
			if(newSJjumpCodeVecSize < oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return -1;
			}
			else if(newSJjumpCodeVecSize == oldSJjumpCodeVecSize)
			{
				if(newSJjumpCodeLength_lastShared <= oldSJjumpCodeLength_lastShared)
					return 0;
				else
					return 1;				
			}
			else // newSJjumpCodeVecSize > oldSJjumpCodeVecSize
			{
				if(newSJjumpCodeLength_lastShared < oldSJjumpCodeLength_lastShared)
					return -1;
				else
					return 1;				
			}
		}
		else // last shared JumpCode inconsistent
		{
			return -1;
		}
	}

	/*
	void updateWithNewAlignmentForOldSJ(
		vector<Jump_Code>& newSJjumpCodeVec_backward, 
		vector<Jump_Code>& newSJjumpCodeVec_forward,
		int newSJjumpCodeReadSeqBaseNum_backward,
		int newSJjumpCodeReadSeqBaseNum_forward)
	{
		int checkRelationWithOldPath = this->checkWithTheOldPaths(newSJjumpCodeVec_backward, newSJjumpCodeVec_forward);
		int indexInOldPathVec = checkRelationWithOldPath;
		if(indexInOldPathVec == -1) // got a new path
		{

		}
		else
		{
			int compliedOldPathReadBaseNum_backward = supportNumVec_backward[indexInOldPathVec];
			int compliedOldPathReadBaseNum_forward = supportNumVec_forward[indexInOldPathVec];
		}
	}*/

	void updateAnchorSize(int tmpSJanchorSizeDonerEnd, int tmpSJanchorSizeAcceptorStart)
	{
		if(tmpSJanchorSizeDonerEnd > anchorSizeMax_doner)
			anchorSizeMax_doner = tmpSJanchorSizeDonerEnd;
		if(tmpSJanchorSizeAcceptorStart > anchorSizeMax_acceptor)
			anchorSizeMax_acceptor = tmpSJanchorSizeAcceptorStart;
	}

	void updateXM(int tmpXMmin, int tmpXMmax)
	{
		if(tmpXMmin < XM_min)
			XM_min = tmpXMmin;
		if(tmpXMmax > XM_max)
			XM_max = tmpXMmax;
	}

	void updateXM(int tmpXM)
	{
		if(tmpXM < XM_min)
			XM_min = tmpXM;
		if(tmpXM > XM_max)
			XM_max = tmpXM;
	}

	void updateSupportNum()
	{
		supportNum ++;
	}

	void updateEncompassingNum()
	{
		encompassingNum ++;
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

	int pushBackNewSupportRefineSAMindex_SE(RefineSAM_SE_Vec_Info* refineSamSEvecInfo)
	{
		int tmpIndex = refineSamSEvecInfo->returnCurrentIndex();
		indexInRefineSAMinfoVec.push_back(tmpIndex);
		//refineSamSEvecInfo->toAddRefineSAMindex_incre();
	}

	int pushBackNewSupportRefineSAMindex_PE(RefineSAM_PE_Vec_Info* refineSamPEvecInfo)
	{
		int tmpIndex = refineSamPEvecInfo->returnCurrentIndex();
		indexInRefineSAMinfoVec.push_back(tmpIndex);
		//refineSamSEvecInfo->toAddRefineSAMindex_incre();
	}	

	void initiateAlignInferInfo_withRefineSAMinfo_SE(
		RefineSAM_SE_Info* refineSamSEinfo, int SJindexInRefineSAMinfo, 
		Index_Info* indexInfo)
	{
		//cout << "initiateAlignInferInfo_withRefineSAMinfo_SE starts ..." << endl;
		chrNameInt = refineSamSEinfo->returnChrNameInt();
		donerEndPos = refineSamSEinfo->returnSJ_donerEndPos(SJindexInRefineSAMinfo);
		acceptorStartPos = refineSamSEinfo->returnSJ_acceptorStartPos(SJindexInRefineSAMinfo);
		flankStringCase = getFlankStringCase(indexInfo);
		supportNum = 1;
		//cout << "chrNameInt: " << endl;
		//cout << "donerEndPos: " << donerEndPos << endl;
		//cout << "acceptorStartPos: " << acceptorStartPos << endl;
		//cout << "flankStringCase: " << flankStringCase << endl;
		pathVec_backward.push_back(
			refineSamSEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo));
		pathSupportNumVec_backward.push_back(1);
		//cout << "end of push_back pathVec_backward " << endl;
		int tmpDonerAnchorSize 
			= refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo);
		//cout << "tmpDonerAnchorSize: " << tmpDonerAnchorSize << endl;
		pathAnchorSizeMax_backward.push_back(tmpDonerAnchorSize);
		pathNum_backward = 1;	

		pathVec_forward.push_back(
			refineSamSEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo));
		//cout << "end of push_back pathVec_forward " << endl;
		pathSupportNumVec_forward.push_back(1);
		int tmpAcceptorAnchorSize 
			= refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo);
		//cout << "tmpAcceptorAnchorSize: " << tmpAcceptorAnchorSize << endl;	
		pathAnchorSizeMax_forward.push_back(tmpAcceptorAnchorSize);
		pathNum_forward = 1;

		anchorSizeMax_doner = tmpDonerAnchorSize;
		anchorSizeMax_acceptor = tmpAcceptorAnchorSize;

		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;	
	}

	void initiateAlignInferInfo_withRefineSAMinfo_PE(
		RefineSAM_PE_Info* refineSamPEinfo, int SJindexInRefineSAMinfo, 
		Index_Info* indexInfo)
	{
		//cout << "initiateAlignInferInfo_withRefineSAMinfo_SE starts ..." << endl;
		chrNameInt = refineSamPEinfo->returnChrNameInt();
		donerEndPos = refineSamPEinfo->returnSJ_donerEndPos(SJindexInRefineSAMinfo);
		acceptorStartPos = refineSamPEinfo->returnSJ_acceptorStartPos(SJindexInRefineSAMinfo);
		flankStringCase = getFlankStringCase(indexInfo);
		supportNum = 1;
		//cout << "chrNameInt: " << endl;
		//cout << "donerEndPos: " << donerEndPos << endl;
		//cout << "acceptorStartPos: " << acceptorStartPos << endl;
		//cout << "flankStringCase: " << flankStringCase << endl;
		pathVec_backward.push_back(
			refineSamPEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo));
		pathSupportNumVec_backward.push_back(1);
		//cout << "end of push_back pathVec_backward " << endl;
		int tmpDonerAnchorSize 
			= refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo);
		//cout << "tmpDonerAnchorSize: " << tmpDonerAnchorSize << endl;
		pathAnchorSizeMax_backward.push_back(tmpDonerAnchorSize);
		pathNum_backward = 1;	

		pathVec_forward.push_back(
			refineSamPEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo));
		//cout << "end of push_back pathVec_forward " << endl;
		pathSupportNumVec_forward.push_back(1);
		int tmpAcceptorAnchorSize 
			= refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo);
		//cout << "tmpAcceptorAnchorSize: " << tmpAcceptorAnchorSize << endl;	
		pathAnchorSizeMax_forward.push_back(tmpAcceptorAnchorSize);
		pathNum_forward = 1;

		anchorSizeMax_doner = tmpDonerAnchorSize;
		anchorSizeMax_acceptor = tmpAcceptorAnchorSize;

		donerAnchorNWDPpenalty_max = 1 + anchorSizeMax_doner/8;
		acceptorAnchorNWDPpenalty_max = 1 + anchorSizeMax_acceptor/8;	
	}

	void updateOldAlignInferJunc_withRefineSAMinfo_SE(
		RefineSAM_SE_Info* refineSamSEinfo, int SJindexInRefineSAMinfo, 
		Index_Info* indexInfo)
	{
		pair<int,int> checkWithOldPath_backward_relation_pair 
			= this->checkWithOldPathVec_backward(refineSamSEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo));
		int compliedOldPathIndex_backward = checkWithOldPath_backward_relation_pair.second;
		int checkWithOldPath_backward_relation = checkWithOldPath_backward_relation_pair.first;			
		if(checkWithOldPath_backward_relation == -1) // no need to update path
		{
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
			int tmpPathMaxAnchorSize_backward = pathAnchorSizeMax_backward[compliedOldPathIndex_backward];
			if(refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_backward)
				pathAnchorSizeMax_backward[compliedOldPathIndex_backward] 
					= refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo);
		}
		else if(checkWithOldPath_backward_relation == -2) // new path
		{
			pathVec_backward.push_back(refineSamSEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo));
			pathSupportNumVec_backward.push_back(1);
			pathAnchorSizeMax_backward.push_back(refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo));
			pathNum_backward ++;
		}
		else if(checkWithOldPath_backward_relation == 0)
		{
			pathVec_backward[compliedOldPathIndex_backward] = refineSamSEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo);
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
			int tmpPathMaxAnchorSize_backward = pathAnchorSizeMax_backward[compliedOldPathIndex_backward];
			if(refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_backward)
				pathAnchorSizeMax_backward[compliedOldPathIndex_backward] 
					= refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo);
		}
		else
		{
			cout << "checkWithOldPath_backward_relation: " << checkWithOldPath_backward_relation << endl;
			cout << "error in checkWithOldPath_backward_relation ..." << endl;
			exit(1);
		}

		pair<int,int> checkWithOldPath_forward_relation_pair
			= this->checkWithOldPathVec_forward(refineSamSEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo));
		int compliedOldPathIndex_forward = checkWithOldPath_forward_relation_pair.second;
		int checkWithOldPath_forward_relation = checkWithOldPath_forward_relation_pair.first; 
		if(checkWithOldPath_forward_relation == -1) // no need to update path
		{
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
			int tmpPathMaxAnchorSize_forward = pathAnchorSizeMax_forward[compliedOldPathIndex_forward];
			if(refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_forward)
				pathAnchorSizeMax_forward[compliedOldPathIndex_forward] 
					= refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo);
		}
		else if(checkWithOldPath_forward_relation == -2) // new path
		{
			pathVec_forward.push_back(refineSamSEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo));
			pathSupportNumVec_forward.push_back(1);
			pathAnchorSizeMax_forward.push_back(refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo));
			pathNum_forward ++;
		}
		else if(checkWithOldPath_forward_relation == 0)
		{
			pathVec_forward[compliedOldPathIndex_forward] = refineSamSEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo);
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
			int tmpPathMaxAnchorSize_forward = pathAnchorSizeMax_forward[compliedOldPathIndex_forward];
			if(refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_forward)
				pathAnchorSizeMax_forward[compliedOldPathIndex_forward] 
					= refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo);
		}
		else
		{
			cout << "checkWithOldPath_forward_relation: " << checkWithOldPath_forward_relation << endl;
			cout << "error in checkWithOldPath_forward_relation ..." << endl;
			exit(1);
		}
		this->updateMaxAnchorSize_NWDPpenaltyMax(
			refineSamSEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo), 
			refineSamSEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo));
		this->updateSupportNum();
	}

	void updateOldAlignInferJunc_withRefineSAMinfo_PE(
		RefineSAM_PE_Info* refineSamPEinfo, int SJindexInRefineSAMinfo, 
		Index_Info* indexInfo)
	{
		pair<int,int> checkWithOldPath_backward_relation_pair 
			= this->checkWithOldPathVec_backward(refineSamPEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo));
		int compliedOldPathIndex_backward = checkWithOldPath_backward_relation_pair.second;
		int checkWithOldPath_backward_relation = checkWithOldPath_backward_relation_pair.first;			
		if(checkWithOldPath_backward_relation == -1) // no need to update path
		{
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
			int tmpPathMaxAnchorSize_backward = pathAnchorSizeMax_backward[compliedOldPathIndex_backward];
			if(refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_backward)
				pathAnchorSizeMax_backward[compliedOldPathIndex_backward] 
					= refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo);
		}
		else if(checkWithOldPath_backward_relation == -2) // new path
		{
			pathVec_backward.push_back(refineSamPEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo));
			pathSupportNumVec_backward.push_back(1);
			pathAnchorSizeMax_backward.push_back(refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo));
			pathNum_backward ++;
		}
		else if(checkWithOldPath_backward_relation == 0)
		{
			pathVec_backward[compliedOldPathIndex_backward] = refineSamPEinfo->returnSJjumpCodeVec_backward(SJindexInRefineSAMinfo);
			pathSupportNumVec_backward[compliedOldPathIndex_backward] ++;
			int tmpPathMaxAnchorSize_backward = pathAnchorSizeMax_backward[compliedOldPathIndex_backward];
			if(refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_backward)
				pathAnchorSizeMax_backward[compliedOldPathIndex_backward] 
					= refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo);
		}
		else
		{
			cout << "checkWithOldPath_backward_relation: " << checkWithOldPath_backward_relation << endl;
			cout << "error in checkWithOldPath_backward_relation ..." << endl;
			exit(1);
		}

		pair<int,int> checkWithOldPath_forward_relation_pair
			= this->checkWithOldPathVec_forward(refineSamPEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo));
		int compliedOldPathIndex_forward = checkWithOldPath_forward_relation_pair.second;
		int checkWithOldPath_forward_relation = checkWithOldPath_forward_relation_pair.first; 
		if(checkWithOldPath_forward_relation == -1) // no need to update path
		{
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
			int tmpPathMaxAnchorSize_forward = pathAnchorSizeMax_forward[compliedOldPathIndex_forward];
			if(refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_forward)
				pathAnchorSizeMax_forward[compliedOldPathIndex_forward] 
					= refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo);
		}
		else if(checkWithOldPath_forward_relation == -2) // new path
		{
			pathVec_forward.push_back(refineSamPEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo));
			pathSupportNumVec_forward.push_back(1);
			pathAnchorSizeMax_forward.push_back(refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo));
			pathNum_forward ++;
		}
		else if(checkWithOldPath_forward_relation == 0)
		{
			pathVec_forward[compliedOldPathIndex_forward] = refineSamPEinfo->returnSJjumpCodeVec_forward(SJindexInRefineSAMinfo);
			pathSupportNumVec_forward[compliedOldPathIndex_forward] ++;
			int tmpPathMaxAnchorSize_forward = pathAnchorSizeMax_forward[compliedOldPathIndex_forward];
			if(refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo) > tmpPathMaxAnchorSize_forward)
				pathAnchorSizeMax_forward[compliedOldPathIndex_forward] 
					= refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo);
		}
		else
		{
			cout << "checkWithOldPath_forward_relation: " << checkWithOldPath_forward_relation << endl;
			cout << "error in checkWithOldPath_forward_relation ..." << endl;
			exit(1);
		}
		this->updateMaxAnchorSize_NWDPpenaltyMax(
			refineSamPEinfo->returnSJanchorSize_doner(SJindexInRefineSAMinfo), 
			refineSamPEinfo->returnSJanchorSize_acceptor(SJindexInRefineSAMinfo));
		this->updateSupportNum();
	}
};





#endif