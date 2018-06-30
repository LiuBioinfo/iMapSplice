// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef PAIRALIGNINFO_REFINEINNERSOFTCLIP_H
#define PAIRALIGNINFO_REFINEINNERSOFTCLIP_H

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

#include "otherFunc.h"
#include "index_info.h"
#include "splice_info.h"
#include "nw_DP.h"

using namespace std;

class PairAlignInfo_RefineInnerSoftClip
{
private:
	string originalSamFieldStr_for;
	string originalSamFieldStr_rcm;

	string readName;
	int chrNameInt;

	// before fixing forward_tail, original SAM info
	int flag_for;
	int chrMapPos_for;
	int chrMapPos_end_for;
	vector<Jump_Code> alignJumpCodeVec_for;
	string readSeq_for;
	int readLength_for;
	string qualSeq_for;
	string otherSamFieldStr_for;

	// after fixing forward_tail, original SAM info
	int chrMapPos_for_final;
	int chrMapPos_end_for_final;
	vector<Jump_Code> alignJumpCodeVec_for_final;

	// before fixing reverseComplement_head, original SAM info
	int flag_rcm;
	int chrMapPos_rcm;
	int chrMapPos_end_rcm;
	vector<Jump_Code> alignJumpCodeVec_rcm;
	string readSeq_rcm;
	int readLength_rcm;
	string qualSeq_rcm;
	string otherSamFieldStr_rcm;

	// after fixing reverseComplement_head, original SAM info
	int chrMapPos_rcm_final;
	int chrMapPos_end_rcm_final;
	vector<Jump_Code> alignJumpCodeVec_rcm_final;

	bool strandSJ_exits_bool;
	bool positiveStrandSJ_exists_bool;
	bool negativeStrandSJ_exists_bool;

	// before fixing forward_tail
	int unfixedTailLen_for;
	int mappedPartEndPos_inChr_for;
	int mappedPartEndLoc_InRead_for;
	// after fixing forward_tail
	vector<Jump_Code> fixedTailJumpCodeVec_for;
	vector<int> fixedTailMismatchLocVec_inTailSeq_for; // mismatch pos in tail sequence 
	vector<char> fixedTailMismatchCharVec_for;
	bool tail_for_fixed_bool;

	// before fixing reverseComplement_head
	int unfixedHeadLen_rcm;
	int mappedPartStartPos_inChr_rcm;
	int mappedPartStartLoc_inRead_rcm;
	// after fixing reverseComplement_head
	vector<Jump_Code> fixedHeadJumpCodeVec_rcm;
	vector<int> fixedHeadMismatchLocVec_inHeadSeq_rcm;
	vector<char> fixedHeadMismatchCharVec_rcm;
	bool head_rcm_fixed_bool;


	int toSearchAnchorSeqLen_max;
public:
	PairAlignInfo_RefineInnerSoftClip()
	{
		unfixedTailLen_for = 0;
		unfixedHeadLen_rcm = 0;
		tail_for_fixed_bool = false;
		head_rcm_fixed_bool = false;
		strandSJ_exits_bool = false;
		positiveStrandSJ_exists_bool = false;
		negativeStrandSJ_exists_bool = false;

		toSearchAnchorSeqLen_max = 5;
	}

	void generateFinalAlignInfo_afterFixing_forTail()
	{
		if(tail_for_fixed_bool)
		{
			//cout << "alignJumpCodeVec_for: " << this->jumpCodeVec2cigarString(alignJumpCodeVec_for) << endl;
			//cout << "fixedTailJumpCodeVec_for: " << this->jumpCodeVec2cigarString(fixedTailJumpCodeVec_for) << endl;
			chrMapPos_for_final = chrMapPos_for;
			this->mergeJumpCodeVec(alignJumpCodeVec_for,
				0, alignJumpCodeVec_for.size()-2, 
				fixedTailJumpCodeVec_for,
				0, fixedTailJumpCodeVec_for.size()-1,
				alignJumpCodeVec_for_final);
			//cout << "alignJumpCodeVec_for_final: " << this->jumpCodeVec2cigarString(alignJumpCodeVec_for_final) << endl;
			chrMapPos_end_for_final = this->getEndPosOfSpecificJumpCode(
				chrMapPos_for_final, alignJumpCodeVec_for_final, alignJumpCodeVec_for_final.size()-1);
		}
	}

	void generateFinalAlignInfo_afterFixing_rcmHead()
	{
		if(head_rcm_fixed_bool)
		{
			chrMapPos_end_rcm_final = chrMapPos_end_rcm;
			this->mergeJumpCodeVec(fixedHeadJumpCodeVec_rcm,
				0, fixedHeadJumpCodeVec_rcm.size()-1,
				alignJumpCodeVec_rcm,
				1, alignJumpCodeVec_rcm.size()-1,
				alignJumpCodeVec_rcm_final);
			int tmpFixedHeadJumpCodeVecPosRange = this->jumpCodeVecPosRange(
				fixedHeadJumpCodeVec_rcm);
			chrMapPos_rcm_final = chrMapPos_rcm - tmpFixedHeadJumpCodeVecPosRange;
		}
	}

	void selectBestShortAnchorSJ_completeInnerSoftClip_for_tail(
		vector<int>& foundSJstartLocVec_positiveStrand,
		vector<int>& foundSegMapPosVec_positiveStrand,
		vector<int>& foundSJstartLocVec_negativeStrand,
		vector<int>& foundSegMapPosVec_negativeStrand)
	{
		this->selectBestShortAnchorSJ_completeInnerSoftClip_for_tail_select1stOne(
			foundSJstartLocVec_positiveStrand,
			foundSegMapPosVec_positiveStrand,
			foundSJstartLocVec_negativeStrand,
			foundSegMapPosVec_negativeStrand);
	}

	void selectBestShortAnchorSJ_completeInnerSoftClip_for_tail_select1stOne(
		vector<int>& foundSJstartLocVec_positiveStrand,
		vector<int>& foundSegMapPosVec_positiveStrand,
		vector<int>& foundSJstartLocVec_negativeStrand,
		vector<int>& foundSegMapPosVec_negativeStrand)
	{
		if(foundSJstartLocVec_positiveStrand.size() 
			+ foundSJstartLocVec_negativeStrand.size() > 0)
		{
			int tmpSJstartLocInRead_1st;
			int tmpSJanchorMapPosInChr_1st;
			if(foundSJstartLocVec_positiveStrand.size() > 0)
			{
				tmpSJstartLocInRead_1st = foundSJstartLocVec_positiveStrand[0];
				tmpSJanchorMapPosInChr_1st = foundSegMapPosVec_positiveStrand[0];
			}
			else
			{
				tmpSJstartLocInRead_1st = foundSJstartLocVec_negativeStrand[0];
				tmpSJanchorMapPosInChr_1st = foundSegMapPosVec_negativeStrand[0];
			}

			int tmpSJanchorSeqLength = readLength_for - tmpSJstartLocInRead_1st + 1;
			int tmpExtendedSeqLength_mappedMidPart = unfixedTailLen_for - tmpSJanchorSeqLength;
			int tmpSJsize = tmpSJanchorMapPosInChr_1st 
			    - (chrMapPos_end_for + tmpExtendedSeqLength_mappedMidPart) - 1;

			// extended match 
			if(tmpExtendedSeqLength_mappedMidPart > 0)
			{
				Jump_Code tmpExtendedMatchJumpCode(tmpExtendedSeqLength_mappedMidPart, "M");
				fixedTailJumpCodeVec_for.push_back(tmpExtendedMatchJumpCode);
			}

			// SJ 
			Jump_Code tmpSJjumpCode(tmpSJsize, "N");
			fixedTailJumpCodeVec_for.push_back(tmpSJjumpCode);

			// short SJanchor
			if(tmpSJanchorSeqLength <= toSearchAnchorSeqLen_max)
			{		
				Jump_Code tmpSJshortAnchorMatchJumpCode(tmpSJanchorSeqLength, "M");
				fixedTailJumpCodeVec_for.push_back(tmpSJshortAnchorMatchJumpCode);
			}
			else // long SJanchor
			{
				// check right sequence of toSearchSJanchorSeq -- Match / NWDP / SoftClip
				int tmpToCheckReadSeq_startLocInRead = 
					tmpSJstartLocInRead_1st + toSearchAnchorSeqLen_max;
				string tmpToCheckReadSeq = readSeq_for.substr(tmpToCheckReadSeq_startLocInRead - 1);
				int tmpToCheckChrSeq_startPosInChr = tmpSJanchorMapPosInChr_1st + toSearchAnchorSeqLen_max;
				int tmpToCheckChrSeq_endPosInChr = tmpSJanchorMapPosInChr_1st + tmpSJanchorSeqLength - 1;
				string tmpToCheckChrSeq_match = indexInfo->returnChromStrSubstr(chrNameInt,
					tmpToCheckChrSeq_startPosInChr, 
					tmpToCheckChrSeq_endPosInChr - tmpToCheckChrSeq_startPosInChr + 1);
				bool matchOrNot_bool = 
				if(matchOrNot_bool) // other read tail seq matched 
				{
					Jump_Code tmpSJshortAnchorMatchJumpCode(tmpSJanchorSeqLength, "M");
					fixedTailJumpCodeVec_for.push_back(tmpSJshortAnchorMatchJumpCode);					
				}
				else // unmatched, try nwdp 
				{
					string tmpToCheckChrSeq_nwdp = indexInfo->returnChromStrSubstr(chrNameInt,
						tmpToCheckChrSeq_startPosInChr, 
						tmpToCheckChrSeq_endPosInChr - tmpToCheckChrSeq_startPosInChr + 1 + 3);
					FixSingleAnchor_NWDP_Info* fixSingleAnchorNWDPinfo = new FixSingleAnchor_NWDP_Info();
					fixSingleAnchorNWDPinfo->doNWDP(tmpToCheckReadSeq, tmpToCheckChrSeq_nwdp);
					bool fixNWDP_bool = fixSingleAnchorNWDPinfo->fixedOrNot();
					if(fixNWDP_bool) // nwdp fixed
					{
						vector<Jump_Code> tmpMatchJumpCodeVec;
						Jump_Code onlyToSearchSJanchorMatchJumpCode(toSearchAnchorSeqLen_max, "M");
						tmpMatchJumpCodeVec.push_back(onlyToSearchSJanchorMatchJumpCode);

						vector<Jump_Code> tmpNwdpJumpCodeVec;
						fixSingleAnchorNWDPinfo->copyJumpCodeVec2TargetVec(tmpNwdpJumpCodeVec);

						//vector<Jump_Code> tmpTotalJumpCodeVec;
						this->mergeJumpCodeVec(tmpMatchJumpCodeVec, 0, 0,
							tmpNwdpJumpCodeVec, 0, tmpNwdpJumpCodeVec.size()-1,
							fixedTailJumpCodeVec_for);
					}
					else // nwdp unfixed
					{
						Jump_Code onlyToSearchSJanchorMatchJumpCode(toSearchAnchorSeqLen_max, "M");
						fixedTailJumpCodeVec_for.push_back(onlyToSearchSJanchorMatchJumpCode);
						int tmpSoftClipTailSeqLength = tmpSJanchorSeqLength - toSearchAnchorSeqLen_max;
						Jump_Code tailSoftClipJumpCode(tmpSoftClipTailSeqLength, "S");
						fixedTailJumpCodeVec_for.push_back(tailSoftClipJumpCode);
					}
				}
			}
			tail_for_fixed_bool = true;
		}
		else
		{}
	}

	void selectBestShortAnchorSJ_completeInnerSoftClip_rcm_head(
		vector<int>& foundSJendLocVec_positiveStrand,
		vector<int>& foundSegMapPosVec_positiveStrand,
		vector<int>& foundSJendLocVec_negativeStrand,
		vector<int>& foundSegMapPosVec_negativeStrand)
	{
		this->selectBestShortAnchorSJ_completeInnerSoftClip_rcm_head_select1stOne(
			foundSJendLocVec_positiveStrand,
			foundSegMapPosVec_positiveStrand,
			foundSJendLocVec_negativeStrand,
			foundSegMapPosVec_negativeStrand);
	}

	void selectBestShortAnchorSJ_completeInnerSoftClip_rcm_head_select1stOne(
		vector<int>& foundSJendLocVec_positiveStrand,
		vector<int>& foundSegMapPosVec_positiveStrand,
		vector<int>& foundSJendLocVec_negativeStrand,
		vector<int>& foundSegMapPosVec_negativeStrand)
	{
		if(foundSJendLocVec_positiveStrand.size() 
			+ foundSJendLocVec_negativeStrand.size() > 0)
		{
			int tmpSJendLocInRead_1st;
			int tmpSJanchorMapPosInChr_1st;
			if(foundSJendLocVec_positiveStrand.size() > 0)
			{
				tmpSJendLocInRead_1st = foundSJendLocVec_positiveStrand[0];
				tmpSJanchorMapPosInChr_1st = foundSegMapPosVec_positiveStrand[0];
			}
			else
			{
				tmpSJendLocInRead_1st = foundSJendLocVec_negativeStrand[0];
				tmpSJanchorMapPosInChr_1st = foundSegMapPosVec_negativeStrand[0];
			}
			int tmpSJanchorSeqLength = tmpSJendLocInRead_1st;
			int tmpExtendedSeqLength_mappedMidPart = unfixedHeadLen_rcm - tmpSJanchorSeqLength;
			int tmpSJsize = (chrMapPos_rcm - tmpExtendedSeqLength_mappedMidPart) 
			    - (tmpSJanchorMapPosInChr_1st + tmpSJanchorSeqLength - 1) - 1;
			
			Jump_Code tmpSJshortAnchorMatchJumpCode(tmpSJanchorSeqLength, "M");
			fixedHeadJumpCodeVec_rcm.push_back(tmpSJshortAnchorMatchJumpCode);			
			Jump_Code tmpSJjumpCode(tmpSJsize, "N");
			fixedHeadJumpCodeVec_rcm.push_back(tmpSJjumpCode);
			if(tmpExtendedSeqLength_mappedMidPart > 0)
			{
				Jump_Code tmpExtendedMatchJumpCode(tmpExtendedSeqLength_mappedMidPart, "M");
				fixedHeadJumpCodeVec_rcm.push_back(tmpExtendedMatchJumpCode);
			}
			tail_for_fixed_bool = true;
		}
		else
		{}
	}	

	void determineStrand(Index_Info* indexInfo)
	{
		vector<string> flankStringVec;
		for(int tmp = 0; tmp < alignJumpCodeVec_for.size(); tmp++)
		{
			string tmpJumpCodeType = alignJumpCodeVec_for[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int tmpSJstartPos = this->getEndPosOfSpecificJumpCode(chrMapPos_for,
					alignJumpCodeVec_for, tmp-1) + 1;
				int tmpSJendPos = this->getEndPosOfSpecificJumpCode(chrMapPos_for,
					alignJumpCodeVec_for, tmp);
				string tmpSJflankString = indexInfo->returnFlankString(chrNameInt,
					tmpSJstartPos, tmpSJendPos);
				flankStringVec.push_back(tmpSJflankString);
			}
		}

		for(int tmp = 0; tmp < alignJumpCodeVec_rcm.size(); tmp++)
		{
			string tmpJumpCodeType = alignJumpCodeVec_rcm[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int tmpSJstartPos = this->getEndPosOfSpecificJumpCode(chrMapPos_rcm,
					alignJumpCodeVec_rcm, tmp-1) + 1;
				int tmpSJendPos = this->getEndPosOfSpecificJumpCode(chrMapPos_rcm,
					alignJumpCodeVec_rcm, tmp);
				string tmpSJflankString = indexInfo->returnFlankString(chrNameInt,
					tmpSJstartPos, tmpSJendPos);
				flankStringVec.push_back(tmpSJflankString);
			}
		}

		for(int tmp = 0; tmp < flankStringVec.size(); tmp++)
		{
			string tmpSJflankString = flankStringVec[tmp];
			if((tmpSJflankString == "GTAG")
				||(tmpSJflankString == "GTAT")
				||(tmpSJflankString == "GCAG"))
			{
				positiveStrandSJ_exists_bool = true;
				strandSJ_exits_bool = true;
			}
			else if((tmpSJflankString == "CTAC")
				||(tmpSJflankString == "ATAC")
				||(tmpSJflankString == "CTGC"))
			{
				negativeStrandSJ_exists_bool = true;
				strandSJ_exits_bool = true;
			}
			else
			{}
		}
	}

	int jumpCodeVecPosRange(vector<Jump_Code>& tmpJumpCodeVec)
	{
		int tmpRange = 0;
		for(int tmp = 0; tmp < tmpJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = tmpJumpCodeVec[tmp].type; 
			int tmpJumpCodeLength = tmpJumpCodeVec[tmp].len;
			if((tmpJumpCodeType == "M")||(tmpJumpCodeType == "N")||(tmpJumpCodeType == "D"))
			{
				tmpRange += tmpJumpCodeLength;
			}
		}
		return tmpRange;
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
			
		string candidateJumpCodeType = "SMNIDX";
		while(1)
		{
			jumpCodeEndPosInCigarStr = 
				jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
			if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
				{break;}
			else
			{
				tmpJumpCodeLength = 
					atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
				tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
				cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	string jumpCodeVec2cigarString(vector<Jump_Code>& alignJumpCodeVec)
	{
		string tmpCigarString;
		//cout << "********" << "jumpCodeVecSize: " << jumpCodeVec.size() << endl;
		for(int tmp = 0; tmp < alignJumpCodeVec.size(); tmp++)
		{
			tmpCigarString += alignJumpCodeVec[tmp].toString();
		}
		return tmpCigarString;
	}	

	int getEndLocInReadOfSpecificJumpCode(
		vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		if(jumpCodeIndex < 0)
			return 0;
		int tmpEndLocInRead = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
				tmpEndLocInRead += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "D")
			{
			}
			else if(tmpJumpCodeType == "N")
			{
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return tmpEndLocInRead;
	}	

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
		int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
			if(tmpJumpCodeType == "S")
			{
			}
			else if(tmpJumpCodeType == "M")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "I")
			{
			}
			else if(tmpJumpCodeType == "D")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else if(tmpJumpCodeType == "N")
			{
				tmpEndPos += tmpJumpCodeLength;
			}
			else
			{
				cout << "incorrect jumpCode type" << endl;
				exit(1);
			}								
		}
		return (tmpEndPos + startPos-1);
	}

	void initiateWithSamStr(const string& samStr_for,
		const string& samStr_rcm, Index_Info* indexInfo)
	{
		// original sam strings
		originalSamFieldStr_for = samStr_for;
		originalSamFieldStr_rcm = samStr_rcm;

		vector<string> samFieldVec_for;
		vector<string> samFieldVec_rcm;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr_for.find("\t", startLoc);
			string tmpSamField = samStr_for.substr(startLoc, tabLoc-startLoc);
			samFieldVec_for.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_for.push_back(samStr_for.substr(startLoc));

		startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr_rcm.find("\t", startLoc);
			string tmpSamField = samStr_rcm.substr(startLoc, tabLoc-startLoc);
			samFieldVec_rcm.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}		
		samFieldVec_rcm.push_back(samStr_rcm.substr(startLoc));

		// read name, chr name int //
		readName = samFieldVec_for[0];
		string mapChrNameStr_for = samFieldVec_for[2];
		chrNameInt = indexInfo->convertStringToInt(mapChrNameStr_for);
		
		// sam_info_forward
		flag_for = atoi(samFieldVec_for[1].c_str());
		string mapChrPosStr_for = samFieldVec_for[3];
		chrMapPos_for = atoi(mapChrPosStr_for.c_str()); // chrMapPos_for
		string cigarString_for = samFieldVec_for[5];
		this->cigarString2jumpCodeVec(cigarString_for, alignJumpCodeVec_for); // jumpCodeVec_for		
		int jumpCodeVecSize_for = alignJumpCodeVec_for.size();
		chrMapPos_end_for = this->getEndPosOfSpecificJumpCode(chrMapPos_for,
			alignJumpCodeVec_for, jumpCodeVecSize_for-1);
		readSeq_for = samFieldVec_for[9];
		readLength_for = readSeq_for.length();
		qualSeq_for = samFieldVec_for[10];
		otherSamFieldStr_for = samFieldVec_for[11];

		// sam_info_reverseComplement
		flag_rcm = atoi(samFieldVec_rcm[1].c_str());
		string mapChrPosStr_rcm = samFieldVec_rcm[3];
		chrMapPos_rcm = atoi(mapChrPosStr_rcm.c_str()); // chrMapPos_for
		string cigarString_rcm = samFieldVec_rcm[5];
		this->cigarString2jumpCodeVec(cigarString_rcm, alignJumpCodeVec_rcm); // jumpCodeVec_for
		int jumpCodeVecSize_rcm = alignJumpCodeVec_rcm.size();
		chrMapPos_end_rcm = this->getEndPosOfSpecificJumpCode(chrMapPos_rcm,
			alignJumpCodeVec_rcm, jumpCodeVecSize_rcm-1);
		readSeq_rcm = samFieldVec_rcm[9];
		readLength_rcm = readSeq_rcm.length();
		qualSeq_rcm = samFieldVec_rcm[10];
		otherSamFieldStr_rcm = samFieldVec_rcm[11];

		// unfixedTailLen_for, unfixedHeadLen_rcm
		if(alignJumpCodeVec_for[jumpCodeVecSize_for-1].type == "S")
		{
			unfixedTailLen_for = alignJumpCodeVec_for[jumpCodeVecSize_for-1].len;
			mappedPartEndPos_inChr_for = this->getEndPosOfSpecificJumpCode(chrMapPos_for,
				alignJumpCodeVec_for, jumpCodeVecSize_for-2);
			mappedPartEndLoc_InRead_for = this->getEndLocInReadOfSpecificJumpCode(
				alignJumpCodeVec_for, jumpCodeVecSize_for-2);
		}
		else
			unfixedTailLen_for = 0;

		if(alignJumpCodeVec_rcm[0].type == "S")
		{
			unfixedHeadLen_rcm = alignJumpCodeVec_rcm[0].len;
			mappedPartStartPos_inChr_rcm = chrMapPos_rcm;
			mappedPartStartLoc_inRead_rcm = unfixedHeadLen_rcm + 1;
		}
		else
			unfixedHeadLen_rcm = 0;

		// determine strand
		this->determineStrand(indexInfo);
	}

	bool innerSoftClip_exists_bool()
	{
		if((unfixedTailLen_for > 0)||(unfixedHeadLen_rcm > 0))
			return true;
		else
			return false;
	}

	bool innerSoftClip_for_tail_exists_bool()
	{
		return (unfixedTailLen_for > 0);
	}

	bool innerSoftClip_rcm_head_exists_bool()
	{
		return (unfixedHeadLen_rcm > 0);
	}

	bool innerSoftClip_for_tail_fixed_bool()
	{
		return tail_for_fixed_bool;
	}

	bool innerSoftClip_rcm_head_fixed_bool()
	{
		return head_rcm_fixed_bool;
	}

	int returnMatePairDistance()
	{
		return (chrMapPos_rcm - chrMapPos_end_for);
	}

	void mergeJumpCodeVec(vector<Jump_Code>& jumpCodeVec_1,
		int startIndex_inJumpCodeVec_1, int endIndex_inJumpCodeVec_1, 
		vector<Jump_Code>& jumpCodeVec_2,
		int startIndex_inJumpCodeVec_2, int endIndex_inJumpCodeVec_2,
		vector<Jump_Code>& jumpCodeVec_target)
	{
		for(int tmp = startIndex_inJumpCodeVec_1; tmp <= endIndex_inJumpCodeVec_1;
			tmp++)
		{
			jumpCodeVec_target.push_back(jumpCodeVec_1[tmp]);
		}
		string lastTypeInCurrentJumpCodeVec = jumpCodeVec_target[jumpCodeVec_target.size()-1].type;
		string nextType = jumpCodeVec_2[startIndex_inJumpCodeVec_2].type;
		if(nextType == lastTypeInCurrentJumpCodeVec)
		{
			jumpCodeVec_target[jumpCodeVec_target.size()-1].len += jumpCodeVec_2[startIndex_inJumpCodeVec_2].len;
		}
		else
		{
			jumpCodeVec_target.push_back(jumpCodeVec_2[startIndex_inJumpCodeVec_2]);
		}
		for(int tmp = startIndex_inJumpCodeVec_2 + 1; tmp <= endIndex_inJumpCodeVec_2;
			tmp++)
		{
			jumpCodeVec_target.push_back(jumpCodeVec_2[tmp]);
		}
	}

	string returnSamStr_fixedTail_for(Index_Info* indexInfo)
	{
		//vector<Jump_Code> finalJumpCodeVec;
		//this->mergeJumpCodeVec(alignJumpCodeVec_for, 0, alignJumpCodeVec_for.size()-2,
		//	fixedTailJumpCodeVec_for, 0, fixedTailJumpCodeVec_for.size()-1, finalJumpCodeVec);
		string tmpCigarString = jumpCodeVec2cigarString(alignJumpCodeVec_for_final);
		string tmpStr;
		tmpStr = readName + "\t" + int_to_str(flag_for) + "\t"
			+ indexInfo->returnChrNameStr(chrNameInt) + "\t"
			+ int_to_str(chrMapPos_for_final) + "\t*\t" + tmpCigarString + "\t*\t*\t"
			+ readSeq_for + "\t" + qualSeq_for + "\t" + otherSamFieldStr_for;
		return tmpStr;
	}

	string returnSamStr_fixedHead_rcm(Index_Info* indexInfo)
	{
		// vector<Jump_Code> finalJumpCodeVec;
		// this->mergeJumpCodeVec(fixedHeadJumpCodeVec_rcm, 0, fixedHeadJumpCodeVec_rcm.size()-1,
		// 	alignJumpCodeVec_rcm, 1, alignJumpCodeVec_rcm.size()-1, finalJumpCodeVec);
		string tmpCigarString = jumpCodeVec2cigarString(alignJumpCodeVec_rcm_final);
		// int fixedHeadJumpCodeVecRange = this->jumpCodeVecPosRange(fixedHeadJumpCodeVec_rcm);
		// int finalMapPos = chrMapPos_rcm - fixedHeadJumpCodeVecRange;
		string tmpStr;
		tmpStr = readName + "\t" + int_to_str(flag_rcm) + "\t"
			+ indexInfo->returnChrNameStr(chrNameInt) + "\t"
			+ int_to_str(chrMapPos_rcm_final) + "\t*\t" + tmpCigarString + "\t*\t*\t"
			+ readSeq_rcm + "\t" + qualSeq_rcm + "\t" + otherSamFieldStr_rcm;
		return tmpStr;
	}

	string returnSamStr_original_for()
	{
		return originalSamFieldStr_for;
	}

	string returnSamStr_original_rcm()
	{
		return originalSamFieldStr_rcm;
	}

	/*int return_type_of_method_toApply(Index_Info* indexInfo)
	{
		int type; // 1--extension; 2--remapping with current splicing structures; 3--search new SJ between two ends.

		return type;
	}*/

	void fix_softClip_tail_for(Index_Info* indexInfo)
	{
		bool doExtensionOrNot_bool = this->doExtensionOrNot(indexInfo);
	}

	bool doExtensionOrNot(Index_Info* indexInfo)
	{
		return false;
	}

	bool doRemappingOrNot(Index_Info* indexInfo, 
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{	
		return false;
	}

	bool doNovelSJsearchOrNot(Index_Info* indexInfo,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{
		return false;
	}

	void fix_doNovelSJsearch_tail_for(Index_Info* indexInfo,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		SJhash_Info* SJhashInfo)
	{
		int backCheckBuffer = 7;
		int anchorSeq_min = 3;
		int maxHit_GT_CT = 1;
		int maxHit_anchorMap = 1;
		int exonicSearchSpace = 200;

		//int toSearchAnchorSeqLen_max = 5;

		vector<int> candidateSJstartLocInReadVec_GT; 
		vector<int> candidateSJstartLocInReadVec_CT;
		//cout << "start to search candidate SJstartLocInRead" << endl;
		//cout << "searchStartChrPos: " << chrMapPos_end_for - backCheckBuffer + 1 << endl;
		//cout << "searchEndChrPos: " << chrMapPos_end_for + unfixedTailLen_for - anchorSeq_min << endl;
		//cout << "searchTargetSeqLength: " << chrMapPos_end_for + unfixedTailLen_for - anchorSeq_min
		//	- (chrMapPos_end_for - backCheckBuffer + 1) + 1 << endl;
		//cout << "baseLoc2add: " << readLength_for - unfixedTailLen_for + 1 - backCheckBuffer << endl;

		string GTstr = "GT";
		if((positiveStrandSJ_exists_bool)||(!strandSJ_exits_bool))
		{	
			this->segmentSearch_shortRange(candidateSJstartLocInReadVec_GT, chrNameInt,
				chrMapPos_end_for - backCheckBuffer + 1, 
				chrMapPos_end_for + unfixedTailLen_for - anchorSeq_min,
				GTstr, indexInfo, maxHit_GT_CT, readLength_for - unfixedTailLen_for + 1 - backCheckBuffer);
		}
		string CTstr = "CT";
		if((negativeStrandSJ_exists_bool)||(!strandSJ_exits_bool))
		{
			this->segmentSearch_shortRange(candidateSJstartLocInReadVec_CT, chrNameInt,
				chrMapPos_end_for - backCheckBuffer + 1, 
				chrMapPos_end_for + unfixedTailLen_for - anchorSeq_min,
				CTstr, indexInfo, maxHit_GT_CT, readLength_for - unfixedTailLen_for + 1 - backCheckBuffer);
		}
		/*
		cout << "segmentSearch_shortRange starts ......" << endl;
		cout << "candidateSJstartLocInReadVec_GT: " << endl;
		for(int tmp = 0; tmp < candidateSJstartLocInReadVec_GT.size(); tmp++)
		{
			cout << "candiSJstartLocInRead: " << candidateSJstartLocInReadVec_GT[tmp] << endl;
		}
		cout << "candidateSJstartLocInReadVec_CT: " << endl;
		for(int tmp = 0; tmp < candidateSJstartLocInReadVec_CT.size(); tmp++)
		{
			cout << "candiSJstartLocInRead: " << candidateSJstartLocInReadVec_CT[tmp] << endl;
		}
		*/

		int closestSJstartPos;
		if(candidateSJstartLocInReadVec_GT.size() + candidateSJstartLocInReadVec_CT.size() > 0)
		{
			// closestSJstartPos = SJhashInfo->returnLeftMostSJdonnerEndPos(
			// 	chrNameInt, chrMapPos_end_for, chrMapPos_rcm);

			vector<int> foundSJstartLocVec_positiveStrand;
			vector<int> foundSegMapPosVec_positiveStrand;
			vector<int> foundSJstartLocVec_negativeStrand; // true readAnchorSeq map position
			vector<int>	foundSegMapPosVec_negativeStrand;

			int searchEndPosInChr = closestSJstartPos;
			//if(searchEndPosInChr > chrMapPos_rcm)
				searchEndPosInChr = chrMapPos_rcm;
			int searchStartPosInChr = searchEndPosInChr - exonicSearchSpace;
			if(searchStartPosInChr < chrMapPos_end_for)
				searchStartPosInChr = chrMapPos_end_for;

			// cout << "start to search for short anchor SJs " << endl;
			// cout << "searchStartPosInChr: " << searchStartPosInChr << endl;
			// cout << "searchEndPosInChr: " << searchEndPosInChr << endl; 
			// cout << "searchRange: " << searchEndPosInChr - searchStartPosInChr + 1 << endl;

			for(int tmp = 0; tmp < candidateSJstartLocInReadVec_GT.size(); tmp++)
			{
				int tmpToSearchSJanchorReadSeqLength 
					= readLength_for - candidateSJstartLocInReadVec_GT[tmp] + 1;
				if(tmpToSearchSJanchorReadSeqLength > toSearchAnchorSeqLen_max)
					tmpToSearchSJanchorReadSeqLength = toSearchAnchorSeqLen_max;

				vector<int> tmpFoundSegMapPosVec;
				string tmp2searchReadAnchorSeq = "AG" + readSeq_for.substr(
					candidateSJstartLocInReadVec_GT[tmp]-1, tmpToSearchSJanchorReadSeqLength);
				this->segmentSearch_longRange(tmpFoundSegMapPosVec,
					chrNameInt, searchStartPosInChr, searchEndPosInChr, 
					tmp2searchReadAnchorSeq, indexInfo, maxHit_anchorMap, 
					searchStartPosInChr+2);
				for(int tmp2 = 0; tmp2 < tmpFoundSegMapPosVec.size(); tmp2++)
				{
					foundSJstartLocVec_positiveStrand.push_back(
						candidateSJstartLocInReadVec_GT[tmp]);
					foundSegMapPosVec_positiveStrand.push_back(
						tmpFoundSegMapPosVec[tmp2]);
				}
			}
			for(int tmp = 0; tmp < candidateSJstartLocInReadVec_CT.size(); tmp++)
			{
				int tmpToSearchSJanchorReadSeqLength 
					= readLength_for - candidateSJstartLocInReadVec_CT[tmp] + 1;
				if(tmpToSearchSJanchorReadSeqLength > toSearchAnchorSeqLen_max)
					tmpToSearchSJanchorReadSeqLength = toSearchAnchorSeqLen_max;

				vector<int> tmpFoundSegMapPosVec;
				string tmp2searchReadAnchorSeq = "AC" + readSeq_for.substr(
					candidateSJstartLocInReadVec_CT[tmp]-1, tmpToSearchSJanchorReadSeqLength);
				this->segmentSearch_longRange(tmpFoundSegMapPosVec,
					chrNameInt, searchStartPosInChr, searchEndPosInChr, 
					tmp2searchReadAnchorSeq, indexInfo, maxHit_anchorMap, 
					searchStartPosInChr+2);			
				for(int tmp2 = 0; tmp2 < tmpFoundSegMapPosVec.size(); tmp2++)
				{
					foundSJstartLocVec_negativeStrand.push_back(
						candidateSJstartLocInReadVec_CT[tmp]);
					foundSegMapPosVec_negativeStrand.push_back(
						tmpFoundSegMapPosVec[tmp2]);
				}
			}

			this->selectBestShortAnchorSJ_completeInnerSoftClip_for_tail(
				foundSJstartLocVec_positiveStrand,
				foundSegMapPosVec_positiveStrand,
				foundSJstartLocVec_negativeStrand,
				foundSegMapPosVec_negativeStrand);
		}
	}

	void fix_doNovelSJsearch_head_rcm(Index_Info* indexInfo,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		SJhash_Info* SJhashInfo)
	{
		int forwardCheckBuffer = 7;
		int anchorSeq_min = 3;
		int maxHit_AG_AC = 1;
		int maxHit_anchorMap = 1;
		int exonicSearchSpace = 200;

		vector<int> candidateSJendLocInReadVec_AG; 
		vector<int> candidateSJendLocInReadVec_AC;

		string AGstr = "AG";
		if((positiveStrandSJ_exists_bool)||(!strandSJ_exits_bool))
		{	
			this->segmentSearch_shortRange(candidateSJendLocInReadVec_AG, chrNameInt,
				chrMapPos_rcm - unfixedHeadLen_rcm + anchorSeq_min - 1, 
				chrMapPos_rcm + forwardCheckBuffer - 2,
				AGstr, indexInfo, maxHit_AG_AC, 
				chrMapPos_rcm - unfixedHeadLen_rcm + anchorSeq_min);
		}
		string ACstr = "AC";
		if((negativeStrandSJ_exists_bool)||(!strandSJ_exits_bool))
		{
			this->segmentSearch_shortRange(candidateSJendLocInReadVec_AC, chrNameInt,
				chrMapPos_rcm - unfixedHeadLen_rcm + anchorSeq_min - 1, 
				chrMapPos_rcm + forwardCheckBuffer - 2,
				ACstr, indexInfo, maxHit_AG_AC, 
				chrMapPos_rcm - unfixedHeadLen_rcm + anchorSeq_min);
		}

		int closestSJendPos;
		if(candidateSJendLocInReadVec_AG.size() + candidateSJendLocInReadVec_AC.size() > 0)
		{	
			// closestSJendPos = SJhashInfo->returnRightMostSJacceptorStartPos(
			// 	chrNameInt, chrMapPos_end_for, chrMapPos_rcm);

			vector<int> foundSJendLocVec_positiveStrand;
			vector<int> foundSegMapPosVec_positiveStrand;
			vector<int> foundSJendLocVec_negativeStrand; // true readAnchorSeq map position
			vector<int>	foundSegMapPosVec_negativeStrand;

			int searchStartPosInChr = closestSJendPos;
			//if(searchStartPosInChr < chrMapPos_end_for)
				searchStartPosInChr = chrMapPos_end_for;
			int searchEndPosInChr = searchStartPosInChr + exonicSearchSpace;
			for(int tmp = 0; tmp < candidateSJendLocInReadVec_AG.size(); tmp++)
			{
				int tmpToSearchSJanchorReadSeqLength = candidateSJendLocInReadVec_AG[tmp];
				if(tmpToSearchSJanchorReadSeqLength > toSearchAnchorSeqLen_max)
					tmpToSearchSJanchorReadSeqLength = toSearchAnchorSeqLen_max;

				vector<int> tmpFoundSegMapPosVec;
				// string tmp2searchReadAnchorSeq = readSeq_rcm.substr(0,
				// 	candidateSJendLocInReadVec_AG[tmp]) + "GT";
				string tmp2searchReadAnchorSeq = readSeq_rcm.substr(
					candidateSJendLocInReadVec_AG[tmp] - tmpToSearchSJanchorReadSeqLength, 
					tmpToSearchSJanchorReadSeqLength) + "GT";
				this->segmentSearch_longRange(tmpFoundSegMapPosVec,
					chrNameInt, searchStartPosInChr, searchEndPosInChr, 
					tmp2searchReadAnchorSeq, indexInfo, maxHit_anchorMap, 
					searchStartPosInChr);
				for(int tmp2 = 0; tmp2 < tmpFoundSegMapPosVec.size(); tmp2++)
				{
					foundSJendLocVec_positiveStrand.push_back(
						candidateSJendLocInReadVec_AG[tmp]);
					foundSegMapPosVec_positiveStrand.push_back(
						tmpFoundSegMapPosVec[tmp2]);
				}
			}
			for(int tmp = 0; tmp < candidateSJendLocInReadVec_AC.size(); tmp++)
			{
				int tmpToSearchSJanchorReadSeqLength = candidateSJendLocInReadVec_AC[tmp];
				if(tmpToSearchSJanchorReadSeqLength > toSearchAnchorSeqLen_max)
					tmpToSearchSJanchorReadSeqLength = toSearchAnchorSeqLen_max;

				vector<int> tmpFoundSegMapPosVec;
				// string tmp2searchReadAnchorSeq = readSeq_rcm.substr(0,
				// 	candidateSJendLocInReadVec_AC[tmp]) + "CT";
				string tmp2searchReadAnchorSeq = readSeq_rcm.substr(
					candidateSJendLocInReadVec_AC[tmp] - tmpToSearchSJanchorReadSeqLength, 
					tmpToSearchSJanchorReadSeqLength) + "CT";
				this->segmentSearch_longRange(tmpFoundSegMapPosVec,
					chrNameInt, searchStartPosInChr, searchEndPosInChr, 
					tmp2searchReadAnchorSeq, indexInfo, maxHit_anchorMap, 
					searchStartPosInChr);
				for(int tmp2 = 0; tmp2 < tmpFoundSegMapPosVec.size(); tmp2++)
				{
					foundSJendLocVec_negativeStrand.push_back(
						candidateSJendLocInReadVec_AC[tmp]);
					foundSegMapPosVec_negativeStrand.push_back(
						tmpFoundSegMapPosVec[tmp2]);
				}
			}
			this->selectBestShortAnchorSJ_completeInnerSoftClip_rcm_head(
				foundSJendLocVec_positiveStrand,
				foundSegMapPosVec_positiveStrand,
				foundSJendLocVec_negativeStrand,
				foundSegMapPosVec_negativeStrand);
		}
	}

	void segmentSearch_shortRange(vector<int>& candiSegMapPosVec, 
		int searchChrNameInt, int searchStartPosInChr, int searchEndPosInChr,
		const string& targetSegSeq, Index_Info* indexInfo, int maxHit, int basePos2add)
	{
		//cout << endl << "segmentSearch_shortRange starts ......" << endl;
		int targetSegSeqLength = targetSegSeq.length();
		//cout << "targetSegSeqLength: " << targetSegSeqLength << endl;
		int tmpChrLength = indexInfo->returnChromLength(searchChrNameInt);
		if((searchStartPosInChr < 0)||(searchEndPosInChr < 0)
			||(searchStartPosInChr + targetSegSeqLength - 1 > tmpChrLength)
			||(searchEndPosInChr + targetSegSeqLength -1 > tmpChrLength)
			||(searchStartPosInChr > searchEndPosInChr))
			return;
		//cout << "searchStartPosInChr: " << searchStartPosInChr << endl;
		//cout << "searchEndPosInChr: " << searchEndPosInChr + targetSegSeqLength - 1 << endl;
		string targetChromSeq = indexInfo->returnChromStrSubstr(
			searchChrNameInt, searchStartPosInChr, searchEndPosInChr + targetSegSeqLength - 1 - searchStartPosInChr + 1);
		//cout << "targetChromSeq: " << targetChromSeq << endl;
		int startLoc = 0;
		for(int tmp = 0; tmp < maxHit; tmp++)
		{
			int tabLoc = targetChromSeq.find(targetSegSeq, startLoc);
			if(tabLoc != string::npos)
			{
				candiSegMapPosVec.push_back(tabLoc + basePos2add);
				startLoc+1;
			}
			else
			{
				break;
			}
		}
	}

	// FIX Me: modify, should not extract targetChromSeq and reuse segmentSearch_shortRange
	int segmentSearch_longRange(vector<int>& candiSegMapPosVec, 
		int searchChrNameInt, int searchStartPosInChr, int searchEndPosInChr,
		const string& targetSegSeq, Index_Info* indexInfo, int maxHit, int basePos2add)
	{ 

		this->segmentSearch_shortRange(candiSegMapPosVec, 
			searchChrNameInt, searchStartPosInChr, searchEndPosInChr,
			targetSegSeq, indexInfo, maxHit, basePos2add);
	}	

	void fix_nwdp_extension_tail_for(Index_Info* indexInfo)
	{
		int backBuffer = 4;
		string unfixedTailSeq_for = readSeq_for.substr(
			mappedPartEndLoc_InRead_for - backBuffer, unfixedTailLen_for + backBuffer);
		string toProcessSeq = indexInfo->returnChromStrSubstr(chrNameInt, 
			mappedPartEndPos_inChr_for + 1 - backBuffer, unfixedTailLen_for + 3 + backBuffer);
		FixSingleAnchor_NWDP_Info* fixSingleAnchorNWDPinfo = new FixSingleAnchor_NWDP_Info();
		fixSingleAnchorNWDPinfo->doNWDP(unfixedTailSeq_for, toProcessSeq);
		tail_for_fixed_bool = fixSingleAnchorNWDPinfo->fixedOrNot();
		if(tail_for_fixed_bool)
		{
			vector<Jump_Code> tmpBackBufferMatchJumpCodeVec;
			Jump_Code tmpBackBufferMatchJumpCode(0-backBuffer, "M");
			tmpBackBufferMatchJumpCodeVec.push_back(tmpBackBufferMatchJumpCode);
			//cout << "tmpBackBufferMatchJumpCodeVec: " << this->jumpCodeVec2cigarString(tmpBackBufferMatchJumpCodeVec) << endl;
			vector<Jump_Code> tmpFixedTailJumpCodeVec_for;
			fixSingleAnchorNWDPinfo->copyJumpCodeVec2TargetVec(tmpFixedTailJumpCodeVec_for);
			//cout << "tmpFixedTailJumpCodeVec_for: " << this->jumpCodeVec2cigarString(tmpFixedTailJumpCodeVec_for) << endl;
			this->mergeJumpCodeVec(tmpBackBufferMatchJumpCodeVec, 0, 0,
				tmpFixedTailJumpCodeVec_for, 0, tmpFixedTailJumpCodeVec_for.size()-1,
				fixedTailJumpCodeVec_for);
			//cout << "fixedTailJumpCodeVec_for: " << this->jumpCodeVec2cigarString(fixedTailJumpCodeVec_for) << endl;
			fixSingleAnchorNWDPinfo->copyMismatchPosVec2TargetVec(fixedTailMismatchLocVec_inTailSeq_for);
			fixSingleAnchorNWDPinfo->copyMismatchCharVec2TargetVec(fixedTailMismatchCharVec_for);
		}
		delete fixSingleAnchorNWDPinfo;
	}

	void fix_nwdp_extension_head_rcm(Index_Info* indexInfo)
	{
		int forwardBuffer = 4;
		string unfixedHeadSeq_rcm = readSeq_rcm.substr(
			0, unfixedHeadLen_rcm + forwardBuffer);
		string toProcessSeq = indexInfo->returnChromStrSubstr(chrNameInt,
			mappedPartStartPos_inChr_rcm - unfixedHeadLen_rcm - 3,
			unfixedHeadLen_rcm + 3 + forwardBuffer);
		// nw_oneEndOpen_withMismatchJumpCode_reverse(unfixedHeadSeq_rcm,
		// 	toProcessSeq, fixedHeadJumpCodeVec_rcm,
		// 	fixedHeadMismatchLocVec_inHeadSeq_rcm, fixedHeadMismatchCharVec_rcm);
		FixSingleAnchor_NWDP_Info* fixSingleAnchorNWDPinfo = new FixSingleAnchor_NWDP_Info();
		fixSingleAnchorNWDPinfo->doNWDP_reverse(unfixedHeadSeq_rcm, toProcessSeq);
		head_rcm_fixed_bool = fixSingleAnchorNWDPinfo->fixedOrNot();
		if(head_rcm_fixed_bool)
		{
			vector<Jump_Code> tmpFixedHeadJumpCodeVec_rcm;
			fixSingleAnchorNWDPinfo->copyJumpCodeVec2TargetVec(tmpFixedHeadJumpCodeVec_rcm);
			
			vector<Jump_Code> tmpForwardBufferMatchJumpCodeVec;
			Jump_Code tmpForwardBufferMatchJumpCode(0-forwardBuffer, "M");
			tmpForwardBufferMatchJumpCodeVec.push_back(tmpForwardBufferMatchJumpCode);

			this->mergeJumpCodeVec(tmpFixedHeadJumpCodeVec_rcm, 0, tmpFixedHeadJumpCodeVec_rcm.size()-1,
				tmpForwardBufferMatchJumpCodeVec, 0, 0, fixedHeadJumpCodeVec_rcm);

			fixSingleAnchorNWDPinfo->copyMismatchPosVec2TargetVec(fixedHeadMismatchLocVec_inHeadSeq_rcm);
			fixSingleAnchorNWDPinfo->copyMismatchCharVec2TargetVec(fixedHeadMismatchCharVec_rcm);
		}
		delete fixSingleAnchorNWDPinfo;
	}

	void fix_tail_for(Index_Info* indexInfo, 
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		SJhash_Info* SJhashInfo,
		bool do_nwdp_bool, bool do_novelSJsearch_bool)
	{
		if(do_nwdp_bool)
			this->fix_nwdp_extension_tail_for(indexInfo);
		if(!tail_for_fixed_bool)
		{
			if(do_novelSJsearch_bool)
				this->fix_doNovelSJsearch_tail_for(indexInfo, 
					alignInferJunctionHashInfo, SJhashInfo);
		}
	
		if(tail_for_fixed_bool)
			this->generateFinalAlignInfo_afterFixing_forTail();
	}

	void fix_head_rcm(Index_Info* indexInfo, 
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo, 
		SJhash_Info* SJhashInfo,
		bool do_nwdp_bool, bool do_novelSJsearch_bool)
	{
		if(do_nwdp_bool)
			this->fix_nwdp_extension_head_rcm(indexInfo);
		if(!head_rcm_fixed_bool)
		{
			if(do_novelSJsearch_bool)
				this->fix_doNovelSJsearch_head_rcm(indexInfo, 
					alignInferJunctionHashInfo, SJhashInfo);
		}

		if(head_rcm_fixed_bool)
			this->generateFinalAlignInfo_afterFixing_rcmHead();
	}

};
#endif