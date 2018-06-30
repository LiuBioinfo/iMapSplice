// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef INCOMPLETEUNIQUEPAIRALIGNMENT2DETECTFUSION_INFO_H
#define INCOMPLETEUNIQUEPAIRALIGNMENT2DETECTFUSION_INFO_H

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"
#include "remapAgainstFusionBreakPoint_info.h"

using namespace std;

bool pairedMappedOrNot(int flag)
{
	if((flag & 0x1) && (flag & 0x2))
		return true;
	else
		return false;
}

class IncompleteUniquePairedAlignment2detectFusion_Info
{
private:
	int XM;

	bool Nor1Rcm2_or_Nor2Rcm1_bool;
	int chrNameInt;

	string readName_1;
	int startPos_1;
	vector<Jump_Code> cigarStringJumpCodeVec_1;
	string readSeq_1;
	int readLength_1;
	string qualSeq_1;
	int unfixedHeadLen_1;
	int unfixedTailLen_1;

	string readName_2;
	int startPos_2;
	vector<Jump_Code> cigarStringJumpCodeVec_2;
	string readSeq_2;
	int readLength_2;
	string qualSeq_2;
	int unfixedHeadLen_2;
	int unfixedTailLen_2;

public:	
	IncompleteUniquePairedAlignment2detectFusion_Info()
	{}

	int returnXM()
	{
		return XM;
	}

	string returnReadSeqSubStr_left(int startLocInRead, int endLocInRead)
	{
		return readSeq_1.substr(startLocInRead-1, endLocInRead - startLocInRead + 1);
	}

	string returnReadSeqSubStr_right(int startLocInRead, int endLocInRead)
	{
		return readSeq_2.substr(startLocInRead-1, endLocInRead - startLocInRead + 1);
	}	

	int returnUnfixedHeadLength_left()
	{
		if(cigarStringJumpCodeVec_1[0].type == "S")
			return cigarStringJumpCodeVec_1[0].len;
		else
			return 0;		
	}

	int returnUnfixedTailLength_right()
	{
		int jumpCodeVecSize_2 = cigarStringJumpCodeVec_2.size();
		if(cigarStringJumpCodeVec_1[jumpCodeVecSize_2 - 1].type == "S")
			return cigarStringJumpCodeVec_1[jumpCodeVecSize_2 - 1].len;
		else
			return 0;		
	}	

	int returnFirstMatchLength_left()
	{
		if(cigarStringJumpCodeVec_1[0].type == "M")
			return cigarStringJumpCodeVec_1[0].len;
		else
			return cigarStringJumpCodeVec_1[1].len;
	}

	int returnLastMatchLength_right()
	{
		int cigarStringJumpCodeVec_2_size = cigarStringJumpCodeVec_2.size();
		if(cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2_size - 1].type == "M")
			return cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2_size - 1].len;
		else
			return cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2_size - 2].len;
	}

	void returnFusionSpanningReadMapRange_canonical(int fusionCase, 
		int otherEndFusionAlignment_startPos,
		int otherEndFusionAlignment_endPos,
		int& mapRange_most_1, int& mapRange_most_2)
	{
		if(fusionCase == 1)
		{
			mapRange_most_1 = startPos_1;
			mapRange_most_2 = otherEndFusionAlignment_endPos;
		}
		else if(fusionCase == 2)
		{
			mapRange_most_1 = otherEndFusionAlignment_startPos;
			mapRange_most_2 = this->returnEndPos_2();
		}
		else if(fusionCase == 4)
		{
			mapRange_most_1 = otherEndFusionAlignment_endPos;
			mapRange_most_2 = startPos_1;
		}
		else if(fusionCase == 5)
		{
			mapRange_most_1 = this->returnEndPos_2();
			mapRange_most_2 = otherEndFusionAlignment_startPos;
		}
		else if(fusionCase == 7)
		{
			mapRange_most_1 = otherEndFusionAlignment_startPos;
			mapRange_most_2 = startPos_1;
		}
		else if(fusionCase == 8)
		{
			mapRange_most_1 = startPos_1;
			mapRange_most_2 = otherEndFusionAlignment_startPos;
		}
		else if(fusionCase == 10)
		{
			mapRange_most_1 = otherEndFusionAlignment_endPos;
			mapRange_most_2 = this->returnEndPos_2();
		}
		else if(fusionCase == 11)
		{
			mapRange_most_1 = this->returnEndPos_2();
			mapRange_most_2 = otherEndFusionAlignment_endPos;
		}
		else
		{
			cout << "error in returnFusionSpanningReadMapRange_canonical ..." << endl;
			exit(1);
		}									
	}

	void supportNumIncrement(
		string& tmpChrNameStr_gene1, string& tmpChrNameStr_gene2,
		int tmpBreakPointPos_gene1, int tmpBreakPointPos_gene2,
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		int tmpChrNameInt_gene1 = indexInfo->convertStringToInt(tmpChrNameStr_gene1);
		int tmpChrNameInt_gene2 = indexInfo->convertStringToInt(tmpChrNameStr_gene2);
		this->supportNumIncrement(tmpChrNameInt_gene1, tmpChrNameInt_gene2,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2,
				tmpFusionBreakPointHashInfo);
	}

	void supportNumIncrement(
		int tmpChrNameInt_gene1, int tmpChrNameInt_gene2,
		int tmpBreakPointPos_gene1, int tmpBreakPointPos_gene2,
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo)
	{
		//cout << "function supportNumIncrement starts ......" << endl;
		int tmpIndexInFusionBreakPointInfoVec 
			= tmpFusionBreakPointHashInfo->searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_gene1, tmpChrNameInt_gene2,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2);
		//cout << "tmpIndexInFusionBreakPointInfoVec: " << tmpIndexInFusionBreakPointInfoVec << endl;
		if(tmpIndexInFusionBreakPointInfoVec < 0)
		{
			cout << "error ! this fusion junc is not in the fusionBreakPointHashInfo" << endl;
			exit(1);
			return;
		}
		else
		{
			tmpFusionBreakPointHashInfo->supportNumIncrementWithIndex(
				tmpIndexInFusionBreakPointInfoVec);
		}
	}

	void supportNumIncrement_encompassing(
		string& tmpChrNameStr_gene1, string& tmpChrNameStr_gene2,
		int tmpBreakPointPos_gene1, int tmpBreakPointPos_gene2,
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		int tmpChrNameInt_gene1 = indexInfo->convertStringToInt(tmpChrNameStr_gene1);
		int tmpChrNameInt_gene2 = indexInfo->convertStringToInt(tmpChrNameStr_gene2);
		this->supportNumIncrement_encompassing(tmpChrNameInt_gene1, tmpChrNameInt_gene2,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2,
				tmpFusionBreakPointHashInfo);
	}	

	void supportNumIncrement_encompassing(
		int tmpChrNameInt_gene1, int tmpChrNameInt_gene2,
		int tmpBreakPointPos_gene1, int tmpBreakPointPos_gene2,
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo)
	{
		int tmpIndexInFusionBreakPointInfoVec 
			= tmpFusionBreakPointHashInfo->searchAndReturnIndexInBreakPointInfoVec_for(
				tmpChrNameInt_gene1, tmpChrNameInt_gene2,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2);
		if(tmpIndexInFusionBreakPointInfoVec < 0)
			return;
		else
		{
			tmpFusionBreakPointHashInfo->supportNumIncrementWithIndex_encompassing(
				tmpIndexInFusionBreakPointInfoVec);
		}
	}

	int returnJumpCodeLength_1(int index)
	{
		return cigarStringJumpCodeVec_1[index].len;
	}

	int returnJumpCodeLength_rev_1(int index)
	{
		return cigarStringJumpCodeVec_1[cigarStringJumpCodeVec_1.size()-1-index].len;
	}

	int returnJumpCodeLength_2(int index)
	{
		return cigarStringJumpCodeVec_2[index].len;
	}

	int returnJumpCodeLength_rev_2(int index)
	{
		return cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2.size()-1-index].len;
	}	

	int returnReadLength_1()
	{
		return readLength_1;
	}

	int returnReadLength_2()
	{
		return readLength_2;
	}

	int returnUnfixedHeadLength_1()
	{
		return unfixedHeadLen_1;
	}

	int returnUnfixedHeadLength_2()
	{
		return unfixedHeadLen_2;
	}

	int returnUnfixedTailLength_1()
	{
		return unfixedTailLen_1;
	}

	int returnUnfixedTailLength_2()
	{
		return unfixedTailLen_2;
	}	

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnCigarStringJumpCodeVecSize_1()
	{
		return cigarStringJumpCodeVec_1.size();
	}

	int returnCigarStringJumpCodeVecSize_2()
	{
		return cigarStringJumpCodeVec_2.size();
	}

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
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

	string returnChrName(Index_Info* indexInfo)
	{
		string tmpChrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		return tmpChrNameStr;
	}

	int returnStartPos_1()
	{
		return startPos_1;
	}

	int returnEndPos_2()
	{
		int endPos_2 = this->getEndPosOfSpecificJumpCode(startPos_2, cigarStringJumpCodeVec_2,
			cigarStringJumpCodeVec_2.size()-1);
		return endPos_2;
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

	void remapUnfixedHeadAgainstFusionBreakPoint_leftRead(
		RemapAgainstFusionBreakPoint_Info* remapInfo_leftReadHead_oriGeneAs1stGene,
		RemapAgainstFusionBreakPoint_Info* remapInfo_leftReadHead_oriGeneAs2ndGene,
		FusionBreakPointHash_Info* fusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		this->remapUnfixedHeadAgainstFusionBreakPoint_leftRead_oriGeneAs1stGene(
			remapInfo_leftReadHead_oriGeneAs1stGene, fusionBreakPointHashInfo, indexInfo);
		this->remapUnfixedHeadAgainstFusionBreakPoint_leftRead_oriGeneAs2ndGene(
			remapInfo_leftReadHead_oriGeneAs2ndGene, fusionBreakPointHashInfo, indexInfo);
	}
	
	void remapUnfixedHeadAgainstFusionBreakPoint_leftRead_oriGeneAs1stGene( // case 5, case 11
		RemapAgainstFusionBreakPoint_Info* remapInfo,
		FusionBreakPointHash_Info* fusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		//cout << "function remapUnfixedHeadAgainstFusionBreakPoint_leftRead_oriGeneAs1stGene( starts ......" << endl;
		// this gene -- gene1
		int leftMostPos = startPos_1;
		int rightMostPos = startPos_1 + cigarStringJumpCodeVec_1[1].len - 1;
		//cout << "leftMostPos: " << leftMostPos << endl;
		//cout << "rightMostPos: " << rightMostPos << endl;
		vector<int> candiBreakPointPosVec_thisGene;
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameInt, leftMostPos, rightMostPos, true, candiBreakPointPosVec_thisGene);

		vector<int> candiFusion_startLocInReadVec_oriGene;
		vector<int> candiFusion_breakPointPosVec_oriGene;
		vector<string> candiFusion_strandVec_oriGene;
		vector<int> candiFusion_chrNameIntVec_theOtherGene;
		vector<int> candiFusion_startPosInChrVec_theOtherGene;
		vector<int> candiFusion_breakPointPosVec_theOtherGene;
		vector<string> candiFusion_strandVec_theOtherGene;
		vector< vector<Jump_Code> > candiFusion_jumpCodeVecVec_theOtherGene;

		// generate candiFusion info for leftReadHead
		for(int tmp = 0; tmp < candiBreakPointPosVec_thisGene.size(); tmp++)
		{	
			int tmpCandiThisGeneBreakPointPos = candiBreakPointPosVec_thisGene[tmp];
			int tmpOffset = tmpCandiThisGeneBreakPointPos - startPos_1;
			int tmpUnfixedHeadLen = tmpOffset + unfixedHeadLen_1;
			int tmpMismatchNumMax = tmpUnfixedHeadLen / MATCH_BASE_PER_MISMATCH_BASE;
			string tmpReadHeadSeq = readSeq_1.substr(0, tmpUnfixedHeadLen);

			vector<int> tmpChrNameIntVec_theOtherGene;
			vector<int> tmpBreakPosVec_theOtherGene;
			vector<string> tmpStrandVec_thisGene;
			vector<string> tmpStrandVec_theOtherGene;
			fusionBreakPointHashInfo->returnTheOtherFusionGeneBreakPoint(
				chrNameInt, tmpCandiThisGeneBreakPointPos, true,
				tmpStrandVec_thisGene, tmpChrNameIntVec_theOtherGene,
				tmpBreakPosVec_theOtherGene, tmpStrandVec_theOtherGene);
			for(int tmp2 = 0; tmp2 < tmpChrNameIntVec_theOtherGene.size(); tmp2++)
			{
				int tmpChrNameInt_theOtherGene = tmpChrNameIntVec_theOtherGene[tmp2];
				int tmpBreakPos_theOtherGene = tmpBreakPosVec_theOtherGene[tmp2];
				string tmpStrand_thisGene = tmpStrandVec_thisGene[tmp2];
				string tmpStrand_theOtherGene = tmpStrandVec_theOtherGene[tmp2];
				string tmpStrand_fusion = tmpStrand_thisGene + tmpStrand_theOtherGene;

				if(tmpStrand_fusion == "--") // case 5
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene - tmpUnfixedHeadLen + 1;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedHeadLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(tmpReadHeadSeq, tmpChrSeq,
						tmpMismatchNumMax, 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedHeadLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_startLocInReadVec_oriGene.push_back(tmpUnfixedHeadLen + 1);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("-");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("-");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);						
					}
					else
					{}
					delete fixMatchInfo;
				}
				else if(tmpStrand_fusion == "-+")// tmpStrand_theOtherGene == "+" // case 11
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedHeadLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(
						convertStringToReverseComplement(tmpReadHeadSeq), 
						tmpChrSeq, tmpMismatchNumMax, 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedHeadLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_startLocInReadVec_oriGene.push_back(tmpUnfixedHeadLen + 1);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("-");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("+");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);							
					}
					else
					{}
					delete fixMatchInfo;
				}
				else
				{}
			}
		}

		//cout << "candiFusion_startLocInReadVec_oriGene.size(): " 
		//	<< candiFusion_startLocInReadVec_oriGene.size() << endl;
		// decide final candiFusion for leftReadHead
		remapInfo->setUnfixedHeadOrTail2remapBool(true);
		remapInfo->setOriGene1or2Bool(true);
		remapInfo->generateCandiFusionResultsVec(
			candiFusion_startLocInReadVec_oriGene,
			candiFusion_breakPointPosVec_oriGene,
			candiFusion_strandVec_oriGene,
			candiFusion_chrNameIntVec_theOtherGene,
			candiFusion_startPosInChrVec_theOtherGene,
			candiFusion_breakPointPosVec_theOtherGene,
			candiFusion_strandVec_theOtherGene,
			candiFusion_jumpCodeVecVec_theOtherGene);
	}

	void remapUnfixedHeadAgainstFusionBreakPoint_leftRead_oriGeneAs2ndGene( // case 2, case 10
		RemapAgainstFusionBreakPoint_Info* remapInfo,
		FusionBreakPointHash_Info* fusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		//cout << "function remapUnfixedHeadAgainstFusionBreakPoint_leftRead_oriGeneAs2ndGene( starts ......" << endl;
		// this gene -- gene2
		int leftMostPos = startPos_1;
		int rightMostPos = startPos_1 + cigarStringJumpCodeVec_1[1].len - 1;	
		//cout << "leftMostPos: " << leftMostPos << endl;
		//cout << "rightMostPos: " << rightMostPos << endl;		
		vector<int> candiBreakPointPosVec_thisGene;
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameInt, leftMostPos, rightMostPos, false, candiBreakPointPosVec_thisGene);
		//cout << "candiBreakPointPosVec_thisGene.size(): " << candiBreakPointPosVec_thisGene.size() << endl;
		vector<int> candiFusion_startLocInReadVec_oriGene;
		vector<int> candiFusion_breakPointPosVec_oriGene;
		vector<string> candiFusion_strandVec_oriGene;
		vector<int> candiFusion_chrNameIntVec_theOtherGene;
		vector<int> candiFusion_startPosInChrVec_theOtherGene;
		vector<int> candiFusion_breakPointPosVec_theOtherGene;
		vector<string> candiFusion_strandVec_theOtherGene;
		vector< vector<Jump_Code> > candiFusion_jumpCodeVecVec_theOtherGene;

		// generate candiFusion info for leftReadHead
		for(int tmp = 0; tmp < candiBreakPointPosVec_thisGene.size(); tmp++)
		{	
			int tmpCandiThisGeneBreakPointPos = candiBreakPointPosVec_thisGene[tmp];
			int tmpOffset = tmpCandiThisGeneBreakPointPos - startPos_1;
			int tmpUnfixedHeadLen = tmpOffset + unfixedHeadLen_1;
			int tmpMismatchNumMax = tmpUnfixedHeadLen / MATCH_BASE_PER_MISMATCH_BASE;
			string tmpReadHeadSeq = readSeq_1.substr(0, tmpUnfixedHeadLen);
			// cout << "tmpCandiThisGeneBreakPointPos: " << tmpCandiThisGeneBreakPointPos << endl;
			// cout << "tmpOffset: " << tmpOffset << endl;
			// cout << "tmpUnfixedHeadLen: " << tmpUnfixedHeadLen << endl;
			// cout << "tmpReadHeadSeq: " << tmpReadHeadSeq << endl;
			vector<int> tmpChrNameIntVec_theOtherGene;
			vector<int> tmpBreakPosVec_theOtherGene;
			vector<string> tmpStrandVec_thisGene;
			vector<string> tmpStrandVec_theOtherGene;
			fusionBreakPointHashInfo->returnTheOtherFusionGeneBreakPoint(
				chrNameInt, tmpCandiThisGeneBreakPointPos, false,
				tmpStrandVec_thisGene, tmpChrNameIntVec_theOtherGene,
				tmpBreakPosVec_theOtherGene, tmpStrandVec_theOtherGene);
			//cout << "tmpBreakPosVec_theOtherGene.size(): " << tmpBreakPosVec_theOtherGene.size() << endl; 
			for(int tmp2 = 0; tmp2 < tmpChrNameIntVec_theOtherGene.size(); tmp2++)
			{
				int tmpChrNameInt_theOtherGene = tmpChrNameIntVec_theOtherGene[tmp2];
				int tmpBreakPos_theOtherGene = tmpBreakPosVec_theOtherGene[tmp2];
				string tmpStrand_thisGene = tmpStrandVec_thisGene[tmp2];
				string tmpStrand_theOtherGene = tmpStrandVec_theOtherGene[tmp2];
				string tmpStrand_fusion = tmpStrand_theOtherGene + tmpStrand_thisGene;
				//cout << "tmpBreakPos_theOtherGene: " << tmpBreakPos_theOtherGene << endl;
				//cout << "tmpStrand_fusion: " << tmpStrand_fusion << endl;
				if(tmpStrand_fusion == "++") // case 2
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene - tmpUnfixedHeadLen + 1;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedHeadLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(tmpReadHeadSeq, tmpChrSeq,
						tmpMismatchNumMax, 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedHeadLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_startLocInReadVec_oriGene.push_back(tmpUnfixedHeadLen + 1);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("+");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("+");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);						
					}
					else
					{}
					delete fixMatchInfo;
				}
				else if(tmpStrand_fusion == "-+") // tmpStrand_theOtherGene == "-", case 10
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedHeadLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(
						convertStringToReverseComplement(tmpReadHeadSeq), 
						tmpChrSeq, tmpMismatchNumMax, 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedHeadLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_startLocInReadVec_oriGene.push_back(tmpUnfixedHeadLen + 1);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("+");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("-");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);							
					}
					else
					{}
					delete fixMatchInfo;
				}
				else
				{}
			}
		}
		//cout << "candiFusion_startLocInReadVec_oriGene.size(): " 
		//	<< candiFusion_startLocInReadVec_oriGene.size() << endl;
		// decide final candiFusion for leftReadHead
		remapInfo->setUnfixedHeadOrTail2remapBool(true);
		remapInfo->setOriGene1or2Bool(false);		
		remapInfo->generateCandiFusionResultsVec(
			candiFusion_startLocInReadVec_oriGene,
			candiFusion_breakPointPosVec_oriGene,
			candiFusion_strandVec_oriGene,
			candiFusion_chrNameIntVec_theOtherGene,
			candiFusion_startPosInChrVec_theOtherGene,
			candiFusion_breakPointPosVec_theOtherGene,
			candiFusion_strandVec_theOtherGene,
			candiFusion_jumpCodeVecVec_theOtherGene);
	}

	void remapUnfixedTailAgainstFusionBreakPoint_rightRead(
		RemapAgainstFusionBreakPoint_Info* remapInfo_rightReadTail_oriGeneAs1stGene,
		RemapAgainstFusionBreakPoint_Info* remapInfo_rightReadTail_oriGeneAs2ndGene,		
		FusionBreakPointHash_Info* fusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		this->remapUnfixedTailAgainstFusionBreakPoint_rightRead_oriGeneAs1stGene(
			remapInfo_rightReadTail_oriGeneAs1stGene,
			fusionBreakPointHashInfo, indexInfo);
		this->remapUnfixedTailAgainstFusionBreakPoint_rightRead_oriGeneAs2ndGene(
			remapInfo_rightReadTail_oriGeneAs2ndGene,
			fusionBreakPointHashInfo, indexInfo);
	}

	void remapUnfixedTailAgainstFusionBreakPoint_rightRead_oriGeneAs1stGene(
		RemapAgainstFusionBreakPoint_Info* remapInfo,
		FusionBreakPointHash_Info* fusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		int endPos_2 = this->returnEndPos_2();
		int leftMostPos = endPos_2 
			- cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2.size()-2].len + 1;
		int rightMostPos = endPos_2;

		vector<int> candiBreakPointPosVec_thisGene;
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameInt, leftMostPos, rightMostPos, true, candiBreakPointPosVec_thisGene);

		vector<int> candiFusion_endLocInReadVec_oriGene;
		vector<int> candiFusion_breakPointPosVec_oriGene;
		vector<string> candiFusion_strandVec_oriGene;
		vector<int> candiFusion_chrNameIntVec_theOtherGene;
		vector<int> candiFusion_startPosInChrVec_theOtherGene;
		vector<int> candiFusion_breakPointPosVec_theOtherGene;
		vector<string> candiFusion_strandVec_theOtherGene;
		vector< vector<Jump_Code> > candiFusion_jumpCodeVecVec_theOtherGene;

		// generate candiFusion info for rightReadTail
		for(int tmp = 0; tmp < candiBreakPointPosVec_thisGene.size(); tmp++)
		{
			int tmpCandiThisGeneBreakPointPos = candiBreakPointPosVec_thisGene[tmp];
			int tmpOffset = endPos_2 - tmpCandiThisGeneBreakPointPos;
			int tmpUnfixedTailLen = tmpOffset + unfixedTailLen_2;
			int tmpMismatchNumMax = tmpUnfixedTailLen / MATCH_BASE_PER_MISMATCH_BASE;
			string tmpReadTailSeq = readSeq_2.substr(readLength_2 - tmpUnfixedTailLen, 
				tmpUnfixedTailLen);

			vector<int> tmpChrNameIntVec_theOtherGene;
			vector<int> tmpBreakPosVec_theOtherGene;
			vector<string> tmpStrandVec_thisGene;
			vector<string> tmpStrandVec_theOtherGene;
			fusionBreakPointHashInfo->returnTheOtherFusionGeneBreakPoint(
				chrNameInt, tmpCandiThisGeneBreakPointPos, true,
				tmpStrandVec_thisGene, tmpChrNameIntVec_theOtherGene,
				tmpBreakPosVec_theOtherGene, tmpStrandVec_theOtherGene);
			for(int tmp2 = 0; tmp2 < tmpChrNameIntVec_theOtherGene.size(); tmp2++)
			{
				int tmpChrNameInt_theOtherGene = tmpChrNameIntVec_theOtherGene[tmp2];
				int tmpBreakPos_theOtherGene = tmpBreakPosVec_theOtherGene[tmp2];
				string tmpStrand_thisGene = tmpStrandVec_thisGene[tmp2];
				string tmpStrand_theOtherGene = tmpStrandVec_theOtherGene[tmp2];
				string tmpStrand_fusion = tmpStrand_thisGene + tmpStrand_theOtherGene;

				if(tmpStrand_fusion == "++") // case 1
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedTailLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(tmpReadTailSeq, tmpChrSeq,
						tmpMismatchNumMax, readLength_2 - tmpUnfixedTailLen + 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedTailLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_endLocInReadVec_oriGene.push_back(readLength_2 - tmpUnfixedTailLen);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("+");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("+");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);
					}
					else
					{}
					delete fixMatchInfo;
				}
				else if(tmpStrand_fusion == "+-") // case 8
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene - tmpUnfixedTailLen + 1;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedTailLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(convertStringToReverseComplement(tmpReadTailSeq), 
						tmpChrSeq, tmpMismatchNumMax, readLength_2 - tmpUnfixedTailLen + 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedTailLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_endLocInReadVec_oriGene.push_back(readLength_2 - tmpUnfixedTailLen);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("+");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("-");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);
					}
					else
					{}
					delete fixMatchInfo;
				}
				else
				{}
			}
		}

		// decide final candiFusion for leftReadHead
		remapInfo->setUnfixedHeadOrTail2remapBool(false);
		remapInfo->setOriGene1or2Bool(true);			
		remapInfo->generateCandiFusionResultsVec(
			candiFusion_endLocInReadVec_oriGene,
			candiFusion_breakPointPosVec_oriGene,
			candiFusion_strandVec_oriGene,
			candiFusion_chrNameIntVec_theOtherGene,
			candiFusion_startPosInChrVec_theOtherGene,
			candiFusion_breakPointPosVec_theOtherGene,
			candiFusion_strandVec_theOtherGene,
			candiFusion_jumpCodeVecVec_theOtherGene);
	}

	void remapUnfixedTailAgainstFusionBreakPoint_rightRead_oriGeneAs2ndGene(
		RemapAgainstFusionBreakPoint_Info* remapInfo,
		FusionBreakPointHash_Info* fusionBreakPointHashInfo, Index_Info* indexInfo)
	{
		int endPos_2 = this->returnEndPos_2();
		int leftMostPos = endPos_2 
			- cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2.size()-2].len + 1;
		int rightMostPos = endPos_2;

		vector<int> candiBreakPointPosVec_thisGene;
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameInt, leftMostPos, rightMostPos, false, candiBreakPointPosVec_thisGene);

		vector<int> candiFusion_endLocInReadVec_oriGene;
		vector<int> candiFusion_breakPointPosVec_oriGene;
		vector<string> candiFusion_strandVec_oriGene;
		vector<int> candiFusion_chrNameIntVec_theOtherGene;
		vector<int> candiFusion_startPosInChrVec_theOtherGene;
		vector<int> candiFusion_breakPointPosVec_theOtherGene;
		vector<string> candiFusion_strandVec_theOtherGene;
		vector< vector<Jump_Code> > candiFusion_jumpCodeVecVec_theOtherGene;

		// generate candiFusion info for rightReadTail
		for(int tmp = 0; tmp < candiBreakPointPosVec_thisGene.size(); tmp++)
		{
			int tmpCandiThisGeneBreakPointPos = candiBreakPointPosVec_thisGene[tmp];
			int tmpOffset = endPos_2 - tmpCandiThisGeneBreakPointPos;
			int tmpUnfixedTailLen = tmpOffset + unfixedTailLen_2;
			int tmpMismatchNumMax = tmpUnfixedTailLen / MATCH_BASE_PER_MISMATCH_BASE;
			string tmpReadTailSeq = readSeq_2.substr(readLength_2 - tmpUnfixedTailLen, 
				tmpUnfixedTailLen);			
		
			vector<int> tmpChrNameIntVec_theOtherGene;
			vector<int> tmpBreakPosVec_theOtherGene;
			vector<string> tmpStrandVec_thisGene;
			vector<string> tmpStrandVec_theOtherGene;
			fusionBreakPointHashInfo->returnTheOtherFusionGeneBreakPoint(
				chrNameInt, tmpCandiThisGeneBreakPointPos, false,
				tmpStrandVec_thisGene, tmpChrNameIntVec_theOtherGene,
				tmpBreakPosVec_theOtherGene, tmpStrandVec_theOtherGene);
			for(int tmp2 = 0; tmp2 < tmpChrNameIntVec_theOtherGene.size(); tmp2++)
			{
				int tmpChrNameInt_theOtherGene = tmpChrNameIntVec_theOtherGene[tmp2];
				int tmpBreakPos_theOtherGene = tmpBreakPosVec_theOtherGene[tmp2];
				string tmpStrand_thisGene = tmpStrandVec_thisGene[tmp2];
				string tmpStrand_theOtherGene = tmpStrandVec_theOtherGene[tmp2];
				string tmpStrand_fusion = tmpStrand_theOtherGene + tmpStrand_thisGene;
				
				if(tmpStrand_fusion == "--") // case 4
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedTailLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(tmpReadTailSeq, tmpChrSeq,
						tmpMismatchNumMax, readLength_2 - tmpUnfixedTailLen + 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedTailLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_endLocInReadVec_oriGene.push_back(readLength_2 - tmpUnfixedTailLen);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("-");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("-");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);						
					}
					else
					{}
					delete fixMatchInfo;					
				}
				else if(tmpStrand_fusion == "+-") // case 7
				{
					int tmpChrSeq_startPos = tmpBreakPos_theOtherGene - tmpUnfixedTailLen + 1;
					string tmpChrSeq = indexInfo->returnChromStrSubstr(
						tmpChrNameInt_theOtherGene, tmpChrSeq_startPos, tmpUnfixedTailLen);
					FixDoubleAnchor_Match_Info* fixMatchInfo = new FixDoubleAnchor_Match_Info();
					bool matchBool = fixMatchInfo->fixMatch(convertStringToReverseComplement(tmpReadTailSeq), 
						tmpChrSeq, tmpMismatchNumMax, readLength_2 - tmpUnfixedTailLen + 1);
					if(matchBool)
					{
						vector<Jump_Code> jumpCodeVec_theOtherGene;
						Jump_Code theOtherGeneMatchJumpCode(tmpUnfixedTailLen, "M");
						jumpCodeVec_theOtherGene.push_back(theOtherGeneMatchJumpCode);

						candiFusion_endLocInReadVec_oriGene.push_back(readLength_2 - tmpUnfixedTailLen);
						candiFusion_breakPointPosVec_oriGene.push_back(tmpCandiThisGeneBreakPointPos);
						candiFusion_strandVec_oriGene.push_back("-");
						candiFusion_chrNameIntVec_theOtherGene.push_back(tmpChrNameInt_theOtherGene);
						candiFusion_startPosInChrVec_theOtherGene.push_back(tmpChrSeq_startPos);
						candiFusion_breakPointPosVec_theOtherGene.push_back(tmpBreakPos_theOtherGene);
						candiFusion_strandVec_theOtherGene.push_back("+");
						candiFusion_jumpCodeVecVec_theOtherGene.push_back(jumpCodeVec_theOtherGene);						
					}
					else
					{}
					delete fixMatchInfo;
				}
				else
				{}
			}
		}

		// decide final candiFusion for leftReadHead
		remapInfo->setUnfixedHeadOrTail2remapBool(false);
		remapInfo->setOriGene1or2Bool(false);
		remapInfo->generateCandiFusionResultsVec(
			candiFusion_endLocInReadVec_oriGene,
			candiFusion_breakPointPosVec_oriGene,
			candiFusion_strandVec_oriGene,
			candiFusion_chrNameIntVec_theOtherGene,
			candiFusion_startPosInChrVec_theOtherGene,
			candiFusion_breakPointPosVec_theOtherGene,
			candiFusion_strandVec_theOtherGene,
			candiFusion_jumpCodeVecVec_theOtherGene);
	}


	bool initiateWith2samStr_globalMapOuterUnfixedEnd2detectFusionBreakPoint(
		const string& samStr_1, const string& samStr_2,
		Index_Info* indexInfo)
	{
		//cout << "start to initiateWith2samStr_outerSoftClipUniquePairedAlignmentOrNot" << endl;
		//cout << "samStr_1: " << samStr_1 << endl;
		//cout << "samStr_2: " << samStr_2 << endl;
		vector<string> samFieldVec_1;
		vector<string> samFieldVec_2;
		int startLoc = 0;
		for(int tmp = 0; tmp < 13; tmp++)
		{
			int tabLoc = samStr_1.find("\t", startLoc);
			string tmpSamField = samStr_1.substr(startLoc, tabLoc-startLoc);
			samFieldVec_1.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_1.push_back(samStr_1.substr(startLoc));	
		startLoc = 0;
		for(int tmp = 0; tmp < 13; tmp++)
		{
			int tabLoc = samStr_2.find("\t", startLoc);
			string tmpSamField = samStr_2.substr(startLoc, tabLoc-startLoc);
			samFieldVec_2.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec_2.push_back(samStr_2.substr(startLoc));	

		string chrNameStr_1 = samFieldVec_1[2];
		string chrNameStr_2 = samFieldVec_2[2];
		//cout << "chrNameStr_1: " << chrNameStr_1 << endl;
		//cout << "chrNameStr_2: " << chrNameStr_2 << endl;

		if((chrNameStr_1 == "*")||(chrNameStr_2 == "*"))
			return false;
		string IHfield_1 = samFieldVec_1[12];
		string IHfield_2 = samFieldVec_2[12];
		string IHintStr_1 = IHfield_1.substr(5);
		string IHintStr_2 = IHfield_2.substr(5);
		//cout << "IHfield_1: " << IHfield_1 << endl;
		//cout << "IHfield_2: " << IHfield_2 << endl;
		int IHint_1 = atoi(IHintStr_1.c_str());
		int IHint_2 = atoi(IHintStr_2.c_str());
		//cout << "IHint_1: " << IHint_1 << endl;
		//cout << "IHint_2: " << IHint_2 << endl;
		if(!((IHint_1 == 1)&&(IHint_2 == 1)))
			return false;

		// string XMfield = samFieldVec_1[14];
		// string XMintStr = XMfield.substr(5);
		// XM = atoi(XMintStr.c_str());
		int XM_loc = samStr_1.find("XM:i:");
		if(XM_loc == string::npos) // XM not found
			XM = 0;
		else // XM found
		{
			int XM_nextTab_loc = samStr_1.find("\t", XM_loc+1);
			int XMintStartLoc = XM_loc + 5;
			int XMintEndLoc;
			if(XM_nextTab_loc == string::npos) // XM found, but it is the last field as no nextTab is found 
				XMintEndLoc = samStr_1.length() - 1;
			else // XM found, but it is not the last field as the nextTab is found
				XMintEndLoc = XM_nextTab_loc - 1;
			string XMfieldStr = samStr_1.substr(XMintStartLoc, XMintStartLoc - XMintEndLoc + 1);
			XM = atoi(XMfieldStr.c_str());
		}
		string flagStr_1 = samFieldVec_1[1];
		string flagStr_2 = samFieldVec_2[1];
		int flagInt_1 = atoi(flagStr_1.c_str());
		int flagInt_2 = atoi(flagStr_2.c_str());
		//cout << "flagInt_1: " << flagInt_1 << endl;
		//cout << "flagInt_2: " << flagInt_2 << endl;
		bool pairedMappedOrNot_1 = pairedMappedOrNot(flagInt_1);
		bool pairedMappedOrNot_2 = pairedMappedOrNot(flagInt_2);
		//cout << "pairedMappedOrNot_1: " << pairedMappedOrNot_1 << endl;
		//cout << "pairedMappedOrNot_2: " << pairedMappedOrNot_2 << endl;
		if(!(pairedMappedOrNot_1 && pairedMappedOrNot_2))
			return false;

		// start to initiate all the private elements
		if((flagInt_1 & 0x40)&&(flagInt_2 & 0x80)
			&&(flagInt_1 & 0x20)&&(flagInt_2 & 0x10))
			Nor1Rcm2_or_Nor2Rcm1_bool = true;
		else if((flagInt_1 & 0x80)&&(flagInt_2 & 0x40)
			&&(flagInt_1 & 0x20)&&(flagInt_2 & 0x10))	
			Nor1Rcm2_or_Nor2Rcm1_bool = false;
		else
			return false;
		//cout << "Nor1Rcm2_or_Nor2Rcm1_bool: " << Nor1Rcm2_or_Nor2Rcm1_bool << endl;
		chrNameInt = indexInfo->convertStringToInt(chrNameStr_1);
		//cout << "chrNameInt: " << chrNameInt << endl;
		readName_1 = samFieldVec_1[0];
		//cout << "readName_1: " << readName_1 << endl;
		string startPosStr_1 = samFieldVec_1[3];
		startPos_1 = atoi(startPosStr_1.c_str());
		//cout << "startPos_1: " << startPos_1 << endl;
		string cigarString_1 = samFieldVec_1[5];
		//cout << "cigarString_1: " << cigarString_1 << endl;
		//vector<Jump_Code> cigarStringJumpCodeVec_1;
		this->cigarString2jumpCodeVec(cigarString_1, cigarStringJumpCodeVec_1);
		//cout << "cigarStringJumpCodeVec_1.size(): " << cigarStringJumpCodeVec_1.size() << endl;
		readSeq_1 = samFieldVec_1[9];
		readLength_1 = readSeq_1.length();
		//cout << "readLength_1: " << readLength_1 << endl;
		qualSeq_1 = samFieldVec_1[10];
		if(cigarStringJumpCodeVec_1[0].type == "S")
			unfixedHeadLen_1 = cigarStringJumpCodeVec_1[0].len;
		else
			unfixedHeadLen_1 = 0;

		if(cigarStringJumpCodeVec_1[cigarStringJumpCodeVec_1.size()-1].type == "S")
			unfixedTailLen_1 = cigarStringJumpCodeVec_1[cigarStringJumpCodeVec_1.size()-1].len;
		else
			unfixedTailLen_1 = 0;


		readName_2 = samFieldVec_2[0];
		//cout << "readName_2: " << readName_2 << endl;
		string startPosStr_2 = samFieldVec_2[3];
		startPos_2 = atoi(startPosStr_2.c_str());
		//cout << "startPos_2: " << startPos_2 << endl;
		string cigarString_2 = samFieldVec_2[5];
		//cout << "cigarString_2: " << cigarString_2 << endl;
		//vector<Jump_Code> cigarStringJumpCodeVec_2;
		this->cigarString2jumpCodeVec(cigarString_2, cigarStringJumpCodeVec_2);
		//cout << "cigarStringJumpCodeVec_2.size(): " << cigarStringJumpCodeVec_2.size() << endl;
		readSeq_2 = samFieldVec_2[9];
		readLength_2 = readSeq_2.length();
		//cout << "readLength_2: " << readLength_2 << endl;
		qualSeq_2 = samFieldVec_2[10];
		if(cigarStringJumpCodeVec_2[0].type == "S")
			unfixedHeadLen_2 = cigarStringJumpCodeVec_2[0].len;
		else
			unfixedHeadLen_2 = 0;
		if(cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2.size()-1].type == "S")
			unfixedTailLen_2 = cigarStringJumpCodeVec_2[cigarStringJumpCodeVec_2.size()-1].len;
		else
			unfixedTailLen_2 = 0;

		if(((unfixedHeadLen_1 > 0)||(unfixedTailLen_2 > 0))
			&&(unfixedHeadLen_2 == 0)&&(unfixedTailLen_1 == 0))
			return true;
		else
			return false;
	}

	bool initiateWith2samStr_remapOuterUnfixedEndAgainstFusionBreakPoint(
		const string& samStr_1, const string& samStr_2,
		Index_Info* indexInfo)
	{
		return this->initiateWith2samStr_globalMapOuterUnfixedEnd2detectFusionBreakPoint(
			samStr_1, samStr_2, indexInfo);
	}


	bool leftReadHeadUnfixedOrNot()
	{
		return (unfixedHeadLen_1 > 0);
	}

	bool leftReadHeadUnfixed_bool()
	{
		return (unfixedHeadLen_1 > 0);
	}	

	bool rightReadHeadUnfixed_bool()
	{
		return (unfixedHeadLen_2 > 0);
	}

	bool leftReadTailUnfixed_bool()
	{
		return (unfixedTailLen_1 > 0);
	}

	bool rightReadTailUnfixed_bool()
	{
		return (unfixedTailLen_2 > 0);
	}

	bool rightReadHeadUnfixedOrNot()
	{
		return (unfixedHeadLen_2 > 0);
	}

	bool leftReadTailUnfixedOrNot()
	{
		return (unfixedTailLen_1 > 0);
	}

	bool rightReadTailUnfixedOrNot()
	{
		return (unfixedTailLen_2 > 0);
	}	

	string returnReadName_1()
	{
		return readName_1;
	}

	string returnReadName_2()
	{
		return readName_2;
	}

	string returnUnfixedLeftReadHead_readSeq()
	{
		return readSeq_1.substr(0, unfixedHeadLen_1);
	}

	string returnUnfixedLeftReadHead_qualSeq(bool fasta_or_fastq)
	{
		string default_qualSeq = "*";
		if(fasta_or_fastq)
			return default_qualSeq;
		else
			return qualSeq_1.substr(0, unfixedHeadLen_1);
	}

	string returnUnfixedRightReadHead_readSeq()
	{
		return readSeq_2.substr(0, unfixedHeadLen_2);
	}

	string returnUnfixedRightReadHead_qualSeq(bool fasta_or_fastq)
	{
		string default_qualSeq = "*";
		if(fasta_or_fastq)
			return default_qualSeq;
		else
			return qualSeq_2.substr(0, unfixedHeadLen_2);
	}

	string returnUnfixedRightReadTail_readSeq()
	{
		return readSeq_2.substr(readLength_2 - unfixedTailLen_2, unfixedTailLen_2);
	}

	string returnUnfixedRightReadTail_qualSeq(bool fasta_or_fastq)
	{
		string default_qualSeq = "*";
		if(fasta_or_fastq)
			return default_qualSeq;
		else
			return qualSeq_2.substr(readLength_2 - unfixedTailLen_2, unfixedTailLen_2);
	}
	
	/*
	string returnFixedFusionSamStr(
		bool leftReadHeadUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd, 
		bool rightReadTailUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs2nd,
		Index_Info* indexInfo)
	{
		string tmpFixedFusionSamStr = "";
		string tmpSAM_leftReadHead; 
		string tmpSAM_rightReadTail;
		string tmpSAM_oriGene;
		
		if(leftReadHeadUnfixed_fixed_bool)
			tmpSAM_leftReadHead = 
		
		if(rightReadTailUnfixed_fixed_bool)
			tmpSAM_rightReadTail = 
		
		if(leftReadHeadUnfixed_fixed_bool || rightReadTailUnfixed_fixed_bool)
			tmpSAM_oriGene = 
		else
			tmpSAM_oriGene = 

		if(leftReadHeadUnfixed_fixed_bool && rightReadTailUnfixed_fixed_bool)
			tmpSAM_oriGene = tmpSAM_leftReadHead + "\t"
				+ tmpSAM_oriGene + "\t" + tmpSAM_rightReadTail;
		else if(leftReadHeadUnfixed_fixed_bool)
			tmpSAM_oriGene = tmpSAM_leftReadHead + "\t" + tmpSAM_oriGene;
		else if(rightReadTailUnfixed_fixed_bool)
			tmpSAM_oriGene = tmpSAM_oriGene + "\t" + tmpSAM_rightReadTail
		else
			tmpFixedFusionSamStr = tmpSAM_oriGene;
		return tmpFixedFusionSamStr;
	}

	void generateFusionInfoFrom2remapInfo_leftHead(
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd,
		string& tmpChrNameStr_leftHead, int& tmpBreakPointPos_leftHead,
		string& tmpChrNameStr_oriGene, int& tmpBreakPointPos_oriGene,
		bool& forOrRevBool_leftHead, string& strand_leftHead,
		bool& forOrRevBool_oriGene, string& strand_oriGene,
		bool& oriGene_1_or_2_bool, string& strand_fusion,
		Index_Info* indexInfo)
	{
		tmpChrNameStr_oriGene = indexInfo->returnChrNameStr(chrNameInt);
		if(tmpRemapInfo_leftHead_oriGeneAs1st->returnResultsSize() == 1)
		{
			//oriGene = gene1, theOtherGene = gene2
			oriGene_1_or_2_bool = true;
			int tmpChrNameInt_leftHead 
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_chrNameInt_theOtherGene(0);
			tmpChrNameStr_leftHead = indexInfo->returnChrNameStr(tmpChrNameInt_leftHead);
			tmpBreakPointPos_leftHead 
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_breakPointPos_theOtherGene(0);
			tmpBreakPointPos_oriGene
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_breakPointPos_oriGene(0);
			strand_leftHead = tmpRemapInfo_leftHead_oriGeneAs1st->return_strand_theOtherGene(0);
			strand_oriGene = tmpRemapInfo_leftHead_oriGeneAs1st->return_strand_oriGene(0);
			strand_fusion = strand_oriGene + strand_leftHead;
			if(strand_theOtherGene == "+") // case 11, 
			{
				forOrRevBool_leftHead = false;
				forOrRevBool_oriGene = true;
			}
			else// case 5 (strand_oriGene == "-")
			{
				forOrRevBool_leftHead = true;
				forOrRevBool_oriGene = true;
			}
		}
		else
		{
			//oriGene = gene2, theOtherGene = gene1
			oriGene_1_or_2_bool = false;
			int tmpChrNameInt_leftHead 
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_chrNameInt_theOtherGene(0);
			tmpChrNameStr_leftHead = indexInfo->returnChrNameStr(tmpChrNameInt_leftHead);
			tmpBreakPointPos_leftHead 
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_breakPointPos_theOtherGene(0);
			tmpBreakPointPos_oriGene
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_breakPointPos_oriGene(0);
			strand_leftHead = tmpRemapInfo_leftHead_oriGeneAs2nd->return_strand_theOtherGene(0);
			strand_oriGene = tmpRemapInfo_leftHead_oriGeneAs2nd->return_strand_oriGene(0);			
			strand_fusion = strand_leftHead + strand_oriGene;
			if(strand_theOtherGene == "+") // case 2
			{
				forOrRevBool_leftHead = true;
				forOrRevBool_oriGene = true;				
			}
			else // case 10
			{
				forOrRevBool_leftHead = false;
				forOrRevBool_oriGene = true;				
			}
		}
	}

	void generateFusionRemapInfo_leftHead_toReportInAdjustedFusionBreakPointStr(
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd,
		string& tmpStrand_gene1, string& tmpStrand_gene2,
		string& tmpChrNameStr_gene1, string& tmpChrNameStr_gene2,
		int& tmpBreakPointPos_gene1, int& tmpBreakPointPos_gene2,
		string& tmpFusionJuncFlankString, int& tmpAnchor_left, int& tmpAnchor_right,
		Index_Info* indexInfo)
	{
		string tmpChrNameStr_leftHead;
		int tmpBreakPointPos_leftHead;
		string tmpChrNameStr_oriGene;
		int tmpBreakPointPos_oriGene;
		bool forOrRevBool_leftHead; 
		string strand_leftHead;
		bool forOrRevBool_oriGene; 
		string strand_oriGene;
		bool oriGene_1_or_2_bool;
		string strand_fusion;
		this->generateFusionInfoFrom2remapInfo_leftHead(
			tmpRemapInfo_leftHead_oriGeneAs1st, tmpRemapInfo_leftHead_oriGeneAs2nd,
			tmpChrNameStr_leftHead, tmpBreakPointPos_leftHead, tmpChrNameStr_oriGene, 
			tmpBreakPointPos_oriGene, forOrRevBool_leftHead, strand_leftHead, 
			forOrRevBool_oriGene, strand_oriGene,
			oriGene_1_or_2_bool, strand_fusion, indexInfo);
		if(oriGene_1_or_2_bool) // case 5, case 11, oriGene == Gene 1
		{
			tmpChrNameStr_gene1 = tmpChrNameStr_oriGene;
			tmpChrNameStr_gene2 = tmpChrNameStr_leftHead;
			tmpBreakPointPos_gene1 = tmpBreakPointPos_oriGene;
			tmpBreakPointPos_gene2 = tmpBreakPointPos_leftHead;			
			if(strand_fusion == "--") // case 5
			{
				tmpStrand_gene1 = "-";
				tmpStrand_gene2 = "-";
				tmpFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					tmpChrNameStr_leftHead, tmpBreakPointPos_leftHead,
					tmpChrNameStr_oriGene, tmpBreakPointPos_oriGene,
					true, true, true);
				tmpAnchor_left 
					= tmpRemapInfo_leftHead_oriGeneAs1st->return_startOrEndLocInReadVec_oriGene(0) - 1;
				tmpAnchor_right = readLength_1 - tmpAnchor_left;
			}
			else if(strand_fusion == "-+") // case 11 
			{
				tmpStrand_gene1 = "-";
				tmpStrand_gene2 = "+";
				tmpFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					tmpChrNameStr_leftHead, tmpBreakPointPos_leftHead,
					tmpChrNameStr_oriGene, tmpBreakPointPos_oriGene,
					false, true, true);	
				tmpAnchor_left 
					= tmpRemapInfo_leftHead_oriGeneAs1st->return_startOrEndLocInReadVec_oriGene(0) - 1;
				tmpAnchor_right = readLength_1 - tmpAnchor_left;
			}
			else
			{
				cout << "error in generateFusionRemapInfo_leftHead ...." << endl;
				exit(1);
			}
		}
		else // case 2, case 10, oriGene == Gene 2
		{
			tmpChrNameStr_gene2 = tmpChrNameStr_oriGene;
			tmpChrNameStr_gene1 = tmpChrNameStr_leftHead;
			tmpBreakPointPos_gene2 = tmpBreakPointPos_oriGene;
			tmpBreakPointPos_gene1 = tmpBreakPointPos_leftHead;	
			if(strand_fusion == "++") // case 2
			{
				tmpStrand_gene1 = "+";
				tmpStrand_gene2 = "+";
				tmpFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					tmpChrNameStr_leftHead, tmpBreakPointPos_leftHead,
					tmpChrNameStr_oriGene, tmpBreakPointPos_oriGene,
					true, true, true);
				tmpAnchor_left 
					= tmpRemapInfo_leftHead_oriGeneAs2nd->return_startOrEndLocInReadVec_oriGene(0) - 1;
				tmpAnchor_right = readLength_1 - tmpAnchor_left;
			}
			else if(strand_fusion == "-+") // case 10 
			{
				tmpStrand_gene1 = "-";
				tmpStrand_gene2 = "+";		
				tmpFusionJuncFlankString = indexInfo->returnRawFusionJuncFlankString(
					tmpChrNameStr_leftHead, tmpBreakPointPos_leftHead,
					tmpChrNameStr_oriGene, tmpBreakPointPos_oriGene,
					false, true, true);
				tmpAnchor_left 
					= tmpRemapInfo_leftHead_oriGeneAs2nd->return_startOrEndLocInReadVec_oriGene(0) - 1;
				tmpAnchor_right = readLength_1 - tmpAnchor_left;				
			}
			else
			{
				cout << "error in generateFusionRemapInfo_leftHead ...." << endl;
				exit(1);
			}			
		}
	}

	string returnFusionBreakPointInfoStr_raw(
		bool leftReadHeadUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd, 
		bool rightReadTailUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs2nd,
		Index_Info* indexInfo)
	{
		string tmpFusionBreakPointInfoStr = "";
		string tmpBreakPointInfoStr_leftHead, tmpBreakPointInfoStr_rightTail;
		if(leftReadHeadUnfixed_fixed_bool)
		{
			string tmpChrNameStr_leftHead, tmpChrNameStr_oriGene, 
				strand_leftHead, strand_oriGene;
			int tmpBreakPointPos_leftHead, tmpBreakPointPos_oriGene;
			bool forOrRevBool_leftHead, forOrRevBool_oriGene, oriGene_1_or_2_bool;
			this->generateFusionInfoFrom2remapInfo_leftHead(
				tmpRemapInfo_leftHead_oriGeneAs1st, tmpRemapInfo_leftHead_oriGeneAs2nd,
				tmpChrNameStr_leftHead, tmpBreakPointPos_leftHead,
				tmpChrNameStr_oriGene, tmpBreakPointPos_oriGene,
				forOrRevBool_leftHead, forOrRevBool_oriGene,
				oriGene_1_or_2_bool, indexInfo);
			if(oriGene_1_or_2_bool)
				tmpBreakPointInfoStr_leftHead = readName_1 + "\t"
					+ tmpChrNameStr_leftHead + "\t" + tmpChrNameStr_oriGene + "\t"
					+ int_to_str(tmpBreakPointPos_leftHead) + "\t" + int_to_str(tmpBreakPointPos_oriGene) 
					+ "\t"
		} 
		if(rightReadTailUnfixed_fixed_bool)
			tmpBreakPointInfoStr_rightTail = 

		if(leftReadHeadUnfixed_fixed_bool && rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead 
				+ "\t" + tmpBreakPointInfoStr_rightTail;
		else if(leftReadHeadUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead;
		else if(rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_rightTail;
		else
		{}
		return tmpFusionBreakPointInfoStr;
	}
	*/
	void generateFusionInfo2reportInBreakPointFile_leftHead(
		RemapAgainstFusionBreakPoint_Info*  tmpRemapInfo_leftHead_oriGeneAs1st, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd,
		string& tmpChrNameStr_gene1, string& tmpChrNameStr_gene2,
		int& tmpBreakPointPos_gene1, int& tmpBreakPointPos_gene2,
		int& tmpFusionAnchorLength_gene1, int& tmpFusionAnchorLength_gene2, Index_Info* indexInfo)
	{
		//cout << "function generateFusionInfo2reportInBreakPointFile_leftHead( starts ......" << endl;
		if(tmpRemapInfo_leftHead_oriGeneAs1st->returnResultsSize()
			+ tmpRemapInfo_leftHead_oriGeneAs2nd->returnResultsSize() != 1)
		{
			cout << "error in generateFusionInfo2reportInBreakPointFile_leftHead ...." << endl;
			cout << "tmpRemapInfo_leftHead_oriGeneAsGene1->returnResultsSize()";
			cout << " + tmpRemapInfo_leftHead_oriGeneAsGene2->returnResultsSize() != 1" << endl;
			exit(1);
		}

		string tmpChrNameStr_leftGene;
		string tmpChrNameStr_rightGene;
		int tmpBreakPointPos_leftGene;
		int tmpBreakPointPos_rightGene;
		int tmpFusionAnchorLength_left;
		int tmpFusionAnchorLength_right;
		if(tmpRemapInfo_leftHead_oriGeneAs1st->returnResultsSize() == 1)
		{
			//cout << "tmpRemapInfo_leftHead_oriGeneAs1st->returnResultsSize() == 1" << endl;
			int tmpChrNameInt_leftGene 
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_chrNameInt_theOtherGene(0);
			tmpChrNameStr_leftGene = indexInfo->returnChrNameStr(tmpChrNameInt_leftGene);
			tmpChrNameStr_rightGene = indexInfo->returnChrNameStr(chrNameInt);
			tmpBreakPointPos_leftGene 
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_breakPointPos_theOtherGene(0);
			tmpBreakPointPos_rightGene
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_breakPointPos_oriGene(0);
			int startLocInRead_oriGene 
				= tmpRemapInfo_leftHead_oriGeneAs1st->return_startOrEndLocInReadVec_oriGene(0);
			tmpFusionAnchorLength_left = startLocInRead_oriGene - 1;
			tmpFusionAnchorLength_right = readLength_1 - tmpFusionAnchorLength_left;
			
			tmpChrNameStr_gene1 = tmpChrNameStr_rightGene;
			tmpChrNameStr_gene2 = tmpChrNameStr_leftGene;
			tmpBreakPointPos_gene1 = tmpBreakPointPos_rightGene;
			tmpBreakPointPos_gene2 = tmpBreakPointPos_leftGene;
			tmpFusionAnchorLength_gene1 = tmpFusionAnchorLength_right;
			tmpFusionAnchorLength_gene2 = tmpFusionAnchorLength_left;			
		}
		else if(tmpRemapInfo_leftHead_oriGeneAs2nd->returnResultsSize() == 1)
		{
			//cout << "tmpRemapInfo_leftHead_oriGeneAs2nd->returnResultsSize() == 1" << endl;
			int tmpChrNameInt_leftGene 
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_chrNameInt_theOtherGene(0);
			//cout << "tmpChrNameInt_leftGene: " << tmpChrNameInt_leftGene << endl;
			tmpChrNameStr_leftGene = indexInfo->returnChrNameStr(tmpChrNameInt_leftGene);
			//cout << "tmpChrNameStr_leftGene: " << tmpChrNameStr_leftGene << endl;
			tmpChrNameStr_rightGene = indexInfo->returnChrNameStr(chrNameInt);
			//cout << "tmpChrNameStr_rightGene: " << tmpChrNameStr_rightGene << endl;
			tmpBreakPointPos_leftGene 
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_breakPointPos_theOtherGene(0);
			//cout << "tmpBreakPointPos_leftGene: " << tmpBreakPointPos_leftGene << endl;
			tmpBreakPointPos_rightGene
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_breakPointPos_oriGene(0);
			//cout << "tmpBreakPointPos_rightGene: " << tmpBreakPointPos_rightGene << endl;
			int startLocInRead_oriGene 
				= tmpRemapInfo_leftHead_oriGeneAs2nd->return_startOrEndLocInReadVec_oriGene(0);
			tmpFusionAnchorLength_left = startLocInRead_oriGene - 1;
			tmpFusionAnchorLength_right = readLength_1 - tmpFusionAnchorLength_left;
			
			tmpChrNameStr_gene1 = tmpChrNameStr_leftGene;
			tmpChrNameStr_gene2 = tmpChrNameStr_rightGene;
			tmpBreakPointPos_gene1 = tmpBreakPointPos_leftGene;
			tmpBreakPointPos_gene2 = tmpBreakPointPos_rightGene;
			tmpFusionAnchorLength_gene1 = tmpFusionAnchorLength_left;
			tmpFusionAnchorLength_gene2 = tmpFusionAnchorLength_right;	
		}
		else
		{
			cout << "error in generateFusionInfo2reportInBreakPointFile_leftHead ...." << endl;
			cout << "tmpRemapInfo_leftHead_oriGeneAsGene1->returnResultsSize()";
			cout << " + tmpRemapInfo_leftHead_oriGeneAsGene2->returnResultsSize() != 1" << endl;
			exit(1);
		}
	}

	void generateFusionInfo2reportInBreakPointFile_rightTail(
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs1st, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs2nd,
		string& tmpChrNameStr_gene1, string& tmpChrNameStr_gene2,
		int& tmpBreakPointPos_gene1, int& tmpBreakPointPos_gene2,
		int& tmpFusionAnchorLength_gene1, int& tmpFusionAnchorLength_gene2, Index_Info* indexInfo)
	{
		if(tmpRemapInfo_rightTail_oriGeneAs1st->returnResultsSize()
			+ tmpRemapInfo_rightTail_oriGeneAs2nd->returnResultsSize() != 1)
		{
			cout << "error in generateFusionInfo2reportInBreakPointFile_leftHead ...." << endl;
			cout << "tmpRemapInfo_rightTail_oriGeneAsGene1->returnResultsSize()";
			cout << " + tmpRemapInfo_rightTail_oriGeneAsGene2->returnResultsSize() != 1" << endl;
			exit(1);
		}

		string tmpChrNameStr_leftGene;
		string tmpChrNameStr_rightGene;
		int tmpBreakPointPos_leftGene;
		int tmpBreakPointPos_rightGene;
		int tmpFusionAnchorLength_left;
		int tmpFusionAnchorLength_right;
		if(tmpRemapInfo_rightTail_oriGeneAs1st->returnResultsSize() == 1)
		{
			int tmpChrNameInt_rightGene 
				= tmpRemapInfo_rightTail_oriGeneAs1st->return_chrNameInt_theOtherGene(0);
			tmpChrNameStr_rightGene = indexInfo->returnChrNameStr(tmpChrNameInt_rightGene);
			tmpChrNameStr_leftGene = indexInfo->returnChrNameStr(chrNameInt);
			tmpBreakPointPos_rightGene 
				= tmpRemapInfo_rightTail_oriGeneAs1st->return_breakPointPos_theOtherGene(0);
			tmpBreakPointPos_leftGene
				= tmpRemapInfo_rightTail_oriGeneAs1st->return_breakPointPos_oriGene(0);
			int endLocInRead_oriGene 
				= tmpRemapInfo_rightTail_oriGeneAs1st->return_startOrEndLocInReadVec_oriGene(0);
			tmpFusionAnchorLength_left = endLocInRead_oriGene;
			tmpFusionAnchorLength_right = readLength_2 - endLocInRead_oriGene;

			tmpChrNameStr_gene1 = tmpChrNameStr_leftGene;
			tmpChrNameStr_gene2 = tmpChrNameStr_rightGene;
			tmpBreakPointPos_gene1 = tmpBreakPointPos_leftGene;
			tmpBreakPointPos_gene2 = tmpBreakPointPos_rightGene;
			tmpFusionAnchorLength_gene1 = tmpFusionAnchorLength_left;
			tmpFusionAnchorLength_gene2 = tmpFusionAnchorLength_right;	
		}
		else if(tmpRemapInfo_rightTail_oriGeneAs2nd->returnResultsSize() == 1)
		{
			int tmpChrNameInt_rightGene 
				= tmpRemapInfo_rightTail_oriGeneAs2nd->return_chrNameInt_theOtherGene(0);
			tmpChrNameStr_rightGene = indexInfo->returnChrNameStr(tmpChrNameInt_rightGene);
			tmpChrNameStr_leftGene = indexInfo->returnChrNameStr(chrNameInt);
			tmpBreakPointPos_rightGene 
				= tmpRemapInfo_rightTail_oriGeneAs2nd->return_breakPointPos_theOtherGene(0);
			tmpBreakPointPos_leftGene
				= tmpRemapInfo_rightTail_oriGeneAs2nd->return_breakPointPos_oriGene(0);
			int endLocInRead_oriGene 
				= tmpRemapInfo_rightTail_oriGeneAs2nd->return_startOrEndLocInReadVec_oriGene(0);
			tmpFusionAnchorLength_left = endLocInRead_oriGene;
			tmpFusionAnchorLength_right = readLength_2 - endLocInRead_oriGene;

			tmpChrNameStr_gene1 = tmpChrNameStr_rightGene;
			tmpChrNameStr_gene2 = tmpChrNameStr_leftGene;
			tmpBreakPointPos_gene1 = tmpBreakPointPos_rightGene;
			tmpBreakPointPos_gene2 = tmpBreakPointPos_leftGene;
			tmpFusionAnchorLength_gene1 = tmpFusionAnchorLength_right;
			tmpFusionAnchorLength_gene2 = tmpFusionAnchorLength_left;
		}
		else//(tmpRemapInfo_rightTail_oriGeneAs1st->returnResultsSize()+ tmpRemapInfo_rightTail_oriGeneAs2nd->returnResultsSize() != 1)
		{
			cout << "error in generateFusionInfo2reportInBreakPointFile_leftHead ...." << endl;
			cout << "tmpRemapInfo_rightTail_oriGeneAsGene1->returnResultsSize()";
			cout << " + tmpRemapInfo_rightTail_oriGeneAsGene2->returnResultsSize() != 1" << endl;
			exit(1);
		}
	}

	string returnFusionBreakPointInfoStr_supportNumIncrement(
		bool leftReadHeadUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd, 
		bool rightReadTailUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs2nd,
		FusionBreakPointHash_Info* tmpFusionBreakPointHashInfo,
		Index_Info* indexInfo)
	{
		string tmpFusionBreakPointInfoStr = "";
		string tmpBreakPointInfoStr_leftHead, tmpBreakPointInfoStr_rightTail;
		//cout << "leftReadHeadUnfixed_fixed_bool: " << leftReadHeadUnfixed_fixed_bool << endl;
		//cout << "rightReadTailUnfixed_fixed_bool: " << rightReadTailUnfixed_fixed_bool << endl;
		if(leftReadHeadUnfixed_fixed_bool)
		{
			string tmpChrNameStr_gene1, tmpChrNameStr_gene2;
				//tmpStrand_gene1, tmpStrand_gene2, tmpFusionJuncFlankString_adjusted;
			int tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, 
				tmpFusionAnchorLength_gene1, tmpFusionAnchorLength_gene2;
			//cout << "start to generateFusionInfo2reportInBreakPointFile_leftHead " << endl;
			this->generateFusionInfo2reportInBreakPointFile_leftHead(
				tmpRemapInfo_leftHead_oriGeneAs1st, tmpRemapInfo_leftHead_oriGeneAs2nd,
				tmpChrNameStr_gene1, tmpChrNameStr_gene2,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2,
				tmpFusionAnchorLength_gene1, tmpFusionAnchorLength_gene2, indexInfo);

			tmpBreakPointInfoStr_leftHead = readName_1 + "\t"
				+ tmpChrNameStr_gene1 + "\t" + tmpChrNameStr_gene2 + "\t"
				+ int_to_str(tmpBreakPointPos_gene1) + "\t"
				+ int_to_str(tmpBreakPointPos_gene2) + "\t"
				+ int_to_str(tmpFusionAnchorLength_gene1) + "\t"
				+ int_to_str(tmpFusionAnchorLength_gene2);
			//cout << "start to do supportNumIncrement" << endl;
			this->supportNumIncrement(
				tmpChrNameStr_gene1, tmpChrNameStr_gene2, 
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, 
				tmpFusionBreakPointHashInfo, indexInfo);
		}
		if(rightReadTailUnfixed_fixed_bool)
		{
			string tmpChrNameStr_gene1, tmpChrNameStr_gene2;
				//tmpStrand_gene1, tmpStrand_gene2, tmpFusionJuncFlankString_adjusted;
			int tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, 
				tmpFusionAnchorLength_gene1, tmpFusionAnchorLength_gene2;
			//cout << "start to generateFusionInfo2reportInBreakPointFile_rightTail " << endl;
			this->generateFusionInfo2reportInBreakPointFile_rightTail(
				tmpRemapInfo_rightTail_oriGeneAs1st, tmpRemapInfo_rightTail_oriGeneAs2nd,
				tmpChrNameStr_gene1, tmpChrNameStr_gene2, 
				//tmpStrand_gene1, tmpStrand_gene2,
				//tmpFusionJuncFlankString_adjusted,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2,
				tmpFusionAnchorLength_gene1, tmpFusionAnchorLength_gene2, indexInfo);

			tmpBreakPointInfoStr_rightTail = readName_2 + "\t"
				+ tmpChrNameStr_gene1 + "\t" + tmpChrNameStr_gene2 + "\t"
				+ int_to_str(tmpBreakPointPos_gene1) + "\t"
				+ int_to_str(tmpBreakPointPos_gene2) + "\t"
				+ int_to_str(tmpFusionAnchorLength_gene1) + "\t"
				+ int_to_str(tmpFusionAnchorLength_gene2);
			//cout << "start to do supportNumIncrement " << endl;
			this->supportNumIncrement(
				tmpChrNameStr_gene1, tmpChrNameStr_gene2, 
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, 
				tmpFusionBreakPointHashInfo, indexInfo);
		}
		//cout << "start to generate final fusionBreakPointHashInfoStr ..." << endl;
		if(leftReadHeadUnfixed_fixed_bool && rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead 
				+ "\t" + tmpBreakPointInfoStr_rightTail;
		else if(leftReadHeadUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead;
		else if(rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_rightTail;
		else
		{}
		//cout << "tmpFusionBreakPointInfoStr: " << tmpFusionBreakPointInfoStr << endl;
		return tmpFusionBreakPointInfoStr;
	}

	/*
	string returnFusionBreakPointInfoStr_raw(
		bool leftReadHeadUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd, 
		bool rightReadTailUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs2nd,
		Index_Info* indexInfo)
	{
		string tmpFusionBreakPointInfoStr = "";
		string tmpBreakPointInfoStr_leftHead, tmpBreakPointInfoStr_rightTail;
		if(leftReadHeadUnfixed_fixed_bool)
		{
			string tmpChrNameStr_leftGene, tmpChrNameStr_rightGene, 
				tmpFusionMatchDirStr, tmpFusionJuncFlankString_raw;
			int tmpBreakPointPos_leftGene, tmpBreakPointPos_rightGene,
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right;
			this->generateFusionInfo2reportInRawBreakPointFile_leftHead(
				tmpRemapInfo_leftHead_oriGeneAs1st, tmpRemapInfo_leftHead_oriGeneAs2nd,
				tmpChrNameStr_leftGene, tmpChrNameStr_rightGene,
				tmpBreakPointPos_leftGene, tmpBreakPointPos_rightGene,
				tmpFusionMatchDirStr, tmpFusionJuncFlankString_raw,
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right, indexInfo);

			tmpBreakPointInfoStr_leftHead = readName_1 + "\t"
				+ tmpChrNameStr_leftGene + "\t" + tmpChrNameStr_rightGene + "\t"
				+ int_to_str(tmpBreakPointPos_leftGene) + "\t"
				+ int_to_str(tmpBreakPointPos_rightGene) + "\t"
				+ tmpFusionMatchDirStr + "\t" + tmpFusionJuncFlankString_raw + "\t"
				+ int_to_str(tmpFusionAnchorLength_left) + "\t"
				+ int_to_str(tmpFusionAnchorLength_right);
		}
		if(rightReadTailUnfixed_fixed_bool)
		{
			string tmpChrNameStr_leftGene, tmpChrNameStr_rightGene, 
				tmpFusionMatchDirStr, tmpFusionJuncFlankString_raw;
			int tmpBreakPointPos_leftGene, tmpBreakPointPos_rightGene,
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right;
			this->generateFusionInfo2reportInRawBreakPointFile_rightTail(
				tmpRemapInfo_rightTail_oriGeneAs1st, tmpRemapInfo_rightTail_oriGeneAs2nd,
				tmpChrNameStr_leftGene, tmpChrNameStr_rightGene,
				tmpBreakPointPos_leftGene, tmpBreakPointPos_rightGene,
				tmpFusionMatchDirStr, tmpFusionJuncFlankString_raw,
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right, indexInfo);

			tmpBreakPointInfoStr_rightTail = readName_2 + "\t"
				+ tmpChrNameStr_leftGene + "\t" + tmpChrNameStr_rightGene + "\t"
				+ int_to_str(tmpBreakPointPos_leftGene) + "\t"
				+ int_to_str(tmpBreakPointPos_rightGene) + "\t"
				+ tmpFusionMatchDirStr + "\t" + tmpFusionJuncFlankString_raw + "\t"
				+ int_to_str(tmpFusionAnchorLength_left) + "\t"
				+ int_to_str(tmpFusionAnchorLength_right);
		}

		if(leftReadHeadUnfixed_fixed_bool && rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead 
				+ "\t" + tmpBreakPointInfoStr_rightTail;
		else if(leftReadHeadUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead;
		else if(rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_rightTail;
		else
		{}
		return tmpFusionBreakPointInfoStr;
	}*/
	/*
	string returnFusionBreakPointInfoStr_adjusted(
		bool leftReadHeadUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_leftHead_oriGeneAs2nd, 
		bool rightReadTailUnfixed_fixed_bool, 
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs1st,
		RemapAgainstFusionBreakPoint_Info* tmpRemapInfo_rightTail_oriGeneAs2nd,
		Index_Info* indexInfo)
	{
		string tmpFusionBreakPointInfoStr = "";
		string tmpBreakPointInfoStr_leftHead, tmpBreakPointInfoStr_rightTail;
		if(leftReadHeadUnfixed_fixed_bool)
		{
			string tmpChrNameStr_gene1, tmpChrNameStr_gene2,
				//tmpStrand_gene1, tmpStrand_gene2, tmpFusionJuncFlankString_adjusted;
			int tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, 
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right;
			this->generateFusionInfo2reportInAdjustedBreakPointFile_leftHead(
				tmpRemapInfo_leftHead_oriGeneAs1st, tmpRemapInfo_leftHead_oriGeneAs2nd,
				tmpChrNameStr_gene1, tmpChrNameStr_gene2, tmpStrand_gene1, tmpStrand_gene2,
				tmpFusionJuncFlankString_adjusted,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2,
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right, indexInfo);

			tmpBreakPointInfoStr_leftHead = readName_1 + "\t"
				+ tmpChrNameStr_gene1 + "\t" + tmpChrNameStr_gene2 + "\t"
				+ int_to_str(tmpBreakPointPos_gene1) + "\t"
				+ int_to_str(tmpBreakPointPos_gene2) + "\t"
				+ tmpStrand_gene1 + "\t" + tmpStrand_gene2 + "\t"
				+ tmpFusionJuncFlankString_adjusted + "\t"
				+ int_to_str(tmpFusionAnchorLength_left) + "\t"
				+ int_to_str(tmpFusionAnchorLength_right);
		}
		if(rightReadTailUnfixed_fixed_bool)
		{
			string tmpChrNameStr_gene1, tmpChrNameStr_gene2, 
				tmpStrand_gene1, tmpStrand_gene2, tmpFusionJuncFlankString_adjusted;
			int tmpBreakPointPos_gene1, tmpBreakPointPos_gene2, 
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right;
			this->generateFusionInfo2reportInAdjustedBreakPointFile_rightTail(
				tmpRemapInfo_rightTail_oriGeneAs1st, tmpRemapInfo_rightTail_oriGeneAs2nd,
				tmpChrNameStr_gene1, tmpChrNameStr_gene2, tmpStrand_gene1, tmpStrand_gene2,
				tmpFusionJuncFlankString_adjusted,
				tmpBreakPointPos_gene1, tmpBreakPointPos_gene2,
				tmpFusionAnchorLength_left, tmpFusionAnchorLength_right, indexInfo);

			tmpBreakPointInfoStr_rightTail = readName_2 + "\t"
				+ tmpChrNameStr_gene1 + "\t" + tmpChrNameStr_gene2 + "\t"
				+ int_to_str(tmpBreakPointPos_gene1) + "\t"
				+ int_to_str(tmpBreakPointPos_gene2) + "\t"
				+ tmpStrand_gene1 + "\t" + tmpStrand_gene2 + "\t"
				+ tmpFusionJuncFlankString_adjusted + "\t"
				+ int_to_str(tmpFusionAnchorLength_left) + "\t"
				+ int_to_str(tmpFusionAnchorLength_right);
		}

		if(leftReadHeadUnfixed_fixed_bool && rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead 
				+ "\t" + tmpBreakPointInfoStr_rightTail;
		else if(leftReadHeadUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_leftHead;
		else if(rightReadTailUnfixed_fixed_bool)
			tmpFusionBreakPointInfoStr = tmpBreakPointInfoStr_rightTail;
		else
		{}
		return tmpFusionBreakPointInfoStr;
	}*/
	/*
	string returnFixedFusionSamStr(
		bool leftReadHeadUnfixed_fixed_bool,
		bool rightReadTailUnfixed_fixed_bool,
		string& tmpSAMstr_1, string& tmpSAMstr_2)
	{
		string tmpFinalSAMstr;
		string tmpSAMstr_fixedFusionInLeftHead,
			tmpSAMStr_fixedFusionInRightTail;
		if(leftReadHeadUnfixed_fixed_bool)
			tmpSAMstr_fixedFusionInLeftHead = 

		if(rightReadTailUnfixed_fixed_bool)
			tmpSAMStr_fixedFusionInRightTail = 

		if(leftReadHeadUnfixed_fixed_bool && rightReadTailUnfixed_fixed_bool)
			tmpFinalSAMstr = tmpSAMstr_fixedFusionInLeftHead + "\n"
				+ tmpSAMstr_1 + "\n" + tmpSAMstr_2 + "\n"
				+ tmpSAMStr_fixedFusionInRightTail;
		else if(leftReadHeadUnfixed_fixed_bool)
			tmpFinalSAMstr = tmpSAMstr_fixedFusionInLeftHead + "\n"
				+ tmpSAMstr_1 + "\n" + tmpSAMstr_2;
		else if(rightReadTailUnfixed_fixed_bool)
			tmpFinalSAMstr = tmpSAMstr_1 + "\n" + tmpSAMstr_2 + "\n"
				+ tmpSAMStr_fixedFusionInRightTail;
		else
		{
			cout << "error inreturn returnFixedFusionSamStr ..." << endl;
			exit(1);
		}
		return tmpFinalSAMstr;
	}*/
};
#endif