// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef INCOMPLETESEALIGNMENT2REMAPAGAINSTFUSIONBREAKPOINT_INFO_H
#define INCOMPLETESEALIGNMENT2REMAPAGAINSTFUSIONBREAKPOINT_INFO_H

#include "../../../general/index_info.h"
#include "../../../general/splice_info.h"

using namespace std;

class IncompleteSeAlignment2remapAgainstFusionBreakPoint_Info
{
private:	
	bool Nor_or_Rcm_bool;
	string readName;
	int chrNameInt;
	int startPos;
	int endPos;
	vector<Jump_Code> cigarStringJumpCodeVec;
	string readSeq;
	int readLength;
	string qualSeq;

	int unfixedHeadLength;
	vector<int> chrNameIntVec_fixedHead;
	vector<int> chrPosVec_fixedHead; 
	vector<bool> NorOrRcmBoolVec_fixedHead;
	vector< vector<Jump_Code> > jumpCodeVecVec_fixedHead;
	vector<int> rightBreakPointLocVecInRead_fixedHead; 

	int unfixedTailLength;
	vector<int> chrNameIntVec_fixedTail;
	vector<int> chrPosVec_fixedTail; 
	vector<bool> NorOrRcmBoolVec_fixedTail;
	vector< vector<Jump_Code> > jumpCodeVecVec_fixedTail;
	vector<int> rightBreakPointLocVecInRead_fixedTail; 

public:
	IncompleteSeAlignment2remapAgainstFusionBreakPoint_Info()
	{}

	void remapAgainstFusionBreakPoint_unfixedHead_oriFor(
		FusionBreakPointHash_Info* fusionBreakPointHashInfo)
	{
		int firstMatchLengthInOriCigarString 
			= cigarStringJumpCodeVec[2].len;
		int startMatchLocInRead = unfixedHeadLength + 1;
		int startMatchPosInChr = startPos;
		int endMatchLocInRead 
			= unfixedHeadLength + 1 + firstMatchLengthInOriCigarString - 1;
		int endMatchPosInChr 
			= startPos + firstMatchLengthInOriCigarString - 1;

		vector<int> fusionBreakPointPosVec_gene1;	
		vector<int> fusionBreakPointPosVec_gene2;	
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameInt, startMatchPosInChr, endMatchPosInChr, true,
			fusionBreakPointPosVec_gene1);
		fusionBreakPointHashInfo->searchFusionBreakPointFromAreaHashWithinRange(
			chrNameInt, startMatchPosInChr, endMatchPosInChr, false,
			fusionBreakPointPosVec_gene2);

		// case 2
		vector< pair<int,int> > candiFusionGeneChrNameIntVec_case2; // chrNameInt_gene1, chrNameInt_gene2
		vector< pair<int,int> > candiFusionGeneChrPosVec_case2; // chrPos_gene1, chrPos_gene2
		for(int tmp = 0; tmp < fusionBreakPointPosVec_gene2.size(); tmp++)
		{
			string tmpThisGeneStrand = "+";
			string tmpTheOtherGeneStrand = "+";
			int tmpFusionBreakPointPos_gene2 = fusionBreakPointPosVec_gene2[tmp];
			vector<int> tmpOtherGeneChrNameIntVec;
			vector<int> tmpOtherGeneBreakPointPosVec;
			fusionBreakPointHashInfo->returnTheOtherFusionBreakPoint_afterCheckingStrand(
				chrNameInt, tmpFusionBreakPointPos_gene2, false,
				tmpThisGeneStrand, tmpTheOtherGeneStrand,
				tmpOtherGeneChrNameIntVec, tmpOtherGeneBreakPointPosVec);
			for(int tmp2 = 0; tmp2 < tmpOtherGeneChrNameIntVec.size(); tmp2 ++)
			{
				int 
			}
		}
		// case 5
		vector< pair<int,int> > candiFusionGeneChrNameIntVec_case5;
		vector< pair<int,int> > candiFusionGeneChrPosVec_case5;
		// case 10
		vector< pair<int,int> > candiFusionGeneChrNameIntVec_case10;
		vector< pair<int,int> > candiFusionGeneChrPosVec_case10;
		// case 11
		vector< pair<int,int> > candiFusionGeneChrNameIntVec_case11;
		vector< pair<int,int> > candiFusionGeneChrPosVec_case11;

	}

	void remapAgainstFusionBreakPoint_unfixedHead_oriRev(
		FusionBreakPointHash_Info* fusionBreakPointHashInfo)
	{
	}	

	void remapAgainstFusionBreakPoint_unfixedTail_oriFor(
		FusionBreakPointHash_Info* fusionBreakPointHashInfo)
	{
		int lastMatchLengthInOriCigarString 
			= cigarStringJumpCodeVec[cigarStringJumpCodeVec.size()-2].len;
		int startMatchLocInRead 
			= readLength - unfixedTailLength - lastMatchLengthInOriCigarString + 1;
		int startMatchPosInChr 
			= endPos - lastMatchLengthInOriCigarString + 1;
		int endMatchLocInRead = readLength - unfixedTailLength;
		int endMatchPosInChr = endPos;

	}

	void remapAgainstFusionBreakPoint_unfixedTail_oriRev(
		FusionBreakPointHash_Info* fusionBreakPointHashInfo)
	{

	}	

	int returnUnfixedHeadLength()
	{
		return unfixedHeadLength;
	}

	int returnUnfixedTailLength()
	{
		return unfixedTailLength;
	}

	bool NorOrRcmBool()
	{
		return Nor_or_Rcm_bool;
	}

	bool headUnfixed()
	{
		return (unfixedHeadLength > 0);
	}

	bool tailUnfixed()
	{
		return (unfixedTailLength > 0);
	}

	bool initiateWithSamStr_remapUnfixedEnd2fusionBreakPoint(
		const string& tmpSamStr, Index_Info* indexInfo)
	{
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

		string chrNameStr_1 = samFieldVec_1[2];
		chrNameInt = indexInfo->convertStringToInt(chrNameStr_1.c_str());
		if(chrNameStr_1 == "*")
			return false;
		string IHfield_1 = samFieldVec_1[12];
		string IHintStr_1 = IHfield_1.substr(5);
		int IHint_1 = atoi(IHintStr_1.c_str());
		if(IHint_1 != 1)
			return false;
		string flagStr_1 = samFieldVec_1[1];
		int flagInt_1 = atoi(flagStr_1.c_str());

		readName = samFieldVec_1[0];
		string startPosStr_1 = samFieldVec_1[3];
		startPos = atoi(startPosStr_1.c_str());
		string cigarString_1 = samFieldVec_1[5];
		this->cigarString2jumpCodeVec(cigarString_1, cigarStringJumpCodeVec);
		
		if(cigarStringJumpCodeVec[0].type == "S")
			unfixedHeadLen = cigarStringJumpCodeVec[0].len;
		else
			unfixedHeadLen = 0;
		
		if(cigarStringJumpCodeVec[cigarString2jumpCodeVec.size()-1].type == "S")
			unfixedTailLen = cigarStringJumpCodeVec[cigarString2jumpCodeVec.size()-1].len;
		else
			unfixedTailLen = 0;
	
		if((unfixedHeadLen == 0)||(unfixedTailLen == 0))
			return false;

		readSeq = samFieldVec_1[9];
		readLength = readSeq.length();
		qualSeq = samFieldVec_1[10];

		if(flagInt_1 * 0x10)
			Nor_or_Rcm_bool = false;
		else
			Nor_or_Rcm_bool = true;

		endPos = getEndPos();

		return true;
	}

	int returnCigarStringJumpCodeVecSize()
	{
		return cigarStringJumpCodeVec.size();
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnStartPos()
	{
		return startPos;
	}

	int returnEndPos()
	{
		return endPos;
	}

	int getEndPos()
	{
		int tmpEndPos = this->getEndPosOfSpecificJumpCode(startPos, cigarStringJumpCodeVec,
			cigarStringJumpCodeVec.size()-1);
		return tmpEndPos;
	}

	int getEndPosOfSpecificJumpCode(int startPos, 
		vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
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

};
#endif