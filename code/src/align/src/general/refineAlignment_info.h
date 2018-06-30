// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REFALIGNMENT_INFO_H
#define REFALIGNMENT_INFO_H

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

using namespace std;


// inline char getCharRevComp(char ch)
// {
// 	int chInt = ch - 'A';
// 	static const char alphatChar[26] = {'T', 'N', 'G', 'N', 'N', 'N', 'C',
// 		'N', 'N', 'N', 'N', 'N', 'N', 'N',
// 		'N', 'N', 'N', 'N', 'N', 'A',
// 		'N', 'N', 'N', 'N', 'N', 'N'};
// 	return alphatChar[chInt];
// }

// string getRcmSeq(const string& readSeq)
// {
// 	int readSeqLength = readSeq.length();

// 	char readRcmSeqChar[readSeqLength];

// 	readRcmSeqChar[0] = getCharRevComp(readSeq.at(readSeqLength-1));

// 	for(int tmp = 1; tmp < readSeqLength; tmp ++)
// 	{
// 		readRcmSeqChar[tmp] = getCharRevComp((readSeq.at(readSeqLength - tmp - 1)));
// 	}
// 	string rcmSeq = readRcmSeqChar;
// 	return rcmSeq.substr(0, readSeqLength);
// }

class SingleAlignment_Info
{
private:
	string readName;
	int chrNameInt;
	int chrMapPos;
	int chrMapPos_end;
	vector<Jump_Code> jumpCodeVec;
	string readSeq;
	string qualSeq;

	int flag;

	//string tagStr;
	//vector<int> mismatchPosVec;
	//vector<char> mismatchCharVec;
public:
	SingleAlignment_Info()
	{}

	string returnReadSeq()
	{
		return readSeq;
	}

	string returnRcmReadSeq()
	{
		string rcmReadSeq = getRcmSeq(readSeq);
		return rcmReadSeq;
	}

	bool wholeMatch_bool()
	{
		if(jumpCodeVec.size() == 1)
		{
			if(jumpCodeVec[0].type == "M")
				return true;
			else
				return false;
		}
		else
			return false;
	}

	bool withoutSJ_bool()
	{
		int jumpCodeVecSize = jumpCodeVec.size();
		for(int tmp = 0; tmp < jumpCodeVecSize; tmp++)
		{
			string tmpJumpCodeType = jumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
				return false;
		}
		return true;
	}

	string returnClippedAlignmentWithInterval(Index_Info* indexInfo, 
		int jumpCodeIndex_start, int jumpCodeIndex_end)
	{
		string tmpStr;
		int tmpNewChrMapPos = this->getEndPosOfSpecificJumpCode(chrMapPos, 
			jumpCodeVec, jumpCodeIndex_start-1) + 1;
		tmpStr = tmpStr + readName + "\t0\t" + indexInfo->returnChrNameStr(chrNameInt)
			+ "\t" + int_to_str(tmpNewChrMapPos) + "\t255\t" 
			+ this->jumpCodeVec2cigarString_withInterval(jumpCodeIndex_start, jumpCodeIndex_end)
			+ "\t*\t0\t0\t" + readSeq + "\t" + qualSeq + "\tHI:i:1\tIH:i:1";
		return tmpStr;
	}

	string returnUnmappedAlignment(Index_Info* indexInfo)
	{
		string tmpStr;
		tmpStr = tmpStr + readName + "\t0\t*\t0\t255\t*\t*\t0\t0\t" + readSeq + "\t" + qualSeq + "\tHI:i:1\tIH:i:1";
		return tmpStr;
	}

	void generateNewJumpCodeVecWithInterval(
		vector<Jump_Code>& targetJumpCodeVec, 
		int jumpCodeIndex_start, int jumpCodeIndex_end)
	{
		int headSoftClippedSeqLen = this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec, jumpCodeIndex_start-1);

		int tmpReadLength = readSeq.length();
		int tailSoftClippedSeqLen = tmpReadLength - this->getEndLocInReadOfSpecificJumpCode(jumpCodeVec, jumpCodeIndex_end);
		if(headSoftClippedSeqLen > 0)
		{
			Jump_Code tmpHeadSoftClipJumpCode(headSoftClippedSeqLen, "S");
			targetJumpCodeVec.push_back(tmpHeadSoftClipJumpCode);
		}
		for(int tmpIndex = jumpCodeIndex_start; tmpIndex <= jumpCodeIndex_end; tmpIndex++)
		{
			targetJumpCodeVec.push_back(jumpCodeVec[tmpIndex]);
		}
		if(tailSoftClippedSeqLen > 0)
		{
			Jump_Code tmpTailSoftClipJumpCode(tailSoftClippedSeqLen, "S");
			targetJumpCodeVec.push_back(tmpTailSoftClipJumpCode);
		}
	}

	string jumpCodeVec2cigarString_withInterval(int jumpCodeIndex_start, int jumpCodeIndex_end)
	{
		vector<Jump_Code> targetJumpCodeVec;
		this->generateNewJumpCodeVecWithInterval(targetJumpCodeVec, 
			jumpCodeIndex_start, jumpCodeIndex_end);
		string tmpStr = this->jumpCodeVec2cigarString(targetJumpCodeVec);
		return tmpStr;
	}

	string jumpCodeVec2cigarString(vector<Jump_Code>& oriJumpCodeVec)
	{
		string tmpCigarString;
		//cout << "********" << "jumpCodeVecSize: " << jumpCodeVec.size() << endl;
		for(int tmp = 0; tmp < oriJumpCodeVec.size(); tmp++)
		{
			tmpCigarString += oriJumpCodeVec[tmp].toString();
		}
		return tmpCigarString;
	}

	int returnFlag()
	{
		return flag;
	}

	string returnReadName()
	{
		return readName;
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	int returnChrMapPos()
	{
		return chrMapPos;
	}

	int returnChrMapPos_end()
	{
		return chrMapPos_end;
	}

	void initiateSingleAlignmentInfo(
		string& samStr, Index_Info* indexInfo)
	{
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = samStr.find("\t", startLoc);
			string tmpSamField = samStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}


		readName = samFieldVec[0];
		string flagStr = samFieldVec[1];
		flag = atoi(flagStr.c_str());
		string chrNameStr = samFieldVec[2];
		chrNameInt = indexInfo->convertStringToInt(chrNameStr);
		string mapChrPosStr = samFieldVec[3];
		chrMapPos = atoi(mapChrPosStr.c_str());
		string cigarString = samFieldVec[5];
		this->cigarString2jumpCodeVec(cigarString, jumpCodeVec);
		chrMapPos_end = this->getEndPosOfSpecificJumpCode(chrMapPos, jumpCodeVec, jumpCodeVec.size()-1);
		// rnext -- [6]
		// pnext -- [7]
		// tlen -- [8]
		readSeq = samFieldVec[9];
		qualSeq = samFieldVec[10];

		//tagStr = samStr.substr(startLoc);
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
	{
		int tmpJumpCodeLength;
		string tmpJumpCodeType;

		int jumpCodeStartPosInCigarStr = 0;
		int jumpCodeEndPosInCigarStr;
		
		string candidateJumpCodeType = "SMNID";
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

	void generateSpliceSitePosVec(vector< pair<int,int> >& spliceSitePosPairVec)
	{
		for(int tmp = 0; tmp < jumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = jumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				int tmpDonerEndPos = getEndPosOfSpecificJumpCode(chrMapPos, jumpCodeVec, lastJumpCodeIndex);
				int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(chrMapPos, jumpCodeVec, currentJumpCodeIndex) + 1;
				//cout << "tmpDonerEndPos: " << tmpDonerEndPos << endl;
				//cout << "tmpAcceptorStartPos: " << tmpAcceptorStartPos << endl;
				//cout << "tmpSJindexInJumpCodeVec: " << tmp << endl;
				spliceSitePosPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				//SJindexVec_cigarStringJumpCodeVec.push_back(tmp);

			}
		}		
	}	

	void generateSJindexVecInJumpCodeVec(int mapChrPos, 
		vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector<int>& SJindexVec_cigarStringJumpCodeVec)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				SJindexVec_cigarStringJumpCodeVec.push_back(tmp);
			}		
		}
	}

	void generateExonIndexVecInJumpCodeVec(
		vector<int>& SJindexVecInJumpCodeVec,
		vector< pair<int,int> >& exonIndexVecInJumpCodeVec)
	{
		int jumpCodeVecSize = jumpCodeVec.size();
		int SJindexVecSizeInJumpCodeVec = SJindexVecInJumpCodeVec.size();
		//cout << "jumpCodeVecSize: " << jumpCodeVecSize << endl;
		//cout << "SJindexVecSizeInJumpCodeVec: " << SJindexVecSizeInJumpCodeVec << endl;
		if(SJindexVecSizeInJumpCodeVec == 0)
		{
			exonIndexVecInJumpCodeVec.push_back(pair<int,int>(0, jumpCodeVecSize-1));
		}
		else if(SJindexVecSizeInJumpCodeVec == 1)
		{
			int singleSJindexInJumpCodeVec = SJindexVecInJumpCodeVec[0];
			exonIndexVecInJumpCodeVec.push_back(
				pair<int,int>(0,singleSJindexInJumpCodeVec-1));
			exonIndexVecInJumpCodeVec.push_back(
				pair<int,int>(singleSJindexInJumpCodeVec+1, jumpCodeVecSize-1));	
		}
		else
		{
			int firstSJindexInJumpCodeVec = SJindexVecInJumpCodeVec[0];
			int lastSJindexInJumpCodeVec = SJindexVecInJumpCodeVec[SJindexVecSizeInJumpCodeVec-1];
			//first exon JumpCodeIndexPair
			exonIndexVecInJumpCodeVec.push_back(
				pair<int,int>(0,firstSJindexInJumpCodeVec-1));
			//other exons JumpCodeIndexPair
			for(int tmp = 0; tmp < SJindexVecSizeInJumpCodeVec-1; tmp++)
			{
				int tmpLastSJindexBefore = SJindexVecInJumpCodeVec[tmp];
				int tmpFirstSJindexBehind = SJindexVecInJumpCodeVec[tmp+1];
				exonIndexVecInJumpCodeVec.push_back(
					pair<int,int>(tmpLastSJindexBefore+1, tmpFirstSJindexBehind-1));				
			}
			//last exon JumpCodeIndexPair
			exonIndexVecInJumpCodeVec.push_back(
				pair<int,int>(lastSJindexInJumpCodeVec+1, jumpCodeVecSize-1));
		}
	}

	int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
	{
		if(jumpCodeIndex < 0)
			return startPos-1;
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

	void generateExonAndSpliceSiteInfoPair(
		vector< pair<int,int> >& exonIndexVec_jumpCodeVec,
		vector<int>& SJindexVec_jumpCodeVec,
		vector< pair<int,int> >& exonLocInRead, // vector< pair<startLocInRead, endLocInRead> > 
		vector< pair<int,int> >& exonPosInChr, // vector < pair<startPosInChr, endPosInChr> >
		vector< pair<int,int> >& spliceSiteLocInRead, // vector < pair<donerSiteInRead, acceptorSiteInRead> >
		vector< pair<int,int> >& spliceSitePosInChr) // vector < pair<donerSiteInChr, acceptorSiteInChr > >
	{
		// exon info pair
		for(int tmp = 0; tmp < exonIndexVec_jumpCodeVec.size(); tmp++)
		{
			int tmpExonFirstJumpCodeIndex = exonIndexVec_jumpCodeVec[tmp].first;
			int tmpExonLastJumpCodeIndex = exonIndexVec_jumpCodeVec[tmp].second;
			//cout << "tmpExonFirstJumpCodeIndex: " << tmpExonFirstJumpCodeIndex << endl;
			//cout << "tmpExonLastJumpCodeIndex: " << tmpExonLastJumpCodeIndex << endl;
			int tmpExonStartLocInRead = this->getEndLocInReadOfSpecificJumpCode(
				jumpCodeVec, tmpExonFirstJumpCodeIndex-1)+1;
			int tmpExonEndLocInRead = this->getEndLocInReadOfSpecificJumpCode(
				jumpCodeVec, tmpExonLastJumpCodeIndex);
			exonLocInRead.push_back(pair<int,int>(tmpExonStartLocInRead,
				tmpExonEndLocInRead));
			int tmpExonStartMapPosInChr = this->getEndPosOfSpecificJumpCode(
				chrMapPos, jumpCodeVec, tmpExonFirstJumpCodeIndex-1)+1;
			int tmpExonEndMapPosInChr = this->getEndPosOfSpecificJumpCode(
				chrMapPos, jumpCodeVec, tmpExonLastJumpCodeIndex);
			//cout << "tmpExonStartPosInChr: " << tmpExonStartMapPosInChr << endl;
			//cout << "tmpExonEndPosInChr: " << tmpExonEndMapPosInChr << endl;
			exonPosInChr.push_back(pair<int,int>(tmpExonStartMapPosInChr,
				tmpExonEndMapPosInChr));
		}
		//SJ info pair
		for(int tmp = 0; tmp < SJindexVec_jumpCodeVec.size(); tmp++)
		{
			int tmpSJindex = SJindexVec_jumpCodeVec[tmp];
			int tmpSJdonerSiteInRead = this->getEndLocInReadOfSpecificJumpCode(
				jumpCodeVec, tmpSJindex-1);
			int tmpSJacceptorSiteInRead = tmpSJdonerSiteInRead + 1;
			spliceSiteLocInRead.push_back(pair<int,int>(tmpSJdonerSiteInRead, tmpSJacceptorSiteInRead));
			int tmpSJdonerSiteInChr = this->getEndPosOfSpecificJumpCode(
				chrMapPos, jumpCodeVec, tmpSJindex-1);
			int tmpSJacceptorSiteInChr = this->getEndPosOfSpecificJumpCode(
				chrMapPos, jumpCodeVec, tmpSJindex)+1;
			spliceSitePosInChr.push_back(pair<int,int>(tmpSJdonerSiteInChr, tmpSJacceptorSiteInChr));
		}
	}

	void generateExonAndSpliceSiteVec(
		vector< pair<int,int> > tmpExonLocInRead,
		vector< pair<int,int> > tmpExonPosInChr,
		vector< pair<int,int> > tmpSpliceSiteLocInRead,
		vector< pair<int,int> > tmpSpliceSitePosInChr
		)
	{
		vector<int> tmpSJindexVec_jumpCodeVec;
		vector< pair<int,int> > tmpExonIndexVec_jumpCodeVec;
		this->generateSJindexVecInJumpCodeVec(
			chrMapPos, jumpCodeVec, tmpSJindexVec_jumpCodeVec);
		this->generateExonIndexVecInJumpCodeVec(
			tmpSJindexVec_jumpCodeVec, tmpExonIndexVec_jumpCodeVec);
		this->generateExonAndSpliceSiteInfoPair(
			tmpExonIndexVec_jumpCodeVec, tmpSJindexVec_jumpCodeVec, 
			tmpExonLocInRead, tmpExonPosInChr, 
			tmpSpliceSiteLocInRead, tmpSpliceSitePosInChr);					
	}

	bool alignmentWithInvalidSJ_bool(AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{
		vector< pair<int,int> > tmpSpliceSitePosPairVec;	
		this->generateSpliceSitePosVec(tmpSpliceSitePosPairVec);
		int tmpSpliceSitePosPairVec_size = tmpSpliceSitePosPairVec.size();
		//cout << "tmpSpliceSitePosPairVec_size: " << tmpSpliceSitePosPairVec_size << endl;
		for(int tmp = 0; tmp < tmpSpliceSitePosPairVec_size; tmp++)
		{
			int tmpSJpos_doner = tmpSpliceSitePosPairVec[tmp].first;
			int tmpSJpos_acceptor = tmpSpliceSitePosPairVec[tmp].second;
			//cout << endl << "tmp: " << tmp << endl;
			//cout << "tmpSJpos_doner: " << tmpSJpos_doner << endl;
			//cout << "tmpSJpos_acceptor: " << tmpSJpos_acceptor << endl;
			int tmpSearch_returnIndex 
				= alignInferJunctionHashInfo->searchAndReturnAlignInferInfoVecIndex(
					chrNameInt, tmpSJpos_doner, tmpSJpos_acceptor);
			//cout << "tmpSearch_returnIndex: " << tmpSearch_returnIndex << endl;
			if(tmpSearch_returnIndex >= 0)
			{
				bool tmpSJinvalidBool = alignInferJunctionHashInfo->SJinvalid(tmpSearch_returnIndex);
				//cout << "tmpSJinvalidBool: " << tmpSJinvalidBool << endl;
				if(tmpSJinvalidBool)
					return true;
				//else
				//	return false;
			}
		}
		return false;
	}

	void refineInvalidSJinAlignment_SE(
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		Index_Info* indexInfo,
		int& leftJumpCodeIndexPairAfterClipping_result,
		pair<int,int>& leftJumpCodeIndexPairAfterClipping)
	{
		//cout << "refineInvalidSJinAlignment_SE starts for SingleAlignment_Info ......" << endl;
		vector<int> tmpSJindexVec_jumpCodeVec;
		vector< pair<int,int> > tmpSJlocInReadVec;
		vector< pair<int,int> > tmpSJposInChrVec;	
		vector< pair<int,int> > tmpExonIndexVec_jumpCodeVec;
		vector< pair<int,int> > tmpExonLocInReadVec;
		vector< pair<int,int> > tmpExonPosInChrVec;		
		this->generateSJindexVecInJumpCodeVec(
			chrMapPos, jumpCodeVec, tmpSJindexVec_jumpCodeVec);
		this->generateExonIndexVecInJumpCodeVec(
			tmpSJindexVec_jumpCodeVec, tmpExonIndexVec_jumpCodeVec);
		this->generateExonAndSpliceSiteInfoPair(
			tmpExonIndexVec_jumpCodeVec, tmpSJindexVec_jumpCodeVec, 
			tmpExonLocInReadVec, tmpExonPosInChrVec, tmpSJlocInReadVec, tmpSJposInChrVec);		
		//cout << "end of generating exon and SJ info pair ......" << endl;
		vector< int > candiAlterPathSJindexVecInJumpCodeVec_doner;
		vector< vector<Jump_Code> > candiAlterPathJumpCodeVecVec_doner;
		vector< vector< int > > candiAlterPathMismatchPosVecVec_doner;
		vector< vector< char > > candiAlterPathMismatchCharVecVec_doner;

		vector< int > candiAlterPathSJindexVecInJumpCodeVec_acceptor;
		vector< vector<Jump_Code> > candiAlterPathJumpCodeVecVec_acceptor;
		vector< vector< int > > candiAlterPathMismatchPosVecVec_acceptor;
		vector< vector< char > > candiAlterPathMismatchCharVecVec_acceptor;
		//cout << "generateCandiAlterPath starts ......" << endl;
		this->generateCandiAlterPath(
			tmpSJindexVec_jumpCodeVec,
			tmpSJposInChrVec,

			candiAlterPathSJindexVecInJumpCodeVec_doner, // output generated
			candiAlterPathJumpCodeVecVec_doner,  // output generated
			candiAlterPathMismatchPosVecVec_doner, // output generated
			candiAlterPathMismatchCharVecVec_doner, // output generated

			candiAlterPathSJindexVecInJumpCodeVec_acceptor, // output generated
			candiAlterPathJumpCodeVecVec_acceptor,// output generated
			candiAlterPathMismatchPosVecVec_acceptor,// output generated
			candiAlterPathMismatchCharVecVec_acceptor,// output generated			
			
			alignInferJunctionHashInfo,
			indexInfo);

		//pair<int,int> leftJumpCodeIndexPairAfterClipping
		// get clipped alignment ....
		leftJumpCodeIndexPairAfterClipping_result 
			= this->returnLeftJumpCodeIndexPairAfterClipping(
				candiAlterPathSJindexVecInJumpCodeVec_doner,
				candiAlterPathSJindexVecInJumpCodeVec_acceptor,
				leftJumpCodeIndexPairAfterClipping);
	}

	int returnLeftJumpCodeIndexPairAfterClipping(
		vector<int>& candiAlterPathSJindexVecInJumpCodeVec_doner,
		vector<int>& candiAlterPathSJindexVecInJumpCodeVec_acceptor,
		pair<int,int>& leftJumpCodeIndexPairAfterClipping) 
	// return -1 -- no need to do refine alignment, 
	//	(candiAlterPathSJindexVecInJumpCodeVec_doner.size() 
	//	 + candiAlterPathSJindexVecInJumpCodeVec_acceptor.size() == 0)
	// return 0 -- unmap, whole sequence can be clipped
	//	max(candiAlterPathSJindexVecInJumpCodeVec_doner) 
	//	 >= min(candiAlterPathSJindexVecInJumpCodeVec_acceptor)
	// return 1 -- parts of read sequence can be clipped, and some parts can be left
	//	max(candiAlterPathSJindexVecInJumpCodeVec_doner) 
	//	 < min(candiAlterPathSJindexVecInJumpCodeVec_acceptor)
	{
		int candiAlterPathSJindexVecInJumpCodeVec_doner_size
			= candiAlterPathSJindexVecInJumpCodeVec_doner.size();
		int candiAlterPathSJindexVecInJumpCodeVec_acceptor_size
			= candiAlterPathSJindexVecInJumpCodeVec_acceptor.size();		
		int candiAlterPathSJindexInJumpCodeVec_doner_max = -1;
		int candiAlterPathSJindexInJumpCodeVec_acceptor_min = 999;
		if(candiAlterPathSJindexVecInJumpCodeVec_doner_size
			+ candiAlterPathSJindexVecInJumpCodeVec_acceptor_size == 0)
		{
			return -1;
		}
		else if(candiAlterPathSJindexVecInJumpCodeVec_doner_size > 0)
		{
			
			for(int tmp = 0; tmp < candiAlterPathSJindexVecInJumpCodeVec_doner_size; tmp++)
			{
				if(candiAlterPathSJindexVecInJumpCodeVec_doner[tmp] > candiAlterPathSJindexInJumpCodeVec_doner_max)
					candiAlterPathSJindexInJumpCodeVec_doner_max = candiAlterPathSJindexVecInJumpCodeVec_doner[tmp];
			}
		}
		else if(candiAlterPathSJindexVecInJumpCodeVec_acceptor_size > 0)
		{
			
			for(int tmp = 0; tmp < candiAlterPathSJindexVecInJumpCodeVec_acceptor_size; tmp++)
			{
				if(candiAlterPathSJindexVecInJumpCodeVec_acceptor[tmp] < candiAlterPathSJindexInJumpCodeVec_acceptor_min)
					candiAlterPathSJindexInJumpCodeVec_acceptor_min = candiAlterPathSJindexVecInJumpCodeVec_acceptor[tmp];
			}
		}
		else
		{}
		int leftJumpCodeIndex_start = candiAlterPathSJindexInJumpCodeVec_doner_max + 1;
		if(candiAlterPathSJindexInJumpCodeVec_doner_max == -1)
			leftJumpCodeIndex_start = 0;
		int leftJumpCodeIndex_end = candiAlterPathSJindexInJumpCodeVec_acceptor_min - 1;
		if(candiAlterPathSJindexInJumpCodeVec_acceptor_min == 999)
			leftJumpCodeIndex_end = jumpCodeVec.size() - 1;
		if(leftJumpCodeIndex_start > leftJumpCodeIndex_end)
			return 0;
		else
		{
			leftJumpCodeIndexPairAfterClipping.first = leftJumpCodeIndex_start;
			leftJumpCodeIndexPairAfterClipping.second = leftJumpCodeIndex_end;
			return 1;
		}
	}

	void generateCandiAlterPath(
		vector<int>& tmpSJindexVec_jumpCodeVec, // input
		vector< pair<int,int> >& tmpSJposInChrVec, // input 		 
		// doner site alter path
		vector< int >& candiAlterPathSJindexVecInJumpCodeVec_doner, // output generated
		vector< vector<Jump_Code> >& candiAlterPathJumpCodeVecVec_doner,// output generated
		vector< vector< int > >& candiAlterPathMismatchPosVecVec_doner,// output generated
		vector< vector< char > >& candiAlterPathMismatchCharVecVec_doner,// output generated
		// acceptor site alter path
		vector< int >& candiAlterPathSJindexVecInJumpCodeVec_acceptor, // output generated
		vector< vector<Jump_Code> >& candiAlterPathJumpCodeVecVec_acceptor,// output generated
		vector< vector< int > >& candiAlterPathMismatchPosVecVec_acceptor,// output generated
		vector< vector< char > >& candiAlterPathMismatchCharVecVec_acceptor,// output generated

		AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		Index_Info* indexInfo
		)
	{
		for(int tmp = 0; tmp < tmpSJposInChrVec.size(); tmp++)
		{
			int tmpSJindex_inJumpCodeVec = tmpSJindexVec_jumpCodeVec[tmp];
			int tmpSJdonerEndPos = tmpSJposInChrVec[tmp].first;
			int tmpSJacceptorStartPos = tmpSJposInChrVec[tmp].second;
			int tmpSJindex_inAlignInferHash 
				= alignInferJunctionHashInfo->searchAndReturnAlignInferInfoVecIndex(
					chrNameInt, tmpSJdonerEndPos, tmpSJacceptorStartPos);
			// cout << "tmp in tmpSJposInChrVec: " << tmp << endl;	
			// cout << "tmpSJindex_inJumpCodeVec: " << tmpSJindex_inJumpCodeVec << endl;
			// cout << "chrNameInt: " << chrNameInt << endl;
			// cout << "tmpSJdonerEndPos: " << tmpSJdonerEndPos << endl;	
			// cout << "tmpSJacceptorStartPos: " << tmpSJacceptorStartPos << endl;	
			// cout << "tmpSJindex_inAlignInferHash: " << tmpSJindex_inAlignInferHash << endl;
			if(tmpSJindex_inAlignInferHash >= 0)
			{
					bool SJcanBeExtendedAtAcceptorSpliceSite_bool 
						= alignInferJunctionHashInfo->SJcanBeExtendedAtAcceptorSpliceSite(tmpSJindex_inAlignInferHash);
					//cout << "SJcanBeExtendedAtAcceptorSpliceSite_bool: " << SJcanBeExtendedAtAcceptorSpliceSite_bool << endl; 
					if(SJcanBeExtendedAtAcceptorSpliceSite_bool)
					{
						// vector<Jump_Code> tmpCandiAlterPathJumpCodeVec_extension_doner;
						// vector< int > tmpCandiAlterPathMismatchPosVec_extension_doner;
						// vector< char > tmpCandiAlterPathMismatchCharVec_extension_doner;	
						// //cout << "extendAtAcceptorSpliceSite starts ......" << endl;					
						// bool extend_success_bool 
						// 	= this->extendAtAcceptorSpliceSite(
						// 		tmpSJindex_inJumpCodeVec, tmpSJacceptorStartPos, indexInfo,
						// 		tmpCandiAlterPathJumpCodeVec_extension_doner,
						// 		tmpCandiAlterPathMismatchPosVec_extension_doner,
						// 		tmpCandiAlterPathMismatchCharVec_extension_doner);
						// //cout << "extend_success_bool : " << extend_success_bool << endl;
						// if(extend_success_bool)
						// {
							candiAlterPathSJindexVecInJumpCodeVec_doner.push_back(
								tmpSJindex_inJumpCodeVec);
						// 	candiAlterPathJumpCodeVecVec_doner.push_back(
						// 		tmpCandiAlterPathJumpCodeVec_extension_doner);
						// 	//cout << "tmpCandiAlterPathJumpCodeVec_extension_doner: " 
						// 	//	<< this->jumpCodeVec2cigarString(tmpCandiAlterPathJumpCodeVec_extension_doner) << endl;
						// 	candiAlterPathMismatchPosVecVec_doner.push_back(
						// 		tmpCandiAlterPathMismatchPosVec_extension_doner);
						// 	candiAlterPathMismatchCharVecVec_doner.push_back(
						// 		tmpCandiAlterPathMismatchCharVec_extension_doner);
						// }
					}
					bool SJcanBeExtendedAtDonerSpliceSite_bool
						= alignInferJunctionHashInfo->SJcanBeExtendedAtDonerSpliceSite(tmpSJindex_inAlignInferHash);
					//cout << "SJcanBeExtendedAtDonerSpliceSite_bool: " << SJcanBeExtendedAtDonerSpliceSite_bool << endl;
					if(SJcanBeExtendedAtDonerSpliceSite_bool)
					{
						// vector<Jump_Code> tmpCandiAlterPathJumpCodeVec_extension_acceptor;
						// vector< int > tmpCandiAlterPathMismatchPosVec_extension_acceptor;
						// vector< char > tmpCandiAlterPathMismatchCharVec_extension_acceptor;	
						// //cout << "extendAtDonerSpliceSite starts ...... " << endl;				
						// bool extend_success_bool 
						// 	= this->extendAtDonerSpliceSite(
						// 		tmpSJindex_inJumpCodeVec, tmpSJdonerEndPos, indexInfo,
						// 		tmpCandiAlterPathJumpCodeVec_extension_acceptor,
						// 		tmpCandiAlterPathMismatchPosVec_extension_acceptor,
						// 		tmpCandiAlterPathMismatchCharVec_extension_acceptor);
						// //cout << "extend_success_bool: " << extend_success_bool << endl;
						// if(extend_success_bool)
						// {
							candiAlterPathSJindexVecInJumpCodeVec_acceptor.push_back(
								tmpSJindex_inJumpCodeVec);
						// 	candiAlterPathJumpCodeVecVec_acceptor.push_back(
						// 		tmpCandiAlterPathJumpCodeVec_extension_acceptor);
						// 	//cout << "tmpCandiAlterPathJumpCodeVec_extension_acceptor: " 
						// 	//	<< this->jumpCodeVec2cigarString(tmpCandiAlterPathJumpCodeVec_extension_acceptor) << endl;
						// 	candiAlterPathMismatchPosVecVec_acceptor.push_back(
						// 		tmpCandiAlterPathMismatchPosVec_extension_acceptor);
						// 	candiAlterPathMismatchCharVecVec_acceptor.push_back(
						// 		tmpCandiAlterPathMismatchCharVec_extension_acceptor);
						// }
					}
					bool exactTheSame2someAlterSpliceSite_doner_bool
						= alignInferJunctionHashInfo->exactTheSame2someAlterSpliceSite_doner(tmpSJindex_inAlignInferHash);
					if(exactTheSame2someAlterSpliceSite_doner_bool)
					{
						candiAlterPathSJindexVecInJumpCodeVec_doner.push_back(
							tmpSJindex_inJumpCodeVec);						
					}
					bool exactTheSame2someAlterSpliceSite_acceptor_bool
						= alignInferJunctionHashInfo->exactTheSame2someAlterSpliceSite_acceptor(tmpSJindex_inAlignInferHash);
					if(exactTheSame2someAlterSpliceSite_acceptor_bool)
					{
						candiAlterPathSJindexVecInJumpCodeVec_acceptor.push_back(
							tmpSJindex_inJumpCodeVec);						
					}
			}
			else
			{}
		}
	}

 	bool extendAtAcceptorSpliceSite(
 		int tmpCandiAlterPathSJindexInJumpCodeVec_extension_doner,
 		int acceptorSpliceSite, Index_Info* indexInfo,
		vector<Jump_Code>& tmpCandiAlterPathJumpCodeVec_extension_doner,
		vector< int >& tmpCandiAlterPathMismatchPosVec_extension_doner,
		vector< char >& tmpCandiAlterPathMismatchCharVec_extension_doner)
 	{
 		//cout << "extendAtAcceptorSpliceSite starts ......" << endl;
 		int readSubSeq_fakeSJdonerSite_length 
 			= this->getEndLocInReadOfSpecificJumpCode(
				jumpCodeVec, tmpCandiAlterPathSJindexInJumpCodeVec_extension_doner-1);
 		string readSubSeq_fakeSJdonerSite
 			= readSeq.substr(0, readSubSeq_fakeSJdonerSite_length);
 		string chromSubSeq_extensionAtAcceptorSpliceSite
 			= indexInfo->returnChromStrSubstr(chrNameInt, 
 				acceptorSpliceSite - readSubSeq_fakeSJdonerSite_length,
 				readSubSeq_fakeSJdonerSite_length);
 		int max_mismatch = 	readSubSeq_fakeSJdonerSite_length/8 + 1;
 		FixDoubleAnchor_Match_Info fixMatchInfo; //= FixDoubleAnchor_Match_Info();
 		bool scoreStringBool = fixMatchInfo.fixMatch(
 			readSubSeq_fakeSJdonerSite, chromSubSeq_extensionAtAcceptorSpliceSite,
 			max_mismatch, acceptorSpliceSite - readSubSeq_fakeSJdonerSite_length);
 		if(scoreStringBool)
 		{
 			Jump_Code tmpMatchJumpCode(readSubSeq_fakeSJdonerSite_length, "M");
 			tmpCandiAlterPathJumpCodeVec_extension_doner.push_back(tmpMatchJumpCode);
 			fixMatchInfo.copyMismatchPos2TargetVec(
 				tmpCandiAlterPathMismatchPosVec_extension_doner);
 			fixMatchInfo.copyMismatchChar2TargetVec(
 				tmpCandiAlterPathMismatchCharVec_extension_doner);
 		}
 		return scoreStringBool;
 	}

 	bool extendAtDonerSpliceSite(
 		int tmpCandiAlterPathSJindexInJumpCodeVec_extension_acceptor,
 		int donerSpliceSite, Index_Info* indexInfo,
		vector<Jump_Code>& tmpCandiAlterPathJumpCodeVec_extension_acceptor,
		vector< int >& tmpCandiAlterPathMismatchPosVec_extension_acceptor,
		vector< char >& tmpCandiAlterPathMismatchCharVec_extension_acceptor)
 	{
 		//cout << "extendAtDoenrSpliceSite starts ......" << endl;
 		int tmpReadLength = readSeq.length();
 		int readSubSeq_fakeSJacceptorSide_startLocInRead
 			= this->getEndLocInReadOfSpecificJumpCode(
 				jumpCodeVec, tmpCandiAlterPathSJindexInJumpCodeVec_extension_acceptor)+1;
 		int readSubSeq_fakeSJacceptorSide_length 
 			= tmpReadLength - readSubSeq_fakeSJacceptorSide_startLocInRead + 1;
 		string readSubSeq_fakeSJacceptorSide
 			= readSeq.substr(readSubSeq_fakeSJacceptorSide_startLocInRead-1, 
 				readSubSeq_fakeSJacceptorSide_length);
 		string chromSubSeq_extensionAtDonerSpliceSite
 			= indexInfo->returnChromStrSubstr(chrNameInt, 
 				donerSpliceSite + 1, readSubSeq_fakeSJacceptorSide_length);
 		int max_mismatch = 	readSubSeq_fakeSJacceptorSide_length/8 + 1;
 		FixDoubleAnchor_Match_Info fixMatchInfo; //= FixDoubleAnchor_Match_Info();
 		bool scoreStringBool = fixMatchInfo.fixMatch(
 			readSubSeq_fakeSJacceptorSide, chromSubSeq_extensionAtDonerSpliceSite,
 			max_mismatch, donerSpliceSite + 1);
 		if(scoreStringBool)
 		{
 			Jump_Code tmpMatchJumpCode(readSubSeq_fakeSJacceptorSide_length, "M");
 			tmpCandiAlterPathJumpCodeVec_extension_acceptor.push_back(tmpMatchJumpCode);
 			fixMatchInfo.copyMismatchPos2TargetVec(
 				tmpCandiAlterPathMismatchPosVec_extension_acceptor);
 			fixMatchInfo.copyMismatchChar2TargetVec(
 				tmpCandiAlterPathMismatchCharVec_extension_acceptor);
 		}
 		return scoreStringBool;
 	}

};

class PairAlignment_Info
{
private:
	SingleAlignment_Info singleAlignmentInfo_1;
	SingleAlignment_Info singleAlignmentInfo_2;
public:
	PairAlignment_Info()
	{

	}

	bool BeersGroundTruth_End1BeforeEnd2_bool()
	{
		int flagInt_Beers_1 = singleAlignmentInfo_1.returnFlag();
		int flagInt_Beers_2 = singleAlignmentInfo_2.returnFlag();
		if((flagInt_Beers_1 == 0)&&(flagInt_Beers_2 == 16))
			return true;
		else if((flagInt_Beers_1 == 16)&&(flagInt_Beers_2 == 0))
			return false;
		else
		{
			cout << "flag error ...." << endl;
			exit(1);
		}
	}

	bool readNamesPaired_bool()
	{
		string readName_1 = singleAlignmentInfo_1.returnReadName();
		string readName_2 = singleAlignmentInfo_2.returnReadName();
		if(readName_1 == readName_2)
		{
			return true;
		}
		else
		{
			int readName_1_length = readName_1.size();
			int readName_2_length = readName_2.size();
			if((readName_1.at(readName_1_length-2) == '/') && (readName_2.at(readName_2_length-2) == '/'))
			{
				if(readName_1.substr(0,readName_1_length-2) == readName_2.substr(0,readName_2_length-2))
					return true;
				else
					return false;
			}
			else
				return false;
		}
	}

	string returnReadName_1()
	{
		return singleAlignmentInfo_1.returnReadName();
	}

	string returnReadName_2()
	{
		return singleAlignmentInfo_2.returnReadName();
	}

	bool bothPairEndReadsWholeMatch_bool()
	{
		if((singleAlignmentInfo_1.wholeMatch_bool())&&(singleAlignmentInfo_2.wholeMatch_bool()))
			return true;
		else
			return false;
	}

	bool bothPairEndReadsWithoutSJ_bool()
	{
		if((singleAlignmentInfo_1.withoutSJ_bool())&&(singleAlignmentInfo_2.withoutSJ_bool()))
			return true;
		else
			return false;
	}


	string returnOriReadSeq_1()
	{
		return singleAlignmentInfo_1.returnReadSeq();
	}

	string returnOriReadSeq_2()
	{
		return singleAlignmentInfo_2.returnRcmReadSeq();
	}

	string returnRcmReadSeq_1()
	{
		return singleAlignmentInfo_1.returnRcmReadSeq();
	}

	string returnRcmReadSeq_2()
	{
		return singleAlignmentInfo_2.returnReadSeq();
	}

	string returnSingleAlignmentInfo_unmap(Index_Info* indexInfo, bool End1orEnd2_bool)
	{
		if(End1orEnd2_bool)
		{
			return singleAlignmentInfo_1.returnUnmappedAlignment(indexInfo);
		}
		else
		{
			return singleAlignmentInfo_2.returnUnmappedAlignment(indexInfo);
		}
	}

	string returnSingleAlignemntInfo_clippedWithInterval(
		Index_Info* indexInfo, int index_start, int index_end, bool End1orEnd2_bool)
	{
		if(End1orEnd2_bool)
		{	
			return singleAlignmentInfo_1.returnClippedAlignmentWithInterval(indexInfo, index_start, index_end);
		}
		else
		{
			return singleAlignmentInfo_2.returnClippedAlignmentWithInterval(indexInfo, index_start, index_end);
		}
	}

	void initiatePairAlignmentInfo(
		string& tmpAlign_1, string& tmpAlign_2, 
		Index_Info* indexInfo)
	{
		singleAlignmentInfo_1.initiateSingleAlignmentInfo(
			tmpAlign_1, indexInfo);
		singleAlignmentInfo_2.initiateSingleAlignmentInfo(
			tmpAlign_2, indexInfo);
	}

	// void initiatePairAlignmentInfo(char* tmpAlignChar_1, char* tmpAlignChar_2, Index_Info* indexInfo)
	// {
	// 	singleAlignmentInfo_1.initiateSingleAlignmentInfo(
	// 		tmpAlignChar_1, indexInfo);
	// 	singleAlignmentInfo_2.initiateSingleAlignmentInfo(
	// 		tmpAlignChar_2, indexInfo);
	// }

	void generateOverlapRegionExonPosVec(
		int overlapRegion_startPos, int overlapRegion_endPos,
		vector< pair<int,int> >& oriExonPosInChr, 
		vector< pair<int,int> >& overlapRegion_exonPosInChr)
	{
		for(int tmp = 0; tmp < oriExonPosInChr.size(); tmp++)
		{
			//cout << "tmpExonInt: " << tmp << endl;
			int tmpExonStartPosInChr = oriExonPosInChr[tmp].first;
			int tmpExonEndPosInChr = oriExonPosInChr[tmp].second;
			//cout << "tmpExonStartPosInChr: " << tmpExonStartPosInChr << endl;
			//cout << "tmpExonEndPosInChr: " << tmpExonEndPosInChr << endl;
			if( (tmpExonStartPosInChr > overlapRegion_endPos)
				||(tmpExonEndPosInChr < overlapRegion_startPos) )
			{}
			else if( (tmpExonStartPosInChr >= overlapRegion_startPos)
				&&(tmpExonEndPosInChr <= overlapRegion_endPos) )
			{
				overlapRegion_exonPosInChr.push_back(
					pair<int,int> (tmpExonStartPosInChr, tmpExonEndPosInChr));
			}
			else if( (tmpExonStartPosInChr <= overlapRegion_startPos)
				&&(tmpExonEndPosInChr >= overlapRegion_endPos) )
			{
				overlapRegion_exonPosInChr.push_back(
					pair<int,int> (overlapRegion_startPos, overlapRegion_endPos));
			}
			else if(tmpExonStartPosInChr >= overlapRegion_startPos)
				//&&(overlapRegion_startPos <= tmpExonEndPosInChr)
				//&&(tmpExonEndPosInChr))
			{
				overlapRegion_exonPosInChr.push_back(
					pair<int,int> (tmpExonStartPosInChr, overlapRegion_endPos));
			}
			else if( tmpExonStartPosInChr <= overlapRegion_startPos)
			{
				overlapRegion_exonPosInChr.push_back(
					pair<int,int> (overlapRegion_startPos, tmpExonEndPosInChr));
			}
			else
			{
				cout << "invalid exon in generateOverlapRegionExonPosVec ..." << endl;
				exit(1);
			}
		}
	}

	bool alignmentWithInvalidSJ_bool(AlignInferJunctionHash_Info* alignInferJunctionHashInfo)
	{
		bool alignmentWithInvalidSJBool_end1 
			= singleAlignmentInfo_1.alignmentWithInvalidSJ_bool(alignInferJunctionHashInfo);
		bool alignmentWithInvalidSJBool_end2
			= singleAlignmentInfo_2.alignmentWithInvalidSJ_bool(alignInferJunctionHashInfo);
		if(alignmentWithInvalidSJBool_end1 || alignmentWithInvalidSJBool_end2)
			return true;
		else
			return false;
	}

	void generateSpliceSitePposPairVec(
		vector< pair<int,int> >& SJposPairVec_end1,
		vector< pair<int,int> >& SJposPairVec_end2)
	{
		singleAlignmentInfo_1.generateSpliceSitePosVec(SJposPairVec_end1);
		singleAlignmentInfo_2.generateSpliceSitePosVec(SJposPairVec_end2);
	}

	bool getOverlapMapRegion(int& overlap_startMapPos, 
		int& overlap_endMapPos)
	{
		int chrNameInt_1 = singleAlignmentInfo_1.returnChrNameInt();
		int chrNameInt_2 = singleAlignmentInfo_2.returnChrNameInt();
		if(chrNameInt_1 != chrNameInt_2)
			return false;
		int mapPos_1 = singleAlignmentInfo_1.returnChrMapPos();
		int mapPos_2 = singleAlignmentInfo_2.returnChrMapPos();
		int mapPos_1_end = singleAlignmentInfo_1.returnChrMapPos_end();
		int mapPos_2_end = singleAlignmentInfo_2.returnChrMapPos_end();

		if((mapPos_1_end < mapPos_2)||(mapPos_2_end < mapPos_1))
			return false;
		else
		{
			if(mapPos_1 < mapPos_2)
				overlap_startMapPos = mapPos_2;
			else
				overlap_startMapPos = mapPos_1;

			if(mapPos_1_end < mapPos_2_end)
				overlap_endMapPos = mapPos_1_end;
			else
				overlap_endMapPos = mapPos_2_end;

			return true;
		}
	}

	void refineInvalidSJinAlignment_PE(AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		Index_Info* indexInfo,
		int& leftJumpCodeIndexPairAfterClipping_result_1,
		pair<int,int>& leftJumpCodeIndexPairAfterClipping_1,
		int& leftJumpCodeIndexPairAfterClipping_result_2,
		pair<int,int>& leftJumpCodeIndexPairAfterClipping_2
		)
	{
		//cout << "refineInvalidSJinAlignment_PE starts ......" << endl;
		singleAlignmentInfo_1.refineInvalidSJinAlignment_SE(alignInferJunctionHashInfo, indexInfo,
			leftJumpCodeIndexPairAfterClipping_result_1, leftJumpCodeIndexPairAfterClipping_1);
		singleAlignmentInfo_2.refineInvalidSJinAlignment_SE(alignInferJunctionHashInfo, indexInfo,
			leftJumpCodeIndexPairAfterClipping_result_2, leftJumpCodeIndexPairAfterClipping_2);
	}

	bool overlapExons_compatible_bool()
	{
		//cout << endl;
		//cout << "readName_1: " << singleAlignmentInfo_1.returnReadName() << endl;
		//cout << "readName_2: " << singleAlignmentInfo_2.returnReadName() << endl;
		int chrNameInt_1 = singleAlignmentInfo_1.returnChrNameInt();
		int chrNameInt_2 = singleAlignmentInfo_2.returnChrNameInt();
		if(chrNameInt_1 != chrNameInt_2)
			return false;
		int overlapRegion_startPos, overlapRegion_endPos;		
		bool overlapOrNot = this->getOverlapMapRegion(overlapRegion_startPos, 
			overlapRegion_endPos);
		if(!overlapOrNot)
			return true;
		vector< pair<int,int> > exonLocInRead_1;
		vector< pair<int,int> > exonPosInChr_1;
		vector< pair<int,int> > SJlocInRead_1;
		vector< pair<int,int> > SJposInChr_1;	
		vector< pair<int,int> > exonLocInRead_2;
		vector< pair<int,int> > exonPosInChr_2;
		vector< pair<int,int> > SJlocInRead_2;
		vector< pair<int,int> > SJposInChr_2;

		singleAlignmentInfo_1.generateExonAndSpliceSiteVec(
			exonLocInRead_1, exonPosInChr_1,
			SJlocInRead_1, SJposInChr_1);	
		singleAlignmentInfo_2.generateExonAndSpliceSiteVec(
			exonLocInRead_2, exonPosInChr_2,
			SJlocInRead_2, SJposInChr_2);

		vector< pair<int,int> > overlapRegion_exonPosInChr_1;
		vector< pair<int,int> > overlapRegion_exonPosInChr_2; 
		//cout << "overlapRegion_startPos: " << overlapRegion_startPos << endl;
		//cout << "overlapRegion_endPos: " << overlapRegion_endPos << endl; 
		this->generateOverlapRegionExonPosVec(
			overlapRegion_startPos, overlapRegion_endPos,
			exonPosInChr_1, overlapRegion_exonPosInChr_1);
		this->generateOverlapRegionExonPosVec(
			overlapRegion_startPos, overlapRegion_endPos,
			exonPosInChr_2, overlapRegion_exonPosInChr_2);

		if(overlapRegion_exonPosInChr_1.size() != 
			overlapRegion_exonPosInChr_2.size())
		{
			return false;
		}
		else
		{
			for(int tmp = 0; tmp < overlapRegion_exonPosInChr_1.size(); 
				tmp++)
			{
				//cout << "Exon_1: " << tmp+1 << " " << overlapRegion_exonPosInChr_1[tmp].first << " " << overlapRegion_exonPosInChr_1[tmp].second << endl;
				//cout << "Exon_2: " << tmp+1 << " " << overlapRegion_exonPosInChr_2[tmp].first << " " << overlapRegion_exonPosInChr_2[tmp].second << endl;
				if((overlapRegion_exonPosInChr_1[tmp].first != overlapRegion_exonPosInChr_2[tmp].first)
					||(overlapRegion_exonPosInChr_1[tmp].second != overlapRegion_exonPosInChr_2[tmp].second))
					return false;
			}
		}
		return true;
	}
};

class ReadPairMap_Info
{

};

class RefineAlignment_Info
{
private:

public:
	RefineAlignment_Info()
	{}

	void getPositiveStrandFastaFromMPS3sam_reissuedReadName(
		ifstream& sam_ifs,
		ofstream& fasta_end1_ofs,
		ofstream& fasta_end2_ofs,
		Index_Info* indexInfo)
	{
		int reissuedReadNameID = 0;
		while(1)
		{	
			if(sam_ifs.eof())
				break;
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			if(tmpAlign_1.at(0) == '@')
				continue;
			if(tmpAlign_1 == "")
				break;
			getline(sam_ifs, tmpAlign_2);
			//if(sam_ifs.eof())
			//	break;
			reissuedReadNameID ++;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);
			string readName_1 = "seq." + int_to_str(reissuedReadNameID) + "/1";
			string readName_2 = "seq." + int_to_str(reissuedReadNameID) + "/2";
			string readOriSeq_1 = pairAlignmentInfo.returnOriReadSeq_1();
			string readOriSeq_2 = pairAlignmentInfo.returnOriReadSeq_2();
			fasta_end1_ofs << readName_1 << endl << readOriSeq_1 << endl;
			fasta_end2_ofs << readName_2 << endl << readOriSeq_2 << endl;
		}		
	}

	void outputWholeMatchPairReadInFasta_groundTruth(
		ifstream& sam_ifs,
		ofstream& fasta_end1_ofs,
		ofstream& fasta_end2_ofs,
		Index_Info* indexInfo)
	{
		//int reissuedReadNameID = 0;
		while(1)
		{	
			if(sam_ifs.eof())
				break;
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			if(tmpAlign_1.at(0) == '@')
				continue;
			if(tmpAlign_1 == "")
				break;
			getline(sam_ifs, tmpAlign_2);
			//if(sam_ifs.eof())
			//	break;
			//reissuedReadNameID ++;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);
			bool tmpReadNamesPairedBool = pairAlignmentInfo.readNamesPaired_bool();
			bool tmpBothPairEndReadsWholeMatchBool = pairAlignmentInfo.bothPairEndReadsWholeMatch_bool();
			if(tmpReadNamesPairedBool && tmpBothPairEndReadsWholeMatchBool)
			{	
				string readName_1 = pairAlignmentInfo.returnReadName_1();//"seq." + int_to_str(reissuedReadNameID) + "/1";
				string readName_2 = pairAlignmentInfo.returnReadName_2();//"seq." + int_to_str(reissuedReadNameID) + "/2";
				string readOriSeq_1 = pairAlignmentInfo.returnOriReadSeq_1();
				string readOriSeq_2 = pairAlignmentInfo.returnRcmReadSeq_2();
				fasta_end1_ofs << readName_1 << endl << readOriSeq_1 << endl;
				fasta_end2_ofs << readName_2 << endl << readOriSeq_2 << endl;
			}
		}		
	}


	void outputWholeMatchPairReadInFasta_forwardDirection_groundTruth(
		ifstream& sam_ifs,
		ofstream& fasta_end1_ofs,
		ofstream& fasta_end2_ofs,
		Index_Info* indexInfo)
	{
		while(1)
		{	
			if(sam_ifs.eof())
				break;
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			if(tmpAlign_1.at(0) == '@')
				continue;
			if(tmpAlign_1 == "")
				break;
			getline(sam_ifs, tmpAlign_2);
			//if(sam_ifs.eof())
			//	break;
			//reissuedReadNameID ++;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);
			bool tmpReadNamesPairedBool = pairAlignmentInfo.readNamesPaired_bool();
			bool tmpBothPairEndReadsWholeMatchBool 
				= pairAlignmentInfo.bothPairEndReadsWholeMatch_bool();
			bool tmpBeersGroundTruthEnd1BeforeEnd2Bool 
				= pairAlignmentInfo.BeersGroundTruth_End1BeforeEnd2_bool();
			if(tmpReadNamesPairedBool && tmpBothPairEndReadsWholeMatchBool)
			{	
				string readName_1 = pairAlignmentInfo.returnReadName_1();//"seq." + int_to_str(reissuedReadNameID) + "/1";
				string readName_2 = pairAlignmentInfo.returnReadName_2();//"seq." + int_to_str(reissuedReadNameID) + "/2";
				string readOriSeq_1; 
				string readOriSeq_2;
				if(tmpBeersGroundTruthEnd1BeforeEnd2Bool)
				{
					readOriSeq_1 = pairAlignmentInfo.returnOriReadSeq_1();
					readOriSeq_2 = pairAlignmentInfo.returnRcmReadSeq_2();
				}
				else
				{
					readOriSeq_1 = pairAlignmentInfo.returnRcmReadSeq_2();
					readOriSeq_2 = pairAlignmentInfo.returnOriReadSeq_1();					
				}
				fasta_end1_ofs << readName_1 << endl << readOriSeq_1 << endl;
				fasta_end2_ofs << readName_2 << endl << readOriSeq_2 << endl;
			}
		}			
	}	

	void outputWithoutSJpairReadInFasta_groundTruth(
		ifstream& sam_ifs,
		ofstream& fasta_end1_ofs,
		ofstream& fasta_end2_ofs,
		Index_Info* indexInfo)
	{
		int reissuedReadNameID = 0;
		while(1)
		{	
			if(sam_ifs.eof())
				break;
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			if(tmpAlign_1.at(0) == '@')
				continue;
			if(tmpAlign_1 == "")
				break;
			getline(sam_ifs, tmpAlign_2);
			//if(sam_ifs.eof())
			//	break;
			//reissuedReadNameID ++;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);
			bool tmpReadNamesPairedBool = pairAlignmentInfo.readNamesPaired_bool();
			bool tmpBothPairEndReadsWitoutSJbool = pairAlignmentInfo.bothPairEndReadsWithoutSJ_bool();
			//bool tmpBeersGroundTruth_End1BeforeEnd2_bool = pairAlignmentInfo.BeersGroundTruth_End1BeforeEnd2();
			if(tmpReadNamesPairedBool && tmpBothPairEndReadsWitoutSJbool)
			{	
				string readName_1 = pairAlignmentInfo.returnReadName_1();//"seq." + int_to_str(reissuedReadNameID) + "/1";
				string readName_2 = pairAlignmentInfo.returnReadName_2();//"seq." + int_to_str(reissuedReadNameID) + "/2";
				string readOriSeq_1 = pairAlignmentInfo.returnOriReadSeq_1();
				string readOriSeq_2 = pairAlignmentInfo.returnRcmReadSeq_2();
				fasta_end1_ofs << readName_1 << endl << readOriSeq_1 << endl;
				fasta_end2_ofs << readName_2 << endl << readOriSeq_2 << endl;
			}
		}		
	}	

	void detectIncompatiblePairedAlignment(
		ifstream& sam_ifs,
		ofstream& compatibleSAM_ofs,
		ofstream& incompatibleSAM_ofs,
		Index_Info* indexInfo)
	{
		while(1)
		{	
			if(sam_ifs.eof())
				break;
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			if(tmpAlign_1.at(0) == '@')
				continue;
			if(tmpAlign_1 == "")
				break;
			getline(sam_ifs, tmpAlign_2);
			//if(sam_ifs.eof())
			//	break;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);

			bool pairedAlignment_compatible_bool
				= pairAlignmentInfo.overlapExons_compatible_bool();
			if(pairedAlignment_compatible_bool)
			{
				compatibleSAM_ofs << tmpAlign_1 << endl;
				compatibleSAM_ofs << tmpAlign_2 << endl;
			}
			else
			{
				incompatibleSAM_ofs << tmpAlign_1 << endl;
				incompatibleSAM_ofs << tmpAlign_2 << endl;
			}
		}
	}

	void refineAlignment(
		string& inputSAMpath,
		string& finalOutputSAMpath,
		string& correctedSAMpath,
		string& correctedSAMpath_ori,
		AlignInferJunctionHash_Info* alignInferJunctionHashInfo,
		Index_Info* indexInfo)
	{
		ifstream sam_ifs(inputSAMpath.c_str());
		ofstream finalSAM_ofs(finalOutputSAMpath.c_str());
		ofstream correctedSAM_ofs(correctedSAMpath.c_str());
		ofstream correctedSAM_ori_ofs(correctedSAMpath_ori.c_str());
		while(1)
		{	
			if(sam_ifs.eof())
				break;
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			if(tmpAlign_1 == "")
				break;			
			if(tmpAlign_1.at(0) == '@')
			{
				finalSAM_ofs << tmpAlign_1 << endl;
				continue;
			}
			getline(sam_ifs, tmpAlign_2);
			//if(sam_ifs.eof())
			//	break;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);

			//bool pairedAlignment_compatible_bool
			//	= pairAlignmentInfo.overlapExons_compatible_bool();
			//cout << endl << "***********************" << endl << "tmpAlign_1: " << endl << tmpAlign_1 << endl;
			bool alignmentWithInvalidSJBool
				= pairAlignmentInfo.alignmentWithInvalidSJ_bool(alignInferJunctionHashInfo);

			if(alignmentWithInvalidSJBool)
			{
				int leftJumpCodeIndexPairAfterClipping_result_1;
				pair<int,int> leftJumpCodeIndexPairAfterClipping_1;
				int leftJumpCodeIndexPairAfterClipping_result_2;
				pair<int,int> leftJumpCodeIndexPairAfterClipping_2;
				pairAlignmentInfo.refineInvalidSJinAlignment_PE(alignInferJunctionHashInfo, indexInfo,
					leftJumpCodeIndexPairAfterClipping_result_1, leftJumpCodeIndexPairAfterClipping_1,
					leftJumpCodeIndexPairAfterClipping_result_2, leftJumpCodeIndexPairAfterClipping_2);
				
				if(leftJumpCodeIndexPairAfterClipping_result_1 == -1) // no need to refine
				{
					finalSAM_ofs << tmpAlign_1 << endl;
				}
				else if(leftJumpCodeIndexPairAfterClipping_result_1 == 0) // unmap
				{
					correctedSAM_ofs << pairAlignmentInfo.returnSingleAlignmentInfo_unmap(indexInfo, true) << endl;
					correctedSAM_ori_ofs << tmpAlign_1 << endl;
					finalSAM_ofs << pairAlignmentInfo.returnSingleAlignmentInfo_unmap(indexInfo, true) << endl;
				}
				else if(leftJumpCodeIndexPairAfterClipping_result_1 == 1)
				{
					correctedSAM_ofs << pairAlignmentInfo.returnSingleAlignemntInfo_clippedWithInterval(
						indexInfo, leftJumpCodeIndexPairAfterClipping_1.first,
						leftJumpCodeIndexPairAfterClipping_1.second, true) << endl; 
					correctedSAM_ori_ofs << tmpAlign_1 << endl;
					finalSAM_ofs << pairAlignmentInfo.returnSingleAlignemntInfo_clippedWithInterval(
						indexInfo, leftJumpCodeIndexPairAfterClipping_1.first,
						leftJumpCodeIndexPairAfterClipping_1.second, true) << endl;
				}

				if(leftJumpCodeIndexPairAfterClipping_result_2 == -1) // no need to refine
				{
					finalSAM_ofs << tmpAlign_2 << endl;
				}
				else if(leftJumpCodeIndexPairAfterClipping_result_2 == 0) // unmap
				{
					correctedSAM_ofs << pairAlignmentInfo.returnSingleAlignmentInfo_unmap(indexInfo, false) << endl;
					correctedSAM_ori_ofs << tmpAlign_2 << endl;
					finalSAM_ofs << pairAlignmentInfo.returnSingleAlignmentInfo_unmap(indexInfo, false) << endl;
				}
				else if(leftJumpCodeIndexPairAfterClipping_result_2 == 1)
				{
					correctedSAM_ofs << pairAlignmentInfo.returnSingleAlignemntInfo_clippedWithInterval(
						indexInfo, leftJumpCodeIndexPairAfterClipping_2.first,
						leftJumpCodeIndexPairAfterClipping_2.second, false) << endl; 
					correctedSAM_ori_ofs << tmpAlign_2 << endl;
					finalSAM_ofs << pairAlignmentInfo.returnSingleAlignemntInfo_clippedWithInterval(
						indexInfo, leftJumpCodeIndexPairAfterClipping_2.first,
						leftJumpCodeIndexPairAfterClipping_2.second, false) << endl;
				}
			}
			else
			{
				finalSAM_ofs << tmpAlign_1 << endl;
				finalSAM_ofs << tmpAlign_2 << endl;
			}
		}
		sam_ifs.close();
		finalSAM_ofs.close();
		correctedSAM_ofs.close();
		correctedSAM_ori_ofs.close();	
	}

	void parsePEalignment(
		string& inputSAMpath,
		string& outputSAMpath,
		Index_Info* indexInfo)
	{
		// ifstream sam_ifs(inputSAMpath.c_str());
		// ofstream finalSAM_ofs(outputSAMpath.c_str());
		// while(1)
		// {	
		// 	if(sam_ifs.eof())
		// 		break;
		// 	string tmpAlign_1, tmpAlign_2;
		// 	getline(sam_ifs, tmpAlign_1);
		// 	if(tmpAlign_1 == "")
		// 		break;			
		// 	if(tmpAlign_1.at(0) == '@')
		// 	{
		// 		finalSAM_ofs << tmpAlign_1 << endl;
		// 		continue;
		// 	}
		// 	getline(sam_ifs, tmpAlign_2);
		// 	//if(sam_ifs.eof())
		// 	//	break;
		// 	PairAlignment_Info pairAlignmentInfo;
		// 	pairAlignmentInfo.initiatePairAlignmentInfo(
		// 		tmpAlign_1, tmpAlign_2, indexInfo);
		// 	finalSAM_ofs << tmpAlign_1 << endl;
		// 	finalSAM_ofs << tmpAlign_2 << endl;
		// }
		// finalSAM_ofs.close();
		// sam_ifs.close();
		FILE *fp;
		fp = fopen(inputSAMpath.c_str(), "r");
		FILE *fp_2;
		fp_2 = fopen(outputSAMpath.c_str(), "w");

		char buf[1000];
		char buf_2[1000];
		while(fgets(buf, 1000, fp)!=NULL)
		{
			int len = strlen(buf);
			if(len < 2)
				break;
			if(buf[0] == '@')
				continue;
			// char a[50], b[50], c[50], d[50], e[50];//, a[600],
			// sscanf(buf, "%s\t%s\t%s\t%s\t%s", a,b,c,d,e);
			// string a_str = a;
			// string b_str = b;
			// string c_str = c;
			// string d_str = d;
			// string e_str = e;			
			// fputs(buf, fp_2);
			string tmpAlign_1 = buf;
			fgets(buf_2, 1000, fp);
			string tmpAlign_2 = buf_2;
			PairAlignment_Info pairAlignmentInfo;
			pairAlignmentInfo.initiatePairAlignmentInfo(
				tmpAlign_1, tmpAlign_2, indexInfo);
			// finalSAM_ofs << tmpAlign_1 << endl;
			// finalSAM_ofs << tmpAlign_2 << endl;
			//fputs(buf, fp_2);
			//fputs(buf_2, fp_2);
		}
		fclose(fp);
		fclose(fp_2);
	}
};

#endif