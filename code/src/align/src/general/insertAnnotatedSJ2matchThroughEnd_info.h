// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef INSERTANNOTATEDSJ2ALIGNMENT_INFO_H
#define INSERTANNOTATEDSJ2ALIGNMENT_INFO_H

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



class InsertAnnotatedSJ2matchThroughEnd_Info
{
private:
	string oriSamStr;
	string readName;
	int chrNameInt;
	int chrMapPos_start;
	int chrMapPos_end;
	vector<Jump_Code> alignJumpCodeVec;
	string readSeq;
	int readLength;
	string qualSeq;
	string otherSamFieldStr;

	int chrMapPos_start_afterRefinementInHead;
	vector<Jump_Code> alignJumpCodeVec_afterRefinementInHead;

	int chrMapPos_start_refined;
	vector<Jump_Code> alignJumpCodeVec_refined;
public:
	InsertAnnotatedSJ2matchThroughEnd_Info()
	{}

	bool completeOrNot()
	{
		int alignJumpCodeVec_refined_size = alignJumpCodeVec_refined.size();
		if((alignJumpCodeVec_refined[0].type == "S")
			||(alignJumpCodeVec_refined[alignJumpCodeVec_refined_size-1].type == "S"))
			return false;
		else
			return true;
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

	bool initiateWithSamStr(string& tmpSamStr,	Index_Info* indexInfo)
	{
		oriSamStr = tmpSamStr;
		vector<string> samFieldVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = tmpSamStr.find("\t", startLoc);
			string tmpSamField = tmpSamStr.substr(startLoc, tabLoc-startLoc);
			samFieldVec.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		samFieldVec.push_back(tmpSamStr.substr(startLoc));		
		// read name, chr name int //
		readName = samFieldVec[0];
		string mapChrNameStr = samFieldVec[2];
		if(mapChrNameStr == "*")
			return false;
		chrNameInt = indexInfo->convertStringToInt(mapChrNameStr);		

		// sam_info
		//flag = atoi(samFieldVec[1].c_str());
		string mapChrPosStr = samFieldVec[3];
		chrMapPos_start = atoi(mapChrPosStr.c_str()); // chrMapPos
		string cigarString = samFieldVec[5];
		this->cigarString2jumpCodeVec(cigarString, alignJumpCodeVec); // jumpCodeVec		
		int jumpCodeVecSize = alignJumpCodeVec.size();
		chrMapPos_end = this->getEndPosOfSpecificJumpCode(chrMapPos_start,
			alignJumpCodeVec, jumpCodeVecSize-1);
		readSeq = samFieldVec[9];
		readLength = readSeq.length();
		qualSeq = samFieldVec[10];
		otherSamFieldStr = samFieldVec[11];
		return true;
	}

	bool insertAnnotatedSJ2matchThroughEnd(
		AlignInferJunctionHash_Info* alignInferJuncHashInfo,
		SJhash_Info* SJhashInfo, Index_Info* indexInfo)
	{
		bool insertAnnotatedSJ2alignment_inHead_bool;
		if((alignJumpCodeVec[0].type == "S")
			||(alignJumpCodeVec[0].len == 1))
			insertAnnotatedSJ2alignment_inHead_bool = false;
		else
		{	
			insertAnnotatedSJ2alignment_inHead_bool = 
				this->insertAnnotatedSJ2matchThroughHead(
					alignInferJuncHashInfo, SJhashInfo, indexInfo);
		}
		//cout << "insertAnnotatedSJ2alignment_inHead_bool: " << insertAnnotatedSJ2alignment_inHead_bool << endl;
		if(!insertAnnotatedSJ2alignment_inHead_bool)
		{
			chrMapPos_start_afterRefinementInHead = chrMapPos_start;
			int alignJumpCodeVecSize = alignJumpCodeVec.size();
			for(int tmp = 0; tmp < alignJumpCodeVecSize; tmp++)
			{
				alignJumpCodeVec_afterRefinementInHead.push_back(
					alignJumpCodeVec[tmp]);
			}
		}
		//cout << "tmpJumpCodeStr_afterRefinementInHead: " << this->jumpCodeVec2cigarString(alignJumpCodeVec_afterRefinementInHead) << endl;
		bool insertAnnotatedSJ2alignment_inTail_bool;
		int tmpJumpCodeVecSize_afterRefinementInHead 
			= alignJumpCodeVec_afterRefinementInHead.size();
		if((alignJumpCodeVec_afterRefinementInHead[tmpJumpCodeVecSize_afterRefinementInHead - 1].type == "S")
			||(alignJumpCodeVec_afterRefinementInHead[tmpJumpCodeVecSize_afterRefinementInHead - 1].len == 1))
			insertAnnotatedSJ2alignment_inTail_bool = false;
		else
		{	
			insertAnnotatedSJ2alignment_inTail_bool = 
				this->insertAnnotatedSJ2matchThroughTail(
					alignInferJuncHashInfo, SJhashInfo, indexInfo);
		}
		//cout << "insertAnnotatedSJ2alignment_inTail_bool: " << insertAnnotatedSJ2alignment_inTail_bool << endl;

		if(!insertAnnotatedSJ2alignment_inTail_bool)
		{
			if(!insertAnnotatedSJ2alignment_inHead_bool)
				return false;
			else
			{
				chrMapPos_start_refined = chrMapPos_start_afterRefinementInHead;
				for(int tmp = 0; tmp < alignJumpCodeVec_afterRefinementInHead.size();
					tmp++)
				{
					alignJumpCodeVec_refined.push_back(
						alignJumpCodeVec_afterRefinementInHead[tmp]);
				}
				//cout << "tmpJumpCodeStr_refined: " << this->jumpCodeVec2cigarString(alignJumpCodeVec_refined) << endl;
				return true;
			}
		}
		else
		{	
			//cout << "tmpJumpCodeStr_refined: " << this->jumpCodeVec2cigarString(alignJumpCodeVec_refined) << endl;
			return true;
		}
	}

	bool insertAnnotatedSJ2matchThroughHead(
		AlignInferJunctionHash_Info* alignInferJuncHashInfo,
		SJhash_Info* SJhashInfo, Index_Info* indexInfo)
	{
		int matchedHeadLen = alignJumpCodeVec[0].len;
		int candiAcceptorStartLocInRead_mostLeft = 2;
		int candiAcceptorStartLocInRead_mostRight = matchedHeadLen;
		int candiAcceptorStartPosInChr_mostLeft = chrMapPos_start + 1;
		int candiAcceptorStartPosInChr_mostRight = chrMapPos_start + matchedHeadLen - 1;

		vector< pair<int,int> > candiAnnotatedSJpairVec;
		vector< pair<int,int> > remapWithAnnotatedSJinfovec; // <acceptorStartLocInRead, extendedLenInSJdoner>
		for(int tmpAcceptorStartLocInRead = candiAcceptorStartLocInRead_mostLeft;
			tmpAcceptorStartLocInRead <= candiAcceptorStartLocInRead_mostRight; 
			tmpAcceptorStartLocInRead ++)
		{
			int tmpAcceptorStartPosInChr = candiAcceptorStartPosInChr_mostLeft
				+ tmpAcceptorStartLocInRead - candiAcceptorStartLocInRead_mostLeft;
			vector<int> tmpCandiDonerEndPosInChrVec;
			SJhashInfo->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(
				chrNameInt, tmpAcceptorStartPosInChr, tmpCandiDonerEndPosInChrVec);
			for(int tmp = 0; tmp < tmpCandiDonerEndPosInChrVec.size(); tmp++)
			{
				vector<int> tmpExtensionMismatchPosVec;
				vector<char> tmpExtensionMismatchCharVec;
				int tmpExtensionHead = extensionBackwards_errorTolerated_finalStepForAligner(
					chrNameInt, tmpCandiDonerEndPosInChrVec[tmp], indexInfo,
					readSeq, tmpAcceptorStartLocInRead-1, 1,
					tmpExtensionMismatchPosVec, tmpExtensionMismatchCharVec,
					MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH_WITH_ANNOTATED_SJ, 
					MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH_WITH_ANNOTATED_SJ);
				if(tmpExtensionHead == 0)
				{}
				else
				{
					candiAnnotatedSJpairVec.push_back(pair<int,int>(tmpCandiDonerEndPosInChrVec[tmp],
						tmpAcceptorStartPosInChr));
					remapWithAnnotatedSJinfovec.push_back(pair<int,int>(tmpAcceptorStartLocInRead,
						tmpExtensionHead));
				}
			}
		}

		if(candiAnnotatedSJpairVec.size() == 0)
			return false;
		else
		{
			int tmpLeftMostStartMappingBaseLocInRead = matchedHeadLen + 1;
			int tmpLeftMostStartMappingIndexInCandiVec;
			for(int tmp = 0; tmp < candiAnnotatedSJpairVec.size(); tmp++)
			{
				int tmpAcceptorStartLocInRead = remapWithAnnotatedSJinfovec[tmp].first;
				int tmpExtensionLen = remapWithAnnotatedSJinfovec[tmp].second;
				int tmpStartMappingBaseLocInRead = tmpAcceptorStartLocInRead - tmpExtensionLen;
				if(tmpStartMappingBaseLocInRead < tmpLeftMostStartMappingBaseLocInRead)
				{
					tmpLeftMostStartMappingBaseLocInRead = tmpStartMappingBaseLocInRead;
					tmpLeftMostStartMappingIndexInCandiVec = tmp;
				}
			}
			int finalSoftClipHeadLength = tmpLeftMostStartMappingBaseLocInRead - 1;
			int finalDonerEndPosInChr = candiAnnotatedSJpairVec[tmpLeftMostStartMappingIndexInCandiVec].first;
			int finalAcceptorStartPosInChr = candiAnnotatedSJpairVec[tmpLeftMostStartMappingIndexInCandiVec].second;
			int finalSJlength = finalAcceptorStartPosInChr - finalDonerEndPosInChr - 1;
			int finalAcceptorStartLocInRead 
				= remapWithAnnotatedSJinfovec[tmpLeftMostStartMappingIndexInCandiVec].first;			
			int finalMatchLength_donerSide = finalAcceptorStartLocInRead - 1 - finalSoftClipHeadLength;
			int finalMatchLength_acceptorSide = matchedHeadLen - finalAcceptorStartLocInRead + 1;

			if(finalSoftClipHeadLength > 0)
			{
				#ifdef COMPLETEONLY_WHEN_INSERT_ANNOTATEDSJ
				return false;
				#endif
				Jump_Code softClipHeadJumpCode(finalSoftClipHeadLength, "S");
				alignJumpCodeVec_afterRefinementInHead.push_back(softClipHeadJumpCode);
			}
			Jump_Code finalDonerSideMatchJumpCode(finalMatchLength_donerSide, "M");
			Jump_Code finalSJjumpCode(finalSJlength, "N");
			Jump_Code finalAcceptorSideMatchJumpCode(finalMatchLength_acceptorSide, "M");
			alignJumpCodeVec_afterRefinementInHead.push_back(finalDonerSideMatchJumpCode);
			alignJumpCodeVec_afterRefinementInHead.push_back(finalSJjumpCode);
			alignJumpCodeVec_afterRefinementInHead.push_back(finalAcceptorSideMatchJumpCode);
			for(int tmp = 1; tmp < alignJumpCodeVec.size(); tmp++)
				alignJumpCodeVec_afterRefinementInHead.push_back(alignJumpCodeVec[tmp]);

			chrMapPos_start_afterRefinementInHead = chrMapPos_start + matchedHeadLen 
				- finalMatchLength_acceptorSide - finalSJlength
				- finalMatchLength_donerSide;
			return true;
		}
	}

	bool insertAnnotatedSJ2matchThroughTail(
		AlignInferJunctionHash_Info* alignInferJuncHashInfo,
		SJhash_Info* SJhashInfo, Index_Info* indexInfo)
	{
		int tmpJumpCodeVecSize = alignJumpCodeVec_afterRefinementInHead.size();
		int matchedTailLen 
			= alignJumpCodeVec_afterRefinementInHead[tmpJumpCodeVecSize - 1].len;
		//cout << "matchedTailLen: " << matchedTailLen << endl;
		int candiDonerEndLocInRead_mostLeft = readLength - matchedTailLen + 1;
		int candiDonerEndLocInRead_mostRight = readLength - 1;
		int candiDonerEndPosInChr_mostLeft = chrMapPos_end - matchedTailLen + 1;
		int candiDonerEndPosInChr_mostRight = chrMapPos_end - 1;
		//cout << "candiDonerEndLocInRead_mostLeft: " << candiDonerEndLocInRead_mostLeft << endl;
		//cout << "candiDonerEndLocInRead_mostRight: " << candiDonerEndLocInRead_mostRight << endl;
		//cout << "candiDonerEndPosInChr_mostLeft: " << candiDonerEndPosInChr_mostLeft << endl;
		//cout << "candiDonerEndPosInChr_mostRight: " << candiDonerEndPosInChr_mostRight << endl;
		vector< pair<int,int> > candiAnnotatedSJpairVec;
		vector< pair<int,int> > remapWithAnnotatedSJinfovec;//<donerEndLocInRead, extendedLenInSJacceptorSide>		
		for(int tmpDonerEndLocInRead = candiDonerEndLocInRead_mostLeft;
			tmpDonerEndLocInRead <= candiDonerEndLocInRead_mostRight;
			tmpDonerEndLocInRead ++)
		{
			int tmpDonerEndPosInChr = candiDonerEndPosInChr_mostLeft
				+ tmpDonerEndLocInRead - candiDonerEndLocInRead_mostLeft;
			vector<int> tmpCandiAcceptorStartPosInChrVec;
			SJhashInfo->returnAlterAcceptorSpliceSiteVecWithDonerEndPos(
				chrNameInt, tmpDonerEndPosInChr, tmpCandiAcceptorStartPosInChrVec);
			for(int tmp = 0; tmp < tmpCandiAcceptorStartPosInChrVec.size(); tmp++)
			{
				vector<int> tmpExtensionMismatchPosVec;
				vector<char> tmpExtensionMismatchCharVec;
				int tmpExtensionTail = extensionForwards_errorTolerated_finalStepForAligner(
					chrNameInt, tmpCandiAcceptorStartPosInChrVec[tmp], indexInfo,
					readSeq, tmpDonerEndLocInRead + 1, readLength, 
					tmpExtensionMismatchPosVec, tmpExtensionMismatchCharVec,
					MIN_MATCH_BASE_TO_SUPPORT_PER_MISMATCH_WITH_ANNOTATED_SJ, 
					MIN_MATCH_BASE_TO_SUPPORT_TWO_MISMATCH_WITH_ANNOTATED_SJ);
				if(tmpExtensionTail == 0)
				{}
				else
				{
					candiAnnotatedSJpairVec.push_back(
						pair<int,int>(tmpDonerEndPosInChr,
						tmpCandiAcceptorStartPosInChrVec[tmp]));
					remapWithAnnotatedSJinfovec.push_back(
						pair<int,int>(tmpDonerEndLocInRead, tmpExtensionTail));	
					//cout << "getCandiSJ: " << tmpDonerEndPosInChr << "~" << tmpCandiAcceptorStartPosInChrVec[tmp] << endl;
					//cout<< "candiSJ loc in read: " << tmpDonerEndLocInRead << " extensionLength: " << tmpExtensionTail << endl; 			
				}	
			}
		}

		if(candiAnnotatedSJpairVec.size() == 0)
			return false;
		else
		{
			int tmpRightMostEndMappingBaseLocInRead = readLength - matchedTailLen;
			int tmpRightMostEndMappingIndexInCandiVec;
			for(int tmp = 0; tmp < candiAnnotatedSJpairVec.size(); tmp++)
			{
				int tmpDonerEndLocInRead = remapWithAnnotatedSJinfovec[tmp].first;
				int tmpExtensionLen = remapWithAnnotatedSJinfovec[tmp].second;
				int tmpEndMappingBaseLocInRead = tmpDonerEndLocInRead + tmpExtensionLen;
				if(tmpEndMappingBaseLocInRead > tmpRightMostEndMappingBaseLocInRead)
				{
					tmpRightMostEndMappingBaseLocInRead = tmpEndMappingBaseLocInRead;
					tmpRightMostEndMappingIndexInCandiVec = tmp;
				}
			}

			int finalSoftClipTailLength = readLength - tmpRightMostEndMappingBaseLocInRead;
			int finalDonerEndPosInChr = candiAnnotatedSJpairVec[tmpRightMostEndMappingIndexInCandiVec].first;
			int finalAcceptorStartPosInChr = candiAnnotatedSJpairVec[tmpRightMostEndMappingIndexInCandiVec].second;
			int finalSJlength = finalAcceptorStartPosInChr - finalDonerEndPosInChr - 1;
			int finalDonerEndLocInRead = remapWithAnnotatedSJinfovec[tmpRightMostEndMappingIndexInCandiVec].first;
			int finalMatchLength_donerSide = matchedTailLen - (readLength - finalDonerEndLocInRead);
			int finalMatchLength_acceptorSide = readLength - finalDonerEndLocInRead - finalSoftClipTailLength;
			Jump_Code finalDonerSideMatchJumpCode(finalMatchLength_donerSide, "M");
			Jump_Code finalSJjumpCode(finalSJlength, "N");
			Jump_Code finalAcceptorSideMatchJumpCode(finalMatchLength_acceptorSide, "M");

			//cout << "finalMatchLength_donerSide: " << finalMatchLength_donerSide << endl;
			//cout << "finalSJlength: " << finalSJlength << endl;
			//cout << "finalMatchLength_acceptorSide: " << finalMatchLength_acceptorSide << endl; 
			//cout << "finalSoftClipTailLength: " << finalSoftClipTailLength << endl;
			if(finalSoftClipTailLength > 0)
			{
				//cout << "finalSoftClipTailLength: " << finalSoftClipTailLength << endl;
				#ifdef COMPLETEONLY_WHEN_INSERT_ANNOTATEDSJ
				return false;
				#endif
				for(int tmp = 0; tmp < alignJumpCodeVec_afterRefinementInHead.size()-1; tmp++)
					alignJumpCodeVec_refined.push_back(alignJumpCodeVec_afterRefinementInHead[tmp]);
				alignJumpCodeVec_refined.push_back(finalDonerSideMatchJumpCode);
				alignJumpCodeVec_refined.push_back(finalSJjumpCode);
				alignJumpCodeVec_refined.push_back(finalAcceptorSideMatchJumpCode);
				Jump_Code softClipTailJumpCode(finalSoftClipTailLength, "S");
				alignJumpCodeVec_refined.push_back(softClipTailJumpCode);
			}
			else
			{
				for(int tmp = 0; tmp < alignJumpCodeVec_afterRefinementInHead.size()-1; tmp++)
					alignJumpCodeVec_refined.push_back(alignJumpCodeVec_afterRefinementInHead[tmp]);
				alignJumpCodeVec_refined.push_back(finalDonerSideMatchJumpCode);
				alignJumpCodeVec_refined.push_back(finalSJjumpCode);
				alignJumpCodeVec_refined.push_back(finalAcceptorSideMatchJumpCode);
			}

			chrMapPos_start_refined = chrMapPos_start_afterRefinementInHead;
			return true;
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

	string returnRefinedSamStr(Index_Info* indexInfo)
	{
		string tmpStr;
		string tmpChrName = indexInfo->returnChrNameStr(chrNameInt);
		string tmpChrMapPosStr = int_to_str(chrMapPos_start_refined);
		string tmpCigarString = this->jumpCodeVec2cigarString(alignJumpCodeVec_refined);
		tmpStr = readName + "\t0\t" + tmpChrName + "\t" 
			+ tmpChrMapPosStr + "\t255\t" + tmpCigarString + "\t*\t0\t0\t"
			+ readSeq + "\t" + qualSeq + "\t" + otherSamFieldStr;
		return tmpStr;
	}

};
#endif