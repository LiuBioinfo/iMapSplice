// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef REFINESAM_INFO_H
#define REFINESAM_INFO_H

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

#define UNIQUE_SAM_TYPE 1
#define MULTI_SAM_TYPE 2
#define UNMAP_SAM_TYPE 3

using namespace std;

class RefineSAM_SE_Info
{
private:
	string readName;
	int chrNameInt;
	int chrMapPos;
	int chrMapPos_end;
	int unfixedHeadLen;
	int unfixedTailLen;
	vector<Jump_Code> alignJumpCodeVec;	
	//vector< pair<int,int> > SJposPairVec;
	//vector< int > SJindexInAlignInferJuncVec;
	string strand;
	string readSeq;
	string qualSeq;
	string otherStr;
	string originalSamStr;// store original SAM str

	// SJ info
	vector< int > SJindexVecInJumpCodeVec;
	vector< pair<int,int> > SJposPairVec;
	vector< pair<int,int> > SJanchorSizePairVec;
	vector< vector<Jump_Code> > SJjumpCodeVecVec_backward;
	vector< vector<Jump_Code> > SJjumpCodeVecVec_forward;

	//
	vector< int > SJindexVecInAlignInferJuncHash;
public:
	RefineSAM_SE_Info()
	{}

	string returnReadName()
	{
		return readName;
	}

	string returnReadSeq()
	{
		return readSeq;
	}

	string returnQualSeq()
	{
		return qualSeq;
	}

	string returnUnmapSAMstr()
	{
		string tmpStr;
		tmpStr = tmpStr + readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t";
		if(strand == "+")
			tmpStr = tmpStr + readSeq + "\t" + qualSeq;
		else
			tmpStr = tmpStr + getRcmSeq(readSeq) + "\t" + getRevSeq(qualSeq);
		tmpStr += "\tIH:i:0\tHI:i:0";
		return tmpStr;
	}

	string returnOriSAMstr()
	{
		return originalSamStr;
	}

	string jumpCodeVec2cigarString()
	{
		string tmpCigarString;
		//cout << "********" << "jumpCodeVecSize: " << jumpCodeVec.size() << endl;
		for(int tmp = 0; tmp < alignJumpCodeVec.size(); tmp++)
		{
			tmpCigarString += alignJumpCodeVec[tmp].toString();
		}
		return tmpCigarString;
	}	

	string returnRefineSAMinfoStr()
	{
		string tmpStr;
		tmpStr = "tmpChr:" + int_to_str(chrNameInt) + " chrMapPos:" + int_to_str(chrMapPos) + " cigarStr:" + this->jumpCodeVec2cigarString();  
		return tmpStr;
	}

	string returnSJindexVecInAlignInferJuncHashStr()
	{
		string tmpStr;
		for(int tmp = 0; tmp < SJindexVecInAlignInferJuncHash.size(); tmp++)
		{
			tmpStr = tmpStr + "SJindex: " + int_to_str(SJindexVecInAlignInferJuncHash[tmp]) + " ";
		}
		return tmpStr;
	}

	int returnSJindexInAlignInferJuncHash(int index)
	{
		return SJindexVecInAlignInferJuncHash[index];
	}

	vector<Jump_Code>& returnSJjumpCodeVec_backward(int tmp)
	{
		return SJjumpCodeVecVec_backward[tmp];
	}

	vector<Jump_Code>& returnSJjumpCodeVec_forward(int tmp)
	{
		return SJjumpCodeVecVec_forward[tmp];
	}

	int returnSJanchorSize_doner(int tmp)
	{
		return SJanchorSizePairVec[tmp].first;
	}

	int returnSJanchorSize_acceptor(int tmp)
	{
		return SJanchorSizePairVec[tmp].second;
	}

	int returnSJvecSize()
	{
		return SJindexVecInJumpCodeVec.size();
	}

	int returnSJ_donerEndPos(int tmp)
	{
		return SJposPairVec[tmp].first;
	}

	int returnSJ_acceptorStartPos(int tmp)
	{
		return SJposPairVec[tmp].second;
	}

	void getUnfixedHeadTailLen()
	{
		int alignJumpCodeVecSize = alignJumpCodeVec.size();
		
		string firstJumpCodeType = alignJumpCodeVec[0].type;
		if(firstJumpCodeType == "S")
			unfixedHeadLen = alignJumpCodeVec[0].len;
		else
			unfixedHeadLen = 0;

		string lastJumpCodeType = alignJumpCodeVec[alignJumpCodeVecSize-1].type;
		if(lastJumpCodeType == "S")
			unfixedTailLen = alignJumpCodeVec[alignJumpCodeVecSize-1].len;
		else
			unfixedTailLen = 0;
	}

	string returnRcmReadSeq()
	{
		string rcmReadSeq = getRcmSeq(readSeq);
		return rcmReadSeq;
	}

	void cigarString2jumpCodeVec(string& jumpCodeStr)
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
				alignJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
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

	int checkMultiUnmap_initiate_SE(
		string& samStr, Index_Info* indexInfo, int maxReadBaseNum)
	{
		//cout << "samStr: " << endl << samStr << endl;
		int tmpIHint = checkIH(samStr);
		if(tmpIHint == 0)
			return UNMAP_SAM_TYPE;
		else if(tmpIHint > 1)
			return MULTI_SAM_TYPE;
		else
		{}

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
		int flagInt = atoi(flagStr.c_str());
		bool forwardOrReverseComplement_bool = forwardOrReverseComplement(flagInt);
		if(forwardOrReverseComplement_bool)
			strand = "+";
		else
			strand = "-";
		string chrNameStr = samFieldVec[2];
		chrNameInt = indexInfo->convertStringToInt(chrNameStr);
		string mapChrPosStr = samFieldVec[3];
		chrMapPos = atoi(mapChrPosStr.c_str());
		//cout << "chrMapPos: " << chrMapPos << endl;
		string cigarString = samFieldVec[5];
		this->cigarString2jumpCodeVec(cigarString);
		this->getUnfixedHeadTailLen();
		chrMapPos_end = this->getEndPosOfSpecificJumpCode(alignJumpCodeVec.size()-1);
		//cout << "chrMapPos_end: " << chrMapPos_end << endl;
		readSeq = samFieldVec[9];
		qualSeq = samFieldVec[10];
		otherStr = samStr.substr(startLoc);
		//cout << "start to generate SJ" << endl;
		this->generate_SJposPairVec_SJanchorSizePairVec_jumpCodeVecVec(maxReadBaseNum);
		
		originalSamStr = samStr;
		//cout << "end of generate SJ " << endl;
		return UNIQUE_SAM_TYPE;
	}

	void generate_SJposPairVec_SJanchorSizePairVec_jumpCodeVecVec(
		int maxReadBaseNum)
	{
		//cout << "start to generateSJindexVec_cigarStringJumpCodeVec_SJposPairVec" << endl;
		this->generateSJindexVec_cigarStringJumpCodeVec_SJposPairVec();
		//cout << "start to generateSJanchorSizeVec" << endl;
		this->generateSJanchorSizeVec();
		//cout << "start to generateSJjumpCodeVec_maxReadBaseNumInPathStruct" << endl;
		this->generateSJjumpCodeVec_maxReadBaseNumInPathStruct(maxReadBaseNum);
	}

	void generateSJindexVec_cigarStringJumpCodeVec_SJposPairVec()
	{
		for(int tmp = 0; tmp < alignJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = alignJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				int tmpDonerEndPos = getEndPosOfSpecificJumpCode(lastJumpCodeIndex);
				int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(currentJumpCodeIndex) + 1;
				SJposPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				SJindexVecInJumpCodeVec.push_back(tmp);
			}
		}				
	}

	void generateSJanchorSizeVec()
	{
		int SJindexVecSize = SJindexVecInJumpCodeVec.size();
		int cigarStringJumpCodeVecSize = alignJumpCodeVec.size();
		if(SJindexVecSize == 0)
		{}
		else if(SJindexVecSize == 1)
		{
			int tmpDonerAnchorFirstJumpCodeIndex = 0;
			int tmpDonerAnchorLastJumpCodeIndex = SJindexVecInJumpCodeVec[0]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex = SJindexVecInJumpCodeVec[0]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex = cigarStringJumpCodeVecSize - 1;
			int tmpDonerAnchorSize = this->getEndLocInReadOfSpecificJumpCode(
				tmpDonerAnchorLastJumpCodeIndex) - unfixedHeadLen;
			int tmpAcceptorAnchorSize
				= this->getEndLocInReadOfSpecificJumpCode(
					tmpAcceptorAnchorLastJumpCodeIndex)
					- this->getEndLocInReadOfSpecificJumpCode(
						tmpAcceptorAnchorFirstJumpCodeIndex-1) - unfixedTailLen;
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize, tmpAcceptorAnchorSize));
		}
		else if(SJindexVecSize == 2)
		{
			int tmpDonerAnchorFirstJumpCodeIndex_1stSJ = 0;
			int tmpDonerAnchorLastJumpCodeIndex_1stSJ = SJindexVecInJumpCodeVec[0]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ = SJindexVecInJumpCodeVec[0]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_1stSJ = SJindexVecInJumpCodeVec[1]-1;
			int tmpDonerAnchorSize_1stSJ = this->getEndLocInReadOfSpecificJumpCode(
				tmpDonerAnchorLastJumpCodeIndex_1stSJ) - unfixedHeadLen;
			int tmpAcceptorAnchorSize_1stSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					tmpAcceptorAnchorLastJumpCodeIndex_1stSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ-1);
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_1stSJ, tmpAcceptorAnchorSize_1stSJ));
			int tmpAcceptorAnchorFirstJumpCodeIndex_2ndSJ = SJindexVecInJumpCodeVec[1]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_2ndSJ = cigarStringJumpCodeVecSize - 1;
			int tmpDonerAnchorSize_2ndSJ = tmpAcceptorAnchorSize_1stSJ;
			int tmpAcceptorAnchorSize_2ndSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					tmpAcceptorAnchorLastJumpCodeIndex_2ndSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						tmpAcceptorAnchorFirstJumpCodeIndex_2ndSJ-1) - unfixedTailLen;
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_2ndSJ, tmpAcceptorAnchorSize_2ndSJ));
		}
		else if(SJindexVecSize > 2)
		{
			// 1st SJ
			int tmpDonerAnchorFirstJumpCodeIndex_1stSJ = 0;
			int tmpDonerAnchorLastJumpCodeIndex_1stSJ = SJindexVecInJumpCodeVec[0]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ = SJindexVecInJumpCodeVec[0]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_1stSJ = SJindexVecInJumpCodeVec[1]-1;
			int tmpDonerAnchorSize_1stSJ = this->getEndLocInReadOfSpecificJumpCode(
				tmpDonerAnchorLastJumpCodeIndex_1stSJ) - unfixedHeadLen;
			int tmpAcceptorAnchorSize_1stSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					tmpAcceptorAnchorLastJumpCodeIndex_1stSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						tmpAcceptorAnchorFirstJumpCodeIndex_1stSJ-1);
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_1stSJ, tmpAcceptorAnchorSize_1stSJ));
			// other SJs
			for(int tmpSJ = 1; tmpSJ <= SJindexVecSize-2; tmpSJ++)
			{
				int tmpDonerAnchorFirstJumpCodeIndex_tmpSJ = SJindexVecInJumpCodeVec[tmpSJ-1]+1;
				int tmpDonerAnchorLastJumpCodeIndex_tmpSJ = SJindexVecInJumpCodeVec[tmpSJ]-1;
				int tmpAcceptorAnchorFirstJumpCodeIndex_tmpSJ = SJindexVecInJumpCodeVec[tmpSJ]+1;
				int tmpAcceptorAnchorLastJumpCodeIndex_tmpSJ = SJindexVecInJumpCodeVec[tmpSJ+1]-1;
				int tmpDonerAnchorSize_tmpSJ 
					= this->getEndLocInReadOfSpecificJumpCode(
						tmpDonerAnchorLastJumpCodeIndex_tmpSJ);
						- this->getEndLocInReadOfSpecificJumpCode(
							tmpDonerAnchorFirstJumpCodeIndex_tmpSJ-1);
				int tmpAcceptorAnchorSize_tmpSJ 
					= this->getEndLocInReadOfSpecificJumpCode(
						tmpAcceptorAnchorLastJumpCodeIndex_tmpSJ)
						- this->getEndLocInReadOfSpecificJumpCode(
							tmpAcceptorAnchorFirstJumpCodeIndex_tmpSJ-1);
				SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_tmpSJ, tmpAcceptorAnchorSize_tmpSJ));							
			}
			// last SJ
			int tmpDonerAnchorFirstJumpCodeIndex_lastSJ = SJindexVecInJumpCodeVec[SJindexVecSize-2]+1;
			int tmpDonerAnchorLastJumpCodeIndex_lastSJ = SJindexVecInJumpCodeVec[SJindexVecSize-1]-1;
			int tmpAcceptorAnchorFirstJumpCodeIndex_lastSJ = SJindexVecInJumpCodeVec[SJindexVecSize-1]+1;
			int tmpAcceptorAnchorLastJumpCodeIndex_lastSJ = cigarStringJumpCodeVecSize - 1;
			int tmpDonerAnchorSize_lastSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					tmpDonerAnchorLastJumpCodeIndex_lastSJ);
					- this->getEndLocInReadOfSpecificJumpCode(
						tmpDonerAnchorFirstJumpCodeIndex_lastSJ-1);
			int tmpAcceptorAnchorSize_lastSJ 
				= this->getEndLocInReadOfSpecificJumpCode(
					tmpAcceptorAnchorLastJumpCodeIndex_lastSJ)
					- this->getEndLocInReadOfSpecificJumpCode(
						tmpAcceptorAnchorFirstJumpCodeIndex_lastSJ-1) - unfixedTailLen; 					
			SJanchorSizePairVec.push_back(pair<int,int>(tmpDonerAnchorSize_lastSJ, tmpAcceptorAnchorSize_lastSJ));
		}
		else
		{}
	}

	void generateSJjumpCodeVec_maxReadBaseNumInPathStruct(int maxReadBaseNum)
	{
		int cigarStringJumpCodeVecSize = alignJumpCodeVec.size();
		//cout << "cigarStringJumpCodeVecSize: " << cigarStringJumpCodeVecSize << endl;
		for(int tmp = 0; tmp < SJindexVecInJumpCodeVec.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int tmpSJindex_cigarStringJumpCodeVec = SJindexVecInJumpCodeVec[tmp];
			//cout << "tmpSJindex_cigarStringJumpCodeVec " << tmpSJindex_cigarStringJumpCodeVec << endl;
			vector<Jump_Code> tmpSJjumpCodeVec_backward;
			vector<Jump_Code> tmpSJjumpCodeVec_forward;

			this->generateSJjumpCodeVec_backward_maxReadBaseNumInPathStruct(
				tmpSJindex_cigarStringJumpCodeVec,
				tmpSJjumpCodeVec_backward, maxReadBaseNum);
			this->generateSJjumpCodeVec_forward_maxReadBaseNumInPathStruct(
				tmpSJindex_cigarStringJumpCodeVec,
				tmpSJjumpCodeVec_forward, maxReadBaseNum);

			SJjumpCodeVecVec_backward.push_back(tmpSJjumpCodeVec_backward);
			SJjumpCodeVecVec_forward.push_back(tmpSJjumpCodeVec_forward);
		}
	}	

	void generateSJjumpCodeVec_backward_maxReadBaseNumInPathStruct(
		int tmpSJindex_cigarStringJumpCodeVec,
		vector<Jump_Code>& tmpSJjumpCodeVec_backward, 
		int maxReadBaseNum)
	{
		int tmpReadBaseNum = 0;
		for(int tmpJumpCodeIndex = tmpSJindex_cigarStringJumpCodeVec-1; tmpJumpCodeIndex >= 0; tmpJumpCodeIndex--)
		{
			int tmpJumpCodeLength = alignJumpCodeVec[tmpJumpCodeIndex].len;
			string tmpJumpCodeType = alignJumpCodeVec[tmpJumpCodeIndex].type;
			if(tmpJumpCodeType == "S")
				return;
			if((tmpJumpCodeType == "I")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "S"))
				tmpReadBaseNum = tmpReadBaseNum + tmpJumpCodeLength;
			if(tmpReadBaseNum < maxReadBaseNum)
			{
				tmpSJjumpCodeVec_backward.push_back(alignJumpCodeVec[tmpJumpCodeIndex]);
			}
			else if(tmpReadBaseNum == maxReadBaseNum)
			{
				tmpSJjumpCodeVec_backward.push_back(alignJumpCodeVec[tmpJumpCodeIndex]);
				return;
			}
			else
			{
				int overJumpCodeLength = tmpReadBaseNum - maxReadBaseNum;
				int lastJumpCodeLength = tmpJumpCodeLength - overJumpCodeLength;
				Jump_Code lastJumpCode(lastJumpCodeLength, tmpJumpCodeType);
				tmpSJjumpCodeVec_backward.push_back(lastJumpCode);
				return;
			}
		}
	}

	void generateSJjumpCodeVec_forward_maxReadBaseNumInPathStruct(
		int tmpSJindex_cigarStringJumpCodeVec,
		vector<Jump_Code>& tmpSJjumpCodeVec_forward,
		int maxReadBaseNum)
	{
		int tmpReadBaseNum = 0;
		int cigarStringJumpCodeVecSize = alignJumpCodeVec.size();
		for(int tmpJumpCodeIndex = tmpSJindex_cigarStringJumpCodeVec+1; 
			tmpJumpCodeIndex < cigarStringJumpCodeVecSize; tmpJumpCodeIndex++)
		{
			int tmpJumpCodeLength = alignJumpCodeVec[tmpJumpCodeIndex].len;
			string tmpJumpCodeType = alignJumpCodeVec[tmpJumpCodeIndex].type;
			if(tmpJumpCodeType == "S")
				return;
			if((tmpJumpCodeType == "I")||(tmpJumpCodeType == "M")||(tmpJumpCodeType == "S"))
				tmpReadBaseNum = tmpReadBaseNum + tmpJumpCodeLength;
			if(tmpReadBaseNum < maxReadBaseNum)
			{
				tmpSJjumpCodeVec_forward.push_back(alignJumpCodeVec[tmpJumpCodeIndex]);
			}
			else if(tmpReadBaseNum == maxReadBaseNum)
			{
				tmpSJjumpCodeVec_forward.push_back(alignJumpCodeVec[tmpJumpCodeIndex]);
				return;
			}
			else
			{
				int overJumpCodeLength = tmpReadBaseNum - maxReadBaseNum;
				int lastJumpCodeLength = tmpJumpCodeLength - overJumpCodeLength;
				Jump_Code lastJumpCode(lastJumpCodeLength, tmpJumpCodeType);
				tmpSJjumpCodeVec_forward.push_back(lastJumpCode);
				return;
			}
		}
	}

	int getEndLocInReadOfSpecificJumpCode(int jumpCodeIndex)
	{
		if(jumpCodeIndex < 0)
			return 0;
		int tmpEndLocInRead = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = alignJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = alignJumpCodeVec[tmpIndex].len;
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

	int getEndPosOfSpecificJumpCode(int jumpCodeIndex)
	{
		int tmpEndPos = 0;
		for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
		{
			string tmpJumpCodeType = alignJumpCodeVec[tmpIndex].type;
			int tmpJumpCodeLength = alignJumpCodeVec[tmpIndex].len;
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
		return (tmpEndPos + chrMapPos-1);
	}

	void pushBackAlignInferJuncIndex(int tmpIndexInAlignInferJuncHash)
	{
		SJindexVecInAlignInferJuncHash.push_back(tmpIndexInAlignInferJuncHash);
	}
};

class RefineSAM_SE_Vec_Info
{
private:
	vector<RefineSAM_SE_Info* > refineSAMvec;
	int toAddRefineSAMindex;
public:
	RefineSAM_SE_Vec_Info()
	{
		toAddRefineSAMindex = 0;
	}

	string returnUnmapSAMstr(int tmpRefineSAMindex)
	{
		return refineSAMvec[tmpRefineSAMindex]->returnUnmapSAMstr();
	}

	string returnOriSAMstr(int tmpRefineSAMindex)
	{
		return refineSAMvec[tmpRefineSAMindex]->returnOriSAMstr();
	}

	bool point2NULLbool(int index)
	{
		if(refineSAMvec[index] == NULL)
			return true;
		else
			return false;
	}

	int returnRefineSAMvecSize()
	{
		return refineSAMvec.size();
	}

	int returnCurrentIndex()
	{
		return toAddRefineSAMindex;
	}
	void toAddRefineSAMindex_incre()
	{
		toAddRefineSAMindex++;
	}

	int returnIndexInAlignInferJuncHash(int tmpIndexInRefineSAMinfoVec,
		int tmpSJindexInRefineSAMinfo)	
	{
		(refineSAMvec[tmpIndexInRefineSAMinfoVec])->returnSJindexInAlignInferJuncHash(
			tmpSJindexInRefineSAMinfo);
	}

	int returnRefineSAMinfoSJvecSize(int tmpIndex)
	{
		return refineSAMvec[tmpIndex]->returnSJvecSize();
	}

	void pushBackRefineSAMseInfo(RefineSAM_SE_Info* refineSAMseInfo)
	{
		refineSAMvec.push_back(refineSAMseInfo);
		this->toAddRefineSAMindex_incre();
	}

	void freeMemory()
	{
		for(int tmp = 0; tmp < refineSAMvec.size(); tmp++)
		{
			if(refineSAMvec[tmp] != NULL)
			{
				delete refineSAMvec[tmp];
				refineSAMvec[tmp] = NULL;
			}
		}
	}

	void deleteRefineSAMinfo(int index)
	{
		// if(index == 225)
		// {	
		// 	cout << "deleteRefineSAMinfo:" << index << endl;
		// 	cout << "refineSAMinfoStr: " << endl 
		// 		<< refineSAMvec[index]->returnSJindexVecInAlignInferJuncHashStr() << endl;
		// 	cout << "refineSAMinfo: " << endl << refineSAMvec[index]->returnRefineSAMinfoStr() << endl;
		// }
		delete refineSAMvec[index];
		refineSAMvec[index] = NULL;
	}
};


class RefineSAM_PE_Info
{
private:
	bool Nor1Rcm2_or_Nor2Rcm1_bool;

	RefineSAM_SE_Info refineSAMinfo_SE_Nor;
	RefineSAM_SE_Info refineSAMinfo_SE_Rcm;

	vector< pair<int,int> > SJposPairVec_PE;
	vector< pair<int,int> > SJanchorSizePairVec_PE;
	vector< vector<Jump_Code> > SJjumpCodeVecVec_backward_PE;
	vector< vector<Jump_Code> > SJjumpCodeVecVec_forward_PE;

	vector< int > SJindexVecInAlignInferJuncHash_PE;	
public:
	RefineSAM_PE_Info()
	{

	}

	string returnUnmapSAMstr()
	{
		string tmpStr;
		string readNameStr_Nor = refineSAMinfo_SE_Nor.returnReadName();
		string readNameStr_Rcm = refineSAMinfo_SE_Rcm.returnReadName();
		string readOriSeqStr_Nor = refineSAMinfo_SE_Nor.returnReadSeq();
		string readOriSeqStr_Rcm = getRcmSeq(refineSAMinfo_SE_Rcm.returnReadSeq());
		string readOriQualSeqStr_Nor = refineSAMinfo_SE_Nor.returnQualSeq();
		string readOriQualSeqStr_Rcm = getRevSeq(refineSAMinfo_SE_Rcm.returnQualSeq());
		if(Nor1Rcm2_or_Nor2Rcm1_bool)
		{
			tmpStr = tmpStr + readNameStr_Nor + "\t77\t*\t0\t0\t*\t*\t0\t0\t" 
				+ readOriSeqStr_Nor + "\t" + readOriQualSeqStr_Nor + "\tIH:i:0\tHI:i:0\n";
			tmpStr = tmpStr + readNameStr_Rcm + "\t144\t*\t0\t0\t*\t*\t0\t0\t"
				+ readOriSeqStr_Rcm + "\t" + readOriQualSeqStr_Rcm + "\tIH:i:0\tHI:i:0";
		}
		else
		{
			tmpStr = tmpStr + readNameStr_Rcm + "\t77\t*\t0\t0\t*\t*\t0\t0\t" 
				+ readOriSeqStr_Rcm + "\t" + readOriQualSeqStr_Rcm + "\tIH:i:0\tHI:i:0\n";
			tmpStr = tmpStr + readNameStr_Nor + "\t144\t*\t0\t0\t*\t*\t0\t0\t"
				+ readOriSeqStr_Nor + "\t" + readOriQualSeqStr_Nor + "\tIH:i:0\tHI:i:0";
		}
		return tmpStr;
	}

	string returnOriSAMstr()
	{
		string tmpOriSamStr_1 = refineSAMinfo_SE_Nor.returnOriSAMstr();
		string tmpOriSamStr_2 = refineSAMinfo_SE_Rcm.returnOriSAMstr();
		string tmpStr = tmpOriSamStr_1 + "\n" + tmpOriSamStr_2;
		return tmpStr;
	}

	int returnChrNameInt()
	{
		return refineSAMinfo_SE_Nor.returnChrNameInt();
	}

	int returnSJvecSize()
	{
		return SJposPairVec_PE.size();
	}

	int returnSJ_donerEndPos(int tmp)
	{
		return SJposPairVec_PE[tmp].first;
	}

	int returnSJ_acceptorStartPos(int tmp)
	{
		return SJposPairVec_PE[tmp].second;
	}

	int returnSJanchorSize_doner(int tmp)
	{
		return SJanchorSizePairVec_PE[tmp].first;
	}

	int returnSJanchorSize_acceptor(int tmp)
	{
		return SJanchorSizePairVec_PE[tmp].second;
	}

	vector<Jump_Code>& returnSJjumpCodeVec_backward(int tmp)
	{
		return SJjumpCodeVecVec_backward_PE[tmp];
	}

	vector<Jump_Code>& returnSJjumpCodeVec_forward(int tmp)
	{
		return SJjumpCodeVecVec_forward_PE[tmp];
	}

	void pushBackAlignInferJuncIndex(int tmpIndexInAlignInferJuncHash)
	{
		SJindexVecInAlignInferJuncHash_PE.push_back(tmpIndexInAlignInferJuncHash);
	}

	int checkMultiUnmap_initiate_PE(
		string& samStr_1, string& samStr_2, 
		Index_Info* indexInfo, int maxReadBaseNum)
	{
		//cout << "samStr: " << endl << samStr << endl;
		int tmpIHint = checkIH(samStr_1);
		if(tmpIHint == 0)
			return UNMAP_SAM_TYPE;
		else if(tmpIHint > 1)
			return MULTI_SAM_TYPE;
		else
		{}
		refineSAMinfo_SE_Nor.checkMultiUnmap_initiate_SE(samStr_1, indexInfo, maxReadBaseNum);
		refineSAMinfo_SE_Rcm.checkMultiUnmap_initiate_SE(samStr_2, indexInfo, maxReadBaseNum);
		this->generateSJinfoVec_PE();
		//cout << "end of generate SJ " << endl;
		return UNIQUE_SAM_TYPE;
	}

	void generateSJinfoVec_PE()
	{
		int SE_nor_SJvecSize = refineSAMinfo_SE_Nor.returnSJvecSize();
		for(int tmpSJindex = 0; tmpSJindex < SE_nor_SJvecSize; tmpSJindex++)
		{
			int tmpSJpos_donerEnd = refineSAMinfo_SE_Nor.returnSJ_donerEndPos(tmpSJindex);
			int tmpSJpos_acceptorStart = refineSAMinfo_SE_Nor.returnSJ_acceptorStartPos(tmpSJindex);
			SJposPairVec_PE.push_back(pair<int,int>(tmpSJpos_donerEnd, tmpSJpos_acceptorStart));
			int tmpSJanchorSize_doner = refineSAMinfo_SE_Nor.returnSJanchorSize_doner(tmpSJindex);
			int tmpSJanchorSize_acceptor = refineSAMinfo_SE_Nor.returnSJanchorSize_acceptor(tmpSJindex);
			SJanchorSizePairVec_PE.push_back(pair<int,int>(tmpSJanchorSize_doner, tmpSJanchorSize_acceptor));
			SJjumpCodeVecVec_backward_PE.push_back(refineSAMinfo_SE_Nor.returnSJjumpCodeVec_backward(tmpSJindex));
			SJjumpCodeVecVec_forward_PE.push_back(refineSAMinfo_SE_Nor.returnSJjumpCodeVec_forward(tmpSJindex));
		}
		int SE_rcm_SJvecSize = refineSAMinfo_SE_Rcm.returnSJvecSize();
		for(int tmpSJindex = 0; tmpSJindex < SE_rcm_SJvecSize; tmpSJindex++)
		{
			int tmpSJpos_donerEnd = refineSAMinfo_SE_Rcm.returnSJ_donerEndPos(tmpSJindex);
			int tmpSJpos_acceptorStart = refineSAMinfo_SE_Rcm.returnSJ_acceptorStartPos(tmpSJindex);
			int tmpSJanchorSize_doner = refineSAMinfo_SE_Rcm.returnSJanchorSize_doner(tmpSJindex);
			int tmpSJanchorSize_acceptor = refineSAMinfo_SE_Rcm.returnSJanchorSize_acceptor(tmpSJindex);			
			int tmpSJindex_2 = 0;
			bool SJalreadyStoredBool = false;
			for(tmpSJindex_2 = 0; tmpSJindex_2 < SE_nor_SJvecSize; tmpSJindex_2++)
			{
				int tmpSJpos_donerEnd_2 = SJposPairVec_PE[tmpSJindex_2].first;
				int tmpSJpos_acceptorStart_2 = SJposPairVec_PE[tmpSJindex_2].second;
				if((tmpSJpos_donerEnd == tmpSJpos_donerEnd_2)&&(tmpSJpos_acceptorStart == tmpSJpos_acceptorStart_2))
				{
					SJalreadyStoredBool = true;
					break;
				}
			}
			if(SJalreadyStoredBool)
			{
				int tmpSJanchorSize_doner_2 = SJanchorSizePairVec_PE[tmpSJindex_2].first;
				int tmpSJanchorSize_acceptor_2 = SJanchorSizePairVec_PE[tmpSJindex_2].second;
				if(tmpSJanchorSize_doner > tmpSJanchorSize_doner_2)
				{
					SJanchorSizePairVec_PE[tmpSJindex_2].first = tmpSJanchorSize_doner;
					SJjumpCodeVecVec_backward_PE[tmpSJindex_2] = refineSAMinfo_SE_Rcm.returnSJjumpCodeVec_backward(tmpSJindex); 
				}
				if(tmpSJanchorSize_acceptor > tmpSJanchorSize_acceptor_2)
				{
					SJanchorSizePairVec_PE[tmpSJindex_2].second = tmpSJanchorSize_acceptor;
					SJjumpCodeVecVec_forward_PE[tmpSJindex_2] = refineSAMinfo_SE_Rcm.returnSJjumpCodeVec_forward(tmpSJindex);
				}
			}
			else
			{
				SJposPairVec_PE.push_back(pair<int,int>(tmpSJpos_donerEnd, tmpSJpos_acceptorStart));
				SJanchorSizePairVec_PE.push_back(pair<int,int>(tmpSJanchorSize_doner, tmpSJanchorSize_acceptor));
				SJjumpCodeVecVec_backward_PE.push_back(refineSAMinfo_SE_Rcm.returnSJjumpCodeVec_backward(tmpSJindex));
				SJjumpCodeVecVec_forward_PE.push_back(refineSAMinfo_SE_Rcm.returnSJjumpCodeVec_forward(tmpSJindex));
			}
		}
	}

	int returnSJindexInAlignInferJuncHash(int index)
	{
		return SJindexVecInAlignInferJuncHash_PE[index];
	}
};

class RefineSAM_PE_Vec_Info
{
private:
	vector<RefineSAM_PE_Info* > refineSAMvec;
	int toAddRefineSAMindex;
public:
	RefineSAM_PE_Vec_Info()
	{
		toAddRefineSAMindex = 0;
	}

	string returnUnmapSAMstr(int tmpRefineSAMindex)
	{
		return refineSAMvec[tmpRefineSAMindex]->returnUnmapSAMstr();
	}

	string returnOriSAMstr(int tmpRefineSAMindex)
	{
		return refineSAMvec[tmpRefineSAMindex]->returnOriSAMstr();
	}

	bool point2NULLbool(int index)
	{
		if(refineSAMvec[index] == NULL)
			return true;
		else
			return false;
	}

	int returnRefineSAMvecSize()
	{
		return refineSAMvec.size();
	}

	int returnCurrentIndex()
	{
		return toAddRefineSAMindex;
	}
	void toAddRefineSAMindex_incre()
	{
		toAddRefineSAMindex++;
	}

	int returnRefineSAMinfoSJvecSize(int tmpIndex)
	{
		return refineSAMvec[tmpIndex]->returnSJvecSize();
	}	

	int returnIndexInAlignInferJuncHash(int tmpIndexInRefineSAMinfoVec,
		int tmpSJindexInRefineSAMinfo)	
	{
		(refineSAMvec[tmpIndexInRefineSAMinfoVec])->returnSJindexInAlignInferJuncHash(
			tmpSJindexInRefineSAMinfo);
	}

	void pushBackRefineSAMseInfo(RefineSAM_PE_Info* refineSAMseInfo)
	{
		refineSAMvec.push_back(refineSAMseInfo);
		this->toAddRefineSAMindex_incre();
	}	

	void deleteRefineSAMinfo(int index)
	{
		// if(index == 225)
		// {	
		// 	cout << "deleteRefineSAMinfo:" << index << endl;
		// 	cout << "refineSAMinfoStr: " << endl 
		// 		<< refineSAMvec[index]->returnSJindexVecInAlignInferJuncHashStr() << endl;
		// 	cout << "refineSAMinfo: " << endl << refineSAMvec[index]->returnRefineSAMinfoStr() << endl;
		// }
		delete refineSAMvec[index];
		refineSAMvec[index] = NULL;
	}

	void freeMemory()
	{
		for(int tmp = 0; tmp < refineSAMvec.size(); tmp++)
		{
			if(refineSAMvec[tmp] != NULL)
			{
				delete refineSAMvec[tmp];
				refineSAMvec[tmp] = NULL;
			}
		}
	}	
};

#endif