// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SEG_PATH_INFO_H
#define SEG_PATH_INFO_H

#include <stdlib.h>
#include <stdio.h>
#include "index_info.h"
#include "seg_info.h"
#include "seg_candi_path_info.h"
#include "quality_score.h"
using namespace std;

class Seg_Path_Info
{
private:
	int chrNameInt;
	string chrNameStr;
	int chrMapPos;
	int chrMapPos_end;
	vector<Jump_Code> pathJumpCodeVec;
	vector< int > pathMismatchPosVec;
	vector< char > pathMismatchCharVec;
	int startLocInRead;
	int endLocInRead;
	string forOrRcm; // "+" or "-"

	// other alignment Features
	int mismatchNum;
	int insertionLength;
	int deletionLength;
	int mappedLength;
	vector< pair<int,int> > SJposPairVec;
	vector< string > SJflankStringVec;
	string SJstrand;
public:
	Seg_Path_Info()
	{}

	string jumpCodeVec2Str()
	{
		string str;
		for(int tmp = 0; tmp < pathJumpCodeVec.size(); tmp++)
		{	
			str += pathJumpCodeVec[tmp].toString();
		}
		return str;
	}

	int getFlag(bool pairEndReadsBool, bool bothEndMappedBool, bool readUnmappedBool, bool anotherEndReadUnmappedBool,
		bool mappedAsRcmBool, bool anotherEndReadMappedAsRcmBool, bool pairEnd_1_bool, bool pairEnd_2_bool, 
		bool notPrimaryAlignmentBool, bool notPassQualControlBool, bool PCRorOpticalDuplicateBool, bool suppleAlignmentBool)
	{
		int flagInt = 0;
		if(pairEndReadsBool)
			flagInt += 1;
		if(bothEndMappedBool)
			flagInt += 2;
		if(readUnmappedBool)
			flagInt += 4;
		if(anotherEndReadUnmappedBool)
			flagInt += 8;
		if(mappedAsRcmBool)
			flagInt += 16;
		if(anotherEndReadMappedAsRcmBool)
			flagInt += 32;
		if(pairEnd_1_bool)
			flagInt += 64;
		if(pairEnd_2_bool)
			flagInt += 128;
		if(notPrimaryAlignmentBool)
			flagInt += 256;
		if(notPassQualControlBool)
			flagInt += 512;
		if(PCRorOpticalDuplicateBool)
			flagInt += 1024;
		if(suppleAlignmentBool)
			flagInt += 2048;
	
		return flagInt;
	}

	int getFlag_paired_secondaryOrNot(bool mapDirection, 
		bool End1OrEnd2Bool, bool SecondaryOrNot)
	{
		//int flagInt = 4;

		bool pairEndsReadBool = true;
		bool bothEndMappedBool = true;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = false;
		
		bool mappedAsRcmBool, anotherEndReadMappedAsRcmBool;
		if(mapDirection)
		{
			mappedAsRcmBool = false;
			anotherEndReadMappedAsRcmBool = true;
		}
		else
		{
			mappedAsRcmBool = true;
			anotherEndReadMappedAsRcmBool = false;
		}

		bool pairEnd_1_bool, pairEnd_2_bool;
		if(End1OrEnd2Bool)
		{
			pairEnd_1_bool = true;
			pairEnd_2_bool = false;
		}
		else
		{
			pairEnd_1_bool = false;
			pairEnd_2_bool = true;
		}

		bool notPrimaryAlignmentBool = SecondaryOrNot;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;

		int flagInt = this->getFlag(pairEndsReadBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, 
			notPrimaryAlignmentBool, notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool); 

		return flagInt;
	}

	int getMapppingQuality(const string& qualitySeq, string quality_scale, const string& readName)
	{
		int tmpReadLength = qualitySeq.length();
		vector<int> mismatchPosVec_all;
		for(int tmp = 0; tmp < pathMismatchPosVec.size(); tmp++)
		{
			mismatchPosVec_all.push_back(pathMismatchPosVec[tmp]);
		}

		int jumpCodeSize = pathJumpCodeVec.size();
		if(jumpCodeSize > 0)
		{
			if(pathJumpCodeVec[0].type == "S")
			{
				int headSoftClippingLen = pathJumpCodeVec[0].len;
				for(int tmp = 1; tmp <= headSoftClippingLen; tmp++)
				{
					mismatchPosVec_all.push_back(tmp);
				}
			}
			if(pathJumpCodeVec[jumpCodeSize-1].type == "S")
			{
				int tailSoftClippingLen = pathJumpCodeVec[jumpCodeSize-1].len;
				for(int tmp = tmpReadLength - tailSoftClippingLen + 1; tmp <= tmpReadLength; tmp++)
				{
					mismatchPosVec_all.push_back(tmp);
				}
			}
		}

		int mappingQuality = GetQualityScore_new(mismatchPosVec_all, qualitySeq, quality_scale);

		return mappingQuality;
	}

	string getSamFormatString_paired_secondaryOrNot(
		const string& readName, const string& readSeq, 
		const string& qualitySeq,
		Seg_Path_Info* mateSegPathInfo, bool End1OrEnd2, int IH_num, int HI_num, 
		bool SecondaryOrNot, int templateLength, bool fasta_or_fastq_bool)
	{
		//cout << "getSamFormatString_paired_secondaryOrNot ......" << endl;
		int MAPQ;
		if(fasta_or_fastq_bool)
			MAPQ = 255;
		else
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);
		//cout << "MAPQ: " << MAPQ << endl;
		string samString;
		int FLAG;
		string RNAME;
		int POS;
		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		RNAME = chrNameStr;
		POS = chrMapPos;
		//cout << "RNAME: " << RNAME << endl;
		//cout << "POS: " << POS << endl;
		if(forOrRcm == "+")
		{
			FLAG = getFlag_paired_secondaryOrNot(true, End1OrEnd2, SecondaryOrNot);//0;
			TLEN = templateLength;
		}
		else if(forOrRcm == "-")
		{		
			FLAG = getFlag_paired_secondaryOrNot(false, End1OrEnd2, SecondaryOrNot); //16;
			TLEN = 0 - templateLength;
		}
		else
		{}
		//cout << "FLAG: " << FLAG << endl;
		//cout << "TLEN: " << TLEN << endl;
		////////////////////////////////////////////////////////////////////////
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = mateSegPathInfo->returnChrMapPos(); // 0;
		string strandStr = SJstrand; 
		string FLAGstr; 
		string POSstr; 
		string MAPQstr; 
		string PNEXTstr; 
		string TLENstr; 
		string mismatchNumStr;
		string IHstr;
		string HIstr;
		FLAGstr = int_to_str(FLAG);
		POSstr = int_to_str(POS);
		MAPQstr = int_to_str(MAPQ);
		PNEXTstr = int_to_str(PNEXT);
		TLENstr = int_to_str(TLEN);
		mismatchNumStr = int_to_str(mismatchNum);
		IHstr = int_to_str(IH_num);
		HIstr = int_to_str(HI_num);

		samString = readName + "\t" 
			+ FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ qualitySeq 
			+ "\tNM:i:" + mismatchNumStr 
			+ "\tIH:i:" + IHstr
			+ "\tHI:i:" + HIstr
			+ "\tXS:A:" + strandStr 
			+ "\tXF:Z:";
		//cout << "samString_1: " << samString << endl;
		for(int tmp = 0; tmp < SJflankStringVec.size(); tmp++)
			samString = samString + SJflankStringVec[tmp] + ","; 
		//cout << "samString_2: " << samString << endl;
		samString = samString + "\tMD:Z:";
		for(int tmp = 0; tmp < pathMismatchPosVec.size(); tmp++)
			samString = samString + int_to_str(pathMismatchPosVec[tmp]) + ",";
		//cout << "samString_3: " << samString << endl;
		samString = samString + "\tMP:Z:";
		for(int tmp = 0; tmp < pathMismatchCharVec.size(); tmp++)
			samString = samString + pathMismatchCharVec[tmp] + ",";
		//cout << "samString_4: " << samString << endl;
		return samString;
	}

	bool completeOrNot(int readLength)
	{
		if((startLocInRead != 1)||(endLocInRead != readLength))
			return false;
		else
			return true;
	}

	int returnMismatchNum()
	{
		return mismatchNum;
	}

	int returnInsertionLength()
	{
		return insertionLength;
	}

	int returnDeletionLength()
	{
		return deletionLength;
	}

	int returnMappedLength()
	{
		return mappedLength;
	}

	void generateMismatchNum()
	{
		mismatchNum = pathMismatchPosVec.size();
	}

	void generateMappedLength()
	{
		mappedLength = endLocInRead - startLocInRead + 1;
	}

	void generateInDelLen()
	{
		int tmp_ins_len = 0;
		int tmp_del_len = 0;
		for(int tmp = 0; tmp < pathJumpCodeVec.size(); tmp++)
		{
			if(pathJumpCodeVec[tmp].type == "I")
			{
				tmp_ins_len += pathJumpCodeVec[tmp].len;
			}
			if(pathJumpCodeVec[tmp].type == "D")
			{
				tmp_del_len += pathJumpCodeVec[tmp].len;
			}
		}
		insertionLength = tmp_ins_len;
		deletionLength = tmp_del_len;		
	}

	string getStrandFromFlankStringVec()
	{	
		string tmpFlankString;
		string currentStrand = "N";
		for(int tmp = 0; tmp < SJflankStringVec.size(); tmp++)
		{
			tmpFlankString = SJflankStringVec[tmp];
			if((tmpFlankString == "ATAC") || (tmpFlankString == "CTGC") || (tmpFlankString == "CTAC"))
			{
				if(currentStrand == "+")
					return "X";
				else
					currentStrand = "-";
			}
			else
			{
				if(currentStrand == "-")
					return "X";
				else
					return "+";
			}
		}
		return currentStrand;
	}

	void generateSJinfo(Index_Info* indexInfo)
	{
		this->generateSJposVecFromJumpCodeVec(chrMapPos, pathJumpCodeVec, 
			SJposPairVec);
		for(int tmp = 0; tmp < SJposPairVec.size(); tmp++)
		{
			int tmpSJpos_doner = SJposPairVec[tmp].first;
			int tmpSJpos_acceptor = SJposPairVec[tmp].second;
			string tmpSJflankString = indexInfo->returnFlankString(
				chrNameInt, tmpSJpos_doner, tmpSJpos_acceptor);
			SJflankStringVec.push_back(tmpSJflankString);
			SJstrand = getStrandFromFlankStringVec();
		}
	}

	void generateOtherAlignmentFeatures(Index_Info* indexInfo)
	{
		this->generateMismatchNum();
		this->generateMappedLength();
		this->generateInDelLen();
		this->generateSJinfo(indexInfo);
	}

	int returnSpliceJunctionConfidenceLevel()
	{
		if(SJflankStringVec.size() == 0)
			return 0;
		string tmpFlankString;
		int spliceJunctionConfidenceLevel = SPLICE_JUNCTION_CANONICAL_ONLY;
		for(int tmp = 0; tmp < SJflankStringVec.size(); tmp++)
		{
			tmpFlankString = SJflankStringVec[tmp];
			if((tmpFlankString == "GTAG") || (tmpFlankString == "CTAC"))
			{

			}
			else if((tmpFlankString == "ATAC")||(tmpFlankString == "GTAT")
				||(tmpFlankString == "CTGC")||(tmpFlankString == "GCAG"))
			{
				if(spliceJunctionConfidenceLevel == SPLICE_JUNCTION_CANONICAL_ONLY)
					spliceJunctionConfidenceLevel = SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL;
			}
			else
			{
				if(spliceJunctionConfidenceLevel != SPLICE_JUNCTION_NONCANONICAL)
				{
					spliceJunctionConfidenceLevel = SPLICE_JUNCTION_NONCANONICAL;
					return spliceJunctionConfidenceLevel;
				}
			}
		}
		return spliceJunctionConfidenceLevel;		
	}

	int returnChrNameInt()
	{
		return chrNameInt;
	}

	const string& returnChrNameStr()
	{
		return chrNameStr;
	}

	int returnChrMapPos()
	{
		return chrMapPos;
	}

	int returnChrMapPos_end()
	{
		return chrMapPos_end;
	}

	// int returnMismatchNum()
	// {
	// 	return pathMismatchPosVec.size();
	// }

	int returnMismatchPosVecValue(int index)
	{
		return pathMismatchPosVec[index];
	}

	char returnMismatchCharVecValue(int index)
	{
		return pathMismatchCharVec[index];
	}

	vector<Jump_Code>& returnPathJumpCodeVec()
	{
		return pathJumpCodeVec;
	}

	vector<int>& returnPathMismatchPosVec()
	{
		return pathMismatchPosVec;
	}

	vector<char>& returnPathMismatchCharVec()
	{
		return pathMismatchCharVec;
	}

	void finalizeCandiPath(const string& alignDirection,
		Seg_Candi_Path_Info* segCandiPathInfo,
		const string& readSeqInProcess, Index_Info* indexInfo)
	{
		//cout << "finalizeCandiPath starts ...." << endl;
		segCandiPathInfo->fixSegGap_candiPath(
			readSeqInProcess, indexInfo);
		this->getFinalStitchedSegSubCandiPath(
			alignDirection, segCandiPathInfo,
			readSeqInProcess, indexInfo);
		this->generateOtherAlignmentFeatures(indexInfo);
	}

	void getFinalStitchedSegSubCandiPath(const string& alignDirection,
		Seg_Candi_Path_Info* segCandiPathInfo,
		const string& readSeqInProcess, Index_Info* indexInfo)
	{
		int bestSubPathSegIndex_first;
		int bestSubPathSegIndex_last;
		segCandiPathInfo->selectBestSubPath(bestSubPathSegIndex_first, 
			bestSubPathSegIndex_last);
		segCandiPathInfo->doExtensionOnBothDirections(readSeqInProcess,
			indexInfo);
		forOrRcm = alignDirection;
		chrNameInt = segCandiPathInfo->returnChrNameInt();
		chrNameStr = indexInfo->returnChrNameStr(chrNameInt);
		chrMapPos = segCandiPathInfo->returnSegPosInChr(bestSubPathSegIndex_first)
			- segCandiPathInfo->returnExtensionLen_head();
		chrMapPos_end = segCandiPathInfo->returnSegPosInChr(bestSubPathSegIndex_last)
			+ segCandiPathInfo->returnSegLen(bestSubPathSegIndex_last)
			+ segCandiPathInfo->returnExtensionLen_tail() - 1;

		startLocInRead = segCandiPathInfo->returnSegLocInRead(bestSubPathSegIndex_first)
			- segCandiPathInfo->returnExtensionLen_head();
		endLocInRead = segCandiPathInfo->returnSegLocInRead(bestSubPathSegIndex_last)
			+ segCandiPathInfo->returnSegLen(bestSubPathSegIndex_last)
			+ segCandiPathInfo->returnExtensionLen_tail() - 1;
	
		int readLength = readSeqInProcess.length();
		segCandiPathInfo->copySubPathJumpCodeVec2TargetJumpCodeVec_withSoftClip(
			startLocInRead, endLocInRead, readLength,
			bestSubPathSegIndex_first, bestSubPathSegIndex_last,
			pathJumpCodeVec, pathMismatchPosVec, pathMismatchCharVec);
	}

	void generateSJposVecFromJumpCodeVec(int mapChrPos, vector<Jump_Code>& cigarStringJumpCodeVec, 
		vector< pair<int,int> >& SJposPairVec//, vector<int>& SJindexVec_cigarStringJumpCodeVec
		)
	{
		for(int tmp = 0; tmp < cigarStringJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = cigarStringJumpCodeVec[tmp].type;
			if(tmpJumpCodeType == "N")
			{
				int lastJumpCodeIndex = tmp-1;
				int currentJumpCodeIndex = tmp;
				int tmpDonerEndPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, lastJumpCodeIndex);
				int tmpAcceptorStartPos = getEndPosOfSpecificJumpCode(mapChrPos, cigarStringJumpCodeVec, currentJumpCodeIndex) + 1;
				SJposPairVec.push_back(pair<int,int>(tmpDonerEndPos, tmpAcceptorStartPos));
				//SJindexVec_cigarStringJumpCodeVec.push_back(tmp);
			}
		}		
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
};

class Seg_Path_Vec_Info
{
private:
	vector<Seg_Path_Info*> segPathInfoVec_Nor1;
	vector<Seg_Path_Info*> segPathInfoVec_Rcm1;
	vector<Seg_Path_Info*> segPathInfoVec_Nor2;
	vector<Seg_Path_Info*> segPathInfoVec_Rcm2;

	int repeatRegion_index_Nor1;
	int repeatRegion_index_Rcm1;
	int repeatRegion_index_Nor2;
	int repeatRegion_index_Rcm2;

	vector< pair<int,int> > segCandiPathPairVec_Nor1Rcm2;
	vector< pair<int,int> > segCandiPathPairVec_Nor2Rcm1;

	vector< int > bestSegPathInfoVecIndexVec_Nor1Rcm2;
	vector< int > bestSegPathInfoVecIndexVec_Nor2Rcm1;
	double highestPairAlignmentScore;
public:
	Seg_Path_Vec_Info()
	{
		repeatRegion_index_Nor1 = -1;
		repeatRegion_index_Rcm1 = -1;
		repeatRegion_index_Nor2 = -1;
		repeatRegion_index_Rcm2 = -1;	
	}

	Seg_Path_Info* returnSegPathInfo(int type, int index)
	{
		if(type == 1)
		{
			return segPathInfoVec_Nor1[index];
		}
		else if(type == 2)
		{
			return segPathInfoVec_Rcm1[index];
		}
		else if(type == 3)
		{
			return segPathInfoVec_Nor2[index];
		}
		else if(type == 4)
		{
			return segPathInfoVec_Rcm2[index];
		}
		else
		{}	
	}

	int returnSegPathInfoVec(int type)
	{
		if(type == 1)
		{
			return segPathInfoVec_Nor1.size();
		}
		else if(type == 2)
		{
			return segPathInfoVec_Rcm1.size();
		}
		else if(type == 3)
		{
			return segPathInfoVec_Nor2.size();
		}
		else if(type == 4)
		{
			return segPathInfoVec_Rcm2.size();
		}
		else
		{}				
	}

	void pairSegCandiPath(
		Seg_Candi_Path_Vec_Info& segCandiPathVecInfo)
	{
		//cout << "pairSegCandiPath starts ......" << endl;
		// pair Nor1Rcm2
		int segCandPathVecSize_Nor1 
			= segCandiPathVecInfo.returnSegCandiPathInfoVecSize(1);
		int segCandPathVecSize_Rcm2 
			= segCandiPathVecInfo.returnSegCandiPathInfoVecSize(4);
		//cout << "segCandPathVecSize_Nor1: " << segCandPathVecSize_Nor1 << endl;
		//cout << "segCandPathVecSize_Rcm2: " << segCandPathVecSize_Rcm2 << endl;
		vector<int> pairedNor1CandiPathIndex_Rcm2;
		vector<int> pairedNor1CandiPathChrMapPos_Rcm2;
		for(int tmp = 0; tmp < segCandPathVecSize_Rcm2; tmp++)
		{
			pairedNor1CandiPathIndex_Rcm2.push_back(-1);
			pairedNor1CandiPathChrMapPos_Rcm2.push_back(0);
		}
		for(int tmp_Nor1 = 0; tmp_Nor1 < segCandPathVecSize_Nor1; tmp_Nor1++)
		{
			//cout << "tmp_Nor1: " << tmp_Nor1 << endl;
			int tmpBestPairableSegCandiPath_index_Rcm2 = -1;
			int tmpBestPairableSegCandiPath_chrMapPosEnd_Rcm2 = 2147483646;
			for(int tmp_Rcm2 = 0; tmp_Rcm2 < segCandPathVecSize_Rcm2; tmp_Rcm2++)
			{
				int tmpSegCandiPath_chrMapEndPos_Rcm; 
				bool twoSegCandiPathPairableBool 
					= segCandiPathVecInfo.twoSegCandiPathCanBePaired_Nor1Rcm2(
						tmp_Nor1, tmp_Rcm2, tmpSegCandiPath_chrMapEndPos_Rcm);
				if(twoSegCandiPathPairableBool)
				{
					if(tmpSegCandiPath_chrMapEndPos_Rcm < tmpBestPairableSegCandiPath_chrMapPosEnd_Rcm2)
					{
						tmpBestPairableSegCandiPath_index_Rcm2 = tmp_Rcm2;
						tmpBestPairableSegCandiPath_chrMapPosEnd_Rcm2 = tmpSegCandiPath_chrMapEndPos_Rcm;
					}
				}
			}
			//cout << "tmpBestPairableSegCandiPath_index_Rcm2: " << tmpBestPairableSegCandiPath_index_Rcm2 << endl;
			if(tmpBestPairableSegCandiPath_index_Rcm2 >= 0) // found one can be paired Rcm2 for 
			{
				int tmpNor1Path_chrMapPos = segCandiPathVecInfo.returnCandiPathChrMapPos(1, tmp_Nor1);
				int tmpPairedNor1CandiPathIndex 
					= pairedNor1CandiPathIndex_Rcm2[tmpBestPairableSegCandiPath_index_Rcm2];
				if(tmpPairedNor1CandiPathIndex < 0)
				{	
					pairedNor1CandiPathIndex_Rcm2[tmpBestPairableSegCandiPath_index_Rcm2] = tmp_Nor1;
					pairedNor1CandiPathChrMapPos_Rcm2[tmpBestPairableSegCandiPath_index_Rcm2] = tmpNor1Path_chrMapPos;
				}
				else
				{
					int tmpPairedNor1CandiPath_chrMapPos 
						= pairedNor1CandiPathChrMapPos_Rcm2[tmpBestPairableSegCandiPath_index_Rcm2];
					if(tmpNor1Path_chrMapPos > tmpPairedNor1CandiPath_chrMapPos)
					{
						pairedNor1CandiPathIndex_Rcm2[tmpBestPairableSegCandiPath_index_Rcm2] = tmp_Nor1;
						pairedNor1CandiPathChrMapPos_Rcm2[tmpBestPairableSegCandiPath_index_Rcm2] = tmpNor1Path_chrMapPos;						
					}
					else
					{}
				}
				//segCandiPathPairVec_Nor1Rcm2.push_back(pair<int,int>(tmp_Nor1, 
				//	tmpBestPairableSegCandiPath_index_Rcm2));
			}		
		}
		for(int tmp = 0; tmp < segCandPathVecSize_Rcm2; tmp++)
		{
			int tmpPairedNor1CandiPathIndex = pairedNor1CandiPathIndex_Rcm2[tmp];
			if(tmpPairedNor1CandiPathIndex >= 0)
			{
				segCandiPathPairVec_Nor1Rcm2.push_back(pair<int,int>(tmpPairedNor1CandiPathIndex, tmp));
				//return;
			}
		}

		// pair Nor2Rcm1
		int segCandPathVecSize_Nor2 
			= segCandiPathVecInfo.returnSegCandiPathInfoVecSize(3);
		int segCandPathVecSize_Rcm1
			= segCandiPathVecInfo.returnSegCandiPathInfoVecSize(2);
		//cout << "segCandPathVecSize_Nor2: " << segCandPathVecSize_Nor2 << endl;
		//cout << "segCandPathVecSize_Rcm1: " << segCandPathVecSize_Rcm1 << endl;
		vector<int> pairedNor2CandiPathIndex_Rcm1;
		vector<int> pairedNor2CandiPathChrMapPos_Rcm1;
		for(int tmp = 0; tmp < segCandPathVecSize_Rcm1; tmp++)
		{
			pairedNor2CandiPathIndex_Rcm1.push_back(-1);
			pairedNor2CandiPathChrMapPos_Rcm1.push_back(0);
		}
		for(int tmp_Nor2 = 0; tmp_Nor2 < segCandPathVecSize_Nor2; tmp_Nor2++)
		{
			//cout << "tmp_Nor2: " << tmp_Nor2 << endl;
			int tmpBestPairableSegCandiPath_index_Rcm1 = -1;
			int tmpBestPairableSegCandiPath_chrMapPosEnd_Rcm1 = 2147483646;
			for(int tmp_Rcm1 = 0; tmp_Rcm1 < segCandPathVecSize_Rcm1; tmp_Rcm1++)
			{
				int tmpSegCandiPath_chrMapEndPos_Rcm; 
				bool twoSegCandiPathPairableBool 
					= segCandiPathVecInfo.twoSegCandiPathCanBePaired_Nor2Rcm1(
						tmp_Nor2, tmp_Rcm1, tmpSegCandiPath_chrMapEndPos_Rcm);
				if(twoSegCandiPathPairableBool)
				{
					if(tmpSegCandiPath_chrMapEndPos_Rcm < tmpBestPairableSegCandiPath_chrMapPosEnd_Rcm1)
					{
						tmpBestPairableSegCandiPath_index_Rcm1 = tmp_Rcm1;
						tmpBestPairableSegCandiPath_chrMapPosEnd_Rcm1 = tmpSegCandiPath_chrMapEndPos_Rcm;
					}
				}
			}
			//cout << "tmpBestPairableSegCandiPath_index_Rcm1: " << tmpBestPairableSegCandiPath_index_Rcm1 << endl;
			if(tmpBestPairableSegCandiPath_index_Rcm1 >= 0)
			{
				int tmpNor2Path_chrMapPos = segCandiPathVecInfo.returnCandiPathChrMapPos(3, tmp_Nor2);
				int tmpPairedNor2CandiPathIndex 
					= pairedNor2CandiPathIndex_Rcm1[tmpBestPairableSegCandiPath_index_Rcm1];
				if(tmpPairedNor2CandiPathIndex < 0)
				{	
					pairedNor2CandiPathIndex_Rcm1[tmpBestPairableSegCandiPath_index_Rcm1] = tmp_Nor2;
					pairedNor2CandiPathChrMapPos_Rcm1[tmpBestPairableSegCandiPath_index_Rcm1] = tmpNor2Path_chrMapPos;
				}
				else
				{
					int tmpPairedNor2CandiPath_chrMapPos 
						= pairedNor2CandiPathChrMapPos_Rcm1[tmpBestPairableSegCandiPath_index_Rcm1];
					if(tmpNor2Path_chrMapPos > tmpPairedNor2CandiPath_chrMapPos)
					{
						pairedNor2CandiPathIndex_Rcm1[tmpBestPairableSegCandiPath_index_Rcm1] = tmp_Nor2;
						pairedNor2CandiPathChrMapPos_Rcm1[tmpBestPairableSegCandiPath_index_Rcm1] = tmpNor2Path_chrMapPos;						
					}
					else
					{}
				}
				//segCandiPathPairVec_Nor2Rcm1.push_back(pair<int,int>(tmp_Nor2, 
				//	tmpBestPairableSegCandiPath_index_Rcm1));			
			}
		}
		for(int tmp = 0; tmp < segCandPathVecSize_Rcm1; tmp++)
		{
			int tmpPairedNor2CandiPathIndex = pairedNor2CandiPathIndex_Rcm1[tmp];
			if(tmpPairedNor2CandiPathIndex >= 0)
			{	
				segCandiPathPairVec_Nor2Rcm1.push_back(pair<int,int>(tmpPairedNor2CandiPathIndex, tmp));
				//return;
			}
		}
	}

	void finalizeCandiPathVec_Paired(Seg_Candi_Path_Vec_Info& segCandiPathVecInfo,
		PE_Read_Info& peReadInfo, Index_Info* indexInfo)
	{
		//cout << "finalizeCandiPathVec_Paired starts ....." << endl;
		int segCandiPathPairVec_Nor1Rcm2_size = segCandiPathPairVec_Nor1Rcm2.size();
		int segCandiPathPairVec_Nor2Rcm1_size = segCandiPathPairVec_Nor2Rcm1.size();
		//cout << "segCandiPathPairVec_Nor1Rcm2_size: " << segCandiPathPairVec_Nor1Rcm2_size << endl;
		//cout << "segCandiPathPairVec_Nor2Rcm1_size: " << segCandiPathPairVec_Nor2Rcm1_size << endl;
		for(int tmp = 0; tmp < segCandiPathPairVec_Nor1Rcm2_size; tmp++)
		{
			//cout << "segCandiPathPairVec_Nor1Rcm2 tmp: " << tmp << endl; 
			int segCandiPathPairVec_Nor1_index = segCandiPathPairVec_Nor1Rcm2[tmp].first;
			int segCandiPathPairVec_Rcm2_index = segCandiPathPairVec_Nor1Rcm2[tmp].second;
			//cout << "segCandiPathPairVec_Nor1_index: " << segCandiPathPairVec_Nor1_index << endl;
			//cout << "segCandiPathPairVec_Rcm2_index: " << segCandiPathPairVec_Rcm2_index << endl;
			Seg_Path_Info* tmpSegPathInfo_Nor1 = new Seg_Path_Info();
			tmpSegPathInfo_Nor1->finalizeCandiPath("+",
				segCandiPathVecInfo.returnCandiPath(segCandiPathPairVec_Nor1_index, 1),
				peReadInfo.returnReadSeq_1(), indexInfo);
			Seg_Path_Info* tmpSegPathInfo_Rcm2 = new Seg_Path_Info();
			tmpSegPathInfo_Rcm2->finalizeCandiPath("-",
				segCandiPathVecInfo.returnCandiPath(segCandiPathPairVec_Rcm2_index, 4),
				peReadInfo.returnRcmReadSeq_2(), indexInfo);
			segPathInfoVec_Nor1.push_back(tmpSegPathInfo_Nor1);
			segPathInfoVec_Rcm2.push_back(tmpSegPathInfo_Rcm2);
		}
		for(int tmp = 0; tmp < segCandiPathPairVec_Nor2Rcm1_size; tmp++)
		{
			//cout << "segCandiPathPairVec_Nor2Rcm1 tmp: " << tmp << endl; 
			int segCandiPathPairVec_Nor2_index = segCandiPathPairVec_Nor2Rcm1[tmp].first;
			int segCandiPathPairVec_Rcm1_index = segCandiPathPairVec_Nor2Rcm1[tmp].second;
			//cout << "segCandiPathPairVec_Nor2_index: " << segCandiPathPairVec_Nor2_index << endl;
			//cout << "segCandiPathPairVec_Rcm1_index: " << segCandiPathPairVec_Rcm1_index << endl;
			Seg_Path_Info* tmpSegPathInfo_Nor2 = new Seg_Path_Info();
			tmpSegPathInfo_Nor2->finalizeCandiPath("+",
				segCandiPathVecInfo.returnCandiPath(segCandiPathPairVec_Nor2_index, 3),
				peReadInfo.returnReadSeq_2(), indexInfo);
			Seg_Path_Info* tmpSegPathInfo_Rcm1 = new Seg_Path_Info();
			tmpSegPathInfo_Rcm1->finalizeCandiPath("-",
				segCandiPathVecInfo.returnCandiPath(segCandiPathPairVec_Rcm1_index, 2),
				peReadInfo.returnRcmReadSeq_1(), indexInfo);
			segPathInfoVec_Nor2.push_back(tmpSegPathInfo_Nor2);
			segPathInfoVec_Rcm1.push_back(tmpSegPathInfo_Rcm1);
		}		
	}

	void finalizeCandiPathVec_Unpaired(Seg_Candi_Path_Vec_Info& segCandiPathVecInfo,
		PE_Read_Info& peReadInfo, Index_Info* indexInfo)
	{
		int segCandiPathVec_Nor1_size = segCandiPathVecInfo.returnSegCandiPathInfoVecSize(1);
		int segCandiPathVec_Rcm1_size = segCandiPathVecInfo.returnSegCandiPathInfoVecSize(2);
		int segCandiPathVec_Nor2_size = segCandiPathVecInfo.returnSegCandiPathInfoVecSize(3);
		int segCandiPathVec_Rcm2_size = segCandiPathVecInfo.returnSegCandiPathInfoVecSize(4);
		for(int tmp = 0; tmp < segCandiPathVec_Nor1_size; tmp++)
		{
			Seg_Path_Info* tmpSegPathInfo_Nor1 = new Seg_Path_Info();
			tmpSegPathInfo_Nor1->finalizeCandiPath("+",
				segCandiPathVecInfo.returnCandiPath(tmp,1),
				peReadInfo.returnReadSeq_1(), indexInfo);
			segPathInfoVec_Nor1.push_back(tmpSegPathInfo_Nor1);
		}
		for(int tmp = 0; tmp < segCandiPathVec_Rcm1_size; tmp++)
		{
			Seg_Path_Info* tmpSegPathInfo_Rcm1 = new Seg_Path_Info();
			tmpSegPathInfo_Rcm1->finalizeCandiPath("-",
				segCandiPathVecInfo.returnCandiPath(tmp,2),
				peReadInfo.returnRcmReadSeq_1(), indexInfo);
			segPathInfoVec_Rcm1.push_back(tmpSegPathInfo_Rcm1);
		}
		for(int tmp = 0; tmp < segCandiPathVec_Nor2_size; tmp++)
		{
			Seg_Path_Info* tmpSegPathInfo_Nor2 = new Seg_Path_Info();
			tmpSegPathInfo_Nor2->finalizeCandiPath("+",
				segCandiPathVecInfo.returnCandiPath(tmp,3),
				peReadInfo.returnReadSeq_2(), indexInfo);
			segPathInfoVec_Nor2.push_back(tmpSegPathInfo_Nor2);
		}
		for(int tmp = 0; tmp < segCandiPathVec_Rcm2_size; tmp++)
		{
			Seg_Path_Info* tmpSegPathInfo_Rcm2 = new Seg_Path_Info();
			tmpSegPathInfo_Rcm2->finalizeCandiPath("-",
				segCandiPathVecInfo.returnCandiPath(tmp,4),
				peReadInfo.returnRcmReadSeq_2(), indexInfo);
			segPathInfoVec_Rcm2.push_back(tmpSegPathInfo_Rcm2);
		}				
	}

	void finalizeCandiPathVec(Seg_Candi_Path_Vec_Info& segCandiPathVecInfo,
		PE_Read_Info& peReadInfo, Index_Info* indexInfo)
	{
		//cout << "finalizeCandiPathPathVec starts ......" << endl;
		//cout << "start to do pairSegCandiPath ......" << endl;
		repeatRegion_index_Nor1 = segCandiPathVecInfo.returnRepeatRegion_index_Nor1();
		repeatRegion_index_Rcm1 = segCandiPathVecInfo.returnRepeatRegion_index_Rcm1();
		repeatRegion_index_Nor2 = segCandiPathVecInfo.returnRepeatRegion_index_Nor2();
		repeatRegion_index_Rcm2 = segCandiPathVecInfo.returnRepeatRegion_index_Rcm2();
		#ifdef CAL_TIME
		pairSegCandiPath_begin = clock();
		#endif
		this->pairSegCandiPath(segCandiPathVecInfo);
		#ifdef CAL_TIME
		pairSegCandiPath_end = clock();
		pairSegCandiPath_cost = pairSegCandiPath_cost + pairSegCandiPath_end - pairSegCandiPath_begin;
		finalizeCandiPathVec_begin = clock();
		#endif
		if(segCandiPathPairVec_Nor1Rcm2.size() + segCandiPathPairVec_Nor2Rcm1.size() > 0)
		{
			//cout << "pairs exist ......" << endl;
			this->finalizeCandiPathVec_Paired(segCandiPathVecInfo,
				peReadInfo, indexInfo);
		}
		else
		{
			//cout << "no pair exists ....." << endl;
			this->finalizeCandiPathVec_Unpaired(segCandiPathVecInfo,
				peReadInfo, indexInfo);
		}
		#ifdef CAL_TIME
		finalizeCandiPathVec_end = clock();
		finalizeCandiPathVec_cost = finalizeCandiPathVec_cost + finalizeCandiPathVec_end - finalizeCandiPathVec_begin;
		#endif	
	}

	void alignmentFilter_fixPhase1(int readLength_1, int readLength_2)
	{
		//cout << "start to do alignmentFilter_fixPhase1 ...." << endl;
		if(segCandiPathPairVec_Nor1Rcm2.size() + segCandiPathPairVec_Nor2Rcm1.size() == 0)
		{
			//cout << "one End unmapped ...." << endl;
			return;
		}
		//cout << "paired ...." << endl;
		double min_map_score_keptAlignmentPair 
			= ((readLength_1 + readLength_2)/100) * Min_Map_Score_keptAlignmentPair_PerHundredBases;
		//cout << "min_map_score_keptAlignmentPair: " << min_map_score_keptAlignmentPair << endl;
		vector< double > alignmentPairScore_Nor1Rcm2;
		vector< double > alignmentPairScore_Nor2Rcm1;
		//cout << "start to do getAlignmentPairScoreForEachAlignmentPair ...." << endl;
		this->getAlignmentPairScoreForEachAlignmentPair(
			alignmentPairScore_Nor1Rcm2, alignmentPairScore_Nor2Rcm1);
		//cout << "end of getAlignmentPairScoreForEachAlignmentPair ....." << endl;
		//cout << "alignmentPairScore_Nor1Rcm2.size(): " << alignmentPairScore_Nor1Rcm2.size() << endl;
		for(int tmp = 0; tmp < alignmentPairScore_Nor1Rcm2.size(); tmp++)
		{
			double tmpScore = alignmentPairScore_Nor1Rcm2[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				//cout << "push_back: tmp: " << tmp << endl;
				bestSegPathInfoVecIndexVec_Nor1Rcm2.push_back(tmp);
			}
		}
		//cout << "alignmentPairScore_Nor2Rcm1.size(): " << alignmentPairScore_Nor2Rcm1.size() << endl;
		for(int tmp = 0; tmp < alignmentPairScore_Nor2Rcm1.size(); tmp++)
		{
			double tmpScore = alignmentPairScore_Nor2Rcm1[tmp];
			if(tmpScore > min_map_score_keptAlignmentPair)
			{
				if(tmpScore > highestPairAlignmentScore)
					highestPairAlignmentScore = tmpScore;
				//cout << "push_back: tmp: " << tmp << endl;
				bestSegPathInfoVecIndexVec_Nor2Rcm1.push_back(tmp);
			}
		}		
	}

	void getAlignmentPairScoreForEachAlignmentPair(
		vector< double >& alignmentPairScore_Nor1Rcm2,
		vector< double >& alignmentPairScore_Nor2Rcm1)
	{
		int mismatch_weight = MISMATCH_PENALTY;
		int deletion_weight = DELETION_PENALTY;
		int insertion_weight = INSERTION_PENALTY;

		int semiCanonicalSJ_penalty = SEMICANONICALSJ_PENALTY;
		int nonCanonicalSJ_penalty = NONCANONICALSJ_PENALTY;
	
		int distant_insert_size_min = DISTANT_INSERT_SIZE_MIN;
		int distant_insert_size_penalty = DISTANT_INSERT_SIZE_PENALTY;

		for(int tmp = 0; tmp < segCandiPathPairVec_Nor1Rcm2.size(); tmp++)
		{
			int tmpMappedLength = segPathInfoVec_Nor1[tmp]->returnMappedLength()
				+ segPathInfoVec_Rcm2[tmp]->returnMappedLength();
			int tmpMismatchNum = segPathInfoVec_Nor1[tmp]->returnMismatchNum()
				+ segPathInfoVec_Rcm2[tmp]->returnMismatchNum();			
			int tmpInsertionLength = segPathInfoVec_Nor1[tmp]->returnInsertionLength()
				+ segPathInfoVec_Rcm2[tmp]->returnInsertionLength();
			int tmpDeletionLength = segPathInfoVec_Nor1[tmp]->returnDeletionLength()
				+ segPathInfoVec_Rcm2[tmp]->returnDeletionLength();
			int tmpPairDistance = segPathInfoVec_Rcm2[tmp]->returnChrMapPos()
				- segPathInfoVec_Nor1[tmp]->returnChrMapPos_end();
			double tmpPairDistance_score;
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			int tmpPairAlignmentSJconfidenceLevel = getPairedAlignmentSJconfidenceLevel(
				segPathInfoVec_Nor1[tmp]->returnSpliceJunctionConfidenceLevel(), 
				segPathInfoVec_Rcm2[tmp]->returnSpliceJunctionConfidenceLevel());
			int SJ_penalty = this->getSJpenalty(tmpPairAlignmentSJconfidenceLevel,
				semiCanonicalSJ_penalty, nonCanonicalSJ_penalty);
			int mappingAreaSizePenalty = this->getMappingAreaSizePenalty(
				segPathInfoVec_Rcm2[tmp]->returnChrMapPos_end() 
					- segPathInfoVec_Nor1[tmp]->returnChrMapPos(), 
				distant_insert_size_min, distant_insert_size_penalty);

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score - SJ_penalty
				- mappingAreaSizePenalty;
			alignmentPairScore_Nor1Rcm2.push_back(tmpScore);	
		}	
		
		for(int tmp = 0; tmp < segCandiPathPairVec_Nor2Rcm1.size(); tmp++)
		{
			int tmpMappedLength = segPathInfoVec_Nor2[tmp]->returnMappedLength()
				+ segPathInfoVec_Rcm1[tmp]->returnMappedLength();
			int tmpMismatchNum = segPathInfoVec_Nor2[tmp]->returnMismatchNum()
				+ segPathInfoVec_Rcm1[tmp]->returnMismatchNum();			
			int tmpInsertionLength = segPathInfoVec_Nor2[tmp]->returnInsertionLength()
				+ segPathInfoVec_Rcm1[tmp]->returnInsertionLength();
			int tmpDeletionLength = segPathInfoVec_Nor2[tmp]->returnDeletionLength()
				+ segPathInfoVec_Rcm1[tmp]->returnDeletionLength();
			int tmpPairDistance = segPathInfoVec_Rcm1[tmp]->returnChrMapPos()
				- segPathInfoVec_Nor2[tmp]->returnChrMapPos_end();
			double tmpPairDistance_score;
			if(tmpPairDistance >= 500)
				tmpPairDistance_score = 0.1;
			else if(tmpPairDistance <= 0)
				tmpPairDistance_score = 0;
			else
				tmpPairDistance_score = (double)tmpPairDistance/5000;
			int tmpPairAlignmentSJconfidenceLevel = getPairedAlignmentSJconfidenceLevel(
				segPathInfoVec_Nor2[tmp]->returnSpliceJunctionConfidenceLevel(), 
				segPathInfoVec_Rcm1[tmp]->returnSpliceJunctionConfidenceLevel());
			int SJ_penalty = this->getSJpenalty(tmpPairAlignmentSJconfidenceLevel,
				semiCanonicalSJ_penalty, nonCanonicalSJ_penalty);
			int mappingAreaSizePenalty = this->getMappingAreaSizePenalty(
				segPathInfoVec_Rcm1[tmp]->returnChrMapPos_end() 
					- segPathInfoVec_Nor2[tmp]->returnChrMapPos(), 
				distant_insert_size_min, distant_insert_size_penalty);

			double tmpScore = tmpMappedLength 
		        - tmpMismatchNum * mismatch_weight 
				- tmpInsertionLength * insertion_weight 
				- tmpDeletionLength * deletion_weight
				- tmpPairDistance_score - SJ_penalty
				- mappingAreaSizePenalty;
			alignmentPairScore_Nor2Rcm1.push_back(tmpScore);	
		}	
	}

	int getMappingAreaSizePenalty(int mappingAreaSize,
		int distant_insert_size_min, int distant_insert_size_penalty)
	{
		if(mappingAreaSize >= distant_insert_size_min)
			return distant_insert_size_penalty;
		else
			return 0;
	}

	int getSJpenalty(int tmpPairAlignmentSJconfidenceLevel, int semiCanonicalSJ_penalty,
		int nonCanonicalSJ_penalty)
	{
		int SJ_penalty = 0;
		if(tmpPairAlignmentSJconfidenceLevel <= SPLICE_JUNCTION_CANONICAL_ONLY) // canonical or no SJ
		{}
		else if(tmpPairAlignmentSJconfidenceLevel == SPLICE_JUNCTION_SEMICANONICAL_NO_NONCANONICAL) // semiCanonical
			SJ_penalty = semiCanonicalSJ_penalty;
		else  // noncanonical
			SJ_penalty = nonCanonicalSJ_penalty;

		return SJ_penalty;				
	}

	int getPairedAlignmentSJconfidenceLevel(int confidenceLevel_1, int confidenceLevel_2)
	{
		if(confidenceLevel_1 > confidenceLevel_2)
		{
			return confidenceLevel_1;
		}
		else 
			return confidenceLevel_2;
	}	

	bool mappedToRepeatRegionBool()
	{
		if( (repeatRegion_index_Nor1 > 0) || (repeatRegion_index_Rcm1 > 0) 
			|| (repeatRegion_index_Nor2 > 0) || (repeatRegion_index_Rcm2 > 0) )
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	int returnAlignmentScoreMinOutput_PE(PE_Read_Info& peReadInfo)
	{
		int peReadLength = peReadInfo.returnReadLength_end1() + peReadInfo.returnReadLength_end2();
		if(peReadLength < 100)
			return (int)(0.5 * peReadLength);
		else
			return 50;
	}

	bool alignPairScoreTooLow_bool(PE_Read_Info& peReadInfo)
	{
		//return false;
		int alignment_score_min_output = this->returnAlignmentScoreMinOutput_PE(peReadInfo);
		if(highestPairAlignmentScore < alignment_score_min_output)
			return true;
		else
			return false;
	}	

	bool alignInfoExistsBool()
	{
		int alignInfoNum = segPathInfoVec_Nor1.size() + segPathInfoVec_Rcm1.size()
			+ segPathInfoVec_Nor2.size() + segPathInfoVec_Rcm2.size();
		if(alignInfoNum == 0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}



	bool checkUniqueOrMulti()
	{
		if(bestSegPathInfoVecIndexVec_Nor1Rcm2.size() 
			+ bestSegPathInfoVecIndexVec_Nor2Rcm1.size() > 1)
			return false;
		else
			return true;
	}



	bool allAlignmentInFinalPairCompleted(int readLength_1, int readLength_2)
	{
		for(int tmp = 0; tmp < bestSegPathInfoVecIndexVec_Nor1Rcm2.size();
			tmp ++)
		{
			int indexInSegPathInfoVec_Nor1Rcm2 = bestSegPathInfoVecIndexVec_Nor1Rcm2[tmp];
			bool tmpNor1SegPathComplete_bool 
				= segPathInfoVec_Nor1[indexInSegPathInfoVec_Nor1Rcm2]->completeOrNot(readLength_1);
			if(!tmpNor1SegPathComplete_bool)
				return false;
			bool tmpRcm2SegPathComplete_bool 
				= segPathInfoVec_Rcm2[indexInSegPathInfoVec_Nor1Rcm2]->completeOrNot(readLength_2);			
			if(!tmpRcm2SegPathComplete_bool)
				return false;
		}
		for(int tmp = 0; tmp < bestSegPathInfoVecIndexVec_Nor2Rcm1.size();
			tmp ++)
		{
			int indexInSegPathInfoVec_Nor2Rcm1 = bestSegPathInfoVecIndexVec_Nor2Rcm1[tmp];
			bool tmpNor2SegPathComplete_bool 
				= segPathInfoVec_Nor2[indexInSegPathInfoVec_Nor2Rcm1]->completeOrNot(readLength_2);
			if(!tmpNor2SegPathComplete_bool)
				return false;
			bool tmpRcm1SegPathComplete_bool 
				= segPathInfoVec_Rcm1[indexInSegPathInfoVec_Nor2Rcm1]->completeOrNot(readLength_1);			
			if(!tmpRcm1SegPathComplete_bool)
				return false;			
		}
		return true;
	}

	bool finalPairExistsBool()
	{
		if(bestSegPathInfoVecIndexVec_Nor1Rcm2.size() + bestSegPathInfoVecIndexVec_Nor2Rcm1.size() > 0)
			return true;
		else
			return false;
	}

	string getSAMformatForUnpairedAlignments_secondaryOrNot(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
	{
		return this->getSAMformatForBothEndsUnmapped(peReadInfo, FastaOrFastq);
	}

	string getSAMformatForBothEndsUnmapped(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
	{
		string qualitySeq_1; //= peReadInfo.returnReadQual_1();
		string qualitySeq_2;
		if(FastaOrFastq)
		{
			qualitySeq_1 = "*";
			qualitySeq_2 = "*";
		}
		else
		{
			qualitySeq_1 = peReadInfo.returnReadQual_1();
			qualitySeq_2 = peReadInfo.returnReadQual_2();
		}
 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();		

		string peAlignSamStr;
		peAlignSamStr = readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_1() 
			+ "\t" + qualitySeq_1 + "\tIH:i:0\tHI:i:0\n" +
			readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + peReadInfo.returnReadSeq_2() 
			+ "\t" + qualitySeq_2 + "\tIH:i:0\tHI:i:0";
		return peAlignSamStr;
	}	

	string getSAMformatForBothEndsUnmapped_mappedToRepeatRegionReads(
		PE_Read_Info& peReadInfo, RepeatRegion_Info* repeatRegionInfo, 
		Index_Info* indexInfo, bool FastaOrFastq)
	{
		string qualitySeq_1; //= peReadInfo.returnReadQual_1();
		string qualitySeq_2;
		if(FastaOrFastq)
		{
			qualitySeq_1 = "*";
			qualitySeq_2 = "*";
		}
		else
		{
			qualitySeq_1 = peReadInfo.returnReadQual_1();
			qualitySeq_2 = peReadInfo.returnReadQual_2();
		}
		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();

		bool mappedToRepeatRegionBool_Nor1 = (repeatRegion_index_Nor1 >= 0);
		bool mappedToRepeatRegionBool_Rcm1 = (repeatRegion_index_Rcm1 >= 0);
		bool mappedToRepeatRegionBool_Nor2 = (repeatRegion_index_Nor2 >= 0);
		bool mappedToRepeatRegionBool_Rcm2 = (repeatRegion_index_Rcm2 >= 0);

		string peAlignSamStr;

		if((!mappedToRepeatRegionBool_Nor1)&&(!mappedToRepeatRegionBool_Rcm1))
		{
			peAlignSamStr = peAlignSamStr + readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_1() + "\t" + qualitySeq_1 + "\tIH:i:0\tHI:i:0\n";			
		}
		else
		{
			peAlignSamStr = peAlignSamStr + readName_1 + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_1() + "\t" + qualitySeq_1 + "\tIH:i:0\tHI:i:0";
			if(repeatRegion_index_Nor1 >= 0)
			{
				peAlignSamStr =  peAlignSamStr + "\tRf:i:" + int_to_str(repeatRegion_index_Nor1);
			}
			if(repeatRegion_index_Rcm1 >= 0)
			{
				peAlignSamStr = peAlignSamStr + "\tRr:i:" + int_to_str(repeatRegion_index_Rcm1);
			}
			peAlignSamStr += "\n";
		}

		if((repeatRegion_index_Nor2 < 0)&&(repeatRegion_index_Rcm2 < 0))
		{
			peAlignSamStr = peAlignSamStr + readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_2() + "\t" + qualitySeq_2 + "\tIH:i:0\tHI:i:0\n";		
		}
		else
		{
			peAlignSamStr = peAlignSamStr + readName_2 + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + 
				peReadInfo.returnReadSeq_2() + "\t" + qualitySeq_2 + "\tIH:i:0\tHI:i:0";

			if(repeatRegion_index_Nor2 >= 0)
			{
				peAlignSamStr = peAlignSamStr + "\tRf:i:" + int_to_str(repeatRegion_index_Nor2);
			}
			if(repeatRegion_index_Rcm2 >= 0)
			{
				peAlignSamStr = peAlignSamStr + "\tRr:i:" + int_to_str(repeatRegion_index_Rcm2);
			}
			peAlignSamStr += "\n";
		}		

		return peAlignSamStr.substr(0,peAlignSamStr.length()-1);
	}


	string getSAMformatForFinalPair_secondaryOrNot(PE_Read_Info& peReadInfo,
		bool FastaOrFastq)
	{
		//cout << "getSAMformatForFinalPair_secondaryOrNot( starts ......." << endl;
		//return "paired";
 		string readName_1 = peReadInfo.returnReadName_beforeSlash_1();
 		string readName_2 = peReadInfo.returnReadName_beforeSlash_2();
		string tmpSamStr;

		int IH_Nor1Rcm2 = bestSegPathInfoVecIndexVec_Nor1Rcm2.size();
		int IH_Nor2Rcm1 = bestSegPathInfoVecIndexVec_Nor2Rcm1.size();

		int IH_allPair = IH_Nor1Rcm2 + IH_Nor2Rcm1;
		//cout << "IH_Nor1Rcm2: " << IH_Nor1Rcm2 << endl;
		for(int tmp = 0; tmp < bestSegPathInfoVecIndexVec_Nor1Rcm2.size(); tmp++)
		{
			int HI_Nor1Rcm2_tmp = tmp + 1;
			int tmpSegPath_index = bestSegPathInfoVecIndexVec_Nor1Rcm2[tmp];
			int template_start = segPathInfoVec_Nor1[tmpSegPath_index]->returnChrMapPos();
			int template_end = segPathInfoVec_Rcm2[tmpSegPath_index]->returnChrMapPos_end();

			int template_length = template_end - template_start + 1;

			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_1, readSeq_1);
				= tmpSamStr 
					+ segPathInfoVec_Nor1[tmpSegPath_index]->getSamFormatString_paired_secondaryOrNot(
						readName_1, peReadInfo.returnReadSeq_1(), peReadInfo.returnQualitySeq_1(), 
						segPathInfoVec_Rcm2[tmpSegPath_index], true, IH_allPair, HI_Nor1Rcm2_tmp, 
						(tmp != 0), template_length, FastaOrFastq);

			//outputFile << tmpSamStr << endl; 
			tmpSamStr += "\n";

			tmpSamStr 
				= tmpSamStr 
					+ segPathInfoVec_Rcm2[tmpSegPath_index]->getSamFormatString_paired_secondaryOrNot(
						readName_2, peReadInfo.returnRcmReadSeq_2(), peReadInfo.returnRcmQualitySeq_2(),
						segPathInfoVec_Nor1[tmpSegPath_index], false, IH_allPair, HI_Nor1Rcm2_tmp, (tmp != 0), template_length, FastaOrFastq);

			tmpSamStr += "\n";
		}
		//cout << "IH_Nor2Rcm1: " << IH_Nor2Rcm1 << endl;
		for(int tmp = 0; tmp < bestSegPathInfoVecIndexVec_Nor2Rcm1.size(); tmp++)
		{
			//cout << "tmp: " << tmp << endl;
			int HI_Nor2Rcm1_tmp = tmp + 1 + IH_Nor1Rcm2;
			int tmpSegPath_index = bestSegPathInfoVecIndexVec_Nor2Rcm1[tmp];
			//cout << "tmpSegPath_index: " << tmpSegPath_index << endl;
			// int tmpNor2NO = finalAlignPair_Nor2Rcm1[tmp].first;
			// int tmpRcm1NO = finalAlignPair_Nor2Rcm1[tmp].second;
			// Alignment_Info* tmpAlignInfo_1 = norAlignmentInfo_PE_2[tmpNor2NO];		
			// Alignment_Info* tmpAlignInfo_2 = rcmAlignmentInfo_PE_1[tmpRcm1NO];
			
			int template_start = segPathInfoVec_Nor2[tmpSegPath_index]->returnChrMapPos();
			int template_end = segPathInfoVec_Rcm1[tmpSegPath_index]->returnChrMapPos_end();
			//cout << "template_start: " << template_start << endl;
			//cout << "template_end: " << template_end << endl;
			int template_length = template_end - template_start + 1;		
			//cout << "template_length: " << template_length << endl;
			tmpSamStr 
				//= tmpAlignInfo_1 -> getSamFormatString(readName_2, readSeq_2);
				= tmpSamStr 
					+ segPathInfoVec_Nor2[tmpSegPath_index]->getSamFormatString_paired_secondaryOrNot(
						readName_2, peReadInfo.returnReadSeq_2(), peReadInfo.returnQualitySeq_2(),
						segPathInfoVec_Rcm1[tmpSegPath_index], false, IH_allPair, HI_Nor2Rcm1_tmp, 
						((tmp != 0)||(IH_Nor1Rcm2 > 0)), 
						template_length, FastaOrFastq);
			tmpSamStr += "\n";
			tmpSamStr 
				= tmpSamStr 
					+ segPathInfoVec_Rcm1[tmpSegPath_index]->getSamFormatString_paired_secondaryOrNot(
						readName_1, peReadInfo.returnRcmReadSeq_1(), peReadInfo.returnRcmQualitySeq_1(),
						segPathInfoVec_Nor2[tmpSegPath_index], true, IH_allPair, HI_Nor2Rcm1_tmp, 
						((tmp != 0)||(IH_Nor1Rcm2 > 0)), 
						template_length, FastaOrFastq);
			tmpSamStr += "\n";
		}
		string returnStr = tmpSamStr.substr(0, tmpSamStr.length()-1);
		//cout << "returnStr: " << returnStr << endl;
		return returnStr;
	}	

	string getTmpAlignInfo(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, bool FastaOrFastq)
	{

		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(segPathInfoVec_Nor1.size()) + "\t"
			+ int_to_str(segPathInfoVec_Rcm1.size()) + "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(segPathInfoVec_Nor2.size()) + "\t"
			+ int_to_str(segPathInfoVec_Rcm2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		string tmpNorAlignInfo_1;
		for(int tmp = 0; tmp < segPathInfoVec_Nor1.size(); tmp++)
		{
			if(segPathInfoVec_Nor1.size() >= 40)
				break;
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ segPathInfoVec_Nor1[tmp]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Nor1[tmp]->returnChrMapPos()) + ","
				+ segPathInfoVec_Nor1[tmp]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Nor1[tmp]->returnMismatchNum()) + ","; //+ "\t";

			int tmpMismatchNum = segPathInfoVec_Nor1[tmp]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + int_to_str(segPathInfoVec_Nor1[tmp]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + int_to_str(segPathInfoVec_Nor1[tmp]->returnMismatchPosVecValue(tmp2));
			tmpNorAlignInfo_1 += ",";

			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + (segPathInfoVec_Nor1[tmp]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" + (segPathInfoVec_Nor1[tmp]->returnMismatchCharVecValue(tmp3));
			tmpNorAlignInfo_1 += ",";

			tmpNorAlignInfo_1 += "\t";
		}
		string tmpRcmAlignInfo_1;
		for(int tmp = 0; tmp < segPathInfoVec_Rcm1.size(); tmp++)
		{
			if(segPathInfoVec_Rcm1.size() >= 40)
				break;
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ segPathInfoVec_Rcm1[tmp]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Rcm1[tmp]->returnChrMapPos()) + ","
				+ segPathInfoVec_Rcm1[tmp]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Rcm1[tmp]->returnMismatchNum()) + ",";//\t";

			int tmpMismatchNum = segPathInfoVec_Rcm1[tmp]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(segPathInfoVec_Rcm1[tmp]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(segPathInfoVec_Rcm1[tmp]->returnMismatchPosVecValue(tmp2));
			tmpRcmAlignInfo_1 += ",";

			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + (segPathInfoVec_Rcm1[tmp]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + (segPathInfoVec_Rcm1[tmp]->returnMismatchCharVecValue(tmp3));
			tmpRcmAlignInfo_1 += ",";

			tmpRcmAlignInfo_1 += "\t";
		}
		string tmpNorAlignInfo_2;
		for(int tmp = 0; tmp < segPathInfoVec_Nor2.size(); tmp++)
		{
			if(segPathInfoVec_Nor2.size() >= 40)
				break;
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ segPathInfoVec_Nor2[tmp]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Nor2[tmp]->returnChrMapPos()) + ","
				+ segPathInfoVec_Nor2[tmp]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Nor2[tmp]->returnMismatchNum()) + ",";//\t";

			int tmpMismatchNum = segPathInfoVec_Nor2[tmp]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + int_to_str(segPathInfoVec_Nor2[tmp]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + int_to_str(segPathInfoVec_Nor2[tmp]->returnMismatchPosVecValue(tmp2));
			tmpNorAlignInfo_2 += ",";

			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + (segPathInfoVec_Nor2[tmp]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" + (segPathInfoVec_Nor2[tmp]->returnMismatchCharVecValue(tmp3));
			tmpNorAlignInfo_2 += ",";

			tmpNorAlignInfo_2 += "\t";				
		}
		string tmpRcmAlignInfo_2;
		for(int tmp = 0; tmp < segPathInfoVec_Rcm2.size(); tmp++)
		{
			if(segPathInfoVec_Rcm2.size() >= 40)
				break;
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ segPathInfoVec_Rcm2[tmp]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Rcm2[tmp]->returnChrMapPos()) + ","
				+ segPathInfoVec_Rcm2[tmp]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Rcm2[tmp]->returnMismatchNum()) + ",";//\t";

			int tmpMismatchNum = segPathInfoVec_Rcm2[tmp]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + int_to_str(segPathInfoVec_Rcm2[tmp]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + int_to_str(segPathInfoVec_Rcm2[tmp]->returnMismatchPosVecValue(tmp2));
			tmpRcmAlignInfo_2 += ",";

			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + (segPathInfoVec_Rcm2[tmp]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" + (segPathInfoVec_Rcm2[tmp]->returnMismatchCharVecValue(tmp3));
			tmpRcmAlignInfo_2 += ",";

			tmpRcmAlignInfo_2 += "\t";	
		}	

		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		return tmpAlignInfoStr;
	}

	string getTmpAlignInfoForFinalPair(const string& readName_1,//_ori,
		const string& readName_2,//_ori, 
		const string& readOriSeq_1, 
		const string& readOriSeq_2, const string& readOriQualSeq_1,
		const string& readOriQualSeq_2, bool FastaOrFastq)
	{
		string tmpAlignInfoStr = "\n" + readName_1 + "\t" 
			+ int_to_str(bestSegPathInfoVecIndexVec_Nor1Rcm2.size()) + "\t"
			+ int_to_str(bestSegPathInfoVecIndexVec_Nor2Rcm1.size()) + "\n"
			+ readOriSeq_1 + "\n" + readOriQualSeq_1 + "\n"
			+ readName_2 + "\t"
			+ int_to_str(bestSegPathInfoVecIndexVec_Nor2Rcm1.size()) + "\t"
			+ int_to_str(bestSegPathInfoVecIndexVec_Nor1Rcm2.size()) + "\n"
			+ readOriSeq_2 + "\n" + readOriQualSeq_2; 

		string tmpNorAlignInfo_1;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < bestSegPathInfoVecIndexVec_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(bestSegPathInfoVecIndexVec_Nor1Rcm2.size() >= 40)
				{break;}
			int tmpSegPath_index = bestSegPathInfoVecIndexVec_Nor1Rcm2[tmpNor1Rcm2];
			tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
				+ segPathInfoVec_Nor1[tmpSegPath_index]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Nor1[tmpSegPath_index]->returnChrMapPos()) + ","
				+ segPathInfoVec_Nor1[tmpSegPath_index]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Nor1[tmpSegPath_index]->returnMismatchNum()) + ",";//\t";

			int tmpMismatchNum = segPathInfoVec_Nor1[tmpSegPath_index]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
					+ int_to_str(segPathInfoVec_Nor1[tmpSegPath_index]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" 
					+ int_to_str(segPathInfoVec_Nor1[tmpSegPath_index]->returnMismatchPosVecValue(tmp2));
			tmpNorAlignInfo_1 += ",";

			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 
					+ (segPathInfoVec_Nor1[tmpSegPath_index]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpNorAlignInfo_1 = tmpNorAlignInfo_1 + "-" 
					+ (segPathInfoVec_Nor1[tmpSegPath_index]->returnMismatchCharVecValue(tmp3));
			tmpNorAlignInfo_1 += ",";
			tmpNorAlignInfo_1 += "\t";				
		}
		string tmpRcmAlignInfo_1;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < bestSegPathInfoVecIndexVec_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(bestSegPathInfoVecIndexVec_Nor2Rcm1.size() >= 40)
				{break;}
			int tmpSegPath_index = bestSegPathInfoVecIndexVec_Nor2Rcm1[tmpNor2Rcm1];
			tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1
				+ segPathInfoVec_Rcm1[tmpSegPath_index]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Rcm1[tmpSegPath_index]->returnChrMapPos()) + ","
				+ segPathInfoVec_Rcm1[tmpSegPath_index]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Rcm1[tmpSegPath_index]->returnMismatchNum()) + ",";//\t";			

			int tmpMismatchNum = segPathInfoVec_Rcm1[tmpSegPath_index]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + int_to_str(
					segPathInfoVec_Rcm1[tmpSegPath_index]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" + int_to_str(
					segPathInfoVec_Rcm1[tmpSegPath_index]->returnMismatchPosVecValue(tmp2));
			tmpRcmAlignInfo_1 += ",";

			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 
					+ (segPathInfoVec_Rcm1[tmpSegPath_index]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpRcmAlignInfo_1 = tmpRcmAlignInfo_1 + "-" 
					+ (segPathInfoVec_Rcm1[tmpSegPath_index]->returnMismatchCharVecValue(tmp3));
			tmpRcmAlignInfo_1 += ",";
			tmpRcmAlignInfo_1 += "\t";
		}		
		string tmpNorAlignInfo_2;
		for(int tmpNor2Rcm1 = 0; tmpNor2Rcm1 < bestSegPathInfoVecIndexVec_Nor2Rcm1.size(); tmpNor2Rcm1++)
		{
			if(bestSegPathInfoVecIndexVec_Nor2Rcm1.size() >= 40)
				{break;}
			int tmpSegPath_index = bestSegPathInfoVecIndexVec_Nor2Rcm1[tmpNor2Rcm1];
			tmpNorAlignInfo_2 = tmpNorAlignInfo_2
				+ segPathInfoVec_Nor2[tmpSegPath_index]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Nor2[tmpSegPath_index]->returnChrMapPos()) + ","
				+ segPathInfoVec_Nor2[tmpSegPath_index]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Nor2[tmpSegPath_index]->returnMismatchNum()) + ",";//\t";

			int tmpMismatchNum = segPathInfoVec_Nor2[tmpSegPath_index]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 
					+ int_to_str(segPathInfoVec_Nor2[tmpSegPath_index]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" 
					+ int_to_str(segPathInfoVec_Nor2[tmpSegPath_index]->returnMismatchPosVecValue(tmp2));
			tmpNorAlignInfo_2 += ",";

			if(tmpMismatchNum >= 1)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 
					+ (segPathInfoVec_Nor2[tmpSegPath_index]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpNorAlignInfo_2 = tmpNorAlignInfo_2 + "-" 
					+ (segPathInfoVec_Nor2[tmpSegPath_index]->returnMismatchCharVecValue(tmp3));
			tmpNorAlignInfo_2 += ",";

			tmpNorAlignInfo_2 += "\t";		
		}
		string tmpRcmAlignInfo_2;
		for(int tmpNor1Rcm2 = 0; tmpNor1Rcm2 < bestSegPathInfoVecIndexVec_Nor1Rcm2.size(); tmpNor1Rcm2++)
		{
			if(bestSegPathInfoVecIndexVec_Nor1Rcm2.size() >= 40)
				{break;}
			int tmpSegPath_index = bestSegPathInfoVecIndexVec_Nor1Rcm2[tmpNor1Rcm2];
			tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2
				+ segPathInfoVec_Rcm2[tmpSegPath_index]->returnChrNameStr() + ","
				+ int_to_str(segPathInfoVec_Rcm2[tmpSegPath_index]->returnChrMapPos()) + ","
				+ segPathInfoVec_Rcm2[tmpSegPath_index]->jumpCodeVec2Str() + ","
				+ int_to_str(segPathInfoVec_Rcm2[tmpSegPath_index]->returnMismatchNum()) + ",";//\t";

			int tmpMismatchNum = segPathInfoVec_Rcm2[tmpSegPath_index]->returnMismatchNum();
			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 
					+ int_to_str(segPathInfoVec_Rcm2[tmpSegPath_index]->returnMismatchPosVecValue(0));
			for(int tmp2 = 1; tmp2 < tmpMismatchNum; tmp2++)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" 
					+ int_to_str(segPathInfoVec_Rcm2[tmpSegPath_index]->returnMismatchPosVecValue(tmp2));
			tmpRcmAlignInfo_2 += ",";

			if(tmpMismatchNum >= 1)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 
					+ (segPathInfoVec_Rcm2[tmpSegPath_index]->returnMismatchCharVecValue(0));
			for(int tmp3 = 1; tmp3 < tmpMismatchNum; tmp3++)
				tmpRcmAlignInfo_2 = tmpRcmAlignInfo_2 + "-" 
					+ (segPathInfoVec_Rcm2[tmpSegPath_index]->returnMismatchCharVecValue(tmp3));
			tmpRcmAlignInfo_2 += ",";

			tmpRcmAlignInfo_2 += "\t";	
		}		
		tmpAlignInfoStr = tmpAlignInfoStr + "\n" 
			+ "Nor_1:\t" + tmpNorAlignInfo_1 + "\n"
			+ "Rcm_1:\t" + tmpRcmAlignInfo_1 + "\n"
			+ "Nor_2:\t" + tmpNorAlignInfo_2 + "\n"
			+ "Rcm_2:\t" + tmpRcmAlignInfo_2; 
		return tmpAlignInfoStr;
	}

	void output_phase1(
		vector<string>& PeAlignSamStrVec_complete,
		vector<string>& PeAlignInfoStrVec_inCompletePair,
		vector<string>& PeAlignInfoStrVec_oneEndUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_lowScore,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
		vector<string>& PeAlignSamStrVec_inCompletePair,
		//vector<string>& PeAlignInfoStrVec_completePaired,
		vector<RepeatRegion_Info*>& repeatRegionInfoVec,
		PE_Read_Info& readInfo, Index_Info* indexInfo, int tmpOpenMP, int threadNO,
		Stats_Info* statsInfo, bool FastaOrFastq
		)
	{
		//cout << "output_phase1 starts ......" << endl;
		bool pairExistsBool = this->finalPairExistsBool();
		//cout << "pairExistsBool: " << pairExistsBool << endl;
		bool mappedToRepeatRegionBool = this->mappedToRepeatRegionBool();
		//cout << "mappedToRepeatRegionBool: " << mappedToRepeatRegionBool << endl;
		if(pairExistsBool) // some pair exits
		{
			bool allAlignmentCompleteBool = this->allAlignmentInFinalPairCompleted(
				readInfo.returnReadLength_end1(), readInfo.returnReadLength_end2());
			//cout << "allAlignmentCompleteBool: " << allAlignmentCompleteBool << endl;
			bool unique_bool = this->checkUniqueOrMulti();
			//cout << "unique_bool: " << unique_bool << endl;
			statsInfo->increPairedNum_phase1(threadNO, allAlignmentCompleteBool, unique_bool);
			if(allAlignmentCompleteBool)
			{
				bool completeAlignmentPairScore_tooLow
				 	= this->alignPairScoreTooLow_bool(readInfo);
				//cout << "completeAlignmentPairScore_tooLow: " << completeAlignmentPairScore_tooLow << endl;
				if(completeAlignmentPairScore_tooLow)
				{
					statsInfo->increLowScoreComplete_phase1(threadNO, unique_bool);
					PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmpOpenMP] = 
						this->getSAMformatForBothEndsUnmapped(readInfo, FastaOrFastq);	
				}
				else
				{
					PeAlignSamStrVec_complete[tmpOpenMP] = 
						this->getSAMformatForFinalPair_secondaryOrNot(readInfo, FastaOrFastq);
				}
			}
			else
			{
				//cout << "to getTmpAlignInfoForFinalPair ......." << endl;
				PeAlignInfoStrVec_inCompletePair[tmpOpenMP] = 
					this->getTmpAlignInfoForFinalPair(
				 		readInfo.returnReadName_1(), readInfo.returnReadName_2(), 
						readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(),
						readInfo.returnReadQual_1(), readInfo.returnReadQual_2(),
						FastaOrFastq);
				//cout << "to getSAMformatForFinalPair_secondaryOrNot ..." << endl;
				PeAlignSamStrVec_inCompletePair[tmpOpenMP] = 
					this->getSAMformatForFinalPair_secondaryOrNot(readInfo, FastaOrFastq);
			}
		}
		else // no pair exists: 1.one end unmapped; 2. both ends unmapped
		{
			bool alignmentExistsBool = this->alignInfoExistsBool();
			statsInfo->increUnpairedNum_phase1(threadNO, alignmentExistsBool, mappedToRepeatRegionBool);
			if(alignmentExistsBool) // one end unmapped
			{	
				PeAlignInfoStrVec_oneEndUnmapped[tmpOpenMP] = 
					this->getTmpAlignInfo(
						readInfo.returnReadName_1(), readInfo.returnReadName_2(), 
						readInfo.returnReadSeq_1(), readInfo.returnReadSeq_2(),
						readInfo.returnReadQual_1(), readInfo.returnReadQual_2(),
						FastaOrFastq);
			}
			else if(mappedToRepeatRegionBool)
			{
				PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmpOpenMP] = 
					this->getSAMformatForBothEndsUnmapped_mappedToRepeatRegionReads(
						readInfo, repeatRegionInfoVec[threadNO], indexInfo,
						FastaOrFastq);
			}
			else // both ends unmapped
			{
				PeAlignSamStrVec_bothEndsUnmapped[tmpOpenMP] = 
					this->getSAMformatForBothEndsUnmapped(
						readInfo, FastaOrFastq);						
			}
		}			
	}

	void memoryFree()
	{
		int segPathInfoVec_Nor1_size = segPathInfoVec_Nor1.size();
		int segPathInfoVec_Rcm1_size = segPathInfoVec_Rcm1.size();
		int segPathInfoVec_Nor2_size = segPathInfoVec_Nor2.size();
		int segPathInfoVec_Rcm2_size = segPathInfoVec_Rcm2.size();
		for(int tmp = 0; tmp < segPathInfoVec_Nor1_size; tmp++)
		{
			delete segPathInfoVec_Nor1[tmp];
		}
		for(int tmp = 0; tmp < segPathInfoVec_Rcm1_size; tmp++)
		{
			delete segPathInfoVec_Rcm1[tmp];
		}
		for(int tmp = 0; tmp < segPathInfoVec_Nor2_size; tmp++)
		{
			delete segPathInfoVec_Nor2[tmp];
		}
		for(int tmp = 0; tmp < segPathInfoVec_Rcm2_size; tmp++)
		{
			delete segPathInfoVec_Rcm2[tmp];
		}				
	}
};
#endif