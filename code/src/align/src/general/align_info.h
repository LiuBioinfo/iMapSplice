// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGN_INFO_H
#define ALIGN_INFO_H

#include "SJ_info.h"
#include "transcript_set.h"
#include "quality_score.h"
//#include "peAlign_info.h"

using namespace std;

class Alignment_Info
{
//public:
private:

	string alignDirection;

	unsigned int mapPosInWholeGenome;

	string alignChromName;
	int alignChromPos;
	int mismatchNum;
	string SJstrand;

	int mappedBaseNum;
	//int startPosInRead;
	//int endPosInRead;
	int endMatchedPosInChr;

	int insertionLength;
	int deletionLength;

	vector<int> mismatchPosVec;
	vector<char> mismatchCharVec;

public:
	vector<Jump_Code> cigarStringJumpCode;
	vector< pair< int, SpliceJunction_Alignment > > spliceJunctionVec; // <posInRead, SpliceJunction>
	
	int returnFirstMatchLen()
	{
		if(cigarStringJumpCode[0].type == "M")
			return cigarStringJumpCode[0].len;
		else
			return cigarStringJumpCode[1].len;
	}

	int returnLastMatchLen()
	{
		int jumpCodeVecSize = cigarStringJumpCode.size();
		if(cigarStringJumpCode[jumpCodeVecSize-1].type == "M")
			return cigarStringJumpCode[jumpCodeVecSize-1].len;
		else
			return cigarStringJumpCode[jumpCodeVecSize-2].len;
	}

	int returnHeadSoftClipLen()
	{
		if(cigarStringJumpCode[0].type == "S")
			return cigarStringJumpCode[0].len;
		else
			return 0;				
	}

	int returnTailSoftClipLen()
	{
		int jumpCodeVecSize = cigarStringJumpCode.size();
		if(cigarStringJumpCode[jumpCodeVecSize-1].type == "S")
			return cigarStringJumpCode[jumpCodeVecSize-1].len;
		else
			return 0;		
	}

	void memoryFree()
	{
		vector<int>().swap(mismatchPosVec);
		vector<char>().swap(mismatchCharVec);
		vector<Jump_Code>().swap(cigarStringJumpCode);
		vector< pair< int, SpliceJunction_Alignment > >().swap(spliceJunctionVec);
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

	void generateCoveredBaseVec(int tmpReadLength, vector<bool>& coveredBaseBoolVec , bool for_or_rev_bool)
	{
		for(int tmp = 0; tmp < tmpReadLength; tmp++)
			coveredBaseBoolVec.push_back(false);
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if(cigarStringJumpCode[tmp].type == "M")
			{
				int tmpStartLocInRead = this->getEndLocInReadOfSpecificJumpCode(cigarStringJumpCode, tmp - 1) + 1;
				int tmpEndLocInRead = this->getEndLocInReadOfSpecificJumpCode(cigarStringJumpCode, tmp);
				if(for_or_rev_bool)
				{
					for(int tmpLoc = tmpStartLocInRead; tmpLoc <= tmpEndLocInRead; tmpLoc++)
						coveredBaseBoolVec[tmpLoc-1] = true;
				}
				else
				{
					for(int tmpLoc = tmpStartLocInRead; tmpLoc <= tmpEndLocInRead; tmpLoc++)
						coveredBaseBoolVec[tmpReadLength - tmpLoc] = true;					
				}
			}
		}
	}

	string returnSJstrand()
	{
		return SJstrand;
	}

	bool backSpliceExists_bool()
	{
		if(spliceJunctionVec.size() == 0)
			return false;
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			int tmpSJ_startPos = (spliceJunctionVec[tmp].second).returnSJdonerEnd();
			int tmpSJ_endPos = (spliceJunctionVec[tmp].second).returnSJacceptorStart();
			if(tmpSJ_startPos > tmpSJ_endPos)
				return true;
		}
		return false;
	}

	int backSplicceDistanceTotal()
	{
		int tmpBackSpliceDistanceTotal = 0;
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			int tmpSJ_startPos = (spliceJunctionVec[tmp].second).returnSJdonerEnd();
			int tmpSJ_endPos = (spliceJunctionVec[tmp].second).returnSJacceptorStart();
			if(tmpSJ_startPos > tmpSJ_endPos)
			{
				tmpBackSpliceDistanceTotal += (tmpSJ_startPos - tmpSJ_endPos - 1);
			}			
		}
		return tmpBackSpliceDistanceTotal;
	}

	void copyAlignInfo(Alignment_Info* oriAlignInfo)
	{
		alignDirection = oriAlignInfo->returnAlignDirection();
		mapPosInWholeGenome = oriAlignInfo->returnMapPosInWholeGenome();
		alignChromName = oriAlignInfo->returnAlignChromName();
		alignChromPos = oriAlignInfo->returnAlignChromPos();
		mismatchNum = oriAlignInfo->returnMismatchNum();
		SJstrand = oriAlignInfo->returnSJstrand();
		mappedBaseNum = oriAlignInfo->returnMappedBaseNum();
		endMatchedPosInChr = oriAlignInfo->returnEndMatchedPosInChr();
		insertionLength = oriAlignInfo->returnInsertionLength();
		deletionLength = oriAlignInfo->returnDeletionLength();
		mismatchPosVec.clear();
		for(int tmp = 0; tmp < oriAlignInfo->returnMismatchPosVecSize(); tmp++)
		{
			mismatchPosVec.push_back(oriAlignInfo->returnMismatchPosVecValue(tmp));
		}
		mismatchCharVec.clear();
		for(int tmp = 0; tmp < oriAlignInfo->returnMismatchCharVecSize(); tmp++)
		{
			mismatchCharVec.push_back(oriAlignInfo->returnMismatchCharVecValue(tmp));
		}
		cigarStringJumpCode.clear();
		for(int tmp = 0; tmp < (oriAlignInfo->cigarStringJumpCode).size(); tmp++)
		{
			cigarStringJumpCode.push_back((oriAlignInfo->cigarStringJumpCode)[tmp]);
		}
		spliceJunctionVec.clear();
		for(int tmp = 0; tmp < (oriAlignInfo->spliceJunctionVec).size(); tmp++)
		{
			spliceJunctionVec.push_back(oriAlignInfo->spliceJunctionVec[tmp]);
		}
	}	

	void removeDuplicateMismatch()
	{
		set<int> uniqueMismatchSet;
		
		vector<int> uniqueMismatchIndexVec;		
		vector<int> uniqueMismatchPosVec;
		vector<char> uniqueMismatchCharVec;

		for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = mismatchPosVec[tmp];
			if(uniqueMismatchSet.find(tmpMismatchPos) == uniqueMismatchSet.end())
			{
				uniqueMismatchSet.insert(tmpMismatchPos);
				uniqueMismatchPosVec.push_back(tmpMismatchPos);
				uniqueMismatchIndexVec.push_back(tmp);
			}
			else
			{
			}
		}

		int uniqueMismatchNum = uniqueMismatchPosVec.size();
		if(mismatchNum == uniqueMismatchNum)
		{}
		else
		{
			mismatchNum = uniqueMismatchNum;
			if(STORE_MISMATCH_POS)
			{
				mismatchPosVec.clear();
				for(int tmp = 0; tmp < uniqueMismatchNum; tmp++)
				{
					mismatchPosVec.push_back(uniqueMismatchPosVec[tmp]);
				}
				if(STORE_MISMATCH_CHA)
				{
					//mismatchCharVec.clear();
					for(int tmp = 0; tmp < uniqueMismatchNum; tmp++)
					{
						int tmpIndex = uniqueMismatchIndexVec[tmp];
						uniqueMismatchCharVec.push_back(mismatchCharVec[tmpIndex]);
					}
					mismatchCharVec.clear();
					for(int tmp = 0; tmp < uniqueMismatchNum; tmp++)
					{
						mismatchCharVec.push_back(uniqueMismatchCharVec[tmp]);
					}
				}
			}
		}
	}

	int getMapppingQuality(const string& qualitySeq, string quality_scale, const string& readName)
	{
		int tmpReadLength = qualitySeq.length();
		vector<int> mismatchPosVec_all;
		for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
		{
			mismatchPosVec_all.push_back(mismatchPosVec[tmp]);
		}

		int jumpCodeSize = cigarStringJumpCode.size();
		if(jumpCodeSize > 0)
		{
			if(cigarStringJumpCode[0].type == "S")
			{
				int headSoftClippingLen = cigarStringJumpCode[0].len;
				for(int tmp = 1; tmp <= headSoftClippingLen; tmp++)
				{
					mismatchPosVec_all.push_back(tmp);
				}
			}
			if(cigarStringJumpCode[jumpCodeSize-1].type == "S")
			{
				int tailSoftClippingLen = cigarStringJumpCode[jumpCodeSize-1].len;
				for(int tmp = tmpReadLength - tailSoftClippingLen + 1; tmp <= tmpReadLength; tmp++)
				{
					mismatchPosVec_all.push_back(tmp);
				}
			}
		}

		int mappingQuality = GetQualityScore_new(mismatchPosVec_all, qualitySeq, quality_scale);

		return mappingQuality;
	}

	string returnAlignInfoStr()
	{
		string tmpStr;
		tmpStr = alignChromName + " " + int_to_str(alignChromPos) + " " + alignDirection + " ";
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(cigarStringJumpCode[tmp].len) + cigarStringJumpCode[tmp].type;
		}
		return tmpStr;
	}

	bool containSJ()
	{
		bool containSJ_bool = false;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if(cigarStringJumpCode[tmp].type == "N")
			{
				return true;
			}
		}
		return containSJ_bool;
	}

	void generateMismatchNum()
	{
		mismatchNum = mismatchPosVec.size();
	}

	int returnJumpCodeSize()
	{
		return cigarStringJumpCode.size();
	}
	int returnCigarStringJumpCodeSize()
	{
		return cigarStringJumpCode.size();
	}
	string returnCigarStringJumpCodeType(int tmp)
	{
		return cigarStringJumpCode[tmp].type;
	}
	int returnCigarStringJumpCodeLen(int tmp)
	{
		return cigarStringJumpCode[tmp].len;
	}
	Jump_Code returnCigarStringJumpCodeElement(int index)
	{
		return cigarStringJumpCode[index];
	}

	void insertSJ2transcriptAlignInfoExonArea(int mapPosInTranscriptExon_start,
		int tmpJumpCodeLen, vector<pair<int,int> > SJvec)
	{
		vector< pair<int, int> > SJvecInExon;
		for(int tmp = 0; tmp < SJvec.size(); tmp++)
		{
			int spliceSite = SJvec[tmp].first;
			if((spliceSite >= mapPosInTranscriptExon_start)&&(spliceSite < mapPosInTranscriptExon_start + tmpJumpCodeLen - 1))
				SJvecInExon.push_back(SJvec[tmp]);
		}
		vector<int> exonLength;
		if(SJvecInExon.size() == 0)
		{
			exonLength.push_back(tmpJumpCodeLen);
		}
		else
		{
			int firstExonLength = SJvecInExon[0].first - mapPosInTranscriptExon_start;
			exonLength.push_back(firstExonLength);

			int SJvecInExon_size = SJvecInExon.size();
			for(int tmp = 0; tmp < SJvecInExon_size-1; tmp ++)
			{
				int spliceSite_1 = SJvecInExon[tmp].first;
				int spliceSite_2 = SJvecInExon[tmp + 1].first;
				int tmpExonLength = spliceSite_2 - spliceSite_1;
				exonLength.push_back(tmpExonLength);
			}

			int lastSpliceSite = SJvecInExon[SJvecInExon_size-1].first;
			int lastExonLength = mapPosInTranscriptExon_start + tmpJumpCodeLen - 1 - lastSpliceSite + 1;
			exonLength.push_back(lastExonLength);
		}

		if(exonLength[0] == 0)
		{
			int cigarStringJumpCodeSize = cigarStringJumpCode.size();
			if(cigarStringJumpCodeSize == 0)
			{}
			else if(cigarStringJumpCode[cigarStringJumpCodeSize - 1].type == "D")
			{
				cigarStringJumpCode[cigarStringJumpCodeSize - 1].type = "N";
				cigarStringJumpCode[cigarStringJumpCodeSize - 1].len += SJvecInExon[0].second;
			}
			else 
			{
				//int tmp_jumpCodeType = "N";
				int tmp_jumpCodeLen = SJvecInExon[0].second;
				Jump_Code newJumpCode(tmp_jumpCodeLen, "N");
				cigarStringJumpCode.push_back(newJumpCode);
			}	

			if(exonLength.size() >= 2)
			{
				Jump_Code newMatchJumpCode(exonLength[1], "M");
				cigarStringJumpCode.push_back(newMatchJumpCode);
				for(int tmp = 2; tmp < exonLength.size(); tmp++)
				{
					Jump_Code SJjumpCode(SJvecInExon[tmp-1].second, "N");
					Jump_Code matchJumpCode(exonLength[tmp], "M");
					cigarStringJumpCode.push_back(SJjumpCode);
					cigarStringJumpCode.push_back(matchJumpCode);
				}
			}
		}
		else
		{
			Jump_Code firstMatchJumpCode(exonLength[0], "M");
			cigarStringJumpCode.push_back(firstMatchJumpCode);
			for(int tmp = 1; tmp < exonLength.size(); tmp++)
			{
				Jump_Code SJjumpCode(SJvecInExon[tmp-1].second, "N");
				Jump_Code matchJumpCode(exonLength[tmp], "M");
				cigarStringJumpCode.push_back(SJjumpCode);
				cigarStringJumpCode.push_back(matchJumpCode);
			}
		}
		/*int left_exon_length = tmpJumpCodeLen;
		int tmp_mapPosInTranscript_start = mapPosInTranscriptExon_start;
		for(int tmp = 0; tmp < SJvec.size(); tmp++)
		{
			cout << "tmp: " << tmp << endl;
			cout << "tmp_mapPosInTranscript_start: " << tmp_mapPosInTranscript_start << endl;
			cout << "tmp_SJ_site: " << tmp_SJ_site << endl;
			int tmp_SJ_site = SJvec[tmp].first;
			if(tmp_SJ_site < tmp_mapPosInTranscript_start)
			{}
			else if((tmp_SJ_site > tmp_mapPosInTranscript_start + left_exon_length - 1))
			{
				Jump_Code newJumpCode(tmpJumpCodeLen, "M");
				cigarStringJumpCode.push_back(newJumpCode);
			}
			else if(tmp_SJ_site == tmp_mapPosInTranscript_start)
			{
				int cigarStringJumpCodeSize = cigarStringJumpCode.size();
				if(cigarStringJumpCodeSize == 0)
				{}
				else if(cigarStringJumpCode[cigarStringJumpCodeSize - 1].type == "D")
				{
					cigarStringJumpCode[cigarStringJumpCodeSize - 1].type = "N";
					cigarStringJumpCode[cigarStringJumpCodeSize - 1].len += SJvec[tmp].second;
				}
				else 
				{
					//int tmp_jumpCodeType = "N";
					int tmp_jumpCodeLen = SJvec[tmp].second;
					Jump_Code newJumpCode(tmp_jumpCodeLen, "N");
					cigarStringJumpCode.push_back(newJumpCode);
				}
			}
			else
			{
				int prefixMatch_length = tmp_SJ_site - tmp_mapPosInTranscript_start;
				int tmp_jumpCodeLen = SJvec[tmp].second;
				Jump_Code prefixMatchJumpCode(prefixMatch_length, "M");
				Jump_Code sjJumpCode(tmp_jumpCodeLen, "N");
				cigarStringJumpCode.push_back(prefixMatchJumpCode);
				cigarStringJumpCode.push_back(sjJumpCode);
			}
		}*/
	}

	void insertSJ2SNPseqAlignInfoExonArea(int mapPosInSNPseqExon_start,
		int tmpJumpCodeLen, vector<pair<int,int> > SJvec)
	{
		vector< pair<int, int> > SJvecInExon;
		for(int tmp = 0; tmp < SJvec.size(); tmp++)
		{
			int spliceSite = SJvec[tmp].first;
			if((spliceSite >= mapPosInSNPseqExon_start)&&(spliceSite < mapPosInSNPseqExon_start + tmpJumpCodeLen - 1))
				SJvecInExon.push_back(SJvec[tmp]);
		}
		vector<int> exonLength;
		if(SJvecInExon.size() == 0)
			exonLength.push_back(tmpJumpCodeLen);
		else
		{
			int firstExonLength = SJvecInExon[0].first - mapPosInSNPseqExon_start;
			exonLength.push_back(firstExonLength);

			int SJvecInExon_size = SJvecInExon.size();
			for(int tmp = 0; tmp < SJvecInExon_size-1; tmp ++)
			{
				int spliceSite_1 = SJvecInExon[tmp].first;
				int spliceSite_2 = SJvecInExon[tmp + 1].first;
				int tmpExonLength = spliceSite_2 - spliceSite_1;
				exonLength.push_back(tmpExonLength);
			}

			int lastSpliceSite = SJvecInExon[SJvecInExon_size-1].first;
			int lastExonLength = mapPosInSNPseqExon_start + tmpJumpCodeLen - 1 - lastSpliceSite + 1;
			exonLength.push_back(lastExonLength);
		}

		if(exonLength[0] == 0)
		{
			int cigarStringJumpCodeSize = cigarStringJumpCode.size();
			if(cigarStringJumpCodeSize == 0)
			{}
			else if(cigarStringJumpCode[cigarStringJumpCodeSize - 1].type == "D")
			{
				cigarStringJumpCode[cigarStringJumpCodeSize - 1].type = "N";
				cigarStringJumpCode[cigarStringJumpCodeSize - 1].len += SJvecInExon[0].second;
			}
			else 
			{
				//int tmp_jumpCodeType = "N";
				int tmp_jumpCodeLen = SJvecInExon[0].second;
				Jump_Code newJumpCode(tmp_jumpCodeLen, "N");
				cigarStringJumpCode.push_back(newJumpCode);
			}	

			if(exonLength.size() >= 2)
			{
				Jump_Code newMatchJumpCode(exonLength[1], "M");
				cigarStringJumpCode.push_back(newMatchJumpCode);
				for(int tmp = 2; tmp < exonLength.size(); tmp++)
				{
					Jump_Code SJjumpCode(SJvecInExon[tmp-1].second, "N");
					Jump_Code matchJumpCode(exonLength[tmp], "M");
					cigarStringJumpCode.push_back(SJjumpCode);
					cigarStringJumpCode.push_back(matchJumpCode);
				}
			}
		}
		else
		{
			Jump_Code firstMatchJumpCode(exonLength[0], "M");
			cigarStringJumpCode.push_back(firstMatchJumpCode);
			for(int tmp = 1; tmp < exonLength.size(); tmp++)
			{
				Jump_Code SJjumpCode(SJvecInExon[tmp-1].second, "N");
				Jump_Code matchJumpCode(exonLength[tmp], "M");
				cigarStringJumpCode.push_back(SJjumpCode);
				cigarStringJumpCode.push_back(matchJumpCode);
			}
		}
	}

	Jump_Code returnJumpCode(int tmp)
	{
		return cigarStringJumpCode[tmp];
	}

	void getNewJumpCodeVec(Alignment_Info* transcriptAlignInfo,
		Transcript_Set* transcriptInfo, vector< pair<int,int> >& SJvec)
	{
		int mapPosInTranscript = transcriptAlignInfo->returnAlignChromPos();
		int tmp_mapPosInTranscript = mapPosInTranscript;
		int tmp_locInRead = 0;
		//int tmp_exon_length = 0;
		for(int tmp = 0; tmp < transcriptAlignInfo->returnJumpCodeSize(); tmp++)
		{
			int tmp_mapPosInTranscript_old = tmp_mapPosInTranscript;
			int tmp_locInRead_old = tmp_locInRead;
			
			int tmp_mapPosInTranscript_new;
			int tmp_locInRead_new;

			int tmp_jumpCodeLen = transcriptAlignInfo->returnCigarStringJumpCodeLen(tmp);
			string tmp_jumpCodeType = transcriptAlignInfo->returnCigarStringJumpCodeType(tmp);

			//cout << "tmp_jumpCodeLen: " << tmp_jumpCodeLen << endl << "tmp_jumpCodeType: " << tmp_jumpCodeType << endl;

			if(tmp_jumpCodeType == "S")
			{
				tmp_mapPosInTranscript += 0;
				tmp_locInRead += tmp_jumpCodeLen;
				cigarStringJumpCode.push_back(
					transcriptAlignInfo->returnJumpCode(tmp));
			}
			else if(tmp_jumpCodeType == "M")
			{
				tmp_mapPosInTranscript += tmp_jumpCodeLen;
				tmp_locInRead += tmp_jumpCodeLen;
				this->insertSJ2transcriptAlignInfoExonArea(tmp_mapPosInTranscript_old,
					tmp_jumpCodeLen, SJvec);
			}
			else if(tmp_jumpCodeType == "I")
			{
				tmp_mapPosInTranscript += 0;
				tmp_locInRead += tmp_jumpCodeLen;
				cigarStringJumpCode.push_back(
					transcriptAlignInfo->returnJumpCode(tmp));
			}
			else if(tmp_jumpCodeType == "D")
			{
				tmp_mapPosInTranscript += tmp_jumpCodeLen;
				tmp_locInRead += 0;				
				cigarStringJumpCode.push_back(
					transcriptAlignInfo->returnJumpCode(tmp));			
			}
			else if(tmp_jumpCodeType == "N")
			{
				tmp_mapPosInTranscript += tmp_jumpCodeLen;
				tmp_locInRead += 0;
				cigarStringJumpCode.push_back(
					transcriptAlignInfo->returnJumpCode(tmp));			
			}
			else
			{
				cout << "invalid jumpCodeType ..." << endl;
			}
		}
	}

	void getNewJumpCodeVec_fromSNPseqAlignInfo(Alignment_Info* SNPseqAlignInfo,
		SyntheticSNPtransSeq_Info& tmpSNPseqInfo, vector< pair<int,int> >& SJvec)
	{
		int mapPosInSNPseq = SNPseqAlignInfo->returnAlignChromPos();
		int tmp_mapPosInSNPseq = mapPosInSNPseq;
		int tmp_locInRead = 0;
		//int tmp_exon_length = 0;
		for(int tmp = 0; tmp < SNPseqAlignInfo->returnJumpCodeSize(); tmp++)
		{
			int tmp_mapPosInSNPseq_old = tmp_mapPosInSNPseq;
			int tmp_locInRead_old = tmp_locInRead;
			
			int tmp_mapPosInSNPseq_new;
			int tmp_locInRead_new;

			int tmp_jumpCodeLen = SNPseqAlignInfo->returnCigarStringJumpCodeLen(tmp);
			string tmp_jumpCodeType = SNPseqAlignInfo->returnCigarStringJumpCodeType(tmp);

			if(tmp_jumpCodeType == "S")
			{
				tmp_mapPosInSNPseq += 0;
				tmp_locInRead += tmp_jumpCodeLen;
				cigarStringJumpCode.push_back(SNPseqAlignInfo->returnJumpCode(tmp));
			}
			else if(tmp_jumpCodeType == "M")
			{
				tmp_mapPosInSNPseq += tmp_jumpCodeLen;
				tmp_locInRead += tmp_jumpCodeLen;
				this->insertSJ2SNPseqAlignInfoExonArea(tmp_mapPosInSNPseq_old,
					tmp_jumpCodeLen, SJvec);
			}
			else if(tmp_jumpCodeType == "I")
			{
				tmp_mapPosInSNPseq += 0;
				tmp_locInRead += tmp_jumpCodeLen;
				cigarStringJumpCode.push_back(SNPseqAlignInfo->returnJumpCode(tmp));
			}
			else if(tmp_jumpCodeType == "D")
			{
				tmp_mapPosInSNPseq += tmp_jumpCodeLen;
				tmp_locInRead += 0;				
				cigarStringJumpCode.push_back(SNPseqAlignInfo->returnJumpCode(tmp));			
			}
			else if(tmp_jumpCodeType == "N")
			{
				tmp_mapPosInSNPseq += tmp_jumpCodeLen;
				tmp_locInRead += 0;
				cigarStringJumpCode.push_back(SNPseqAlignInfo->returnJumpCode(tmp));			
			}
			else
			{
				cout << "invalid jumpCodeType ..." << endl;
			}
		}
	}

	bool syntheticSNPtransSeqAlignInfo_coverSNPorNot_bool(int tmpReadLength)
	{
		string tmpSyntheticSNPtransSeq_name = this->returnAlignChromName();
		int tmpSyntheticSNPtransSeq_startPos = this->returnAlignChromPos();
		int tmpSyntheticSNPtransSeq_endPos = this->returnEndMatchedPosInChr();
		int tabLoc_1 = tmpSyntheticSNPtransSeq_name.find(":");
		int tabLoc_2 = tmpSyntheticSNPtransSeq_name.find(":", tabLoc_1 + 1);
		int tabLoc_3 = tmpSyntheticSNPtransSeq_name.find(":", tabLoc_2 + 1);
		int tabLoc_4 = tmpSyntheticSNPtransSeq_name.find(":", tabLoc_3 + 1);
		string tmpSyntheticSNPtransSeq_SNPlocInSeqStr = tmpSyntheticSNPtransSeq_name.substr(
			tabLoc_3 + 1, tabLoc_4 - tabLoc_3 - 1);
		int tmpSyntheticSNPtransSeq_SNPlocInSeqInt 
			= atoi(tmpSyntheticSNPtransSeq_SNPlocInSeqStr.c_str());
		if((tmpSyntheticSNPtransSeq_SNPlocInSeqInt >= tmpSyntheticSNPtransSeq_startPos)
			&&(tmpSyntheticSNPtransSeq_SNPlocInSeqInt <= tmpSyntheticSNPtransSeq_endPos))
			return true;
		else
			return false;
	}

	void convertSyntheticSNPtransSeqAlignInfo2GenomeAlignInfo(Alignment_Info* syntheticSNPtransSeqAlignInfo,
		Index_Info* syntheticSNPtransSeqIndexInfo, Index_Info* genomeIndexInfo)
	{
		string tmpSyntheticSNPtransSeq_name = syntheticSNPtransSeqAlignInfo->returnAlignChromName();		
		SyntheticSNPtransSeq_Info tmpSNPseqInfo;
		tmpSNPseqInfo.initiateWithSyntheticSNPtransSeqName(tmpSyntheticSNPtransSeq_name, genomeIndexInfo);
		vector< pair<int, int> > SJvec; // <spliceSite, spliceSize>
		tmpSNPseqInfo.getSJsiteAndSize(syntheticSNPtransSeqAlignInfo->returnAlignChromPos(),
			syntheticSNPtransSeqAlignInfo->returnEndMatchedPosInChr(), SJvec);
		// alignDirection //
		alignDirection = syntheticSNPtransSeqAlignInfo->returnAlignDirection();
		// alignChromPos //
		alignChromName = tmpSNPseqInfo.return_chrName();
		// alignChromPos //
		alignChromPos = tmpSNPseqInfo.returnMapPosInChrWithMapPosInSNPseq(syntheticSNPtransSeqAlignInfo->returnAlignChromPos());
		// mismatchNum //
		mismatchNum = syntheticSNPtransSeqAlignInfo->returnMismatchNum();
		// mappedBaseNum //
		mappedBaseNum = syntheticSNPtransSeqAlignInfo->returnMappedBaseNum();
		// endMatchedPosInChr //
		endMatchedPosInChr = tmpSNPseqInfo.returnMapPosInChrWithMapPosInSNPseq(syntheticSNPtransSeqAlignInfo->returnEndMatchedPosInChr());
		// insertionLength //
		insertionLength = syntheticSNPtransSeqAlignInfo->returnInsertionLength();
		// deletionLength //
		deletionLength = syntheticSNPtransSeqAlignInfo->returnDeletionLength();
		// generate new jump code vec //
		this->getNewJumpCodeVec_fromSNPseqAlignInfo(syntheticSNPtransSeqAlignInfo, tmpSNPseqInfo, SJvec);
		// spliceJunctionVec //
		this->jumpCodeVec2spliceJunctionVec(genomeIndexInfo);
		// SJstrand //
		SJstrand = this->getStrandFromSJ();
		// mismatch pos/cha vec
		int mismatchPosVecSize = syntheticSNPtransSeqAlignInfo->returnMismatchPosVecSize();
		for(int tmp = 0; tmp < mismatchPosVecSize; tmp++)
			mismatchPosVec.push_back(syntheticSNPtransSeqAlignInfo->returnMismatchPosVecValue(tmp));
		int mismatchCharVecSize = syntheticSNPtransSeqAlignInfo->returnMismatchCharVecSize();
		for(int tmp2 = 0; tmp2 < mismatchCharVecSize; tmp2++)
			mismatchCharVec.push_back(syntheticSNPtransSeqAlignInfo->returnMismatchCharVecValue(tmp2));
	}

	void convertTranscriptAlignInfo2GenomeAlignInfo(Alignment_Info* transcriptAlignInfo, 
		Transcript_Set* transcriptInfo, Index_Info* genomeIndexInfo)
	{
		//cout << " ------ start to convertTranscriptAlignInfo2GenomeAlignInfo " << endl;
		string transcript_name = transcriptAlignInfo->returnAlignChromName();
		//cout << "transcript_name: " << transcript_name << endl;
		int index_transcriptSet = transcriptInfo->returnTranscriptSetIndex(transcript_name);
		//cout << "index_transcriptSet: " << index_transcriptSet << endl;
		int align_pos_start_transcript = transcriptAlignInfo->returnAlignChromPos();
		int align_pos_end_transcript = transcriptAlignInfo->returnEndMatchedPosInChr();
		//cout << "align_pos_start_transcript: " << align_pos_start_transcript << endl;
		//cout << "align_pos_end_transcript: " << align_pos_end_transcript << endl;
		vector< pair<int, int> > SJvec; // <spliceSite, spliceSize>
		//cout << "start to getSJsiteAndSize " << endl;
		transcriptInfo->getSJsiteAndSize(index_transcriptSet, align_pos_start_transcript,
			align_pos_end_transcript, SJvec);
		/*cout << endl << "SJvec info: " << endl;
		for(int tmp = 0; tmp < SJvec.size(); tmp++)
		{
			cout << SJvec[tmp].first << " " << SJvec[tmp].second << endl;
		}
		cout << endl << "finish getting SJsiteAndSize " << endl;*/
		// alignDirection //
		alignDirection = transcriptAlignInfo->returnAlignDirection();
		//cout << "alignDirection: " << alignDirection << endl;
		// mapPosInWholeGenome //

		// alignChromName //
		alignChromName = transcriptInfo->returnChromName(index_transcriptSet);
		//cout << "alignChromName: " << alignChromName << endl;
		// alignChromPos //
		alignChromPos = transcriptInfo->returnMapPosInChrWithMapPosInTranscript(
			index_transcriptSet, align_pos_start_transcript);
		//cout << "alignChromPos: " << alignChromPos << endl;
		// mismatchNum //
		mismatchNum = transcriptAlignInfo->returnMismatchNum();
		//cout << "mismatchNum: " << mismatchNum << endl;
		// mappedBaseNum //
		mappedBaseNum = transcriptAlignInfo->returnMappedBaseNum();
		//cout << "mappedBaseNum: " << mappedBaseNum << endl;
		// endMatchedPosInChr //
		endMatchedPosInChr = transcriptInfo->returnMapPosInChrWithMapPosInTranscript(
			index_transcriptSet, align_pos_end_transcript);
		// insertionLength //
		insertionLength = transcriptAlignInfo->returnInsertionLength();
		// deletionLength //
		deletionLength = transcriptAlignInfo->returnDeletionLength();
		// cigarStringJumpCode //
		//cout << "start to get cigarString" << endl;
		this->getNewJumpCodeVec(transcriptAlignInfo, transcriptInfo, SJvec);

		// spliceJunctionVec //
		this->jumpCodeVec2spliceJunctionVec(genomeIndexInfo);
		// SJstrand //
		SJstrand = this->getStrandFromSJ();

		// mismatch pos/cha vec
		if(STORE_MISMATCH_POS)
		{
			int mismatchPosVecSize = transcriptAlignInfo->returnMismatchPosVecSize();
			for(int tmp = 0; tmp < mismatchPosVecSize; tmp++)
			{
				mismatchPosVec.push_back(transcriptAlignInfo->returnMismatchPosVecValue(tmp));
			}
			if(STORE_MISMATCH_CHA)
			{
				int mismatchCharVecSize = transcriptAlignInfo->returnMismatchCharVecSize();
				for(int tmp2 = 0; tmp2 < mismatchCharVecSize; tmp2++)
				{
					mismatchCharVec.push_back(transcriptAlignInfo->returnMismatchCharVecValue(tmp2));
				}
			}
		}
	}

	int returnMismatchPosVecSize()
	{
		return mismatchPosVec.size();
	}

	int returnMismatchPosVecValue(int index)
	{
		return mismatchPosVec[index];
	}

	int returnMismatchCharVecSize()
	{
		return mismatchCharVec.size();
	}

	char returnMismatchCharVecValue(int index)
	{
		return mismatchCharVec[index];
	}
	int returnMappedBaseNum()
	{
		return mappedBaseNum;
	}

	void generateIndelLength()
	{
		int tmp_ins_len = 0;
		int tmp_del_len = 0;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if(cigarStringJumpCode[tmp].type == "I")
			{
				tmp_ins_len += cigarStringJumpCode[tmp].len;
			}
			if(cigarStringJumpCode[tmp].type == "D")
			{
				tmp_del_len += cigarStringJumpCode[tmp].len;
			}
		}
		insertionLength = tmp_ins_len;
		deletionLength = tmp_del_len;
	}
	int returnInsertionLength()
	{
		return insertionLength;
	}
	int returnDeletionLength()
	{
		return deletionLength;
	}
	int returnMismatchNum()
	{
		return mismatchNum;
	}
	int returnEndMatchedPosInChr()
	{
		return endMatchedPosInChr;
	}
	vector<Jump_Code>& returnCigarStringJumpCode()
	{
		return cigarStringJumpCode;
	}
	string returnCigarString()
	{
		string str;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}
	int returnAlignChromPos()
	{
		return alignChromPos;
	}
	string returnAlignChromName()
	{
		return alignChromName;
	}
	string returnAlignDirection()
	{
		return alignDirection;
	}
	unsigned int returnMapPosInWholeGenome()
	{
		return mapPosInWholeGenome;
	}
	bool checkAllJumpCodeValid()
	{
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if((cigarStringJumpCode[tmp].len < 0) || (cigarStringJumpCode[tmp].len > 500000))
				return false;
		}
		return true;
	}

	Alignment_Info()
	{}

	int unfixedHeadLength()
	{
		return cigarStringJumpCode[0].len;
	}

	int unfixedTailLength()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		return cigarStringJumpCode[jumpCodeNum - 1].len;
	}
	bool unfixedHeadExistsBool()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		if(jumpCodeNum < 0)
		{
			return true;
		}

		bool unfixedHeadBool = (cigarStringJumpCode[0].type == "S");
		return unfixedHeadBool;
	}

	bool unfixedTailExistsBool()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		if(jumpCodeNum < 0)
		{
			return true;
		}

		bool unfixedHeadBool = (cigarStringJumpCode[jumpCodeNum-1].type == "S");
		return unfixedHeadBool;		
	}

	bool alignmentIncompleteBool()
	{
		bool unfixedHeadExists_bool = this->unfixedHeadExistsBool();
		bool unfixedTailExists_bool = this->unfixedTailExistsBool();
		if(unfixedHeadExists_bool || unfixedTailExists_bool)
			return true;
		else
			return false;
	}

	bool noUnfixedHeadTailBool()
	{
		int jumpCodeNum = cigarStringJumpCode.size();
		if(jumpCodeNum < 0)
		{
			return false;
		}

		bool unfixedHeadBool = (cigarStringJumpCode[0].type == "S");
		bool unfixedTailBool = (cigarStringJumpCode[jumpCodeNum - 1].type == "S");

		if((!unfixedHeadBool)&&(!unfixedTailBool))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	bool checkJumpCodeValid()
	{
		#ifdef DETECT_CIRCULAR_RNA
		bool jumpCodeValid = true;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if(cigarStringJumpCode[tmp].len <= 0)
			{
				if(cigarStringJumpCode[tmp].type != "N")	
					return false;	
			}
		}
		return jumpCodeValid;
		#else
		bool jumpCodeValid = true;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if(cigarStringJumpCode[tmp].len <= 0)
				return false;	
		}
		return jumpCodeValid;
		#endif
	}

	Alignment_Info(string alignDir, const string& mapChromName, 
		int mapChromPos, vector<Jump_Code>& cigarString, int mismatch, Index_Info* indexInfo)
	{
		alignDirection = alignDir;
		alignChromName = mapChromName;		
		alignChromPos = mapChromPos;
		cigarStringJumpCode = cigarString;
		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		mismatchNum = mismatch;
		SJstrand = this->getStrandFromSJ();
		//SJstrand
	}	

	void pushBack2MismatchPosCharVec(vector<int>& tmpMismatchPosVec, vector<char>& tmpMismatchCharVec)
	{
		if(STORE_MISMATCH_POS)
		{
			for(int tmp = 0; tmp < tmpMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = tmpMismatchPosVec[tmp];
				mismatchPosVec.push_back(tmpMismatchPos);
			}
			if(STORE_MISMATCH_CHA)
			{
				for(int tmp = 0; tmp < tmpMismatchCharVec.size(); tmp++)
				{
					char tmpMismatchChar = tmpMismatchCharVec[tmp];
					mismatchCharVec.push_back(tmpMismatchChar);
				}
			}
		}
	}

	Alignment_Info(string alignDir, const string& mapChromName, 
		int mapChromPos, vector<Jump_Code>& cigarString, int mismatch, Index_Info* indexInfo,
		vector<int>& tmpMismatchPosVec, vector<char>& tmpMismatchCharVec)
	{
		alignDirection = alignDir;
		alignChromName = mapChromName;		
		alignChromPos = mapChromPos;
		cigarStringJumpCode = cigarString;
		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		mismatchNum = mismatch;
		SJstrand = this->getStrandFromSJ();
		this->pushBack2MismatchPosCharVec(tmpMismatchPosVec, tmpMismatchCharVec);
		//SJstrand
	}

	Alignment_Info(const string& alignInfoStr, const string& alignDir, Index_Info* indexInfo)
	{
		alignDirection = alignDir;

		int tmpFieldStartPosInStr = 0;
		int tmpFieldEndPosInStr;
		string tmpFieldStr;

		//chrom name
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);
		alignChromName = tmpFieldStr;

		//chrom mapping pos
		tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);
		alignChromPos = atoi(tmpFieldStr.c_str());

		//jumpCodeVec
		tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);		
		this->jumpCodeStr2jumpCodeVec(tmpFieldStr);

		//mismatch number
		tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
		tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
		tmpFieldStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
				tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);
		mismatchNum = atoi(tmpFieldStr.c_str());

		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		SJstrand = this->getStrandFromSJ();

		// mismatch pos vec
		if(STORE_MISMATCH_POS)
		{
			tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
			tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
			if(tmpFieldEndPosInStr >= tmpFieldStartPosInStr)
			{
				string tmpMismatchPosVecStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
					tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);			
				this->mismatchPosVecStr2Vec(tmpMismatchPosVecStr);
			}
			// mismatch char vec
			if(STORE_MISMATCH_CHA)
			{
				tmpFieldStartPosInStr = tmpFieldEndPosInStr + 2;
				tmpFieldEndPosInStr = alignInfoStr.find(",",tmpFieldStartPosInStr) - 1;
				if(tmpFieldEndPosInStr >= tmpFieldStartPosInStr)
				{
					string tmpMismatchCharVecStr = alignInfoStr.substr(tmpFieldStartPosInStr, 
						tmpFieldEndPosInStr - tmpFieldStartPosInStr + 1);					
					this->mismatchCharVecStr2Vec(tmpMismatchCharVecStr);
				}
			}
		}
	}

	void mismatchPosVecStr2Vec(const string& tmpMismatchPosVecStr)
	{
		// example: "" / "12" / "1-14-29-42"
		/*int tmpMismatchPosVecStrLen = tmpMismatchPosVecStr.length();
		if(tmpMismatchPosVecStrLen == 1) // only 1 mismatch
		{
			int tmpMismatchPos = atoi(tmpMismatchPosVecStr.c_str());
			mismatchPosVec.push_back(tmpMismatchPos);
		}
		else
		{
			int firstMismatchPos = atoi((tmpMismatchPosVecStr.substr(0,1)).c_str());
			mismatchPosVec.push_back(firstMismatchPos);*/
			int startSearchLoc = 0;
			int foundCommaLoc = 0;
			while(1)
			{
				foundCommaLoc = tmpMismatchPosVecStr.find("-", startSearchLoc);
				if(foundCommaLoc == string::npos)
				{
					string tmpMismatchPosStr = tmpMismatchPosVecStr.substr(startSearchLoc);
					int tmpMismatchPosInt = atoi(tmpMismatchPosStr.c_str());
					mismatchPosVec.push_back(tmpMismatchPosInt);
					break;
				}
				else
				{
					string tmpMismatchPosStr = tmpMismatchPosVecStr.substr(startSearchLoc,
						foundCommaLoc - 1 - startSearchLoc + 1);
					int tmpMismatchPosInt = atoi(tmpMismatchPosStr.c_str());
					mismatchPosVec.push_back(tmpMismatchPosInt);
					startSearchLoc = foundCommaLoc + 1;
				}
			}
		//}
	}

	void mismatchCharVecStr2Vec(const string& tmpMismatchCharVecStr)
	{
		// example: "" / "A" / "C-T-G-A-N"
		int startSearchLoc = 0;
		int foundCommaLoc = 0;
		while(1)
		{
			foundCommaLoc = tmpMismatchCharVecStr.find("-", startSearchLoc);
			if(foundCommaLoc == string::npos)
			{
				char tmpMismatchChar = tmpMismatchCharVecStr.at(startSearchLoc);
				mismatchCharVec.push_back(tmpMismatchChar);
				break;
			}
			else
			{
				char tmpMismatchChar = tmpMismatchCharVecStr.at(startSearchLoc);
				mismatchCharVec.push_back(tmpMismatchChar);
				startSearchLoc = foundCommaLoc + 1;
			}
		}
	}

	void addNewMismatchPosVec(vector<int>& toAddMismatchPosVec, vector<int>& targetNewMismatchPosVec)
	{
		// add original mismatchPosVec
		for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = mismatchPosVec[tmp];
			targetNewMismatchPosVec.push_back(tmpMismatchPos);
		}

		// add other mismatchPosVec
		for(int tmp = 0; tmp < toAddMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = toAddMismatchPosVec[tmp];
			targetNewMismatchPosVec.push_back(tmpMismatchPos);
		}
	}

	void addNewMismatchPosVec(vector<int>& toAddMismatchPosVec, vector<int>& targetNewMismatchPosVec, int toAddPartStartLocInRead)
	{
		// add original mismatchPosVec
		for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = mismatchPosVec[tmp];
			targetNewMismatchPosVec.push_back(tmpMismatchPos);
		}

		// add other mismatchPosVec
		for(int tmp = 0; tmp < toAddMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = toAddMismatchPosVec[tmp] + toAddPartStartLocInRead - 1;
			targetNewMismatchPosVec.push_back(tmpMismatchPos);
		}
	}

	void addNewMismatchCharVec(vector<char>& toAddMismatchCharVec, vector<char>& targetNewMismatchCharVec)
	{
		// add original mismatchCharVec
		for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
		{
			char tmpMismatchChar = mismatchCharVec[tmp];
			targetNewMismatchCharVec.push_back(tmpMismatchChar);
		}

		// add other mismatchCharVec
		for(int tmp = 0; tmp < toAddMismatchCharVec.size(); tmp++)
		{
			char tmpMismatchChar = toAddMismatchCharVec[tmp];
			targetNewMismatchCharVec.push_back(tmpMismatchChar);
		}
	}	

	/*
	Alignment_Info* newAlignInfo_addIncompleteHead(
		Path_Info* pathInfo, int tmpPath, int mapChromPosInt, 
		int mismatchNum, Index_Info* indexInfo)
	{
		//((pathInfo->finalPathVec)[tmpPath].second)
		(pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->addOriginalJumpCode_incompleteHead(cigarStringJumpCode);

		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		if(STORE_MISMATCH_POS)
		{
			this->addNewMismatchPosVec(pathInfo->fixedPathMismatchPosVec[tmpPath], newMismatchPosVec);
			if(STORE_MISMATCH_CHA)
			{
				this->addNewMismatchCharVec(pathInfo->fixedPathMismatchCharVec[tmpPath], newMismatchCharVec);
			}
		}

		Alignment_Info* newAlignInfo = new Alignment_Info(alignDirection, alignChromName,
			mapChromPosInt, (pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
			mismatchNum, indexInfo, newMismatchPosVec, newMismatchCharVec);
		return newAlignInfo;
	}

	Alignment_Info* newAlignInfo_addIncompleteTail(
		Path_Info* pathInfo, int tmpPath, int mapChromPosInt, 
		int mismatchNum, Index_Info* indexInfo, int midPartStartLocInRead)
	{
		(pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->addOriginalJumpCode_incompleteTail(
			cigarStringJumpCode);

		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		if(STORE_MISMATCH_POS)
		{
			this->addNewMismatchPosVec(pathInfo->fixedPathMismatchPosVec[tmpPath], newMismatchPosVec, midPartStartLocInRead);
			if(STORE_MISMATCH_CHA)
			{
				this->addNewMismatchCharVec(pathInfo->fixedPathMismatchCharVec[tmpPath], newMismatchCharVec);
			}
		}

		Alignment_Info* newAlignInfo = new Alignment_Info(alignDirection, alignChromName,
			mapChromPosInt, (pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
			mismatchNum, indexInfo, newMismatchPosVec, newMismatchCharVec);
		return newAlignInfo;
	}*/	

	void addJumpCodeVecToCigarStringJumpCodeVec_startIndex_endIndex_jumpCodeVec(
		int startIndex, int endIndex, vector<Jump_Code>& newJumpCodeVec)
	{
		for(int tmp = startIndex; tmp <= endIndex; tmp++) 
		{
			cigarStringJumpCode.push_back(newJumpCodeVec[tmp]);
		}		
	}	
	void addJumpCodeVecToCigarStringJumpCodeVec_startIndex_endIndex_alignInfo(
		int startIndex, int endIndex, Alignment_Info* alignInfo)
	{
		for(int tmp = startIndex; tmp <= endIndex; tmp++) 
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(tmp));
		}		
	}	

	void newAlignInfo_addIncompleteHead_new(
		Path_Info* pathInfo, int tmpPath, int mapChromPosInt, 
		int mismatchNum, Index_Info* indexInfo, Alignment_Info* alignInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = mapChromPosInt;		

		this->addJumpCodeVecToCigarStringJumpCodeVec_startIndex_endIndex_jumpCodeVec(
			0, ((pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code).size()-1,
			(pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code);
		this->addJumpCodeVecToCigarStringJumpCodeVec_startIndex_endIndex_alignInfo(
			2, alignInfo->returnCigarStringJumpCodeSize()-1,
			alignInfo);// 1st Jump_code: Soft clipping; 2nd jump_code: midPartSeg Match
		

		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		//(pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->addOriginalJumpCode_incompleteHead(cigarStringJumpCode);

		//if(STORE_MISMATCH_POS)
		//{	
			this->generateNewMismatchPosVec_new(pathInfo->fixedPathMismatchPosVec[tmpPath], 
				alignInfo);
			//if(STORE_MISMATCH_CHA)
			//{
				this->generateNewMismatchCharVec_new(pathInfo->fixedPathMismatchCharVec[tmpPath], 
					alignInfo);
			//}
		//}		
		SJstrand = this->getStrandFromSJ();
	}

	void newAlignInfo_addIncompleteTail_new(
		Path_Info* pathInfo, int tmpPath, int mapChromPosInt, 
		int mismatchNum, Index_Info* indexInfo, int midPartStartLocInRead,
		Alignment_Info* alignInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = mapChromPosInt;

		this->addJumpCodeVecToCigarStringJumpCodeVec_startIndex_endIndex_alignInfo(
			0, alignInfo->returnCigarStringJumpCodeSize()-3, alignInfo);

		this->addJumpCodeVecToCigarStringJumpCodeVec_startIndex_endIndex_jumpCodeVec(
			0, ((pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code).size()-1,
			((pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code));

		if(STORE_MISMATCH_POS)
		{
			this->generateNewMismatchPosVec_new(pathInfo->fixedPathMismatchPosVec[tmpPath], alignInfo, midPartStartLocInRead);
			if(STORE_MISMATCH_CHA)
			{
				this->generateNewMismatchCharVec_new(pathInfo->fixedPathMismatchCharVec[tmpPath], alignInfo);
			}
		}
		SJstrand = this->getStrandFromSJ();
	}	

	/*	
	Alignment_Info* newAlignInfo_addIncompleteTail(
		Path_Info* pathInfo, int tmpPath, int mapChromPosInt, 
		int mismatchNum, Index_Info* indexInfo)
	{
		(pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->addOriginalJumpCode_incompleteTail(
			cigarStringJumpCode);

		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		if(STORE_MISMATCH_POS)
		{
			this->addNewMismatchPosVec(pathInfo->fixedPathMismatchPosVec[tmpPath], newMismatchPosVec);
			if(STORE_MISMATCH_CHA)
			{
				this->addNewMismatchCharVec(pathInfo->fixedPathMismatchCharVec[tmpPath], newMismatchCharVec);
			}
		}

		Alignment_Info* newAlignInfo = new Alignment_Info(alignDirection, alignChromName,
			mapChromPosInt, (pathInfo->returnSpliceInfoInFinalPathVec(tmpPath))->final_jump_code, 
			mismatchNum, indexInfo, newMismatchPosVec, newMismatchCharVec);
		return newAlignInfo;
	}*/	

	Alignment_Info* newAlignInfoAfterExtension2HeadSoftClippingAlignment(Index_Info* indexInfo, 
		int unfixedHeadLen, int newMismatchToAddInHead, vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd)
	{
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos - unfixedHeadLen;
		
		int oldfirstMatchLen = cigarStringJumpCode[1].len;
		int newFirstMatchLen = oldfirstMatchLen + unfixedHeadLen;
		Jump_Code newJumpCode(newFirstMatchLen, "M");
		vector<Jump_Code> newCigarStringJumpCode;		
		newCigarStringJumpCode.push_back(newJumpCode);
		for(int tmp = 2; tmp < cigarStringJumpCode.size(); tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}		

		int newMismatch = mismatchNum + newMismatchToAddInHead;
		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		if(STORE_MISMATCH_POS)
		{	
			this->generateNewMismatchPosVec(newMismatchPosVecToAdd, newMismatchPosVec);
			if(STORE_MISMATCH_CHA)
			{
				this->generateNewMismatchCharVec(newMismatchCharVecToAdd, newMismatchCharVec);
			}
		}

		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo,
			newMismatchPosVec, newMismatchCharVec);

		return newAlignInfo;		
	}
	/*
	Alignment_Info* newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated(
		Index_Info* indexInfo, int unfixedHeadLen, 
		int extensionLength, vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd)
	{
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos - extensionLength;

		vector<Jump_Code> newCigarStringJumpCode;
		
		int softClippingLength_new = unfixedHeadLen - extensionLength;
		if(softClippingLength_new > 0)
		{
			Jump_Code newSoftClippingJumpCode(softClippingLength_new, "S");
			newCigarStringJumpCode.push_back(newSoftClippingJumpCode);
		}

		int old2ndJumpCodeLen = cigarStringJumpCode[1].len;
		string old2ndJumpCodeType = cigarStringJumpCode[1].type;

		if(old2ndJumpCodeType == "M") // 2nd JumpCode of old align is "match"
		{
			int newFirstMatchLength = extensionLength + old2ndJumpCodeLen;
			Jump_Code newFirstMatchJumpCode(newFirstMatchLength, "M");			
			newCigarStringJumpCode.push_back(newFirstMatchJumpCode);	
		}
		else
		{
			Jump_Code newFirstMatchJumpCode(extensionLength, "M");			
			newCigarStringJumpCode.push_back(newFirstMatchJumpCode);
			newCigarStringJumpCode.push_back(cigarStringJumpCode[1]);
		}
	
		for(int tmp = 2; tmp < cigarStringJumpCode.size(); tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}

		int newMismatch = mismatchNum + newMismatchPosVecToAdd.size();
		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		//if(STORE_MISMATCH_POS)
		//{	
			this->generateNewMismatchPosVec(newMismatchPosVecToAdd, newMismatchPosVec);
			if(STORE_MISMATCH_CHA)
			{
				this->generateNewMismatchCharVec(newMismatchCharVecToAdd, newMismatchCharVec);
			}
		//}

		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo,
			newMismatchPosVec, newMismatchCharVec);

		return newAlignInfo;		
	}*/
	
	void newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated_new(
		Index_Info* indexInfo, int unfixedHeadLen, 
		int extensionLength, vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd, Alignment_Info* alignInfo)
	{
		//cout << "newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated_new starts " << endl;
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos()  - extensionLength;

		//vector<Jump_Code> newCigarStringJumpCode;
		
		int softClippingLength_new = unfixedHeadLen - extensionLength;
		//cout << "unfixedHeadLen: " << 
		if(softClippingLength_new > 0)
		{
			Jump_Code newSoftClippingJumpCode(softClippingLength_new, "S");
			cigarStringJumpCode.push_back(newSoftClippingJumpCode);
		}

		int old2ndJumpCodeLen = alignInfo->returnCigarStringJumpCodeLen(1);
		string old2ndJumpCodeType = alignInfo->returnCigarStringJumpCodeType(1);

		if(old2ndJumpCodeType == "M") // 2nd JumpCode of old align is "match"
		{
			int newFirstMatchLength = extensionLength + old2ndJumpCodeLen;
			Jump_Code newFirstMatchJumpCode(newFirstMatchLength, "M");			
			cigarStringJumpCode.push_back(newFirstMatchJumpCode);	
		}
		else
		{
			Jump_Code newFirstMatchJumpCode(extensionLength, "M");			
			cigarStringJumpCode.push_back(newFirstMatchJumpCode);
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(1));
		}
	
		for(int tmp = 2; tmp < alignInfo->returnCigarStringJumpCodeSize(); tmp++)
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(tmp));
		}
		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		mismatchNum = alignInfo->returnMismatchNum() + newMismatchPosVecToAdd.size();

		//if(STORE_MISMATCH_POS)
		//{	
			this->generateNewMismatchPosVec_new(newMismatchPosVecToAdd, 
				alignInfo);
		//	if(STORE_MISMATCH_CHA)
		//	{
				this->generateNewMismatchCharVec_new(newMismatchCharVecToAdd, 
					alignInfo);
		//	}
		//}		
		SJstrand = this->getStrandFromSJ();
	}	


	void newAlignInfoAfterExtension2Head_fixIndelBesideReadEnd(
		Index_Info* indexInfo, vector<Jump_Code>& fixedIndelJumpCodeVec, 
		int fixedIndelChromMapPos, Alignment_Info* alignInfo)
	{
		//cout << "newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated_new starts " << endl;
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = fixedIndelChromMapPos;

		//vector<Jump_Code> newCigarStringJumpCode;
		
		// int softClippingLength_new = unfixedHeadLen - extensionLength;
		// //cout << "unfixedHeadLen: " << 
		// if(softClippingLength_new > 0)
		// {
		// 	Jump_Code newSoftClippingJumpCode(softClippingLength_new, "S");
		// 	cigarStringJumpCode.push_back(newSoftClippingJumpCode);
		// }
		for(int tmp = 0; tmp < fixedIndelJumpCodeVec.size()-1; tmp++)
		{
			cigarStringJumpCode.push_back(fixedIndelJumpCodeVec[tmp]);
		}
		int lastMatchLen_extendedHead_len = fixedIndelJumpCodeVec[fixedIndelJumpCodeVec.size()-1].len;

		int old2ndJumpCodeLen = alignInfo->returnCigarStringJumpCodeLen(1);
		string old2ndJumpCodeType = alignInfo->returnCigarStringJumpCodeType(1);

		if(old2ndJumpCodeType == "M") // 2nd JumpCode of old align is "match"
		{
			int newFirstMatchLength = lastMatchLen_extendedHead_len + old2ndJumpCodeLen;
			Jump_Code newFirstMatchJumpCode(newFirstMatchLength, "M");			
			cigarStringJumpCode.push_back(newFirstMatchJumpCode);	
		}
		else
		{
			Jump_Code newFirstMatchJumpCode_midPart(lastMatchLen_extendedHead_len, "M");			
			cigarStringJumpCode.push_back(newFirstMatchJumpCode_midPart);
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(1));
		}
	
		for(int tmp = 2; tmp < alignInfo->returnCigarStringJumpCodeSize(); tmp++)
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(tmp));
		}
		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		mismatchNum = alignInfo->returnMismatchNum(); //+ newMismatchPosVecToAdd.size();

		this->copyNewMismatchPosVec(
			alignInfo);
		this->copyNewMismatchCharVec( 
			alignInfo);
	
		SJstrand = this->getStrandFromSJ();
	}	

	void newAlignInfoAfterExtension2Tail_fixIndelBesideReadEnd(
		Index_Info* indexInfo, vector<Jump_Code>& fixedIndelJumpCodeVec, 
		Alignment_Info* alignInfo)
	{
		//cout << "newAlignInfoAfterExtension2HeadSoftClippingAlignment_errorTolerated_new starts " << endl;
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos();

		int jumpCodeSize = alignInfo->returnCigarStringJumpCodeSize();
		for(int tmp = 0; tmp < jumpCodeSize - 2; tmp++)
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(tmp));
		}

		int penultimateJumpCodeLength = alignInfo->returnCigarStringJumpCodeLen(jumpCodeSize-2);
		string penultimateJumpCodeType = alignInfo->returnCigarStringJumpCodeType(jumpCodeSize-2);

		int firstMatchLen_extendedTail_len = fixedIndelJumpCodeVec[0].len;
		if(penultimateJumpCodeType == "M")
		{
			int newLastMatchLength_midPart = penultimateJumpCodeLength + firstMatchLen_extendedTail_len;
			Jump_Code newLastMatchJumpCode_midPart(newLastMatchLength_midPart, "M");
			cigarStringJumpCode.push_back(newLastMatchJumpCode_midPart);
		}
		else
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(jumpCodeSize-2));
			Jump_Code newFirstMatchJumpCode_extendedTail(firstMatchLen_extendedTail_len, "M");
			cigarStringJumpCode.push_back(newFirstMatchJumpCode_extendedTail);
		}


		for(int tmp = 1; tmp < fixedIndelJumpCodeVec.size(); tmp++)
		{
			cigarStringJumpCode.push_back(fixedIndelJumpCodeVec[tmp]);
		}

		this->jumpCodeVec2spliceJunctionVec(indexInfo);
		mismatchNum = alignInfo->returnMismatchNum(); //+ newMismatchPosVecToAdd.size();

		this->copyNewMismatchPosVec(
			alignInfo);
		this->copyNewMismatchCharVec( 
			alignInfo);
	
		SJstrand = this->getStrandFromSJ();
	}	

	int getChrPosLength(vector<Jump_Code>& tmpJumpCodeVec)
	{
		int tmpChrPosLength = 0;
		int tmpJumpCodeVecSize = tmpJumpCodeVec.size();
		for(int tmp = 0; tmp < tmpJumpCodeVec.size(); tmp++)
		{
			string tmpJumpCodeType = tmpJumpCodeVec[tmp].type;
			int tmpJumpCodeLength = tmpJumpCodeVec[tmp].len;
			if((tmpJumpCodeType == "M")
				||(tmpJumpCodeType == "N")
				||(tmpJumpCodeType == "D"))
				tmpChrPosLength += tmpJumpCodeLength;
		}
		return tmpChrPosLength;
	}

	void newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_withAlignInfer(
		vector<Jump_Code>& inferedPathJumpCodeVec,
		vector<int>& inferedPathMismatchPosVec, 
		vector<char>& inferedPathMismatchCharVec,
		Index_Info* indexInfo, Alignment_Info* alignInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		// FIX ME; fix the alignChromPos
		int extendedChrPos = this->getChrPosLength(inferedPathJumpCodeVec);
		alignChromPos = alignInfo->returnAlignChromPos() - extendedChrPos;

		for(int tmp = 0; tmp < inferedPathJumpCodeVec.size()-1; tmp++)
		{
			cigarStringJumpCode.push_back(inferedPathJumpCodeVec[tmp]);
		}
		int new1stMatchLengthInOriPartialAlignment
			= alignInfo->returnCigarStringJumpCodeLen(1) 
				+ inferedPathJumpCodeVec[inferedPathJumpCodeVec.size()-1].len;
		Jump_Code new1stMatchJumpCodeInOriPartialAlignment(new1stMatchLengthInOriPartialAlignment, "M");  
		cigarStringJumpCode.push_back(new1stMatchJumpCodeInOriPartialAlignment);
		int tmpOriJumpCodeVecSize 
			= alignInfo->returnCigarStringJumpCodeSize();
		for(int tmp = 2; tmp < tmpOriJumpCodeVecSize; tmp++)
		{
			cigarStringJumpCode.push_back(
				alignInfo->returnCigarStringJumpCodeElement(tmp)
				);
		}		
		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		int tmpMismatchNum = inferedPathMismatchPosVec.size();
		mismatchNum = alignInfo->returnMismatchNum() + tmpMismatchNum;	
		this->generateNewMismatchPosVec_new(inferedPathMismatchPosVec,
			alignInfo);
		this->generateNewMismatchCharVec_new(inferedPathMismatchCharVec,
			alignInfo);
		SJstrand = this->getStrandFromSJ();
	}	


	void newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new(
		int firstMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, 
		int newMismatchToAddInHead,
		vector<int>& newMismatchPosVecToAdd, vector<char>& newMismatchCharVecToAdd,
		Alignment_Info* alignInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos() 
			- spliceJunctionDistance - alignInfo->returnCigarStringJumpCodeLen(0);

		if(spliceJunctionDistance > MAX_DELETION_LENGTH)
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
		
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			cigarStringJumpCode.push_back(tmpSjJumpCode);

			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}	
		else if (spliceJunctionDistance > 0)
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "D");
			cigarStringJumpCode.push_back(tmpSjJumpCode);
		
			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}
		else if (spliceJunctionDistance == 0)
		{
			Jump_Code tmpMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				, "M");
			cigarStringJumpCode.push_back(tmpMatchJumpCode);
		}
		else if (spliceJunctionDistance < 0)
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpInsJumpCode(0 - spliceJunctionDistance, "I");
			cigarStringJumpCode.push_back(tmpInsJumpCode);

			int tmpInsertionLength = 0 - spliceJunctionDistance;
			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength - tmpInsertionLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}			
		else
		{}

		int tmpJumpCodeVecSize = alignInfo->returnCigarStringJumpCodeSize();

		for(int tmp = 2; tmp < tmpJumpCodeVecSize; tmp++)
		{
			cigarStringJumpCode.push_back(
				alignInfo->returnCigarStringJumpCodeElement(tmp));
		}
		
		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		int tmpMismatchNum;

		if(spliceJunctionDistance < 0)
		{
			int insertionStartPosInRead = firstMatchLength + 1; 
			int insertionEndPosInRead = firstMatchLength - spliceJunctionDistance;
			tmpMismatchNum = 
				this->generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos_new(
					newMismatchPosVecToAdd, newMismatchCharVecToAdd,
					alignInfo,
					insertionStartPosInRead, insertionEndPosInRead);
		}
		else
		{
			//newMismatch = mismatchNum + newMismatchToAddInHead;
			tmpMismatchNum = alignInfo->returnMismatchPosVecSize();
				+ newMismatchPosVecToAdd.size();	
			this->generateNewMismatchPosVec_new(newMismatchPosVecToAdd, 
				alignInfo);
			this->generateNewMismatchCharVec_new(newMismatchCharVecToAdd, 
				alignInfo);
		}
		mismatchNum = tmpMismatchNum;	
		SJstrand = this->getStrandFromSJ();
	}

	void newAlignInfoAfterAddHeadInfo2SoftClippingAlignment_new_backSpliceIncluded(
		int firstMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, 
		int newMismatchToAddInHead,
		vector<int>& newMismatchPosVecToAdd, vector<char>& newMismatchCharVecToAdd,
		Alignment_Info* alignInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos() 
			- spliceJunctionDistance - alignInfo->returnCigarStringJumpCodeLen(0);

		if(spliceJunctionDistance > MAX_DELETION_LENGTH) // normal splice junction
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
		
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			cigarStringJumpCode.push_back(tmpSjJumpCode);

			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}	
		else if (spliceJunctionDistance > 0) // deletion
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "D");
			cigarStringJumpCode.push_back(tmpSjJumpCode);
		
			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}
		else if (spliceJunctionDistance == 0) // match
		{
			Jump_Code tmpMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				, "M");
			cigarStringJumpCode.push_back(tmpMatchJumpCode);
		}
		else if (0 - spliceJunctionDistance <= MAX_INSERTION_LENGTH) 
		// insertion
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpInsJumpCode(0 - spliceJunctionDistance, "I");
			cigarStringJumpCode.push_back(tmpInsJumpCode);

			int tmpInsertionLength = 0 - spliceJunctionDistance;
			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength - tmpInsertionLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}			
		else if(0 - spliceJunctionDistance > MAX_INSERTION_LENGTH)
		// back splice junction, similar to normal splice junction
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
		
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			cigarStringJumpCode.push_back(tmpSjJumpCode);

			Jump_Code tmp2ndMatchJumpCode(
				alignInfo->returnCigarStringJumpCodeLen(0)
				+ alignInfo->returnCigarStringJumpCodeLen(1)
				- firstMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);			
		}
		else
		{}

		int tmpJumpCodeVecSize = alignInfo->returnCigarStringJumpCodeSize();

		for(int tmp = 2; tmp < tmpJumpCodeVecSize; tmp++)
		{
			cigarStringJumpCode.push_back(
				alignInfo->returnCigarStringJumpCodeElement(tmp));
		}
		
		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		int tmpMismatchNum;

		if((spliceJunctionDistance < 0)&&(0 - spliceJunctionDistance <= MAX_INSERTION_LENGTH)) // INSERTION
		{
			int insertionStartPosInRead = firstMatchLength + 1; 
			int insertionEndPosInRead = firstMatchLength - spliceJunctionDistance;
			tmpMismatchNum = 
				this->generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos_new(
					newMismatchPosVecToAdd, newMismatchCharVecToAdd,
					alignInfo,
					insertionStartPosInRead, insertionEndPosInRead);
		}
		else
		{
			//newMismatch = mismatchNum + newMismatchToAddInHead;
			tmpMismatchNum = alignInfo->returnMismatchPosVecSize();
				+ newMismatchPosVecToAdd.size();	
			this->generateNewMismatchPosVec_new(newMismatchPosVecToAdd, 
				alignInfo);
			this->generateNewMismatchCharVec_new(newMismatchCharVecToAdd, 
				alignInfo);
		}
		mismatchNum = tmpMismatchNum;
		SJstrand = this->getStrandFromSJ();
	}
	
	int generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos_new(
		vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd,
		Alignment_Info* alignInfo,
		int insertionStartPosInRead, int insertionEndPosInRead)
	{
		if(STORE_MISMATCH_POS && STORE_MISMATCH_CHA)
		{
			for(int tmp = 0; 
					tmp < alignInfo->returnMismatchPosVecSize(); 
					tmp++)
			{
				int tmpMismatchPos 
					= alignInfo->returnMismatchPosVecValue(tmp);
				char tmpMismatchChar 
					= alignInfo->returnMismatchCharVecValue(tmp);

				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					mismatchPosVec.push_back(tmpMismatchPos);
					mismatchCharVec.push_back(tmpMismatchChar);
				}
			}
			for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
			{
				int tmpMismatchPos = newMismatchPosVecToAdd[tmp];
				char tmpMismatchChar = newMismatchCharVecToAdd[tmp];
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					mismatchPosVec.push_back(tmpMismatchPos);
					mismatchCharVec.push_back(tmpMismatchChar);
				}
			}			
		}
		else if(STORE_MISMATCH_POS)
		{
			for(int tmp = 0; 
					tmp < alignInfo->returnMismatchPosVecSize(); 
					tmp++)
			{
				int tmpMismatchPos 
					= alignInfo->returnMismatchPosVecValue(tmp);
					
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					mismatchPosVec.push_back(tmpMismatchPos);
				}
			}
			for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
			{
				int tmpMismatchPos = newMismatchPosVecToAdd[tmp];
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					mismatchPosVec.push_back(tmpMismatchPos);
				}
			}			
		}
		else
		{
			cout << "STORE_MISMATCH_POS must be set as true to get the correct mismatch number" << endl;
			exit(1);
		}		
		return mismatchPosVec.size();
	}

	Alignment_Info* newAlignInfoAfterAddHeadInfo2SoftClippingAlignment(
		int firstMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, 
		int newMismatchToAddInHead,
		vector<int>& newMismatchPosVecToAdd, vector<char>& newMismatchCharVecToAdd)
	{
		//cout << "start to get new alignInfo" << endl;
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos - spliceJunctionDistance - cigarStringJumpCode[0].len;

		vector<Jump_Code> newCigarStringJumpCode;

		if(spliceJunctionDistance > MAX_DELETION_LENGTH)
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp1stMatchJumpCode);
		
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			newCigarStringJumpCode.push_back(tmpSjJumpCode);

			Jump_Code tmp2ndMatchJumpCode(cigarStringJumpCode[0].len + 
				cigarStringJumpCode[1].len - firstMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}	
		else if (spliceJunctionDistance > 0)
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "D");
			newCigarStringJumpCode.push_back(tmpSjJumpCode);
		
			Jump_Code tmp2ndMatchJumpCode(cigarStringJumpCode[0].len + 
				cigarStringJumpCode[1].len - firstMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}
		else if (spliceJunctionDistance == 0)
		{
			Jump_Code tmpMatchJumpCode(cigarStringJumpCode[0].len + 
				cigarStringJumpCode[1].len, "M");
			newCigarStringJumpCode.push_back(tmpMatchJumpCode);
		}
		else if (spliceJunctionDistance < 0)
		{
			Jump_Code tmp1stMatchJumpCode(firstMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpInsJumpCode(0 - spliceJunctionDistance, "I");
			newCigarStringJumpCode.push_back(tmpInsJumpCode);

			int tmpInsertionLength = 0 - spliceJunctionDistance;
			Jump_Code tmp2ndMatchJumpCode(cigarStringJumpCode[0].len + 
				cigarStringJumpCode[1].len - firstMatchLength - tmpInsertionLength, "M");
			newCigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}			
		else
		{}

		for(int tmp = 2; tmp < cigarStringJumpCode.size(); tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}
		
		int newMismatch;

		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		if(spliceJunctionDistance < 0)
		{
			int insertionStartPosInRead = firstMatchLength + 1; 
			int insertionEndPosInRead = firstMatchLength - spliceJunctionDistance;
			newMismatch = this->generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos(
				newMismatchPosVecToAdd, newMismatchPosVec,
				newMismatchCharVecToAdd, newMismatchCharVec,
				insertionStartPosInRead, insertionEndPosInRead);
		}
		else
		{
			newMismatch = mismatchNum + newMismatchToAddInHead;
			//if(STORE_MISMATCH_POS)
			//{	
				this->generateNewMismatchPosVec(newMismatchPosVecToAdd, newMismatchPosVec);
				//if(STORE_MISMATCH_CHA)
				//{
					this->generateNewMismatchCharVec(newMismatchCharVecToAdd, newMismatchCharVec);
				//}
			//}
		}
		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo,
			newMismatchPosVec, newMismatchCharVec);
		return newAlignInfo;	
	}

	int generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos(
		vector<int>& newMismatchPosVecToAdd, vector<int>& newMismatchPosVec,
		vector<char>& newMismatchCharVecToAdd, vector<char>& newMismatchCharVec,
		int insertionStartPosInRead, int insertionEndPosInRead)
	{
		if(STORE_MISMATCH_POS && STORE_MISMATCH_CHA)
		{
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = mismatchPosVec[tmp];
				char tmpMismatchChar = mismatchCharVec[tmp];
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					newMismatchPosVec.push_back(tmpMismatchPos);
					newMismatchCharVec.push_back(tmpMismatchChar);
				}
			}
			for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
			{
				int tmpMismatchPos = newMismatchPosVecToAdd[tmp];
				char tmpMismatchChar = newMismatchCharVecToAdd[tmp];
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					newMismatchPosVec.push_back(tmpMismatchPos);
					newMismatchCharVec.push_back(tmpMismatchChar);
				}
			}			
		}
		else if(STORE_MISMATCH_POS)
		{
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = mismatchPosVec[tmp];
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					newMismatchPosVec.push_back(tmpMismatchPos);
				}
			}
			for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
			{
				int tmpMismatchPos = newMismatchPosVecToAdd[tmp];
				if((tmpMismatchPos < insertionStartPosInRead)||(tmpMismatchPos > insertionEndPosInRead))
				{
					newMismatchPosVec.push_back(tmpMismatchPos);
				}
			}				
		}
		else
		{
			cout << "STORE_MISMATCH_POS must be set as true to get the correct mismatch number" << endl;
			exit(1);
		}		
		return newMismatchPosVec.size();
	}

	void generateNewMismatchPosVec_new(
		vector<int>& newMismatchPosVecToAdd, 
		Alignment_Info* alignInfo)
	{
		for(int tmp = 0; tmp < alignInfo->returnMismatchPosVecSize(); 
			tmp++)
		{
			int tmpMismatchPos = alignInfo->returnMismatchPosVecValue(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);
		}

		for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
		{
			int tmpMismatchPos = newMismatchPosVecToAdd[tmp];
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}

	void copyNewMismatchPosVec(Alignment_Info* alignInfo)
	{
		for(int tmp = 0; tmp < alignInfo->returnMismatchPosVecSize(); 
			tmp++)
		{
			int tmpMismatchPos = alignInfo->returnMismatchPosVecValue(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);
		}	
	}

	void generateNewMismatchPosVec_new(
		vector<int>& newMismatchPosVecToAdd, 
		Alignment_Info* alignInfo, int toAddPartStartLocInRead)
	{
		for(int tmp = 0; tmp < alignInfo->returnMismatchPosVecSize(); 
			tmp++)
		{
			int tmpMismatchPos = alignInfo->returnMismatchPosVecValue(tmp);
			mismatchPosVec.push_back(tmpMismatchPos);
		}

		for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
		{
			int tmpMismatchPos = newMismatchPosVecToAdd[tmp] + toAddPartStartLocInRead - 1;
			mismatchPosVec.push_back(tmpMismatchPos);
		}
	}

	void generateNewMismatchCharVec_new(
		vector<char>& newMismatchCharVecToAdd, 
		Alignment_Info* alignInfo)
	{
		for(int tmp = 0; tmp < alignInfo->returnMismatchPosVecSize(); 
			tmp++)
		{
			char tmpMismatchChar = alignInfo->returnMismatchCharVecValue(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}
		
		for(int tmp = 0; tmp < newMismatchCharVecToAdd.size(); tmp++)
		{
			char tmpMismatchChar = newMismatchCharVecToAdd[tmp];
			mismatchCharVec.push_back(tmpMismatchChar);
		}
	}	

	void copyNewMismatchCharVec(Alignment_Info* alignInfo)
	{
		for(int tmp = 0; tmp < alignInfo->returnMismatchPosVecSize(); 
			tmp++)
		{
			char tmpMismatchChar = alignInfo->returnMismatchCharVecValue(tmp);
			mismatchCharVec.push_back(tmpMismatchChar);
		}		
	}


	void generateNewMismatchPosVec(vector<int>& newMismatchPosVecToAdd, vector<int>& newMismatchPosVec)
	{
		for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = mismatchPosVec[tmp];
			newMismatchPosVec.push_back(tmpMismatchPos);
		}
		for(int tmp = 0; tmp < newMismatchPosVecToAdd.size(); tmp++)
		{
			int tmpMismatchPos = newMismatchPosVecToAdd[tmp];
			newMismatchPosVec.push_back(tmpMismatchPos);
		}
	}

	void generateNewMismatchCharVec(vector<char>& newMismatchCharVecToAdd, vector<char>& newMismatchCharVec)
	{
		for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
		{
			char tmpMismatchChar = mismatchCharVec[tmp];
			newMismatchCharVec.push_back(tmpMismatchChar);
		}
		for(int tmp = 0; tmp < newMismatchCharVecToAdd.size(); tmp++)
		{
			char tmpMismatchChar = newMismatchCharVecToAdd[tmp];
			newMismatchCharVec.push_back(tmpMismatchChar);
		}
	}

	/*
	Alignment_Info* newAlignInfoAfterExtension2TailSoftClippingAlignment_errorTolerated(Index_Info* indexInfo, 
		int unfixedTailLen, int extensionLength, vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd)
	{
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos;// - spliceJunctionDistance - cigarStringJumpCode[0].len;
		
		//cout << "start to creat new jumpCode " << endl;
		vector<Jump_Code> newCigarStringJumpCode;

		int jumpCodeSize = cigarStringJumpCode.size();
		for(int tmp = 0; tmp < jumpCodeSize - 2; tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}

		int penultimateJumpCodeLength = cigarStringJumpCode[jumpCodeSize-2].len;
		string penultimateJumpCodeType = cigarStringJumpCode[jumpCodeSize-2].type;

		if(penultimateJumpCodeType == "M")
		{
			int newLastMatchLength = penultimateJumpCodeLength + extensionLength;
			Jump_Code newLastMatchJumpCode(newLastMatchLength, "M");
			newCigarStringJumpCode.push_back(newLastMatchJumpCode);
		}
		else
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[jumpCodeSize-2]);
			Jump_Code newExtensionMatchJumpCode(extensionLength, "M");
			newCigarStringJumpCode.push_back(newExtensionMatchJumpCode);
		}

		int softClippingLength_new = unfixedTailLen - extensionLength;
		if(softClippingLength_new > 0)
		{
			Jump_Code newSoftClippingJumpCode(softClippingLength_new, "S");
			newCigarStringJumpCode.push_back(newSoftClippingJumpCode);
		}

		int newMismatch = mismatchNum + newMismatchPosVecToAdd.size();
		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		//if(STORE_MISMATCH_POS)
		//{	
			this->generateNewMismatchPosVec(newMismatchPosVecToAdd, newMismatchPosVec);
			if(STORE_MISMATCH_CHA)
			{
				this->generateNewMismatchCharVec(newMismatchCharVecToAdd, newMismatchCharVec);
			}
		//}
		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo,
			newMismatchPosVec, newMismatchCharVec);

		return newAlignInfo;		
	}*/

	void newAlignInfoAfterExtension2TailSoftClippingAlignment_errorTolerated_new(Index_Info* indexInfo, 
		int unfixedTailLen, int extensionLength, vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd, Alignment_Info* alignInfo)
	{
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos();
		
		//cout << "start to creat new jumpCode " << endl;
		//vector<Jump_Code> newCigarStringJumpCode;

		int jumpCodeSize = alignInfo->returnCigarStringJumpCodeSize();
		for(int tmp = 0; tmp < jumpCodeSize - 2; tmp++)
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(tmp));
		}

		int penultimateJumpCodeLength = alignInfo->returnCigarStringJumpCodeLen(jumpCodeSize-2);
		string penultimateJumpCodeType = alignInfo->returnCigarStringJumpCodeType(jumpCodeSize-2);

		if(penultimateJumpCodeType == "M")
		{
			int newLastMatchLength = penultimateJumpCodeLength + extensionLength;
			Jump_Code newLastMatchJumpCode(newLastMatchLength, "M");
			cigarStringJumpCode.push_back(newLastMatchJumpCode);
		}
		else
		{
			cigarStringJumpCode.push_back(alignInfo->returnCigarStringJumpCodeElement(jumpCodeSize-2));
			Jump_Code newExtensionMatchJumpCode(extensionLength, "M");
			cigarStringJumpCode.push_back(newExtensionMatchJumpCode);
		}

		int softClippingLength_new = unfixedTailLen - extensionLength;
		if(softClippingLength_new > 0)
		{
			Jump_Code newSoftClippingJumpCode(softClippingLength_new, "S");
			cigarStringJumpCode.push_back(newSoftClippingJumpCode);
		}

		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		mismatchNum = alignInfo->returnMismatchNum() + newMismatchPosVecToAdd.size();
		//if(STORE_MISMATCH_POS)
		//{	
			this->generateNewMismatchPosVec_new(newMismatchPosVecToAdd, alignInfo);
			//if(STORE_MISMATCH_CHA)
			//{
				this->generateNewMismatchCharVec_new(newMismatchCharVecToAdd, alignInfo);
			//}
		//}	
		SJstrand = this->getStrandFromSJ();
	}



	Alignment_Info* newAlignInfoAfterExtension2TailSoftClippingAlignment(Index_Info* indexInfo, 
		int unfixedTailLen, int newMismatchToAddInTail, vector<int>& newMismatchPosVecToAdd,
		vector<char>& newMismatchCharVecToAdd)
	{
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos;// - spliceJunctionDistance - cigarStringJumpCode[0].len;
		
		//cout << "start to creat new jumpCode " << endl;
		vector<Jump_Code> newCigarStringJumpCode;

		for(int tmp = 0; tmp < cigarStringJumpCode.size()-2; tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}		
		int jumpCodeSize = cigarStringJumpCode.size();
		int penultimateMatchLength = cigarStringJumpCode[jumpCodeSize-2].len;
		int lastMatchLength = penultimateMatchLength + unfixedTailLen;
		Jump_Code newJumpCode(lastMatchLength, "M");
		newCigarStringJumpCode.push_back(newJumpCode);
	
		int newMismatch = mismatchNum + newMismatchToAddInTail;
		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;

		if(STORE_MISMATCH_POS)
		{	
			this->generateNewMismatchPosVec(newMismatchPosVecToAdd, newMismatchPosVec);
			if(STORE_MISMATCH_CHA)
			{
				this->generateNewMismatchCharVec(newMismatchCharVecToAdd, newMismatchCharVec);
			}
		}
		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo,
			newMismatchPosVec, newMismatchCharVec);

		return newAlignInfo;		
	}

	Alignment_Info* newAlignInfoAfterAddTailInfo2SoftClippingAlignment(
		int lastMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, int newMismatchToAddInTail,
		vector<int>& newMismatchPosVecToAdd, vector<char>& newMismatchCharVecToAdd, int readLength)
	{
		//cout << "start to get new alignInfo" << endl;
		string newAlignDir = alignDirection;
		string newAlignChromName = alignChromName;
		int newAlignChromPos = alignChromPos;// - spliceJunctionDistance - cigarStringJumpCode[0].len;
		
		//cout << "start to creat new jumpCode " << endl;
		vector<Jump_Code> newCigarStringJumpCode;

		for(int tmp = 0; tmp < cigarStringJumpCode.size()-2; tmp++)
		{
			newCigarStringJumpCode.push_back(cigarStringJumpCode[tmp]);
		}

		int tmpInsertionLength = 0;
		if(spliceJunctionDistance < 0)
		{
			tmpInsertionLength = 0 - spliceJunctionDistance;
		}

		int jumpCodeSize = cigarStringJumpCode.size();
		int penultimateMatchLength 
			= cigarStringJumpCode[jumpCodeSize-2].len + cigarStringJumpCode[jumpCodeSize-1].len - lastMatchLength - tmpInsertionLength; 

		if(spliceJunctionDistance > MAX_DELETION_LENGTH)
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			newCigarStringJumpCode.push_back(tmpSjJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}	
		else if (spliceJunctionDistance > 0)
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpDelJumpCode(spliceJunctionDistance, "D");
			newCigarStringJumpCode.push_back(tmpDelJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp2ndMatchJumpCode);			
		}
		else if (spliceJunctionDistance == 0)
		{
			Jump_Code tmpJumpCode(penultimateMatchLength + lastMatchLength, "M");
			newCigarStringJumpCode.push_back(tmpJumpCode);
		}
		else if (spliceJunctionDistance < 0)
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpInsJumpCode(0 - spliceJunctionDistance, "I");
			newCigarStringJumpCode.push_back(tmpInsJumpCode);

			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			newCigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
			//tmpInsertionLength = 0 - spliceJunctionDistance;
		}
		else
		{

		}

		int newMismatch;// = mismatchNum + newMismatchToAddInTail;
		vector<int> newMismatchPosVec;
		vector<char> newMismatchCharVec;
		if(spliceJunctionDistance < 0)
		{
			int insertionStartPosInRead = readLength - lastMatchLength + spliceJunctionDistance + 1;
			int insertionEndPosInRead = readLength - lastMatchLength;
			newMismatch = this->generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos(
				newMismatchPosVecToAdd, newMismatchPosVec,
				newMismatchCharVecToAdd, newMismatchCharVec,
				insertionStartPosInRead, insertionEndPosInRead);
		}
		else
		{
			newMismatch = mismatchNum + newMismatchToAddInTail;
			if(STORE_MISMATCH_POS)
			{	
				this->generateNewMismatchPosVec(newMismatchPosVecToAdd, newMismatchPosVec);
				if(STORE_MISMATCH_CHA)
				{
					this->generateNewMismatchCharVec(newMismatchCharVecToAdd, newMismatchCharVec);
				}
			}	
		}
		Alignment_Info* newAlignInfo = new Alignment_Info(newAlignDir, newAlignChromName,
			newAlignChromPos, newCigarStringJumpCode, newMismatch, indexInfo,
			newMismatchPosVec, newMismatchCharVec);
		return newAlignInfo;		
	}


	void newAlignInfoAfterAddTailInfo2SoftClippingAlignment_withAlignInfer(
		vector<Jump_Code>& inferedPathJumpCodeVec,
		vector<int>& inferedPathMismatchPosVec, 
		vector<char>& inferedPathMismatchCharVec,		
		Index_Info* indexInfo, Alignment_Info* alignInfo)
	{
		//cout << "start to get new alignInfo" << endl;
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos();

		//cout << "start to creat new jumpCode " << endl;
		int tmpJumpCodeVecSize 
			= alignInfo->returnCigarStringJumpCodeSize();		
		for(int tmp = 0; tmp < tmpJumpCodeVecSize-2; tmp++)
		{
			cigarStringJumpCode.push_back(
				alignInfo->returnCigarStringJumpCodeElement(tmp)
				);
		}
		int newLastMatchLengthInOriPartialAlignment
			= alignInfo->returnCigarStringJumpCodeLen(tmpJumpCodeVecSize-2)
				+ inferedPathJumpCodeVec[0].len;
		Jump_Code newLastMatchJumpCodeInOriPartialAlignment(
			newLastMatchLengthInOriPartialAlignment, "M");
		cigarStringJumpCode.push_back(newLastMatchJumpCodeInOriPartialAlignment);
		for(int tmp = 1; tmp < inferedPathJumpCodeVec.size(); tmp++)
		{
			cigarStringJumpCode.push_back(inferedPathJumpCodeVec[tmp]);
		}

		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		int tmpMismatchNum = inferedPathMismatchPosVec.size();
		mismatchNum = alignInfo->returnMismatchNum() + tmpMismatchNum;	
		this->generateNewMismatchPosVec_new(inferedPathMismatchPosVec,
			alignInfo);
		this->generateNewMismatchCharVec_new(inferedPathMismatchCharVec,
			alignInfo);
		SJstrand = this->getStrandFromSJ();
	}


	void newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new(
		int lastMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, int newMismatchToAddInTail,
		vector<int>& newMismatchPosVecToAdd, vector<char>& newMismatchCharVecToAdd, int readLength,
		Alignment_Info* alignInfo)
	{
		//cout << "start to get new alignInfo" << endl;
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos();

		//cout << "start to creat new jumpCode " << endl;
		int tmpJumpCodeVecSize 
			= alignInfo->returnCigarStringJumpCodeSize();		
		for(int tmp = 0; tmp < tmpJumpCodeVecSize-2; tmp++)
		{
			cigarStringJumpCode.push_back(
				alignInfo->returnCigarStringJumpCodeElement(tmp)
				);
		}

		int tmpInsertionLength = 0;
		if(spliceJunctionDistance < 0)
		{
			tmpInsertionLength = 0 - spliceJunctionDistance;
		}

		int jumpCodeSize = tmpJumpCodeVecSize;
		int penultimateMatchLength 
			= alignInfo->returnCigarStringJumpCodeLen(jumpCodeSize-2)
				+ alignInfo->returnCigarStringJumpCodeLen(jumpCodeSize-1)
				- lastMatchLength - tmpInsertionLength; 

		if(spliceJunctionDistance > MAX_DELETION_LENGTH)
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			cigarStringJumpCode.push_back(tmpSjJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}	
		else if (spliceJunctionDistance > 0)
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpDelJumpCode(spliceJunctionDistance, "D");
			cigarStringJumpCode.push_back(tmpDelJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);			
		}
		else if (spliceJunctionDistance == 0)
		{
			Jump_Code tmpJumpCode(penultimateMatchLength + lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmpJumpCode);
		}
		else if (spliceJunctionDistance < 0)
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpInsJumpCode(0 - spliceJunctionDistance, "I");
			cigarStringJumpCode.push_back(tmpInsJumpCode);

			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
			//tmpInsertionLength = 0 - spliceJunctionDistance;
		}
		else
		{

		}

		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		int tmpMismatchNum;
		if(spliceJunctionDistance < 0)
		{
			int insertionStartPosInRead = readLength - lastMatchLength + spliceJunctionDistance + 1;
			int insertionEndPosInRead = readLength - lastMatchLength;
			tmpMismatchNum = this->generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos_new(
				newMismatchPosVecToAdd, newMismatchCharVecToAdd,
				alignInfo,
				insertionStartPosInRead, insertionEndPosInRead);
		}
		else
		{
			tmpMismatchNum = alignInfo->returnMismatchPosVecSize();
				+ newMismatchPosVecToAdd.size();
			//if(STORE_MISMATCH_POS)
			//{	
				this->generateNewMismatchPosVec_new(newMismatchPosVecToAdd, alignInfo);
				//if(STORE_MISMATCH_CHA)
				//{
					this->generateNewMismatchCharVec_new(newMismatchCharVecToAdd, alignInfo);
				//}
			//}	
		}	
		mismatchNum = tmpMismatchNum;
		SJstrand = this->getStrandFromSJ();
	}

	void newAlignInfoAfterAddTailInfo2SoftClippingAlignment_new_backSpliceIncluded(
		int lastMatchLength, int spliceJunctionDistance, Index_Info* indexInfo, int newMismatchToAddInTail,
		vector<int>& newMismatchPosVecToAdd, vector<char>& newMismatchCharVecToAdd, int readLength,
		Alignment_Info* alignInfo)
	{
		//cout << "start to get new alignInfo" << endl;
		alignDirection = alignInfo->returnAlignDirection();
		alignChromName = alignInfo->returnAlignChromName();
		alignChromPos = alignInfo->returnAlignChromPos();

		//cout << "start to creat new jumpCode " << endl;
		int tmpJumpCodeVecSize 
			= alignInfo->returnCigarStringJumpCodeSize();		
		for(int tmp = 0; tmp < tmpJumpCodeVecSize-2; tmp++)
		{
			cigarStringJumpCode.push_back(
				alignInfo->returnCigarStringJumpCodeElement(tmp)
				);
		}

		int tmpInsertionLength = 0;
		if((spliceJunctionDistance < 0)&&(0 - spliceJunctionDistance <= MAX_INSERTION_LENGTH))
		{
			tmpInsertionLength = 0 - spliceJunctionDistance;
		}

		int jumpCodeSize = tmpJumpCodeVecSize;
		int penultimateMatchLength 
			= alignInfo->returnCigarStringJumpCodeLen(jumpCodeSize-2)
				+ alignInfo->returnCigarStringJumpCodeLen(jumpCodeSize-1)
				- lastMatchLength - tmpInsertionLength; 

		if(spliceJunctionDistance > MAX_DELETION_LENGTH) // normal splice junction
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			cigarStringJumpCode.push_back(tmpSjJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}	
		else if (spliceJunctionDistance > 0) // deletion
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpDelJumpCode(spliceJunctionDistance, "D");
			cigarStringJumpCode.push_back(tmpDelJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);			
		}
		else if (spliceJunctionDistance == 0) // match
		{
			Jump_Code tmpJumpCode(penultimateMatchLength + lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmpJumpCode);
		}
		else if (0 - spliceJunctionDistance <= MAX_INSERTION_LENGTH) // insertion
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);

			Jump_Code tmpInsJumpCode(0 - spliceJunctionDistance, "I");
			cigarStringJumpCode.push_back(tmpInsJumpCode);

			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
			//tmpInsertionLength = 0 - spliceJunctionDistance;
		}
		else if(0 - spliceJunctionDistance > MAX_INSERTION_LENGTH) // back splice junction
		{
			Jump_Code tmp1stMatchJumpCode(penultimateMatchLength, "M");
			cigarStringJumpCode.push_back(tmp1stMatchJumpCode);
			Jump_Code tmpSjJumpCode(spliceJunctionDistance, "N");
			cigarStringJumpCode.push_back(tmpSjJumpCode);
			Jump_Code tmp2ndMatchJumpCode(lastMatchLength, "M");
			cigarStringJumpCode.push_back(tmp2ndMatchJumpCode);
		}
		else
		{}

		this->jumpCodeVec2spliceJunctionVec(indexInfo);

		int tmpMismatchNum;
		if((spliceJunctionDistance < 0)&&(0 - spliceJunctionDistance <= MAX_INSERTION_LENGTH))
		{
			int insertionStartPosInRead = readLength - lastMatchLength + spliceJunctionDistance + 1;
			int insertionEndPosInRead = readLength - lastMatchLength;
			tmpMismatchNum = this->generateNewMismatchPosAndCharVec_filterMismatchInInsertionPos_new(
				newMismatchPosVecToAdd, newMismatchCharVecToAdd,
				alignInfo,
				insertionStartPosInRead, insertionEndPosInRead);
		}
		else
		{
			tmpMismatchNum = alignInfo->returnMismatchPosVecSize();
				+ newMismatchPosVecToAdd.size();
			//if(STORE_MISMATCH_POS)
			//{	
				this->generateNewMismatchPosVec_new(newMismatchPosVecToAdd, alignInfo);
				//if(STORE_MISMATCH_CHA)
				//{
					this->generateNewMismatchCharVec_new(newMismatchCharVecToAdd, alignInfo);
				//}
			//}	
		}	
		mismatchNum = tmpMismatchNum;
		SJstrand = this->getStrandFromSJ();
	}

	void jumpCodeStr2jumpCodeVec(const string& jumpCodeStr)
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
				cigarStringJumpCode.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
				jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
			}
		}
	}

	void jumpCodeVec2spliceJunctionVec(Index_Info* indexInfo) 
	{
		int tmpPosInRead = 0; 
		int tmpDonerPosInChr = alignChromPos - 1;
		int tmpAcceptorPosInChr = 0;
		//cout << "start to get sjVec ..." << endl;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			//cout << "JumpCode[" << tmp << "]: " << cigarStringJumpCode[tmp].len << " " << cigarStringJumpCode[tmp].type << endl;
			if(cigarStringJumpCode[tmp].type == "S")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "M")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;
			}
			else if (cigarStringJumpCode[tmp].type == "N")
			{
				//cout << "tmpDonerPosInChr: " << tmpDonerPosInChr << endl;

				tmpAcceptorPosInChr = tmpDonerPosInChr + cigarStringJumpCode[tmp].len + 1;
				//cout << "tmpAcceptorPosInChr: " << tmpAcceptorPosInChr << endl;
				//SpliceJunction_Alignment* tmpSJ = 
				//	new SpliceJunction_Alignment(alignChromName, tmpDonerPosInChr, tmpAcceptorPosInChr, indexInfo);

				SpliceJunction_Alignment tmpSJ(alignChromName, tmpDonerPosInChr, tmpAcceptorPosInChr, indexInfo);
				//cout << "tmpPosInRead: " << tmpPosInRead << endl;
				//tmpSJ
				spliceJunctionVec.push_back(pair <int, SpliceJunction_Alignment> (tmpPosInRead+1, (tmpSJ)));
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "I")
			{
				tmpPosInRead += cigarStringJumpCode[tmp].len;
				//tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else if (cigarStringJumpCode[tmp].type == "D")
			{
				tmpDonerPosInChr += cigarStringJumpCode[tmp].len;
			}
			else
			{
				cout << "other jumpCodeType" << endl;
			}
		}
	}

	int spliceJunctionConfidenceLevel()
	{
		if(spliceJunctionVec.size() == 0)
			return 0;
		string tmpFlankString;
		int spliceJunctionConfidenceLevel = SPLICE_JUNCTION_CANONICAL_ONLY;
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			tmpFlankString = (spliceJunctionVec[tmp].second).returnFlankString();
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

	string getStrandFromSJ()
	{	
		string tmpFlankString;
		string currentStrand = "N";
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			tmpFlankString = (spliceJunctionVec[tmp].second).returnFlankString();
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

	string jumpCodeVec2Str()
	{
		string str;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}

	string otherJumpCodeVec2Str() // only for Head
	{
		string str;
		for(int tmp = 1; tmp < cigarStringJumpCode.size(); tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}

	string otherJumpCodeVec2StrForTail() // only for Tail
	{
		string str;
		for(int tmp = 0; tmp < cigarStringJumpCode.size()-1; tmp++)
		{	
			str += cigarStringJumpCode[tmp].toString();
		}
		return str;
	}

	int mappedLength()
	{
		int pos = 0;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if ((cigarStringJumpCode[tmp].type == "M")||(cigarStringJumpCode[tmp].type == "I"))
			{
				pos += cigarStringJumpCode[tmp].len;
			}
			else
			{
				continue;
				//cout << "error, unexpected jumpCodeType: " << cigarStringJumpCode[tmp].type << endl;
			}
		}
		mappedBaseNum = pos;
		return pos;
	}

	int getEndMatchedPosInChr()
	{
		int pos = 0;//alignChromPos;
		//int tmpJumpCodeLength = 0;
		for(int tmp = 0; tmp < cigarStringJumpCode.size(); tmp++)
		{
			if ((cigarStringJumpCode[tmp].type == "S")||(cigarStringJumpCode[tmp].type == "I"))
			{
				continue;
			}
			else if ((cigarStringJumpCode[tmp].type == "M")||(cigarStringJumpCode[tmp].type == "N")||(cigarStringJumpCode[tmp].type == "D"))
			{
				pos += cigarStringJumpCode[tmp].len;
			}
			else
			{
				cout << "error, unexpected jumpCodeType: " << cigarStringJumpCode[tmp].type << endl;
			}
		}

		endMatchedPosInChr = (alignChromPos + pos - 1);
		return endMatchedPosInChr;
		//return (alignChromPos + pos - 1);
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

	int getFlag_unpaired_secondaryOrNot(bool mapDirection,
		bool End1OrEnd2Bool, bool SecondaryOrNot )
	{
		bool pairEndsReadBool = true;
		bool bothEndMappedBool = false;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = true;

		bool mappedAsRcmBool, anotherEndReadMappedAsRcmBool;
		if(mapDirection)
		{
			mappedAsRcmBool = false;
			anotherEndReadMappedAsRcmBool = false;
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

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	/*
	string getSamFormatString_unpaired_secondaryOrNot(const string& readName, const string& readSeq, 
		//Alignment_Info* mateAlignInfo, 
		bool End1OrEnd2, int IH_num, int HI_num, bool secondaryOrNot)
	{
		if(!this->checkJumpCodeValid())
		{
			return readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + (readSeq) + "\t*\tIH:i:0\tHI:i:0";;
		}

		string samString;

		int FLAG;
		string RNAME;
		int POS;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_unpaired_secondaryOrNot(true, End1OrEnd2, secondaryOrNot);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_unpaired_secondaryOrNot(false, End1OrEnd2, secondaryOrNot); //16;
		}
		else
		{
			RNAME = "*";
			POS = 0;
			FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////
		int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*"; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = 0;//mateAlignInfo->alignChromPos; // 0;
		int TLEN = 0;
		//string SEQ;
		string QUAL = "*";
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
			+ 
			FLAGstr + "\t" 
			+ RNAME + "\t"
			+ POSstr + "\t" 
			+ MAPQstr + "\t" 
			+ CIGAR + "\t" 
			+ RNEXT + "\t" 
			+ PNEXTstr + "\t" 
			+ TLENstr + "\t" 
			+ readSeq + "\t" 
			+ QUAL 
			+ "\tNM:i:" + mismatchNumStr 
			+ "\tIH:i:" + IHstr
			+ "\tHI:i:" + HIstr
			+ "\tXS:A:" + strandStr 
			+ "\tXF:Z:"
			;


		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			samString = samString + (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
		}

		return samString;

	}*/
	string getSamFormatString_SE(const string& readName, const string& readSeq, 
		const string& qualitySeq, int IH_num, int HI_num, 
		bool SecondaryOrNot, bool fasta_or_fastq_bool)
	{
		if(!this->checkJumpCodeValid())
		{
			return readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + (readSeq) + "\t" + qualitySeq + "\tIH:i:0\tHI:i:0";
		}
		int MAPQ;
		if(fasta_or_fastq_bool)
		{
			MAPQ = 255;
		}
		else
		{
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);
		}

		string samString;

		int FLAG;
		string RNAME;
		int POS;

		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_SE(true, SecondaryOrNot);//0;
			//TLEN = 0 when only single-segment exists in a template;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_SE(false, SecondaryOrNot); //16;
			//TLEN = 0 when only single-segment exists in a template;
		}
		else
		{
		}
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*";
		int PNEXT = 0;
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
			//+ "\tXF:Z:"
			;	
		if(spliceJunctionVec.size() > 0)
		{	
			samString += "\tXF:Z:";
			for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
				samString = samString + (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
		}
		if(mismatchPosVec.size() > 0)
		{
			samString = samString + "\tMD:Z:";
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
				samString = samString + int_to_str(mismatchPosVec[tmp]) + ",";
			samString = samString + "\tMP:Z:";
			for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
				samString = samString + mismatchCharVec[tmp] + ",";
		}
		return samString;
	}

	string getSamFormatString_unpaired_secondaryOrNot_PEasSE(
		const string& readName, const string& readSeq, const string& qualitySeq,
		bool End1OrEnd2Bool, int IH_num, int HI_num, bool otherEndReadMappedOrNotBool,
		bool secondaryOrNot, bool fasta_or_fastq_bool, int multiMapSeg_maxLength)
	{
		if(!this->checkJumpCodeValid())
		{
			//cout << "error in checkJumpCodeValid_PEasSE, getSamFormatString_unpaired_secondary_PEasSE, from align_info.h" << endl;
			//exit(1);
			return readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + (readSeq) + "\t" + qualitySeq + "\tIH:i:0\tHI:i:0";
		}
		int MAPQ;
		if(fasta_or_fastq_bool)
			MAPQ = 255;
		else
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);

		string samString;
		int FLAG;
		string RNAME;
		int POS;
		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_unpaired_secondaryOrNot_PEasSE(true, secondaryOrNot,
				otherEndReadMappedOrNotBool, End1OrEnd2Bool);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_unpaired_secondaryOrNot_PEasSE(false, secondaryOrNot,
				otherEndReadMappedOrNotBool, End1OrEnd2Bool); //16;
		}
		else
		{}

		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*";
		int PNEXT = 0;
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
		string XMstr = int_to_str(multiMapSeg_maxLength);

		samString = readName + "\t" + FLAGstr + "\t" + RNAME + "\t" + POSstr + "\t" + MAPQstr + "\t" 
			+ CIGAR + "\t" + RNEXT + "\t" + PNEXTstr + "\t" + TLENstr + "\t" + readSeq + "\t" + qualitySeq 
			+ "\tNM:i:" + mismatchNumStr + "\tIH:i:" + IHstr + "\tHI:i:" + HIstr 
			+ "\tXM:i:" + XMstr
			+ "\tXS:A:" + strandStr;
			//+ "\tXF:Z:";	
		if(spliceJunctionVec.size() > 0)
		{	
			samString += "\tXF:Z:";
			for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
				samString = samString + (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
		}
		if(mismatchPosVec.size() > 0)
		{	
			samString = samString + "\tMD:Z:";
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
				samString = samString + int_to_str(mismatchPosVec[tmp]) + ",";
			samString = samString + "\tMP:Z:";
			for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
				samString = samString + mismatchCharVec[tmp] + ",";
		}
		return samString;

	}

	string getSamFormatString_unpaired_secondaryOrNot_PEasSE(
		const string& readName, const string& readSeq, const string& qualitySeq,
		bool End1OrEnd2Bool, int IH_num, int HI_num, bool otherEndReadMappedOrNotBool,
		bool secondaryOrNot, bool fasta_or_fastq_bool)
	{
		int multiMapSeg_maxLength = 0;
		if(!this->checkJumpCodeValid())
		{
			//cout << "error in checkJumpCodeValid_PEasSE, getSamFormatString_unpaired_secondary_PEasSE, from align_info.h" << endl;
			//exit(1);
			return readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + (readSeq) + "\t" + qualitySeq + "\tIH:i:0\tHI:i:0";
		}
		int MAPQ;
		if(fasta_or_fastq_bool)
			MAPQ = 255;
		else
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);

		string samString;
		int FLAG;
		string RNAME;
		int POS;
		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_unpaired_secondaryOrNot_PEasSE(true, secondaryOrNot,
				otherEndReadMappedOrNotBool, End1OrEnd2Bool);//0;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_unpaired_secondaryOrNot_PEasSE(false, secondaryOrNot,
				otherEndReadMappedOrNotBool, End1OrEnd2Bool); //16;
		}
		else
		{}

		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "*";
		int PNEXT = 0;
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
		string XMstr = int_to_str(multiMapSeg_maxLength);

		samString = readName + "\t" + FLAGstr + "\t" + RNAME + "\t" + POSstr + "\t" + MAPQstr + "\t" 
			+ CIGAR + "\t" + RNEXT + "\t" + PNEXTstr + "\t" + TLENstr + "\t" + readSeq + "\t" + qualitySeq 
			+ "\tNM:i:" + mismatchNumStr + "\tIH:i:" + IHstr + "\tHI:i:" + HIstr 
			+ "\tXM:i:" + XMstr
			+ "\tXS:A:" + strandStr;
			//+ "\tXF:Z:";	
		if(spliceJunctionVec.size() > 0)	
		{
			samString = samString + "\tXF:Z:";
			for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
				samString = samString + (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
		}
		if(mismatchPosVec.size() > 0)
		{
			samString = samString + "\tMD:Z:";
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
				samString = samString + int_to_str(mismatchPosVec[tmp]) + ",";
			samString = samString + "\tMP:Z:";
			for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
				samString = samString + mismatchCharVec[tmp] + ",";
		}
		return samString;

	}

	int getFlag_unpaired_secondaryOrNot_PEasSE(bool mapDirection, bool secondaryOrNot, 
		bool otherEndReadMappedOrNotBool, bool End1OrEnd2Bool)
	{
		bool pairEndReadsBool = true;
		bool bothEndMappedBool = false;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = otherEndReadMappedOrNotBool;
		bool mappedAsRcmBool = (!mapDirection);
		bool anotherEndReadMappedAsRcmBool = false;
		bool pairEnd_1_bool = End1OrEnd2Bool;
		bool pairEnd_2_bool = (!End1OrEnd2Bool);
		bool notPrimaryAlignmentBool = secondaryOrNot;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;	
		int flag_SE = this->getFlag(pairEndReadsBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, notPrimaryAlignmentBool,
			notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool);
		return flag_SE;
	}

	int getFlag_SE(bool mapDirection, bool secondaryOrNot)
	{
		bool pairEndReadsBool = false;
		bool bothEndMappedBool = false;
		bool readUnmappedBool = false;
		bool anotherEndReadUnmappedBool = false;
		bool mappedAsRcmBool = (!mapDirection);
		bool anotherEndReadMappedAsRcmBool = false;
		bool pairEnd_1_bool = false;
		bool pairEnd_2_bool = false;
		bool notPrimaryAlignmentBool = secondaryOrNot;
		bool notPassQualControlBool = false;
		bool PCRorOpticalDuplicateBool = false;
		bool suppleAlignmentBool = false;	
		int flag_SE = this->getFlag(pairEndReadsBool, bothEndMappedBool, readUnmappedBool, anotherEndReadUnmappedBool,
			mappedAsRcmBool, anotherEndReadMappedAsRcmBool, pairEnd_1_bool, pairEnd_2_bool, notPrimaryAlignmentBool,
			notPassQualControlBool, PCRorOpticalDuplicateBool, suppleAlignmentBool);
		return flag_SE;
	}

	string getSamFormatString_paired_secondaryOrNot(const string& readName, const string& readSeq, 
		const string& qualitySeq,
		Alignment_Info* mateAlignInfo, bool End1OrEnd2, int IH_num, int HI_num, 
		bool SecondaryOrNot, int templateLength, bool fasta_or_fastq_bool)
	{
		int multiMapSeg_maxLength = 0;
		if(!this->checkJumpCodeValid())
		{
			return readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + (readSeq) + "\t" + qualitySeq + "\tIH:i:0\tHI:i:0";;
		}

		int MAPQ;
		if(fasta_or_fastq_bool)
		{
			MAPQ = 255;
		}
		else
		{
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);
		}

		string samString;

		int FLAG;
		string RNAME;
		int POS;

		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_paired_secondaryOrNot(true, End1OrEnd2, SecondaryOrNot);//0;
			TLEN = templateLength;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_paired_secondaryOrNot(false, End1OrEnd2, SecondaryOrNot); //16;
			TLEN = 0 - templateLength;
		}
		else
		{
			//RNAME = "*";
			//POS = 0;
			//FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		//int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = mateAlignInfo->alignChromPos; // 0;

		//string QUAL = "*";
		string strandStr = SJstrand;

		/*cout << readName << "\t" << FLAG << "\t" << alignChromName << "\t" << alignChromPos << "\t255\t" << CIGAR << "\t=\t" << PNEXT 
			<< "\t" << TLEN << "\t" << readSeq << "\t*" << "\tNM:i:"
			<< mismatchNum << "\tIH:i:"  << IH_num << "\tHI:i:" << HI_num << "\tXS:A:" << SJstrand << "\tXF:Z:";
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			cout << (spliceJunctionVec[tmp].second).flankString + ","; 
		}
		cout << endl;*/		 
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
		string XMstr = int_to_str(multiMapSeg_maxLength);

		samString = readName + "\t" 
			+ 
				FLAGstr + "\t" 
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
			+ "\tXM:i:" + XMstr
			+ "\tXS:A:" + strandStr
			//+ "\tXF:Z:"
			;
		if(spliceJunctionVec.size() > 0)
		{	
			samString += "\tXF:Z:";
			for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
			{
				samString = samString + (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
			}
		}
		if(mismatchPosVec.size() > 0)
		{	
			samString = samString + "\tMD:Z:";
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
				samString = samString + int_to_str(mismatchPosVec[tmp]) + ",";

			samString = samString + "\tMP:Z:";
			for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
				samString = samString + mismatchCharVec[tmp] + ",";
		}
		return samString;
	}

	string getSamFormatString_paired_secondaryOrNot(const string& readName, const string& readSeq, 
		const string& qualitySeq,
		Alignment_Info* mateAlignInfo, bool End1OrEnd2, int IH_num, int HI_num, 
		bool SecondaryOrNot, int templateLength, bool fasta_or_fastq_bool, int multiMapSeg_maxLength)
	{
		if(!this->checkJumpCodeValid())
		{
			return readName + "\t4\t*\t0\t0\t*\t*\t0\t0\t" + (readSeq) + "\t" + qualitySeq + "\tIH:i:0\tHI:i:0";;
		}

		int MAPQ;
		if(fasta_or_fastq_bool)
		{
			MAPQ = 255;
		}
		else
		{
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);
		}

		string samString;

		int FLAG;
		string RNAME;
		int POS;

		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			RNAME = alignChromName;
			POS = alignChromPos;
			FLAG = getFlag_paired_secondaryOrNot(true, End1OrEnd2, SecondaryOrNot);//0;
			TLEN = templateLength;
		}
		else if(alignDirection == "-")
		{
			RNAME = alignChromName;
			POS = alignChromPos;			
			FLAG = getFlag_paired_secondaryOrNot(false, End1OrEnd2, SecondaryOrNot); //16;
			TLEN = 0 - templateLength;
		}
		else
		{
			//RNAME = "*";
			//POS = 0;
			//FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		//int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		int PNEXT = mateAlignInfo->alignChromPos; // 0;

		//string QUAL = "*";
		string strandStr = SJstrand;

		/*cout << readName << "\t" << FLAG << "\t" << alignChromName << "\t" << alignChromPos << "\t255\t" << CIGAR << "\t=\t" << PNEXT 
			<< "\t" << TLEN << "\t" << readSeq << "\t*" << "\tNM:i:"
			<< mismatchNum << "\tIH:i:"  << IH_num << "\tHI:i:" << HI_num << "\tXS:A:" << SJstrand << "\tXF:Z:";
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			cout << (spliceJunctionVec[tmp].second).flankString + ","; 
		}
		cout << endl;*/		 
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
		string XMstr = int_to_str(multiMapSeg_maxLength);

		samString = readName + "\t" 
			+ 
				FLAGstr + "\t" 
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
			+ "\tXM:i:" + XMstr
			+ "\tXS:A:" + strandStr;
		if(spliceJunctionVec.size() > 0)
		{	
			samString += "\tXF:Z:";
			for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
				samString = samString + (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
		}
		if(mismatchPosVec.size() > 0)
		{
			samString = samString + "\tMD:Z:";
			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
				samString = samString + int_to_str(mismatchPosVec[tmp]) + ",";
			samString = samString + "\tMP:Z:";
			for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
				samString = samString + mismatchCharVec[tmp] + ",";
		}
		return samString;
	}

	/*
	void outputSamFormatString_paired_secondaryOrNot(const string& readName, const string& readSeq,
		const string& qualitySeq, 
		Alignment_Info* mateAlignInfo, bool End1OrEnd2, int IH_num, int HI_num, bool SecondaryOrNot, 
		int templateLength, ofstream* output, bool fasta_or_fastq_bool)
	{
		int MAPQ;
		if(fasta_or_fastq_bool)
		{
			MAPQ = 255;
		}
		else
		{
			MAPQ = this->getMapppingQuality(qualitySeq, "phred33", readName);
		}

		int FLAG;
		//string RNAME;
		//int POS;

		int TLEN = 0;
		///////////////////// get FLAG, RNAME and POS ///////////////////////////
		if(alignDirection == "+")
		{
			//RNAME = alignChromName;
			//POS = alignChromPos;
			FLAG = getFlag_paired_secondaryOrNot(true, End1OrEnd2, SecondaryOrNot);//0;
			TLEN = templateLength;
		}
		else if(alignDirection == "-")
		{
			//RNAME = alignChromName;
			//POS = alignChromPos;			
			FLAG = getFlag_paired_secondaryOrNot(false, End1OrEnd2, SecondaryOrNot); //16;
			TLEN = 0 - templateLength;
		}
		else
		{
			//RNAME = "*";
			//POS = 0;
			//FLAG = 4;
		}
		////////////////////////////////////////////////////////////////////////

		//int MAPQ = 255;
		string CIGAR = this->jumpCodeVec2Str();
		//string RNEXT = "="; //mateAlignInfo->alignChromName; //"*";
		//int PNEXT = mateAlignInfo->alignChromPos; // 0;

		//string QUAL = "*";
		string strandStr = SJstrand;


		((*output)) << readName << "\t" << FLAG << "\t" << alignChromName << "\t" << alignChromPos << "\t" << MAPQ <<"\t" << CIGAR << "\t=\t" << mateAlignInfo->alignChromPos
			<< "\t" << TLEN << "\t" << readSeq << "\t" + qualitySeq + "\tNM:i:"
			<< mismatchNum << "\tIH:i:"  << IH_num << "\tHI:i:" << HI_num << "\tXS:A:" << SJstrand << "\tXF:Z:";
		//cout << "here 2" << endl;
		for(int tmp = 0; tmp < spliceJunctionVec.size(); tmp++)
		{
			((*output)) << (spliceJunctionVec[tmp].second).returnFlankString() + ","; 
		}

		//if(STORE_MISMATCH_POS)
		//{
			(*output) << "\tMD:Z:";

			for(int tmp = 0; tmp < mismatchPosVec.size(); tmp++)
			{
				(*output) << int_to_str(mismatchPosVec[tmp]) << ",";
			}

			//if(STORE_MISMATCH_CHA)
			//{
				(*output) << "\tMP:Z:";
				for(int tmp = 0; tmp < mismatchCharVec.size(); tmp++)
				{
					(*output) << (mismatchCharVec[tmp]) << ",";
				}
			//}
		//}

		((*output)) << endl;		 
	}*/
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	bool checkOverlapPairAlignment_new(Alignment_Info* secondAlignmentInfo) // FIX ME: some bugs contained  // parameter is another alignment 
	{
		//cout << "start to check overlap SJ" << endl;
		bool SJinOverlapAreaNoConflict = true;

		return SJinOverlapAreaNoConflict;


		int overlapStartPos_1 = this->alignChromPos;
		int overlapStartPos_2 = secondAlignmentInfo->alignChromPos;
		int overlapEndPos_1 = this->endMatchedPosInChr;
		int overlapEndPos_2 = secondAlignmentInfo->endMatchedPosInChr;

		int overlapStartPos;
		int overlapEndPos;
		
		if(overlapStartPos_1 > overlapStartPos_2)
			overlapStartPos = overlapStartPos_1;
		else
			overlapStartPos = overlapStartPos_2;

		if(overlapEndPos_1 < overlapEndPos_2)
			overlapEndPos = overlapEndPos_1;
		else
			overlapEndPos = overlapEndPos_2;

		//cout << "overlapStartPos: " << overlapStartPos << endl;
		//cout << "overlapEndPos: " << overlapEndPos << endl;

		int tmpDonerPosInChr, tmpAcceptorPosInChr;
		vector< pair<int, int> > norAlignSJvecInOverlapArea;
		vector< pair<int, int> > rcmAlignSJvecInOverlapArea;
		for(int tmp = 0; tmp < (this->spliceJunctionVec).size(); tmp++)
		{
			if((((this->spliceJunctionVec)[tmp].second).returnSJdonerEnd() >= overlapStartPos)
				&&(((this->spliceJunctionVec)[tmp].second).returnSJacceptorStart() <= overlapEndPos))
			{
				tmpDonerPosInChr = ((this->spliceJunctionVec)[tmp].second).returnSJdonerEnd();
				tmpAcceptorPosInChr = ((this->spliceJunctionVec)[tmp].second).returnSJacceptorStart();
				norAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 

		for(int tmp = 0; tmp < (secondAlignmentInfo->spliceJunctionVec).size(); tmp++)
		{
			if((((secondAlignmentInfo->spliceJunctionVec)[tmp].second).returnSJacceptorStart() <= overlapEndPos)
				&&(((secondAlignmentInfo->spliceJunctionVec)[tmp].second).returnSJdonerEnd() >= overlapStartPos))
			{
				tmpDonerPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).returnSJdonerEnd();
				tmpAcceptorPosInChr = ((secondAlignmentInfo->spliceJunctionVec)[tmp].second).returnSJacceptorStart();
				rcmAlignSJvecInOverlapArea.push_back(pair<int,int> (tmpDonerPosInChr, tmpAcceptorPosInChr));
			}
			else
			{}
		} 
		//cout << "norAlignSJvecInOverlapArea.size(): " << norAlignSJvecInOverlapArea.size() << endl;
		//cout << "rcmAlignSJvecInOverlapArea.size(): " << rcmAlignSJvecInOverlapArea.size() << endl;
		if(norAlignSJvecInOverlapArea.size() == rcmAlignSJvecInOverlapArea.size())
		{
			for(int tmp = 0; tmp < norAlignSJvecInOverlapArea.size(); tmp++)
			{
				if((norAlignSJvecInOverlapArea[tmp].first == rcmAlignSJvecInOverlapArea[tmp].first)
					&&(norAlignSJvecInOverlapArea[tmp].second == rcmAlignSJvecInOverlapArea[tmp].second))
				{
					SJinOverlapAreaNoConflict = true;
				}
				else
				{
					return false;
				}
			}
		}
		else
		{
			return false;
		}

		return SJinOverlapAreaNoConflict;
	}
};

#endif


