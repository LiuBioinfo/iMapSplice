// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SEG_INFO_H
#define SEG_INFO_H

#include <stdlib.h>
#include <stdio.h>
#include "enhanced_suffix_array_info.h"
#include "local_seg_info.h"
#include "missingLongSeg_info.h"

using namespace std;

class Seg_Info
{
//public:
private:
	int repeatRegion_index; // if -1 or < 0, is not mapped to repeatRegion 

	unsigned int segmentNum;
	unsigned int norSegmentLength[SEGMENTNUM];
	unsigned int norSegmentLocInRead[SEGMENTNUM];
	unsigned int norSegmentAlignNum[SEGMENTNUM];
	unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC];

	int longSegMinLength;

	#ifdef PERSONALIZED_CHR_SEQ
	bool norSegmentUpdatedWithSNPorNotBool[SEGMENTNUM];
	#endif
	//vector<Seg_Sub_Info> subSegInfoVec;
public:
	Seg_Info()
	{
		longSegMinLength = LONG_SEG_LENGTH_THRESHOLD_PHASE1;
		segmentNum = 0;
		repeatRegion_index = -1;
	}

	#ifdef PERSONALIZED_CHR_SEQ
	void generateSNPlocatedSegInfoVec(vector<int>& snpLocatedSeg_startLocVec,
		vector<int>& snpLocatedSeg_lengthVec, vector< vector <unsigned int> >& snpLocatedSeg_mapPosVecVecInWholeGenome,
		int SNPlocInSyntheticSNPseq, Index_Info* snpSeqIndexInfo, Index_Info* indexInfo)
	{
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			set<unsigned int> tmpSet;
			int tmpSeg_startLocInRead = norSegmentLocInRead[tmpSeg];
			int tmpSeg_length = norSegmentLength[tmpSeg];
			int tmpSeg_alignNum = norSegmentAlignNum[tmpSeg];
			if((tmpSeg_alignNum <= 0)||(tmpSeg_alignNum > CANDALILOC))
				continue;
			for(int tmpLocIndex = 0; tmpLocIndex < tmpSeg_alignNum; tmpLocIndex++)
			{
				unsigned int snpLocatedSeg_tmpMapPos = this->returnSegmentMapPos(tmpSeg, tmpLocIndex);
				int snpLocatedSeg_tmpMapPos_snpSeqNameInt = snpSeqIndexInfo->getChr(snpLocatedSeg_tmpMapPos);
				int snpLocatedSeg_tmpMapPos_snpSeqLoc = snpSeqIndexInfo->getChrLocation_withChrNameInt(
					snpLocatedSeg_tmpMapPos, snpLocatedSeg_tmpMapPos_snpSeqNameInt);
				if(!((SNPlocInSyntheticSNPseq >= snpLocatedSeg_tmpMapPos_snpSeqLoc)
					&&(SNPlocInSyntheticSNPseq < snpLocatedSeg_tmpMapPos_snpSeqLoc + tmpSeg_length - 1)))
					continue;
				string snpLocatedSeg_tmpMapPos_snpSeqName = snpSeqIndexInfo->returnChrNameStr(
					snpLocatedSeg_tmpMapPos_snpSeqNameInt);			
				vector<string> tmpFieldVec;
				int startLoc = 0;
				for(int tmp = 0; tmp < 4; tmp++)
				{
					int commaLoc = snpLocatedSeg_tmpMapPos_snpSeqName.find(":", startLoc);
					string tmpFieldStr = snpLocatedSeg_tmpMapPos_snpSeqName.substr(startLoc, commaLoc - startLoc);
					tmpFieldVec.push_back(tmpFieldStr);
					startLoc = commaLoc + 1;
				}
				string chrNameInGenome = tmpFieldVec[0];
				int chrNameInGenome_int = indexInfo->convertStringToInt(chrNameInGenome);
				int chrStartPos_tmpSnpSeq = atoi(tmpFieldVec[1].c_str());
				int mapPos_inTmpSnpSeq = snpLocatedSeg_tmpMapPos_snpSeqLoc; 
				int mapPos_inChr = mapPos_inTmpSnpSeq - 1 + chrStartPos_tmpSnpSeq;
				unsigned int mapPos_inWholeGenome = indexInfo->getWholeGenomeLocation(chrNameInGenome_int, mapPos_inChr);
				tmpSet.insert(mapPos_inWholeGenome);
			}
			if(tmpSet.size() > 0)
			{
				vector<unsigned int> tmpVec;
				snpLocatedSeg_startLocVec.push_back(tmpSeg_startLocInRead);
				snpLocatedSeg_lengthVec.push_back(tmpSeg_length);
				for(set<unsigned int>::iterator tmpIter = tmpSet.begin(); tmpIter != tmpSet.end(); tmpIter ++)
				{
					unsigned int tmpVal = (*tmpIter);
					tmpVec.push_back(tmpVal);
				}
				snpLocatedSeg_mapPosVecVecInWholeGenome.push_back(tmpVec);
			}
		}
	}

	int	returnSNPlocInSyntheticSNPseqFromSNPseqNameStr(string& snpLocatedSeg_tmpMapPos_snpSeqNameStr)
	{
		int comma_loc_1 = snpLocatedSeg_tmpMapPos_snpSeqNameStr.find(":");
		int comma_loc_2 = snpLocatedSeg_tmpMapPos_snpSeqNameStr.find(":", comma_loc_1 + 1);
		int comma_loc_3 = snpLocatedSeg_tmpMapPos_snpSeqNameStr.find(":", comma_loc_2 + 1);
		int comma_loc_4 = snpLocatedSeg_tmpMapPos_snpSeqNameStr.find(":", comma_loc_3 + 1);
		int comma_loc_5 = snpLocatedSeg_tmpMapPos_snpSeqNameStr.find(":", comma_loc_4 + 1);
		string snpLocStr = snpLocatedSeg_tmpMapPos_snpSeqNameStr.substr(comma_loc_3 + 1, comma_loc_4 - comma_loc_3 - 1);
		return atoi(snpLocStr.c_str());
	}

	void generateSNPlocatedSegInfoVec_varySNPmer(vector<int>& snpLocatedSeg_startLocVec,
		vector<int>& snpLocatedSeg_lengthVec, vector< vector <unsigned int> >& snpLocatedSeg_mapPosVecVecInWholeGenome,
		//int SNPlocInSyntheticSNPseq, 
		Index_Info* snpSeqIndexInfo, Index_Info* indexInfo)
	{
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			set<unsigned int> tmpSet;
			int tmpSeg_startLocInRead = norSegmentLocInRead[tmpSeg];
			int tmpSeg_length = norSegmentLength[tmpSeg];
			int tmpSeg_alignNum = norSegmentAlignNum[tmpSeg];
			if((tmpSeg_alignNum <= 0)||(tmpSeg_alignNum > CANDALILOC))
				continue;
			for(int tmpLocIndex = 0; tmpLocIndex < tmpSeg_alignNum; tmpLocIndex++)
			{
				unsigned int snpLocatedSeg_tmpMapPos = this->returnSegmentMapPos(tmpSeg, tmpLocIndex);
				int snpLocatedSeg_tmpMapPos_snpSeqNameInt = snpSeqIndexInfo->getChr(snpLocatedSeg_tmpMapPos);
				string snpLocatedSeg_tmpMapPos_snpSeqNameStr = snpSeqIndexInfo->returnChrNameStr(snpLocatedSeg_tmpMapPos_snpSeqNameInt);
				int snpLocatedSeg_tmpMapPos_snpSeq_SNPlocInSyntheticSNPseq
					= this->returnSNPlocInSyntheticSNPseqFromSNPseqNameStr(snpLocatedSeg_tmpMapPos_snpSeqNameStr);
				int snpLocatedSeg_tmpMapPos_snpSeqLoc = snpSeqIndexInfo->getChrLocation_withChrNameInt(
					snpLocatedSeg_tmpMapPos, snpLocatedSeg_tmpMapPos_snpSeqNameInt);
				// FIX ME ???
				if(!((snpLocatedSeg_tmpMapPos_snpSeq_SNPlocInSyntheticSNPseq >= snpLocatedSeg_tmpMapPos_snpSeqLoc)
					&&(snpLocatedSeg_tmpMapPos_snpSeq_SNPlocInSyntheticSNPseq < snpLocatedSeg_tmpMapPos_snpSeqLoc + tmpSeg_length - 1)))
					continue;
				string snpLocatedSeg_tmpMapPos_snpSeqName = snpSeqIndexInfo->returnChrNameStr(snpLocatedSeg_tmpMapPos_snpSeqNameInt);			
				vector<string> tmpFieldVec;
				int startLoc = 0;
				for(int tmp = 0; tmp < 4; tmp++)
				{
					int commaLoc = snpLocatedSeg_tmpMapPos_snpSeqName.find(":", startLoc);
					string tmpFieldStr = snpLocatedSeg_tmpMapPos_snpSeqName.substr(startLoc, commaLoc - startLoc);
					tmpFieldVec.push_back(tmpFieldStr);
					startLoc = commaLoc + 1;
				}
				string chrNameInGenome = tmpFieldVec[0];
				int chrNameInGenome_int = indexInfo->convertStringToInt(chrNameInGenome);
				int chrStartPos_tmpSnpSeq = atoi(tmpFieldVec[1].c_str());
				int mapPos_inTmpSnpSeq = snpLocatedSeg_tmpMapPos_snpSeqLoc; 
				int mapPos_inChr = mapPos_inTmpSnpSeq - 1 + chrStartPos_tmpSnpSeq;
				unsigned int mapPos_inWholeGenome = indexInfo->getWholeGenomeLocation(chrNameInGenome_int, mapPos_inChr);
				tmpSet.insert(mapPos_inWholeGenome);
			}
			if(tmpSet.size() > 0)
			{
				vector<unsigned int> tmpVec;
				snpLocatedSeg_startLocVec.push_back(tmpSeg_startLocInRead);
				snpLocatedSeg_lengthVec.push_back(tmpSeg_length);
				for(set<unsigned int>::iterator tmpIter = tmpSet.begin(); tmpIter != tmpSet.end(); tmpIter ++)
				{
					unsigned int tmpVal = (*tmpIter);
					tmpVec.push_back(tmpVal);
				}
				snpLocatedSeg_mapPosVecVecInWholeGenome.push_back(tmpVec);
			}
		}
	}	

	void update_segGroupWithSNPlocatedSeg(int tmp_snpLocatedSeg_startLocInRead, 
		vector<unsigned int>& snpLocatedSeg_mapPosVecInWholeGenome, int toUpdateSegIndexInThisSegInfo)
	{
		int tmpSeg_startLocInRead = norSegmentLocInRead[toUpdateSegIndexInThisSegInfo];
		int snpLocatedSeg_mapPosVecInWholeGenome_size = snpLocatedSeg_mapPosVecInWholeGenome.size();
		int segGroup_alignNum_ori = this->returnSegmentAlignNum(toUpdateSegIndexInThisSegInfo);
		if(segGroup_alignNum_ori > CANDALILOC)
			segGroup_alignNum_ori = 0;
		if(segGroup_alignNum_ori + snpLocatedSeg_mapPosVecInWholeGenome_size > CANDALILOC)
		{
			norSegmentAlignNum[toUpdateSegIndexInThisSegInfo] = CANDALILOC;
			for(int tmp = 0; tmp < snpLocatedSeg_mapPosVecInWholeGenome_size; tmp++)
			{
				int tmpAlignLocIndexInThisSegGroup = CANDALILOC - 1 - tmp;
				unsigned int tmpToUpdatedSegInThisSegInfo_mapPos = snpLocatedSeg_mapPosVecInWholeGenome[tmp]
					- tmp_snpLocatedSeg_startLocInRead + tmpSeg_startLocInRead;
				this->assignSegmentAlignLoc(toUpdateSegIndexInThisSegInfo, tmpAlignLocIndexInThisSegGroup,
					tmpToUpdatedSegInThisSegInfo_mapPos);
			}
		}
		else
		{
			norSegmentAlignNum[toUpdateSegIndexInThisSegInfo] = segGroup_alignNum_ori + snpLocatedSeg_mapPosVecInWholeGenome_size;
			for(int tmp = 0; tmp < snpLocatedSeg_mapPosVecInWholeGenome_size; tmp++)
			{
				int tmpAlignLocIndexInThisSegGroup = segGroup_alignNum_ori + tmp;
				unsigned int tmpToUpdatedSegInThisSegInfo_mapPos = snpLocatedSeg_mapPosVecInWholeGenome[tmp]
					- tmp_snpLocatedSeg_startLocInRead + tmpSeg_startLocInRead;
				this->assignSegmentAlignLoc(toUpdateSegIndexInThisSegInfo, tmpAlignLocIndexInThisSegGroup,
					tmpToUpdatedSegInThisSegInfo_mapPos);
			}
		}
		norSegmentUpdatedWithSNPorNotBool[toUpdateSegIndexInThisSegInfo] = true;
	}

	void update_includeSNPseqMapSegInfo(Seg_Info& snpSeqMapSegInfo, int SNPlocInSyntheticSNPseq, 
		Index_Info* snpSeqIndexInfo, Index_Info* indexInfo)
	{
		for(int tmpSeg = 0; tmpSeg < SEGMENTNUM; tmpSeg ++)
			norSegmentUpdatedWithSNPorNotBool[tmpSeg] = false;
		vector<int> snpLocatedSeg_startLocVec;
		vector<int> snpLocatedSeg_lengthVec; 
		vector< vector <unsigned int> > snpLocatedSeg_mapPosVecVecInWholeGenome;
		snpSeqMapSegInfo.generateSNPlocatedSegInfoVec(snpLocatedSeg_startLocVec, snpLocatedSeg_lengthVec,
			snpLocatedSeg_mapPosVecVecInWholeGenome, SNPlocInSyntheticSNPseq, snpSeqIndexInfo, indexInfo);
		if(snpLocatedSeg_startLocVec.size() <= 0)
			return;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			int tmpSeg_startLocInRead = norSegmentLocInRead[tmpSeg];
			int tmpSeg_length = norSegmentLength[tmpSeg];
			int tmpSeg_endLocInRead = tmpSeg_startLocInRead + tmpSeg_length - 1;
			for(int tmpSnpLocatedSegIndex = 0; tmpSnpLocatedSegIndex < snpLocatedSeg_startLocVec.size(); tmpSnpLocatedSegIndex ++)
			{	
				int tmp_snpLocatedSeg_startLocInRead = snpLocatedSeg_startLocVec[tmpSnpLocatedSegIndex];
				int tmp_snpLocatedSeg_length = snpLocatedSeg_lengthVec[tmpSnpLocatedSegIndex];
				int tmp_snpLocatedSeg_endLocInRead = tmp_snpLocatedSeg_startLocInRead + tmp_snpLocatedSeg_length - 1;
				if((tmpSeg_startLocInRead >= tmp_snpLocatedSeg_startLocInRead)&&(tmpSeg_endLocInRead <= tmp_snpLocatedSeg_endLocInRead))
				{
					this->update_segGroupWithSNPlocatedSeg(tmp_snpLocatedSeg_startLocInRead, 
						snpLocatedSeg_mapPosVecVecInWholeGenome[tmpSnpLocatedSegIndex], tmpSeg);
					break;
				}
			}
		}
	}

	void update_includeSNPseqMapSegInfo_varySNPmer(Seg_Info& snpSeqMapSegInfo, //int SNPlocInSyntheticSNPseq, 
		Index_Info* snpSeqIndexInfo, Index_Info* indexInfo)
	{
		for(int tmpSeg = 0; tmpSeg < SEGMENTNUM; tmpSeg ++)
			norSegmentUpdatedWithSNPorNotBool[tmpSeg] = false;
		vector<int> snpLocatedSeg_startLocVec;
		vector<int> snpLocatedSeg_lengthVec; 
		vector< vector <unsigned int> > snpLocatedSeg_mapPosVecVecInWholeGenome;
		snpSeqMapSegInfo.generateSNPlocatedSegInfoVec_varySNPmer(snpLocatedSeg_startLocVec, snpLocatedSeg_lengthVec,
			snpLocatedSeg_mapPosVecVecInWholeGenome, //SNPlocInSyntheticSNPseq, 
			snpSeqIndexInfo, indexInfo);
		if(snpLocatedSeg_startLocVec.size() <= 0)
			return;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			int tmpSeg_startLocInRead = norSegmentLocInRead[tmpSeg];
			int tmpSeg_length = norSegmentLength[tmpSeg];
			int tmpSeg_endLocInRead = tmpSeg_startLocInRead + tmpSeg_length - 1;
			for(int tmpSnpLocatedSegIndex = 0; tmpSnpLocatedSegIndex < snpLocatedSeg_startLocVec.size(); tmpSnpLocatedSegIndex ++)
			{	
				int tmp_snpLocatedSeg_startLocInRead = snpLocatedSeg_startLocVec[tmpSnpLocatedSegIndex];
				int tmp_snpLocatedSeg_length = snpLocatedSeg_lengthVec[tmpSnpLocatedSegIndex];
				int tmp_snpLocatedSeg_endLocInRead = tmp_snpLocatedSeg_startLocInRead + tmp_snpLocatedSeg_length - 1;
				if((tmpSeg_startLocInRead >= tmp_snpLocatedSeg_startLocInRead)&&(tmpSeg_endLocInRead <= tmp_snpLocatedSeg_endLocInRead))
				{
					this->update_segGroupWithSNPlocatedSeg(tmp_snpLocatedSeg_startLocInRead, 
						snpLocatedSeg_mapPosVecVecInWholeGenome[tmpSnpLocatedSegIndex], tmpSeg);
					break;
				}
			}
		}
	}	

	bool targetMap_generateSNPlocatedSegMapPosVec(string& tmpSegSeq, unsigned int* sa_snp, BYTE* lcpCompress_snp,
		unsigned int* childTab_snp, char* chrom_snp, BYTE* verifyChild_snp, Index_Info* indexInfo_snp, 
		Index_Info* indexInfo, int& targetMapAlignNum, vector<unsigned int>& targetMapAlignLocVec, int SNPlocInSyntheticSNPseq)
	{
		//cout << "start to do targetMap_generateSNPlocatedSegMapPosVec " << endl;
		char* tmpSegChar = const_cast<char*>(tmpSegSeq.c_str());
		int targetMapAlignNum_ori;
		vector<unsigned int> targetMapAlignLocVec_ori;
		bool targetMapMain_bool = targetMapWithoutPreIndexArray(tmpSegChar, sa_snp, lcpCompress_snp, 
			childTab_snp, chrom_snp, verifyChild_snp, tmpSegSeq.length(), indexInfo_snp,
			targetMapAlignNum_ori, targetMapAlignLocVec_ori);
		//cout << "targetMapMain_bool: " << targetMapMain_bool << endl;
		if(!targetMapMain_bool)
			return false; 
		if((targetMapAlignNum_ori <= 0)||(targetMapAlignNum_ori > CANDALILOC))
			return false;
		//cout << "targetMapAlignNum_ori: " << targetMapAlignNum_ori << endl;
		set<unsigned int> tmpSet;
		for(int tmpLocIndex = 0; tmpLocIndex < targetMapAlignNum_ori; tmpLocIndex++)
		{
			unsigned int snpLocatedSeg_tmpMapPos = targetMapAlignLocVec_ori[tmpLocIndex];
			int snpLocatedSeg_tmpMapPos_snpSeqNameInt = indexInfo_snp->getChr(snpLocatedSeg_tmpMapPos);
			int snpLocatedSeg_tmpMapPos_snpSeqLoc = indexInfo_snp->getChrLocation_withChrNameInt(
					snpLocatedSeg_tmpMapPos, snpLocatedSeg_tmpMapPos_snpSeqNameInt);
			if(!((SNPlocInSyntheticSNPseq >= snpLocatedSeg_tmpMapPos_snpSeqLoc)
				&&(SNPlocInSyntheticSNPseq < snpLocatedSeg_tmpMapPos_snpSeqLoc + tmpSegSeq.length() - 1)))
				continue;
			string snpLocatedSeg_tmpMapPos_snpSeqName = indexInfo_snp->returnChrNameStr(
				snpLocatedSeg_tmpMapPos_snpSeqNameInt);				
			vector<string> tmpFieldVec;
			int startLoc = 0;
			for(int tmp = 0; tmp < 4; tmp++)
			{
				int commaLoc = snpLocatedSeg_tmpMapPos_snpSeqName.find(":", startLoc);
				string tmpFieldStr = snpLocatedSeg_tmpMapPos_snpSeqName.substr(startLoc, commaLoc - startLoc);
				tmpFieldVec.push_back(tmpFieldStr);
				startLoc = commaLoc + 1;
			}
			string chrNameInGenome = tmpFieldVec[0];
			int chrNameInGenome_int = indexInfo->convertStringToInt(chrNameInGenome);
			int chrStartPos_tmpSnpSeq = atoi(tmpFieldVec[1].c_str());
			int mapPos_inTmpSnpSeq = snpLocatedSeg_tmpMapPos_snpSeqLoc; 
			int mapPos_inChr = mapPos_inTmpSnpSeq - 1 + chrStartPos_tmpSnpSeq;
			unsigned int mapPos_inWholeGenome = indexInfo->getWholeGenomeLocation(chrNameInGenome_int, mapPos_inChr);
			tmpSet.insert(mapPos_inWholeGenome);			
		}
		int tmpSet_size = tmpSet.size();
		if(tmpSet_size > 0)
		{
			targetMapAlignNum = tmpSet_size;
			for(set<unsigned int>::iterator tmpIter = tmpSet.begin(); tmpIter != tmpSet.end(); tmpIter ++)
			{
				unsigned int tmpVal = (*tmpIter);
				targetMapAlignLocVec.push_back(tmpVal);
			}
			return true;	
		}
		else
			return false;
	}

	bool targetMap_generateSNPlocatedSegMapPosVec_varySNPmer(string& tmpSegSeq, unsigned int* sa_snp, BYTE* lcpCompress_snp,
		unsigned int* childTab_snp, char* chrom_snp, BYTE* verifyChild_snp, Index_Info* indexInfo_snp, 
		Index_Info* indexInfo, int& targetMapAlignNum, vector<unsigned int>& targetMapAlignLocVec)//, int SNPlocInSyntheticSNPseq)
	{
		//cout << "start to do targetMap_generateSNPlocatedSegMapPosVec " << endl;
		char* tmpSegChar = const_cast<char*>(tmpSegSeq.c_str());
		int targetMapAlignNum_ori;
		vector<unsigned int> targetMapAlignLocVec_ori;
		bool targetMapMain_bool = targetMapWithoutPreIndexArray(tmpSegChar, sa_snp, lcpCompress_snp, 
			childTab_snp, chrom_snp, verifyChild_snp, tmpSegSeq.length(), indexInfo_snp,
			targetMapAlignNum_ori, targetMapAlignLocVec_ori);
		//cout << "targetMapMain_bool: " << targetMapMain_bool << endl;
		if(!targetMapMain_bool)
			return false; 
		if((targetMapAlignNum_ori <= 0)||(targetMapAlignNum_ori > CANDALILOC))
			return false;
		//cout << "targetMapAlignNum_ori: " << targetMapAlignNum_ori << endl;
		set<unsigned int> tmpSet;
		for(int tmpLocIndex = 0; tmpLocIndex < targetMapAlignNum_ori; tmpLocIndex++)
		{
			unsigned int snpLocatedSeg_tmpMapPos = targetMapAlignLocVec_ori[tmpLocIndex];
			int snpLocatedSeg_tmpMapPos_snpSeqNameInt = indexInfo_snp->getChr(snpLocatedSeg_tmpMapPos);
			string snpLocatedSeg_tmpMapPos_snpSeqName = indexInfo_snp->returnChrNameStr(snpLocatedSeg_tmpMapPos_snpSeqNameInt);
			int snpLocatedSeg_tmpMapPos_snpSeq_SNPlocInSyntheticSNPseq
				= this->returnSNPlocInSyntheticSNPseqFromSNPseqNameStr(snpLocatedSeg_tmpMapPos_snpSeqName);
			int snpLocatedSeg_tmpMapPos_snpSeqLoc = indexInfo_snp->getChrLocation_withChrNameInt(
					snpLocatedSeg_tmpMapPos, snpLocatedSeg_tmpMapPos_snpSeqNameInt);
			if(!((snpLocatedSeg_tmpMapPos_snpSeq_SNPlocInSyntheticSNPseq >= snpLocatedSeg_tmpMapPos_snpSeqLoc)
				&&(snpLocatedSeg_tmpMapPos_snpSeq_SNPlocInSyntheticSNPseq < snpLocatedSeg_tmpMapPos_snpSeqLoc + tmpSegSeq.length() - 1)))
				continue;				
			vector<string> tmpFieldVec;
			int startLoc = 0;
			for(int tmp = 0; tmp < 4; tmp++)
			{
				int commaLoc = snpLocatedSeg_tmpMapPos_snpSeqName.find(":", startLoc);
				string tmpFieldStr = snpLocatedSeg_tmpMapPos_snpSeqName.substr(startLoc, commaLoc - startLoc);
				tmpFieldVec.push_back(tmpFieldStr);
				startLoc = commaLoc + 1;
			}
			string chrNameInGenome = tmpFieldVec[0];
			int chrNameInGenome_int = indexInfo->convertStringToInt(chrNameInGenome);
			int chrStartPos_tmpSnpSeq = atoi(tmpFieldVec[1].c_str());
			int mapPos_inTmpSnpSeq = snpLocatedSeg_tmpMapPos_snpSeqLoc; 
			int mapPos_inChr = mapPos_inTmpSnpSeq - 1 + chrStartPos_tmpSnpSeq;
			unsigned int mapPos_inWholeGenome = indexInfo->getWholeGenomeLocation(chrNameInGenome_int, mapPos_inChr);
			tmpSet.insert(mapPos_inWholeGenome);			
		}
		int tmpSet_size = tmpSet.size();
		if(tmpSet_size > 0)
		{
			targetMapAlignNum = tmpSet_size;
			for(set<unsigned int>::iterator tmpIter = tmpSet.begin(); tmpIter != tmpSet.end(); tmpIter ++)
			{
				unsigned int tmpVal = (*tmpIter);
				targetMapAlignLocVec.push_back(tmpVal);
			}
			return true;	
		}
		else
			return false;
	}	

	void update_segGroupWithTargetMap2SNPseq(int tmpSegGroupIndexInThisSegInfo, string& tmpSegSeq, 
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp, 
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, Index_Info* indexInfo, int SNPlocInSyntheticSNPseq)
	{
		//char* tmpSegChar = const_cast<char*>(tmpSegSeq.c_str());
		//cout << "start to do update_segGroupWithTargetMap2SNPseq" << endl;
		int targetMapAlignNum;
		vector<unsigned int> targetMapAlignLocVec;
		bool targetMapSuccess_bool = targetMap_generateSNPlocatedSegMapPosVec(tmpSegSeq, sa_snp, lcpCompress_snp, 
			childTab_snp, chrom_snp, verifyChild_snp, indexInfo_snp, indexInfo, targetMapAlignNum, 
			targetMapAlignLocVec, SNPlocInSyntheticSNPseq);
		//cout << "targetMapSuccess_bool: " << targetMapSuccess_bool << endl; 
		if(!targetMapSuccess_bool)
			return;
		int segGroup_alignNum_ori = this->returnSegmentAlignNum(tmpSegGroupIndexInThisSegInfo);
		//cout << "segGroup_alignNum_ori: " << segGroup_alignNum_ori << endl;
		//cout << "targetMapAlignNum: " << targetMapAlignNum << endl;
		if(segGroup_alignNum_ori > CANDALILOC)
			segGroup_alignNum_ori = 0;
		if(segGroup_alignNum_ori + targetMapAlignNum > CANDALILOC)
		{
			norSegmentAlignNum[tmpSegGroupIndexInThisSegInfo] = CANDALILOC;
			// assign map pos
			for(int tmp = 0; tmp < targetMapAlignNum; tmp++) // reversely copy align poses
			{
				int tmpAlignLocIndexInThisSegGroup = CANDALILOC - 1 - tmp;
				this->assignSegmentAlignLoc(tmpSegGroupIndexInThisSegInfo, tmpAlignLocIndexInThisSegGroup,
					targetMapAlignLocVec[tmp]);
			}
		}
		else
		{
			norSegmentAlignNum[tmpSegGroupIndexInThisSegInfo] = segGroup_alignNum_ori + targetMapAlignNum;
			// assign map ps
			for(int tmp = 0; tmp < targetMapAlignNum; tmp++)
			{
				int tmpAlignLocIndexInThisSegGroup = segGroup_alignNum_ori + tmp;
				this->assignSegmentAlignLoc(tmpSegGroupIndexInThisSegInfo, tmpAlignLocIndexInThisSegGroup,
					targetMapAlignLocVec[tmp]);
			}
		}
		norSegmentUpdatedWithSNPorNotBool[tmpSegGroupIndexInThisSegInfo] = true;
	}

	void update_segGroupWithTargetMap2SNPseq_varySNPmer(int tmpSegGroupIndexInThisSegInfo, string& tmpSegSeq, 
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp, 
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, Index_Info* indexInfo)//, int SNPlocInSyntheticSNPseq)
	{
		//char* tmpSegChar = const_cast<char*>(tmpSegSeq.c_str());
		//cout << "start to do update_segGroupWithTargetMap2SNPseq" << endl;
		int targetMapAlignNum;
		vector<unsigned int> targetMapAlignLocVec;
		bool targetMapSuccess_bool = targetMap_generateSNPlocatedSegMapPosVec_varySNPmer(tmpSegSeq, sa_snp, lcpCompress_snp, 
			childTab_snp, chrom_snp, verifyChild_snp, indexInfo_snp, indexInfo, targetMapAlignNum, targetMapAlignLocVec);//, SNPlocInSyntheticSNPseq);
		//cout << "targetMapSuccess_bool: " << targetMapSuccess_bool << endl; 
		if(!targetMapSuccess_bool)
			return;
		int segGroup_alignNum_ori = this->returnSegmentAlignNum(tmpSegGroupIndexInThisSegInfo);
		//cout << "segGroup_alignNum_ori: " << segGroup_alignNum_ori << endl;
		//cout << "targetMapAlignNum: " << targetMapAlignNum << endl;
		if(segGroup_alignNum_ori > CANDALILOC)
			segGroup_alignNum_ori = 0;
		if(segGroup_alignNum_ori + targetMapAlignNum > CANDALILOC)
		{
			norSegmentAlignNum[tmpSegGroupIndexInThisSegInfo] = CANDALILOC;
			// assign map pos
			for(int tmp = 0; tmp < targetMapAlignNum; tmp++) // reversely copy align poses
			{
				int tmpAlignLocIndexInThisSegGroup = CANDALILOC - 1 - tmp;
				this->assignSegmentAlignLoc(tmpSegGroupIndexInThisSegInfo, tmpAlignLocIndexInThisSegGroup,
					targetMapAlignLocVec[tmp]);
			}
		}
		else
		{
			norSegmentAlignNum[tmpSegGroupIndexInThisSegInfo] = segGroup_alignNum_ori + targetMapAlignNum;
			// assign map ps
			for(int tmp = 0; tmp < targetMapAlignNum; tmp++)
			{
				int tmpAlignLocIndexInThisSegGroup = segGroup_alignNum_ori + tmp;
				this->assignSegmentAlignLoc(tmpSegGroupIndexInThisSegInfo, tmpAlignLocIndexInThisSegGroup,
					targetMapAlignLocVec[tmp]);
			}
		}
		norSegmentUpdatedWithSNPorNotBool[tmpSegGroupIndexInThisSegInfo] = true;
	}	

	void update_targetMap2SNPseq(unsigned int* sa_snp, BYTE* lcpCompress_snp, 
		unsigned int* childTab_snp, char* chrom_snp, BYTE* verifyChild_snp, Index_Info* indexInfo_snp, 
		const string& tmpReadSeq, Index_Info* indexInfo, int SNPlocInSyntheticSNPseq)
	{
		//cout << "update_targetMap2SNPseq starts ........" << endl;
		//char& tmpReadChar = const_cast<char*>(tmpReadSeq.c_str());
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{		
			bool alreadyUpdatedWithSNPseqSegMap_bool = norSegmentUpdatedWithSNPorNotBool[tmpSeg];
			if(alreadyUpdatedWithSNPseqSegMap_bool)
				continue;
			int tmpSeg_startLocInRead = norSegmentLocInRead[tmpSeg];
			int tmpSeg_length = norSegmentLength[tmpSeg];
			if(tmpSeg_length <= IGNORED_SEG_LENGTH_WHEN_UPDATE_SEGINFO_WITH_TARGETMAP2SNPSEQINDEX_MAX)
				continue;
			//cout << "start to update seg: " << tmpSeg << endl;
			if(!((tmpSeg_startLocInRead - 1 >= 0)&&(tmpSeg_startLocInRead - 1 < tmpReadSeq.length())
				&&(tmpSeg_startLocInRead - 1 + tmpSeg_length - 1 >= 0)
				&&(tmpSeg_startLocInRead - 1 + tmpSeg_length - 1 < tmpReadSeq.length())))
				continue;
			string tmpSegSeq = tmpReadSeq.substr(tmpSeg_startLocInRead - 1, tmpSeg_length);
			//if((tmpSeg_startLocInRead < 1)||(tmpSeg_startLocInRead - 1 + tmpSeg_length - 1 >= tmpReadSeq.length()))
			//	continue;
			this->update_segGroupWithTargetMap2SNPseq(tmpSeg, tmpSegSeq, sa_snp, lcpCompress_snp, 
				childTab_snp, chrom_snp, verifyChild_snp, indexInfo_snp, indexInfo, SNPlocInSyntheticSNPseq);
		}
	}

	void update_targetMap2SNPseq_varySNPmer(unsigned int* sa_snp, BYTE* lcpCompress_snp, 
		unsigned int* childTab_snp, char* chrom_snp, BYTE* verifyChild_snp, Index_Info* indexInfo_snp, 
		const string& tmpReadSeq, Index_Info* indexInfo)//, int SNPlocInSyntheticSNPseq)
	{
		//cout << "update_targetMap2SNPseq starts ........" << endl;
		//char& tmpReadChar = const_cast<char*>(tmpReadSeq.c_str());
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{		
			bool alreadyUpdatedWithSNPseqSegMap_bool = norSegmentUpdatedWithSNPorNotBool[tmpSeg];
			if(alreadyUpdatedWithSNPseqSegMap_bool)
				continue;
			int tmpSeg_startLocInRead = norSegmentLocInRead[tmpSeg];
			int tmpSeg_length = norSegmentLength[tmpSeg];
			if(tmpSeg_length <= IGNORED_SEG_LENGTH_WHEN_UPDATE_SEGINFO_WITH_TARGETMAP2SNPSEQINDEX_MAX)
				continue;
			//cout << "start to update seg: " << tmpSeg << endl;
			if(!((tmpSeg_startLocInRead - 1 >= 0)&&(tmpSeg_startLocInRead - 1 < tmpReadSeq.length())
				&&(tmpSeg_startLocInRead - 1 + tmpSeg_length - 1 >= 0)
				&&(tmpSeg_startLocInRead - 1 + tmpSeg_length - 1 < tmpReadSeq.length())))
				continue;
			string tmpSegSeq = tmpReadSeq.substr(tmpSeg_startLocInRead - 1, tmpSeg_length);
			//if((tmpSeg_startLocInRead < 1)||(tmpSeg_startLocInRead - 1 + tmpSeg_length - 1 >= tmpReadSeq.length()))
			//	continue;
			this->update_segGroupWithTargetMap2SNPseq_varySNPmer(tmpSeg, tmpSegSeq, sa_snp, lcpCompress_snp, 
				childTab_snp, chrom_snp, verifyChild_snp, indexInfo_snp, indexInfo);//, SNPlocInSyntheticSNPseq);
		}
	}	

	#endif
	void assignSegmentNum(int tmpSegmentNum)
	{
		segmentNum = tmpSegmentNum;
	}

	void copySegmentLength(Seg_Info* tmpSegInfo)
	{
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{
			norSegmentLength[tmp] = tmpSegInfo->returnSegmentLength(tmp);
		}
	}

	void copySegmentLocInRead(Seg_Info* tmpSegInfo)
	{
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{
			norSegmentLocInRead[tmp] = tmpSegInfo->returnSegmentLocInRead(tmp);
		}		
	}

	void assignZeroAlignNumToAllSeg()
	{
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{
			norSegmentAlignNum[tmp] = 0;
		}			
	}

	void assignSegmentAlignNum(int tmpSegGroupIndex, int tmpSegAlignNum)
	{
		norSegmentAlignNum[tmpSegGroupIndex] = tmpSegAlignNum;
	}

	void assignSegmentAlignLoc(int tmpSegGroupIndex, int tmpSegCandiIndex,
		unsigned int tmpSegAlignLoc)
	{
		norSegmentAlignLoc[tmpSegGroupIndex * CANDALILOC + tmpSegCandiIndex]
			= tmpSegAlignLoc;
	}

	void filterLowQualitySeg(PE_Read_Info& peReadInfo, int type)
	{
		//cout << endl << "start to filterLowQualitySeg" << endl;
		for(int tmpSegIndex = 0; tmpSegIndex < segmentNum; tmpSegIndex++)
		{
			//cout << "tmpSegIndex: " << tmpSegIndex << endl;
			int tmpSegLen = norSegmentLength[tmpSegIndex];
			int tmpSegStartLocInRead = norSegmentLocInRead[tmpSegIndex];
			int tmpSegEndLocInRead = tmpSegStartLocInRead + tmpSegLen - 1;
			//cout << "start to checkSegTrustableOrNot" << endl;
			bool tmpSegSeqTrustable = peReadInfo.checkSegTrustableOrNot(type, tmpSegStartLocInRead, tmpSegEndLocInRead);
			//cout << "tmpSegSeqTrustable: " << tmpSegSeqTrustable << endl;
			//cout << "end of checking segTrustableOrNot" << endl;
			if(!tmpSegSeqTrustable)
				norSegmentAlignNum[tmpSegIndex] = 0;
		}
		//cout << "end of filtering lowQualitySeg" << endl << endl;
	}

	bool returnTooManyValSegsOrNot(int readLength, int valLengthToCal, int multiMapRepeatNum_toCal)
	{
		int standard_tooManyValSegs = (readLength/valLengthToCal)*multiMapRepeatNum_toCal;
		int valSegsNum_tmp = 0;
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{	
			if(norSegmentLength[tmp] >= valLengthToCal)
				valSegsNum_tmp += norSegmentAlignNum[tmp];
		}
		if(valSegsNum_tmp >= standard_tooManyValSegs)
			return true;
		else 
			return false;
	}

	bool combineOriSegInfoWithMissingLongSegInfo(Seg_Info& oriSegInfo, MissingLongSeg_Info& missingLongSegInfo)
	{
			int oriSegInfo_repeatRegionIndex = oriSegInfo.returnRepeatRegion_index();
			if(oriSegInfo_repeatRegionIndex >= 0)
			{
				repeatRegion_index = oriSegInfo_repeatRegionIndex;
				segmentNum = 1;//oriSegInfo->returnSegmentNum();
				norSegmentLength[0] = oriSegInfo.returnSegmentLength(0);
				norSegmentLocInRead[0] = 1;
				norSegmentAlignNum[0] = oriSegInfo.returnSegmentAlignNum(0);
				return true;
			}
			else
			{		
				repeatRegion_index = oriSegInfo_repeatRegionIndex;
				unsigned int oriSegInfo_segmentNum = oriSegInfo.returnSegmentNum();
				unsigned int missingLongSegInfo_segmentNum = missingLongSegInfo.returnMissingSegGroupNum();
				
				segmentNum = oriSegInfo_segmentNum + missingLongSegInfo_segmentNum;
				if(segmentNum > SEGMENTNUM)
					return false;
				//cout << "oriSegInfo_segmentNum: " << oriSegInfo_segmentNum << endl;
				//cout << "missingLongSegInfo_segmentNum: " << missingLongSegInfo_segmentNum << endl;

				int tmpSegNum = 0;
				int index_start_tmp = 0;
				for(int tmp = 0; tmp < missingLongSegInfo_segmentNum; tmp++)
				{
					int index_oriSegInfo = missingLongSegInfo.returnMissingSegGroupIndexInOriSegInfo(tmp);
					int index_end_tmp = index_oriSegInfo - 1;

					//copy oriSegInfo to segInfo in the range of (index_start_tmp ~ index_end_tmp) 
					for(int tmp_index = index_start_tmp; tmp_index <= index_end_tmp; tmp_index++)
					{
						norSegmentLength[tmpSegNum] = oriSegInfo.returnSegmentLength(tmp_index);
						norSegmentLocInRead[tmpSegNum] = oriSegInfo.returnSegmentLocInRead(tmp_index);
						norSegmentAlignNum[tmpSegNum] = oriSegInfo.returnSegmentAlignNum(tmp_index);
						if(norSegmentAlignNum[tmpSegNum] <= CANDALILOC)
						{
							for(int tmpAlignNO = 0; tmpAlignNO < norSegmentAlignNum[tmpSegNum]; tmpAlignNO++)
							{
								*(norSegmentAlignLoc + tmpSegNum*CANDALILOC + tmpAlignNO) 
									= oriSegInfo.returnSegmentMapPos(tmp_index, tmpAlignNO);
							}
						}
						tmpSegNum ++;
					}

					//insert missingSegInfo to segInfo 
					norSegmentLength[tmpSegNum] = missingLongSegInfo.returnMissingSegLength(tmp);
					norSegmentLocInRead[tmpSegNum] = oriSegInfo.returnSegmentLocInRead(index_oriSegInfo);
					norSegmentAlignNum[tmpSegNum] = missingLongSegInfo.returnMissingSegAlignNum(tmp);
					if(norSegmentAlignNum[tmpSegNum] <= CANDALILOC)
					{					
						for(int tmpAlignNO = 0; tmpAlignNO < norSegmentAlignNum[tmpSegNum]; tmpAlignNO++)
						{
							*(norSegmentAlignLoc + tmpSegNum*CANDALILOC + tmpAlignNO) 
								= missingLongSegInfo.returnMissingSegAlignLoc(tmp, tmpAlignNO);
						}
					}
					tmpSegNum ++;

					//insert splitted oriSegInfo to segInfo
					norSegmentLength[tmpSegNum] = oriSegInfo.returnSegmentLength(index_oriSegInfo) - norSegmentLength[tmpSegNum-1];
					norSegmentLocInRead[tmpSegNum] = norSegmentLocInRead[tmpSegNum-1] + norSegmentLength[tmpSegNum-1];
					norSegmentAlignNum[tmpSegNum] = oriSegInfo.returnSegmentAlignNum(index_oriSegInfo);
					if(norSegmentAlignNum[tmpSegNum] <= CANDALILOC)
					{										
						for(int tmpAlignNO = 0; tmpAlignNO < norSegmentAlignNum[tmpSegNum]; tmpAlignNO++)
						{
							*(norSegmentAlignLoc + tmpSegNum*CANDALILOC + tmpAlignNO) 
								= oriSegInfo.returnSegmentMapPos(index_oriSegInfo, tmpAlignNO) + norSegmentLength[tmpSegNum-1];
						}
					}
					tmpSegNum ++;
					index_start_tmp = index_oriSegInfo + 1;
				}

				int lastMissingSegIndex = missingLongSegInfo.returnMissingSegGroupIndexInOriSegInfo(
					missingLongSegInfo_segmentNum-1);
				for(int tmp = lastMissingSegIndex + 1; tmp < oriSegInfo_segmentNum; tmp++)
				{
					norSegmentLength[tmpSegNum] = oriSegInfo.returnSegmentLength(tmp);
					norSegmentLocInRead[tmpSegNum] = oriSegInfo.returnSegmentLocInRead(tmp);
					norSegmentAlignNum[tmpSegNum] = oriSegInfo.returnSegmentAlignNum(tmp);
					if(norSegmentAlignNum[tmpSegNum] <= CANDALILOC)
					{
						for(int tmpAlignNO = 0; tmpAlignNO < norSegmentAlignNum[tmpSegNum]; tmpAlignNO++)
						{
							*(norSegmentAlignLoc + tmpSegNum*CANDALILOC + tmpAlignNO)
								= oriSegInfo.returnSegmentMapPos(tmp, tmpAlignNO);
						}
					}
					tmpSegNum++;	
				}
				return true;			
			}		
	}

	void assignLongSegMinLength(int newLongSegMinLength)
	{
		longSegMinLength = newLongSegMinLength;
	}
	int returnRepeatRegion_index()
	{
		return repeatRegion_index;
	}
	unsigned int returnSegmentNum()
	{
		return segmentNum;
	}
	int returnSegmentLength(int segGroupNO)
	{
		return (int)norSegmentLength[segGroupNO];
	}
	int returnSegmentLocInRead(int segGroupNO)
	{
		return (int)norSegmentLocInRead[segGroupNO];
	}
	int returnSegLocInRead(int segGroupNO)
	{
		return (int)norSegmentLocInRead[segGroupNO];
	}
	int returnSegLocInRead_end(int segGroupNO)
	{
		int endLoc = this->returnSegLocInRead(segGroupNO)
			+ this->returnSegmentLength(segGroupNO) - 1;
		return endLoc;
	}		
	int returnSegmentAlignNum(int segGroupNO)
	{
		return (int)norSegmentAlignNum[segGroupNO];
	}
	unsigned int returnSegmentMapPos(int segGroupNO, int segCandiNO)
	{
		unsigned int mapPos = *(norSegmentAlignLoc + segGroupNO*CANDALILOC + segCandiNO);
		return mapPos;
	}
	int returnLongSegMinLength()
	{
		return longSegMinLength;
	}
	unsigned int returnSegmentMapPos_end(int segGroupNO, int segCandiNO)
	{
		unsigned int mapPos = this->returnSegmentMapPos(segGroupNO, segCandiNO) 
			+ this->returnSegmentLength(segGroupNO) - 1;
		//*(norSegmentAlignLoc + segGroupNO*CANDALILOC + segCandiNO) + ;
		return mapPos;		
	}

	void addMidPartSeg_incompleteHead(int midPartSegLength, int midPartSegLocInRead, 
		int midPartMapChrInt, int midPartMapPosInChr, Index_Info* indexInfo)
	{
		segmentNum ++;
		norSegmentLength[segmentNum - 1] = midPartSegLength;
		norSegmentLocInRead[segmentNum - 1] = midPartSegLocInRead;
		norSegmentAlignNum[segmentNum - 1] = 1;
		unsigned int tmpMidPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			midPartMapChrInt, midPartMapPosInChr);
		*( norSegmentAlignLoc + (segmentNum - 1)*CANDALILOC ) = tmpMidPartMapPosInWholeGenome;
	}

	void addMidPartSeg_incompleteTail(int midPartSegLength, int midPartSegLocInRead, 
		int midPartMapChrInt, int midPartMapPosInChr, Index_Info* indexInfo, Seg_Info* segInfo)
	{
		segmentNum = segInfo->segmentNum + 1;
		norSegmentLength[0] = midPartSegLength;
		norSegmentLocInRead[0] = 1;//midPartSegLocInRead;
		norSegmentAlignNum[0] = 1;
		
		//int tmpUnfixedTailLocInRead = midPartSegLocInRead + midPartSegLength;

		unsigned int tmpMidPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			midPartMapChrInt, midPartMapPosInChr);

		norSegmentAlignLoc[0] = tmpMidPartMapPosInWholeGenome; //midPartMapPosInChr;

		for(int tmp = 1; tmp < segmentNum; tmp ++)
		{
			norSegmentLength[tmp] = segInfo->norSegmentLength[tmp - 1];
			norSegmentLocInRead[tmp] = segInfo->norSegmentLocInRead[tmp - 1] + midPartSegLength;
			norSegmentAlignNum[tmp] = segInfo->norSegmentAlignNum[tmp - 1];
			
			for(int tmp2 = 0; tmp2 < norSegmentAlignNum[tmp]; tmp2++)
			{
				*(norSegmentAlignLoc + (tmp*CANDALILOC) + tmp2) 
					= *(segInfo->norSegmentAlignLoc + ((tmp-1) * CANDALILOC) + tmp2); 
			}
		}
	}


	Seg_Info(Seg2ndOri_Info* seg2ndOriInfo, int mapPosIntervalStart, int mapPosIntervalEnd, 
		int chrPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& chromNameStr)
	{
		repeatRegion_index = -1;

		longSegMinLength = LONG_SEG_LENGTH_THRESHOLD_PHASE2;
		segmentNum = seg2ndOriInfo->returnSegmentNum();
		//cout << "segmentNum: " << segmentNum << endl;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			//cout << "tmpSeg: " << tmpSeg << endl;

			norSegmentLength[tmpSeg] = seg2ndOriInfo->returnSegmentLength(tmpSeg);
			norSegmentLocInRead[tmpSeg] = seg2ndOriInfo->returnSegmentLocInRead(tmpSeg);
			int tmpSegCandiNum = 0;

			if(seg2ndOriInfo->returnSegmentAlignNum(tmpSeg) <= CANDALILOC)
			{
				for(int tmpSegCandi = 0; tmpSegCandi < seg2ndOriInfo->returnSegmentAlignNum(tmpSeg);
					tmpSegCandi++)
				{

					int tmpLoc //= *(seg2ndOriInfo->norSegmentAlignLoc + tmpSeg*CANDALILOC + tmpSegCandi) 
						= seg2ndOriInfo->returnSegmentMapPos(tmpSeg, tmpSegCandi)
							+ chrPosStartIn2ndLevelIndex - 1;
					
					if((tmpLoc <= mapPosIntervalEnd)||(tmpLoc >= mapPosIntervalStart))
					{
						*(norSegmentAlignLoc + tmpSeg*CANDALILOC + tmpSegCandi) 
							= indexInfo->getWholeGenomeLocation(chromNameStr, tmpLoc);
						//seg2ndOriInfo->assignSegmentMapPos(
						//	indexInfo->getWholeGenomeLocation(chromNameStr, tmpLoc),
						//	tmpSeg, tmpSegCandi);
						tmpSegCandiNum ++;
					}
				}
			}
			else
			{
				tmpSegCandiNum = 0;
			}
			norSegmentAlignNum[tmpSeg] = tmpSegCandiNum;
		}
	}		

	bool checkSegLongOrNot(int segGroupNO)
	{
		#ifdef PERSONALIZED_CHR_SEQ
		if(norSegmentUpdatedWithSNPorNotBool[segGroupNO])
			return true;
		#endif

		return (norSegmentLength[segGroupNO] >= longSegMinLength);
	}
	
	int checkSegRelation(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO, int spliceJunctionDistanceMax)
	{
		unsigned int segEndNum1 = segNO_1;
		unsigned int segStartNum2 = segNO_2;

		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		if (segEndNum1 < segStartNum2)
		{
			if ((segStartNum2 - segEndNum1) == 1)
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
				{
					unsigned int gapInChr = alignLoc1 - alignLoc2;
					if (gapInChr > MAX_INSERTION_LENGTH)					
					{
						#ifdef DETECT_CIRCULAR_RNA
						if(gapInChr > spliceJunctionDistanceMax)
						{
							return FIX_TOO_CLOSE;
						}
						else
						{
							return FIX_CIRCULAR_RNA_NEIGHBOUR;
						}
						#else
						return FIX_TOO_CLOSE;
						#endif
					}
					else
						return FIX_INSERTION_NEIGHBOUR; 
				}
				else
				{
					unsigned int gapInChr = alignLoc2 - alignLoc1;
					if (gapInChr > spliceJunctionDistanceMax)
						return FIX_TOO_FAR;
					else if(gapInChr > MAX_DELETION_LENGTH)
						return FIX_SPLICE_NEIGHBOUR;
					else
						return FIX_DELETION_NEIGHBOUR;
				}
			}
			else
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
				{
					unsigned int gapInChr = alignLoc1 - alignLoc2;
					if (gapInChr > MAX_INSERTION_LENGTH)
					{
						#ifdef DETECT_CIRCULAR_RNA
						if(gapInChr > spliceJunctionDistanceMax)
						{
							return FIX_TOO_CLOSE;
						}
						else
						{
							return FIX_CIRCULAR_RNA_GAP;
						}						
						#else
						return FIX_TOO_CLOSE;
						#endif
					}
					else
						return FIX_INSERTION_GAP; 
				}
				else
				{
					unsigned int gapInChr = alignLoc2 - alignLoc1;
					if (gapInChr > spliceJunctionDistanceMax)
						return FIX_TOO_FAR;
					else if(gapInChr > MAX_DELETION_LENGTH)
						return FIX_SPLICE_GAP;
					else
						return FIX_DELETION_GAP;
				}
			}		
		}

		else
			return FIX_NO_RELATIONSHIP;		
	}
	int checkSegRelation_phase2(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO, int spliceJunctionDistanceMax)
	{
		unsigned int segEndNum1 = segNO_1;
		unsigned int segStartNum2 = segNO_2;

		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		if (segEndNum1 < segStartNum2)
		{
			if ((segStartNum2 - segEndNum1) == 1)
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
				{
					unsigned int gapInChr = alignLoc1 - alignLoc2;
					if (gapInChr > MAX_INSERTION_LENGTH)					
						return FIX_TOO_CLOSE;
					else
						return FIX_INSERTION_NEIGHBOUR; 
				}
				else
				{
					unsigned int gapInChr = alignLoc2 - alignLoc1;
					if (gapInChr > spliceJunctionDistanceMax)
						return FIX_TOO_FAR;
					else if(gapInChr > MAX_DELETION_LENGTH)
						return FIX_SPLICE_NEIGHBOUR;
					else
						return FIX_DELETION_NEIGHBOUR;
				}
			}
			else
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
				{
					unsigned int gapInChr = alignLoc1 - alignLoc2;
					if (gapInChr > MAX_INSERTION_LENGTH)
						return FIX_TOO_CLOSE;
					else
						return FIX_INSERTION_GAP; 
				}
				else
				{
					unsigned int gapInChr = alignLoc2 - alignLoc1;
					if (gapInChr > spliceJunctionDistanceMax)
						return FIX_TOO_FAR;
					else if(gapInChr > MAX_DELETION_LENGTH)
						return FIX_SPLICE_GAP;
					else
						return FIX_DELETION_GAP;
				}
			}		
		}

		else
			return FIX_NO_RELATIONSHIP;		
	}
	bool checkSegRelationValidOrNotForNormalAlignment(int relation)
	{
		if(
			(relation == FIX_MATCH) || (relation == FIX_INSERTION_NEIGHBOUR) 
			|| (relation == FIX_SPLICE_NEIGHBOUR) || (relation == FIX_DELETION_NEIGHBOUR)
			|| (relation == FIX_INSERTION_GAP) || (relation == FIX_SPLICE_GAP)
			|| (relation == FIX_DELETION_GAP) 
			)
			return true;
		else
			return false;

	}

	int checkSegRelation_cirRNA(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		unsigned int segEndNum1 = segNO_1;
		unsigned int segStartNum2 = segNO_2;

		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		if (segEndNum1 < segStartNum2)
		{
			if ((segStartNum2 - segEndNum1) == 1)
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
				{
					unsigned int gapInChr = alignLoc1 - alignLoc2;
					if (gapInChr > MAX_INSERTION_LENGTH)
					{
						if(gapInChr <= MAX_SPLICE_DISTANCE_PHASE1)
						{
							return FIX_CIRCULAR_RNA;
						}
						else
						{
							return FIX_TOO_CLOSE;						
						}
					}					
					else
					{
						return FIX_INSERTION_NEIGHBOUR; 
					}
				}
				else
				{
					unsigned int gapInChr = alignLoc2 - alignLoc1;
					if (gapInChr > MAX_SPLICE_DISTANCE_PHASE1)
						return FIX_TOO_FAR;
					else if(gapInChr > MAX_DELETION_LENGTH)
						return FIX_SPLICE_NEIGHBOUR;
					else
						return FIX_DELETION_NEIGHBOUR;
				}
			}
			else
			{
				if(alignLoc1 == alignLoc2)
				{
					return FIX_MATCH;
				}
				else if (alignLoc1 > alignLoc2)
				{
					unsigned int gapInChr = alignLoc1 - alignLoc2;
					if (gapInChr > MAX_INSERTION_LENGTH)
						return FIX_TOO_CLOSE;
					else
						return FIX_INSERTION_GAP; 
				}
				else
				{
					unsigned int gapInChr = alignLoc2 - alignLoc1;
					if (gapInChr > MAX_SPLICE_DISTANCE_PHASE1)
						return FIX_TOO_FAR;
					else if(gapInChr > MAX_DELETION_LENGTH)
						return FIX_SPLICE_GAP;
					else
						return FIX_DELETION_GAP;
				}
			}		
		}

		else
			return FIX_NO_RELATIONSHIP;		
	}

	int distanceBetweenSegment(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO, int spliceJunctionDistanceMax)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
			
			if(segmentDistance <= spliceJunctionDistanceMax)
				return segmentDistance;
			else
				return segDist;
		}
		else
		{
			segmentDistance = alignLoc1 - alignLoc2;
			if(segmentDistance <= MAX_INSERTION_LENGTH)
				return segmentDistance;
			else
				return segDist;
		}
	}

	int distanceBetweenSegment_signed(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO, int spliceJunctionDistanceMax)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
			
			if(segmentDistance <= spliceJunctionDistanceMax)
				return segmentDistance;
			else
				return segDist;
		}
		else
		{
			segmentDistance = alignLoc1 - alignLoc2;
			if(segmentDistance <= MAX_INSERTION_LENGTH)
				return 0-segmentDistance;
			else if(segmentDistance <= spliceJunctionDistanceMax)
			{
				return 0-segmentDistance;
			}
			else
				return segDist;
		}
	}	


	int distanceBetweenSegment_cirRNA(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 = *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
			
			if(segmentDistance <= MAX_SPLICE_DISTANCE_PHASE1)
				return segmentDistance;
			else
				return segDist;
		}
		else
		{
			segmentDistance = alignLoc1 - alignLoc2;
			if(segmentDistance <= MAX_INSERTION_LENGTH)
				return segmentDistance;
			else if(segmentDistance <= MAX_SPLICE_DISTANCE_PHASE1)
				return segmentDistance;
			else
				return segDist;
		}
		//unsigned int segmentDistance = alignLoc2 - alignLoc1;
		//if((segmentDistance < 300000))
		//{
		//	segDist = segmentDistance;
		//} 		
		//return segDist;
	}

	bool checkSegRelatedOrNot(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO, int spliceJunctionDistanceMax)
	{
		int segRelation = checkSegRelation(segNO_1, segNO_1_candiNO, segNO_2, segNO_2_candiNO, spliceJunctionDistanceMax);
		//cout << endl << "SegRelation: " << segRelation << endl;
		if((segRelation == FIX_TOO_FAR)||(segRelation == FIX_TOO_CLOSE)||(segRelation == FIX_NO_RELATIONSHIP))
		{
			return false;
		} 
		else
		{
			return true;
		}
	}

	int getFirstLongSegNO()
	{
		bool longSegExists = false;
		for(int tmpSegNO = 0; tmpSegNO < segmentNum; tmpSegNO++)
		{
			#ifdef PERSONALIZED_CHR_SEQ
			if((norSegmentUpdatedWithSNPorNotBool[tmpSegNO])&&(norSegmentAlignNum[tmpSegNO] <= CANDALILOC))
			{
				longSegExists = true;
				return tmpSegNO;				
			}
			#endif
			if((norSegmentLength[tmpSegNO] >= longSegMinLength)&&(norSegmentAlignNum[tmpSegNO] <= CANDALILOC))
			{
				longSegExists = true;
				return tmpSegNO;
			}
		}
		if(!longSegExists)
			return -1;
	}

	bool greedyMapWithoutPreIndexArray(
		char *read, unsigned int* sa, BYTE* lcpCompress, unsigned int* child, char* chrom, 
		BYTE* verifyChild, int readLength, Index_Info* indexInfo)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		
		//cout << "start to mapMainSecondLevel_compressedIndex " << endl;
		//*(read + readLength) = 'Y';
		unsigned int norSegmentNum;

		(norSegmentNum) = 0;
		bool mapMain = false;	
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_num = 0;
		unsigned int segment_length = 0; 
		unsigned int segment_length_max = 0;//used to compare with segment_length for eache segment to get the maximum length segment
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int read_length = readLength; //READ_LENGTH;
		unsigned int interval_begin, interval_end;
		unsigned int n = (indexInfo->returnIndexSize());//size of SA
		unsigned int norAlignLoc;
		unsigned int align_chr_location;
		unsigned int align_chr;
		//*valLength = 0;
		char* read_local = read;
		while (stop_loc_overall < read_length) //- 15)
		{
			segment_num++;
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
	   	 		//align_length[0] ++;
	   	 		stop_loc = 1;
	   	 		segment_align_SArange[0] = 1;
	   	 		segment_align_SArange[1] = 0;
	   	 		segment_align_rangeNum = 0;
	   	 		queryFound = false;   	
	   	 		if(norSegmentNum >= SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	   	 		
	   	 		(norSegmentNum) ++;
	   	 		norSegmentLength[norSegmentNum - 1] = 1;
				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = 0;				
				//cout << "norSegmentNum: " << norSegmentNum << endl;
				//cout << "norSegmentLength: " << norSegmentLength[norSegmentNum - 1] << endl;
				//cout << "norSegmentLocInRead: "<< norSegmentLocInRead[norSegmentNum - 1] << endl;
				//cout << "norSegmentAlignNum: " << norSegmentAlignNum[norSegmentNum - 1] << endl << endl;
				stop_loc = 1;	
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
	   	 		continue;
	   	 	}
	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;
  		 	if(interval_end < interval_begin - 1)
  	 			segment_align_rangeNum = 0;
	   	 	unsigned int iterateNum = 0;//debug;
	   	 	while((c + stop_loc_overall< read_length) && (queryFound == true))
	   	 	{
	   	 		iterateNum++;
	   	 		if(iterateNum>read_length)
	   	 			return false;
	   	 		unsigned int c_old = c;
				
				if(interval_begin != interval_end)
				{ 
	           		//Xinan: COMPRESS INDEX
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild);
					Min = min(lcp_length, read_length - stop_loc_overall);
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

	            	if(!queryFound)
	            	{
	            		stop_loc = c_old + loc_pos;
	            		break;
	            	}
	            	
	            	c = Min;
	            	if(*(read_local+c) == 'N')
	            	{
	            		queryFound = false; 
	            		stop_loc = c;
	            		break;
	            	}
					start = interval_begin; end = interval_end;
					if (c + stop_loc_overall== read_length)			
						break;			
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
			    	//cout << "char: " << *(read_local+c) << " interval: " << interval_begin << " ~ " << interval_end << endl;
			    	if(interval_begin > interval_end)
			    	{
			    		queryFound = false;
			    		stop_loc = c-1;
	          			segment_align_SArange[0] = interval_begin_ori;
	            		segment_align_SArange[1] = interval_end_ori;
	            		segment_align_rangeNum = interval_end_ori - interval_begin_ori + 1;		    			
		   			 	if(interval_end < interval_begin - 1)
		   		 			segment_align_rangeNum = 0;
			    		break;
			    	}
			    	else
			    	{
	          			segment_align_SArange[0] = interval_begin;
	            		segment_align_SArange[1] = interval_end;
	            		segment_align_rangeNum = interval_end - interval_begin + 1;
		   			 	if(interval_end < interval_begin - 1)
		   		 			segment_align_rangeNum = 0;
			    	}
				}//end if
				else 
				{
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < read_length - c- stop_loc_overall; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(queryFound) 
		    		{}
		    		else 
		    			stop_loc = c+loc_pos;	
	          		segment_align_SArange[0] = interval_begin;
	            	segment_align_SArange[1] = interval_end;
	            	segment_align_rangeNum = interval_end - interval_begin + 1;   	
	   			 	if(interval_end < interval_begin - 1)
			 			segment_align_rangeNum = 0;
		    		break;
		    	}
			} //end while
			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////   	 
			/////////////////////////////////////////SEGMENT MAP RESULT////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{
	    		(norSegmentNum) ++;
	   	 		if(norSegmentNum > SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
	   	 		if(norSegmentNum > (int)(read_length/5))
	   	 		{
	   	 			segmentNum = (int)(read_length/5);
	   	 			return false;
	   	 		}	    		
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;
	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;								
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)  			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
				segment_length = read_length-stop_loc_overall;
				break;
			}
			else 
			{    
				(norSegmentNum) ++;
	   	 		if(norSegmentNum > SEGMENTNUM)
	   	 		{
	   	 			segmentNum = SEGMENTNUM;
	   	 			return false;
	   	 		}
				
				if(norSegmentNum > (int)(read_length/5))
				{
					segmentNum = (int)(read_length/5);
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[norSegmentNum - 1] = stop_loc;
				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++) 			
	    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
	    		segment_length = stop_loc;
			}		
	   	}

		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//debugln("mapMain ended!!!");
		// Xinan: fixMe: some potential bugs -- segmentLength > read length, to fix
		//int tmpSegLengthSum = 0;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			//tmpSegLengthSum = tmpSegLengthSum + norSegmentLength[tmpSeg]
			if(norSegmentLength[tmpSeg] > readLength)
			{
				segmentNum = 0;
				norSegmentNum = 0;
				return false;
	 		}
			if((norSegmentAlignNum[tmpSeg] > CANDALILOC)||(norSegmentAlignNum[tmpSeg] < 0))
				norSegmentAlignNum[tmpSeg] = 0;
		}
		
		return mapMain;
	}	

	bool mapMain_SegInfo_preIndex_repeatRegion_keepMissingLongSeg(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray, const string& readStringStr, 
		RepeatRegion_Info* repeatRegionInfo, MissingLongSeg_Info& missingLongSegInfo)
	{
		//cout << endl << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		//cout << "mapMain_SegInfo_preIndex function starts ...... "<< endl; 
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		vector<int> mismatchPosVecInRead;
		vector<char> mismatchCharVec;
		unsigned int norSegmentNum = 0;
		bool mapMain = false;	
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_num = 0;
		//unsigned int segment_length = 0; 
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int read_length = readLength; //READ_LENGTH;
		unsigned int interval_begin, interval_end;
		unsigned int n = (indexInfo->returnIndexSize());//size of SA
		*valLength = 0;
		char* read_local = read;

		bool noKeepingMissingLongSeg_bool = true;
		while (stop_loc_overall < read_length) //- 15)
		{
			//cout << "stop_loc_overall: " << stop_loc_overall << endl;
			segment_num++;
			bool queryFound = true;
			noKeepingMissingLongSeg_bool = true;
			/////////////////////////////////////////////////////////////////////////
			/////////////////////////   K-mer search  ///////////////////////////////
			bool firstCheckForHeadOfSegmentBool = true;
			bool KmerSearchFound = false;
			int KmerMappedLength;
			unsigned int KmerSearchIndexIntervalStart = 0;
			unsigned int KmerSearchIndexIntervalEnd = 0;			
			//#ifdef CAL_TIME
			//searchPrefix_begin = clock();
			//#endif	
			//cout << "readLength - stop_loc_overall: " << read_length - stop_loc_overall << " INDEX_KMER_LENGTH: " << INDEX_KMER_LENGTH << endl;
			if(read_length - stop_loc_overall < INDEX_KMER_LENGTH)
			{
				KmerSearchFound = false;
				firstCheckForHeadOfSegmentBool = false;

				// to fix, to reduce the search time for the last segment
				norSegmentNum ++;
	    		norSegmentLength[norSegmentNum - 1] = read_length - stop_loc_overall;//READ_LENGTH - stop_loc_overall;
	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = 0;		
				break;
			}
			else
			{
				KmerSearchFound = this->getIndexInterval_PreIndex(readStringStr.substr(stop_loc_overall, INDEX_KMER_LENGTH), &KmerMappedLength,
					&KmerSearchIndexIntervalStart, &KmerSearchIndexIntervalEnd, PreIndexMappedLengthArray, 
					PreIndexIntervalStartArray,	PreIndexIntervalEndArray);
				//cout << "KmerSearchFound: " << KmerSearchFound << " KmerMappedLength: " << KmerMappedLength << endl;
			}
			
			///////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////

	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	   	 	unsigned int iterateNum = 0;//debug;

	   	 	if(!KmerSearchFound)
	   	 	{	
	   	 		//cout << "Kmer not found !" << endl; 
		   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
		   	 	{
		   	 		queryFound = false;
		   	 		stop_loc = 0;
		   	 		segment_align_SArange[0] = 1;
		   	 		segment_align_SArange[1] = 0;
		   	 		segment_align_rangeNum = 0;
		   	 		queryFound = false;   	 			
		   	 		
		   	 		if(norSegmentNum >= SEGMENTNUM)
		   	 		{
		   	 			segmentNum = SEGMENTNUM;
		   	 			return false;
		   	 		}
		   	 		
		   	 		(norSegmentNum) ++;

		   	 		norSegmentLength[norSegmentNum - 1] = 1;
					norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
					norSegmentAlignNum[norSegmentNum - 1] = 0;
					
					stop_loc = 1;	
					read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
					stop_loc_overall = stop_loc_overall + stop_loc + 1;   	 		
		   	 		continue;		   	 		
		   	 	}
	   	 		lcp_length = 0;
	   	 		start = 0, end = n-1;
	   	 		c = 0;
		   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
		   	 	segment_align_SArange[0] = interval_begin;
		   	 	segment_align_SArange[1] = interval_end;
		   	 	segment_align_rangeNum = interval_end - interval_begin + 1;
		   	 	if(interval_end < interval_begin - 1)
		   	 		segment_align_rangeNum = 0;
		   	 	iterateNum = 0;
	   	 	}
	   	 	else // K-mer found in preIndex base, then start to extend the found segment (Kmer) in the target SA interval
	   	 	{
	   	 		//cout << "Kmer found !" << endl;
	   	 		firstCheckForHeadOfSegmentBool = true;
	   	 		lcp_length = 0;
	   	 		start = 0, end = n-1;
	   	 		c = 0;
	   	 		interval_begin = KmerSearchIndexIntervalStart;
	   	 		interval_end = KmerSearchIndexIntervalEnd;
	    	 	segment_align_SArange[0] = interval_begin;
	   	 		segment_align_SArange[1] = interval_end;
	   	 		segment_align_rangeNum = interval_end - interval_begin + 1;
		   	 	if(interval_end < interval_begin - 1)
		   	 		segment_align_rangeNum = 0;
	   	 		//cout << "interval_begin: " << interval_begin << endl;
	   	 		//cout << "interval_end: " << interval_end << endl;
	   	 		//cout << "segment_align_rangeNum: " << segment_align_rangeNum << endl;
	   	 		iterateNum = 0;//debug;  	 		
	   	 	}

	   	 	while((c + stop_loc_overall < read_length) && (queryFound == true))
	   	 	{
	   	 		//cout << "c: " << c << " stop_loc_overall: " << stop_loc_overall << endl;
	   	 		firstCheckForHeadOfSegmentBool = false;
	   	 		iterateNum++;
	   	 		if(iterateNum + stop_loc_overall >read_length)
	   	 		{
	   	 			return false;
	   	 		}
	   	 		unsigned int c_old = c;

				if(interval_begin != interval_end)
				{ 
					//cout << "interval_begin != interval_end " << endl;
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
	 				//cout << "lcp_length: " << lcp_length << endl;
					Min = min(lcp_length, read_length - stop_loc_overall);
					//cout << "Min: " << Min << endl;
					unsigned int loc_pos = 0;

					int startToCheckPos = 0;
					if(firstCheckForHeadOfSegmentBool)
						startToCheckPos = KmerMappedLength;
					//cout << "startToCheckPos: " << startToCheckPos << " Min-c_old: " << Min-c_old << endl;
		            for(loc_pos = startToCheckPos; loc_pos < Min - c_old; loc_pos++)
		            {
		            	queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            	//cout << "...queryFound: " << queryFound << endl;
		            	if (!queryFound)
		            	{	
		   					//cout << "queryFound failed ...." << endl << " loc_pos: " << loc_pos << " c_old: " << c_old << endl;
		            		break;
		            	}
		            }
            		//cout << "queryFound: " << queryFound << endl;
	            	if(!queryFound)
	            	{
	            		//cout << "query not found for lcp ! " << endl;
	            		stop_loc = c_old + loc_pos;
	            		//cout << "stop_loc: " << stop_loc << endl;
	            		break;
	            	}
	            	//cout << "query found for lcp !" << endl;
	            	c = Min;

	            	if(*(read_local+c) == 'N')
	            	{
	            		queryFound = false; 
	            		stop_loc = c;
	            		break;
	            	}
					start = interval_begin; end = interval_end;

					if (c + stop_loc_overall == read_length)
					{		
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
					//cout << "interval_begin_ori: " << interval_begin_ori << " interval_end_ori: " << interval_end << endl;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
			    	//cout << "interval_begin: " << interval_begin << " interval_end: " << interval_end << endl;
			    	////////////////////////////////////////////
			    	///////    keep all missing long segments	
			    	////////////////////////////////////////////
			    	// 08/06/2015 fixMe: try speedUp
					// bool missingLongSeg_bool = missingLongSegInfo.checkSeg_missingLong_orNot(
					// 	interval_begin_ori, interval_end_ori, interval_begin, interval_end, c);

					// if(missingLongSeg_bool && noKeepingMissingLongSeg_bool)
					// {
					// 	missingLongSegInfo.pushBackNewMissingLongSegGroup(interval_begin_ori, 
					// 		interval_end_ori, c, norSegmentNum, sa, indexInfo);
					// 	noKeepingMissingLongSeg_bool = false;
					// }
			    	////////////////////////////////////////////
			    	///////    keep all missing long segments	
			    	////////////////////////////////////////////
			    	if(interval_begin > interval_end)
			    	{
			    		//cout << "interval_begin > interval_end" << endl;
			    		queryFound = false;
			    		stop_loc = c;//c-1;// fixed 11/04/2015
	          			segment_align_SArange[0] = interval_begin_ori;
	            		segment_align_SArange[1] = interval_end_ori;       		
	            		segment_align_rangeNum = interval_end_ori - interval_begin_ori + 1;		    			
		   			 	if(interval_end < interval_begin - 1)
		   		 			segment_align_rangeNum = 0;
			    		break;
			    	}
			    	else
			    	{
			    		//cout << "interval_begin <= interval_end" << endl;
	          			segment_align_SArange[0] = interval_begin;
	            		segment_align_SArange[1] = interval_end;
	            		segment_align_rangeNum = interval_end - interval_begin + 1; // == 1
	            		if(interval_end < interval_begin - 1)
		   	 				segment_align_rangeNum = 0;
			    	}
				}//end if
				else // interval_begin == interval_end
				{
					//cout << "interval_begin == interval_end" << endl;
					unsigned int loc_pos = 0;
		            for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
		            {
		            	queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            	if (!queryFound)
		            	{	
		            		//cout << "queryFound failed ...." << endl << " loc_pos: " << loc_pos << " c: " << c << endl;
		            		break;
		            	}
		            }
		            //cout << "queryFound !" << endl;
		            //cout << "stop_loc: " << stop_loc << " c: " << c << " loc_pos: " << loc_pos << endl;
		    		if(queryFound) 
		    		{}
		    		else 
		    			stop_loc = c+loc_pos;	
	          		segment_align_SArange[0] = interval_begin;
	            	segment_align_SArange[1] = interval_end;
	            	segment_align_rangeNum = interval_end - interval_begin + 1;   	
		   		 	if(interval_end < interval_begin - 1)
		   	 			segment_align_rangeNum = 0;
		    		break;
		    	}
			} //end while
			///////////////////////////////////////////////////////////////////////////////////////////////   	 
			/////////////////////////////////////////SEGMENT MAP RESULT////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{
	    		//cout << "\n$$-01\nqueryFound && (interval_end >= interval_begin)" << endl;
	    		(norSegmentNum) ++;

	    		if(norSegmentNum > SEGMENTNUM)
	    		{
	    			segmentNum = (SEGMENTNUM);
	    			return false;
	    		}
	   	 		if(norSegmentNum > (int)(read_length/5))
	   	 		{
	   	 			segmentNum = (int)(read_length/5);
	   	 			return false;
	   	 		}
	    		//cout << "segmentNum: " << norSegmentNum << endl;
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;

	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				
				// cout << endl << "new seg !" << endl;
				// cout << "segNum: " << norSegmentNum << endl;
	   			//cout << "segmentLength: " << tmpSegLength << endl;
	   			//cout << "segmentStartLocInRead: " << stop_loc_overall + 1 << endl;
	   			//cout << "segmentAlignRangeNum: " << segment_align_rangeNum << endl;

				////////////////////  check if mapped to repeatRegion  ////////////////////
	    		if((norSegmentNum == 1)&&(tmpSegLength == read_length)&&(segment_align_rangeNum > CANDALILOC))
	    		{
	    			int tmpIntervalBegin = interval_begin;
	    			int tmpIntervalEnd = interval_end;
	    			repeatRegion_index = 
	    				repeatRegionInfo->push2RepeatRegionInfo(
	    					tmpIntervalBegin, tmpIntervalEnd, sa, indexInfo);
	    		}

	    		if(segment_align_rangeNum <= CANDALILOC)
	    		{
					for (unsigned int alignment_num = 0; alignment_num < segment_align_rangeNum; alignment_num++)
					{    			
		    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
		    			// cout << "tmpSegPosInChr: " << *(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) << endl;
		    		}
	    		}
	    		//cout << "endl of seg !" << endl << endl;

				break;
			}
			else 
			{    
				(norSegmentNum) ++;
				if((norSegmentNum > (int)(read_length/5)) || (norSegmentNum > SEGMENTNUM))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					if((int)(read_length/5) > SEGMENTNUM)
						segmentNum = (SEGMENTNUM);
					else
						segmentNum = (int)(read_length/5);

					return false;
				}
				//cout << "stop_loc: " << stop_loc << endl;
				norSegmentLength[norSegmentNum - 1] = stop_loc;

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				////////////////////  check if mapped to repeatRegion  ////////////////////
	    		if((norSegmentNum == 1)&&(stop_loc == read_length)&&(segment_align_rangeNum > CANDALILOC))
	    		{
	    			int tmpIntervalBegin = segment_align_SArange[0];
	    			int tmpIntervalEnd = segment_align_SArange[0] + segment_align_rangeNum - 1;
	    			repeatRegion_index = 
	    				repeatRegionInfo->push2RepeatRegionInfo(
	    					tmpIntervalBegin, tmpIntervalEnd, sa, indexInfo);
	    		}

	    		if(segment_align_rangeNum <= CANDALILOC)
	    		{
	    			/*if((segment_align_rangeNum == 1)
	    				&&(stop_loc >= MIN_SEGMENTLENGTH_TOEXTEND_DURINGSEGMENTMAPPING)
	    				&&(stop_loc_overall + 1 + (stop_loc + 1) <= readLength)) // unique mapping pos; or enough long segment
	    			{
	    				vector<unsigned int> tmpSegMapPosVec_extendStartPosInWholeGenome;
	    				vector<int> tmpExtensionLengthVec;
	    				int tmpExtensionLength_Max = 0;
						for (unsigned int alignment_num = 0; alignment_num < segment_align_rangeNum; alignment_num++)
					    {    			
			    			//*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) 
			    			unsigned int tmpSegMapPosInWholeGenome_extendStartPos = sa[segment_align_SArange[0] + alignment_num] + 1 + (stop_loc + 1);
			    			tmpSegMapPosVec_extendStartPosInWholeGenome.push_back(tmpSegMapPosInWholeGenome_extendStartPos);
			    			int tmpExtensionLength 
			    				= extensionForwards_errorTolerated_finalStepForAligner_wholeGenome_duringSegmengMapping(
			    					tmpSegMapPosInWholeGenome_extendStartPos, indexInfo, read, stop_loc_overall + 1 + (stop_loc + 1), readLength,
			    					mismatchPosVecInRead, mismatchCharVec, 5, 8);
			    			tmpExtensionLengthVec.push_back(tmpExtensionLength);
			    			if(tmpExtensionLength > tmpExtensionLength_Max)
			    				tmpExtensionLength_Max = tmpExtensionLength;			    				
			    		}
			    		int tmpSegmentMappedPosNum = 0;
			    		for(int tmp = 0; tmp < tmpSegMapPosVec_extendStartPosInWholeGenome.size(); tmp++)
			    		{
			    			int tmpExtensionLength_inVec = tmpExtensionLengthVec[tmp];
			    			if(tmpExtensionLength_inVec == tmpExtensionLength_Max)
			    			{
			    				unsigned int tmpSegMapPosVec_segStart = tmpSegMapPosVec_extendStartPosInWholeGenome[tmp] - (stop_loc + 1);
			    				*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + tmpSegmentMappedPosNum) 
			    					= tmpSegMapPosVec_segStart;
			    				tmpSegmentMappedPosNum ++;
			    			}
			    		}
			    		stop_loc = stop_loc + tmpExtensionLength_Max + 1;
			    		norSegmentLength[norSegmentNum - 1] = stop_loc;
			    		norSegmentAlignNum[norSegmentNum - 1] = tmpSegmentMappedPosNum;			
	    			}
	    			else
	    			{*/	
						for (unsigned int alignment_num = 0; alignment_num < segment_align_rangeNum; alignment_num++)
					    {    			
			    			*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
			    		}
			    	//}
		    	}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				//cout << "stop_loc + 1: " << stop_loc +1 << endl;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
			}		
	   	}
		segmentNum = (norSegmentNum);
		mapMain = true;

		// Xinan: fixMe: some potential bugs -- segmentLength > read length, to fix
		//int tmpSegLengthSum = 0;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			//tmpSegLengthSum = tmpSegLengthSum + norSegmentLength[tmpSeg]
			if(norSegmentLength[tmpSeg] > readLength)
			{
				segmentNum = 0;
				norSegmentNum = 0;
				return false;
	 		}
			if((norSegmentAlignNum[tmpSeg] > CANDALILOC)||(norSegmentAlignNum[tmpSeg] < 0))
				norSegmentAlignNum[tmpSeg] = 0;
		}
		return mapMain;
	}
	string segInfoStr(Index_Info* indexInfo)
	{
		string segInfoStr;
		segInfoStr += "\nsegment Info: \n";
		unsigned int align_chr, align_chr_location;
		//cout << "segmentNum: " << segmentNum << endl;
	   	for(unsigned int k1 = 0; k1 < segmentNum; k1++)
	   	{
	   		//cout << "k1: " << k1 << endl;
	   		//cout << "location in read: " << norSegmentLocInRead[k1] << endl;
	   		//cout << "segment length: " << norSegmentLength[k1] << endl;
	   		//cout << "align num: " << norSegmentAlignNum[k1] << endl;
	  		segInfoStr = segInfoStr + "... segment " + int_to_str(k1+1) + ": " + int_to_str(norSegmentLocInRead[k1]) + "~"   
	  		+ int_to_str(norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) + 
	  		"  Length: " + int_to_str(norSegmentLength[k1]) + " Num: " + int_to_str(norSegmentAlignNum[k1]) + "\n";

	      	if(//(norSegmentLength[k1]>=10)&&
	      		(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		//segInfoStr += "...... Align Location: \n";
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			indexInfo->getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			segInfoStr = segInfoStr + "\t" + int_to_str(k2+1) +
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			+ ". " 
	      			+ (indexInfo->returnChrNameStr(align_chr)) + " " + int_to_str(align_chr_location) + "\n";
	      		}
	      	}
	   	}
		return segInfoStr;//+"\n";
	}

	unsigned int getPreIndexNO(const string& readPreStr)
	{
		int preIndexStrSize = readPreStr.length();

		unsigned int preIndexNO = 0;

		int baseForCount = 1;

		for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
		{
			preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
			baseForCount = baseForCount * 4;
		}		
		return preIndexNO;
	}

	bool getIndexInterval_PreIndex(const string& readPreStr, int* mappedLength, 
		unsigned int* indexIntervalStart, unsigned int* indexIntervalEnd,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray)
	{
		////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////	
		bool getIndexInterval = false;

		if (readPreStr.find("N") != readPreStr.npos)
			return false;

		int preIndexStrSize = readPreStr.length();

		unsigned int preIndexNO = 0;

		int baseForCount = 1;
		for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
		{
			preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
			baseForCount = baseForCount * 4;
		}

			(*mappedLength) = PreIndexMappedLengthArray[preIndexNO];
			(*indexIntervalStart) = PreIndexIntervalStartArray[preIndexNO];
			(*indexIntervalEnd) = PreIndexIntervalEndArray[preIndexNO];

		getIndexInterval = true;
		return getIndexInterval;
	}

};


#endif
