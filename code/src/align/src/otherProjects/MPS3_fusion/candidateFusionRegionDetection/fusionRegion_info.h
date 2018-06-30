// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FUSIONREGION_INFO_H
#define FUSIONREGION_INFO_H

#include "genomicRegion_info.h"

using namespace std;

typedef map<int, int> RegionIndexMap;
typedef map<int, RegionIndexMap> RegionPairIndexMap;

class FusionRegion_Info
{
private:
	int genomicRegionIndex_Nor;
	int genomicRegionIndex_Rcm;

	// before fusion, scanning to generate candidateFusionRegions
	int supportReadNum_unpairCompleteBeforeFusion; // encompassing reads, but not the only encompassing reads 
	int supportReadNum_unpairIncompleteBeforeFusion;
	int supportReadNum_pairIncompleteBeforeFusion;

	// after fusion
public:

	FusionRegion_Info()
	{
		supportReadNum_unpairCompleteBeforeFusion = 0; 
		supportReadNum_unpairIncompleteBeforeFusion = 0;
		supportReadNum_pairIncompleteBeforeFusion = 0;		
	}

	int returnTotalSupportNum()
	{
		return (supportReadNum_unpairCompleteBeforeFusion 
			+ supportReadNum_unpairIncompleteBeforeFusion
			+ supportReadNum_pairIncompleteBeforeFusion);
	}

	int returnGenomicRegionIndex_Nor()
	{
		return genomicRegionIndex_Nor;
	}

	int returnGenomicRegionIndex_Rcm()
	{
		return genomicRegionIndex_Rcm;
	}
	string returnFuionRegionInfoStr(
		GenomicRegionVec_Info* genomicRegionVecInfo,
		Index_Info* indexInfo)
	{
		int fusion_chrNameInt_Nor 
			= genomicRegionVecInfo->returnChrNameInt_genomicRegionInfoVec(
				genomicRegionIndex_Nor);
		int fusion_chrNameInt_Rcm 
			= genomicRegionVecInfo->returnChrNameInt_genomicRegionInfoVec(
				genomicRegionIndex_Rcm);

		string fusion_chrNameStr_Nor 
			= indexInfo->returnChrNameStr(fusion_chrNameInt_Nor);
		string fusion_chrNameStr_Rcm 
			= indexInfo->returnChrNameStr(fusion_chrNameInt_Rcm);

		int fusion_chrStartPos_Nor
			= genomicRegionVecInfo->returnChrStartPos_genomicRegionInfoVec(
				genomicRegionIndex_Nor);
		int fusion_chrEndPos_Nor
			= genomicRegionVecInfo->returnChrEndPos_genomicRegionInfoVec(
				genomicRegionIndex_Nor);
		int fusion_chrStartPos_Rcm
			= genomicRegionVecInfo->returnChrStartPos_genomicRegionInfoVec(
				genomicRegionIndex_Rcm);
		int fusion_chrEndPos_Rcm
			= genomicRegionVecInfo->returnChrEndPos_genomicRegionInfoVec(
				genomicRegionIndex_Rcm);

		string tmpStr = int_to_str(genomicRegionIndex_Nor) + "\t"
			+ int_to_str(genomicRegionIndex_Rcm) + "\t"
			+ fusion_chrNameStr_Nor + "\t"
			+ int_to_str(fusion_chrStartPos_Nor) + "\t"
			+ int_to_str(fusion_chrEndPos_Nor) + "\t"
			+ fusion_chrNameStr_Rcm + "\t"
			+ int_to_str(fusion_chrStartPos_Rcm) + "\t"
			+ int_to_str(fusion_chrEndPos_Rcm) + "\t"
			+ int_to_str(supportReadNum_unpairCompleteBeforeFusion) + "\t"
			+ int_to_str(supportReadNum_unpairIncompleteBeforeFusion) + "\t"
			+ int_to_str(supportReadNum_pairIncompleteBeforeFusion);
		return tmpStr;
	}

	void initiate_withGenomicRegionIndex_pair(
		int tmpIndex_genomicRegion_Nor,
		int tmpIndex_genomicRegion_Rcm)
	{
		genomicRegionIndex_Nor = tmpIndex_genomicRegion_Nor;
		genomicRegionIndex_Rcm = tmpIndex_genomicRegion_Rcm;
	}

	/*void initiate_chrNameInt_chrMapPos_pair(
		GenomicRegionVec_Info* genomicRegionVecInfo,
		int chrNameInt_Nor, int chrNameInt_Rcm,
		int chrMapPos_Nor, int chrMapPos_Rcm)
	{
		genomicRegionIndex_Nor
			= genomicRegionVecInfo->returnGenomicRegionIndex(
				chrNameInt_Nor, chrMapPos_Nor);
		genomicRegionIndex_Rcm
			= genomicRegionVecInfo->returnGenomicRegionIndex(
				chrNameInt_Rcm, chrMapPos_Rcm);
	}*/

	void increment_unpairCompleteBeforeFusion()
	{
		supportReadNum_unpairCompleteBeforeFusion ++;
	}

	void increment_unpairIncompleteBeforeFusion()
	{
		supportReadNum_unpairIncompleteBeforeFusion ++;
	}	

	void increment_pairIncompleteBeforeFusion()
	{
		supportReadNum_pairIncompleteBeforeFusion ++;
	}
};

class FusionRegionVec_Info
{
private:

	vector<FusionRegion_Info*> fusionRegionInfoVec;
	RegionPairIndexMap fusionRegion_genomicRegionIndexPairMap;

public:
	FusionRegionVec_Info()
	{}

	void output_fusionRegionVecInfo(ofstream& candidateFusionRegion_ofs,
		GenomicRegionVec_Info* genomicRegionVecInfo,
		Index_Info* indexInfo, int validFusionRegionPair_supportNumMin)
	{
		// candidateFusionRegion_ofs << "FusionRegionPairIndex\t";
		// candidateFusionRegion_ofs << "GenomicRegionIndex_Nor\t";
		// candidateFusionRegion_ofs << "GenomicRegionIndex_Rcm\t";
		// candidateFusionRegion_ofs<<  "chrName_Nor\tstartPos_Nor\tendPos_Nor\t";
		// candidateFusionRegion_ofs<<  "chrName_Rcm\tstartPos_Rcm\tendPos_Rcm\t";
		// candidateFusionRegion_ofs << "unpairCompleteReadNum\t";
		// candidateFusionRegion_ofs << "unpairIncompleteReadNum\t";
		// candidateFusionRegion_ofs << "pairIncompleteReadNum" << endl;

		int fusionRegionInfoVecSize = fusionRegionInfoVec.size();
		int tmpCandiFusionRegionPairIndex = 1;
		for(int tmp = 0; tmp < fusionRegionInfoVecSize; tmp++)
		{
			int tmpFusionRegionInfo_supportNum = fusionRegionInfoVec[tmp]->returnTotalSupportNum();
			if(tmpFusionRegionInfo_supportNum >= validFusionRegionPair_supportNumMin)
			{	
				candidateFusionRegion_ofs //<< tmpCandiFusionRegionPairIndex << "\t"
					<< fusionRegionInfoVec[tmp]->returnFuionRegionInfoStr(
						genomicRegionVecInfo, indexInfo) << endl;
				tmpCandiFusionRegionPairIndex ++;
			}
		}
	}

	void insertNewFusionRegionInfo2RegionFusionInfoHash(
		int genomicRegionIndex_Nor, int genomicRegionIndex_Rcm,
		int tmpIndex_inFusionRegionInfoVec)
	{
		RegionPairIndexMap::iterator tmpRegionPairIndexMapIter
			= fusionRegion_genomicRegionIndexPairMap.find(genomicRegionIndex_Nor);
		if(tmpRegionPairIndexMapIter == fusionRegion_genomicRegionIndexPairMap.end())
		{
			RegionIndexMap tmpNewRegionIndexMap;
			tmpNewRegionIndexMap.insert(pair<int,int>(
				genomicRegionIndex_Rcm, tmpIndex_inFusionRegionInfoVec));
			fusionRegion_genomicRegionIndexPairMap.insert(
				pair<int, RegionIndexMap>(
					genomicRegionIndex_Nor, tmpNewRegionIndexMap));
		}
		else
		{
			(tmpRegionPairIndexMapIter->second).insert(pair<int,int>(
				genomicRegionIndex_Rcm, tmpIndex_inFusionRegionInfoVec));
		}
	}

	int searchIndexInExistingFusionRegionHash(
		int genomicRegionIndex_Nor, int genomicRegionIndex_Rcm)
	{
		RegionPairIndexMap::iterator tmpRegionPairIndexMapIter
			= fusionRegion_genomicRegionIndexPairMap.find(genomicRegionIndex_Nor);
		if(tmpRegionPairIndexMapIter 
			!= fusionRegion_genomicRegionIndexPairMap.end())
		{
			RegionIndexMap::iterator tmpRegionIndexMapIter
				= (tmpRegionPairIndexMapIter->second).find(
					genomicRegionIndex_Rcm);
			if(tmpRegionIndexMapIter != (tmpRegionPairIndexMapIter->second).end())
				return (tmpRegionIndexMapIter->second);
			else
				return -1;
		}
		else
			return -1;
	}

	void checkAndGenerateCandidateFusionRegionInfo_fromUnpairCompleteReadPair(
		PE_Read_Alignment_Info* peAlignInfo, 
		GenomicRegionVec_Info* genomicRegionVecInfo, Index_Info* indexInfo)
	{
		//cout << "checkAndGenerateCandidateFusionRegionInfo_fromUnpairCompleteReadPair starts ..." << endl;
		vector< pair<int,int> > tmpGenomicRegionIndexPairVec;
		genomicRegionVecInfo->generateFusionRegionIndexPairVec_unpairComplete(
			peAlignInfo, tmpGenomicRegionIndexPairVec, indexInfo);
		//cout << "tmpGenomicRegionIndexPairVec.size(): " << tmpGenomicRegionIndexPairVec.size() << endl;
		for(int tmpIndex_regionIndexPairVec = 0; 
			tmpIndex_regionIndexPairVec < tmpGenomicRegionIndexPairVec.size();
			tmpIndex_regionIndexPairVec ++)
		{
			int tmpIndex_genomicRegion_Nor 
				= tmpGenomicRegionIndexPairVec[tmpIndex_regionIndexPairVec].first;
			int tmpIndex_genomicRegion_Rcm
				= tmpGenomicRegionIndexPairVec[tmpIndex_regionIndexPairVec].second;
			int tmpIndex_inFusionRegionInfoVec = this->searchIndexInExistingFusionRegionHash(
				tmpIndex_genomicRegion_Nor, tmpIndex_genomicRegion_Rcm);
			//cout << "tmpIndex_genomicRegion_Nor: " << tmpIndex_genomicRegion_Nor << endl;
			//cout << "tmpIndex_genomicRegion_Rcm: " << tmpIndex_genomicRegion_Rcm << endl;
			//cout << "tmpIndex_inFusionRegionInfoVec: " << tmpIndex_inFusionRegionInfoVec << endl;
			if(tmpIndex_inFusionRegionInfoVec < 0)
			{
				this->initiateAndAddFusionRegionInfo_withGenomicRegionIndex_pair(
					tmpIndex_genomicRegion_Nor, tmpIndex_genomicRegion_Rcm);
			}
			else // existed, update suppportNum
			{	
				this->increment_unpairCompleteBeforeFusion_fusionRegionInfo(
					tmpIndex_inFusionRegionInfoVec);
			}
		}
	}

	void initiateAndAddFusionRegionInfo_withGenomicRegionIndex_pair(
		int tmpGenomicRegionIndex_Nor_fusion,
		int tmpGenomicRegionIndex_Rcm_fusion)
	{
		// generate fusionRegionInfo
		FusionRegion_Info* tmpFusionRegionInfo
			= new FusionRegion_Info();
		tmpFusionRegionInfo->initiate_withGenomicRegionIndex_pair(
			tmpGenomicRegionIndex_Nor_fusion,
			tmpGenomicRegionIndex_Rcm_fusion);
		tmpFusionRegionInfo->increment_unpairCompleteBeforeFusion();

		// push back to fusionRegionInfoVec
		int tmpIndex_inFusionRegionInfoVec = fusionRegionInfoVec.size();
		fusionRegionInfoVec.push_back(tmpFusionRegionInfo);

		// insert 2 fusionRegionInfoHash
		this->insertNewFusionRegionInfo2RegionFusionInfoHash(
			tmpGenomicRegionIndex_Nor_fusion,
			tmpGenomicRegionIndex_Rcm_fusion,
			tmpIndex_inFusionRegionInfoVec);
	}

	/*
	void initiateAndAddFusionRegionInfo_withPeAlignInfo(
		PE_Read_Alignment_Info* peAlignInfo)
	{
		FusionRegion_Info* tmpFusionRegionInfo
			 = new FusionRegion_Info;
		tmpFusionRegionInfo->initiate_withPeAlignInfo(peAlignInfo);
		fusionRegionInfoVec.push_back(tmpFusionRegionInfo);
	}*/

	void increment_unpairCompleteBeforeFusion_fusionRegionInfo(int index)
	{
		fusionRegionInfoVec[index]->increment_unpairCompleteBeforeFusion();
	}

	void increment_unpairIncompleteBeforeFusion_fusionRegionInfo(int index)
	{
		fusionRegionInfoVec[index]->increment_unpairIncompleteBeforeFusion();
	}	

	void increment_pairIncompleteBeforeFusion_fusionRegionInfo(int index)
	{
		fusionRegionInfoVec[index]->increment_pairIncompleteBeforeFusion();
	}

	void memoryFree()
	{
		for(int tmp = 0; tmp < fusionRegionInfoVec.size(); tmp++)
		{
			delete fusionRegionInfoVec[tmp];
		}
	}
};
#endif