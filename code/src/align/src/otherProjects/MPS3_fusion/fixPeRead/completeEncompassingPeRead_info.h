// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef COMPLETEENCOMPASSINGPEREAD_INFO_H
#define COMPLETEENCOMPASSINGPEREAD_INFO_H

#include "../../../general/index_info.h"
#include "fusionRegion_info.h"

using namespace std;

class CompleteEncompassionPeRead_Info
{
public:
	vector< pair<int,int> > innerMapPosPairVec;

private:	
	CompleteEncompassionPeRead_Info()
	{}

	void addNewInnerMapPosPair(int chrMapPos_end_Nor, int chrMapPos_Rcm)
	{
		innerMapPosPairVec.push_back(pair<int,int>(chrMapPos_end_Nor,
			chrMapPos_Rcm));
	}

	void generateInnerMapPosPairVec(PE_Read_Alignment_Info* peAlignInfo,
		fusionRegion_info* fusionRegionInfo,
		vector< pair<int,int> >& tmpInnerMapPosPairVec, Index_Info* indexInfo)
	{
		int fusion_genomicRegionIndex_Nor
			= fusionRegionInfo->returnGenomicRegionIndex_Nor();
		int fusion_genomicRegionIndex_Rcm
			= fusionRegionInfo->returnGenomicRegionIndex_Rcm();

		vector<int> chrMapPos_end_Nor1;
		vector<int> chrMapPos_Rcm2;

		vector<int> chrMapPos_end_Nor2;
		vector<int> chrMapPos_Rcm1;

		int tmpAlignInfoNum_Nor1 = peAlignInfo->getAlignInfoVecSize(1);
		int tmpAlignInfoNum_Rcm1 = peAlignInfo->getAlignInfoVecSize(2);
		int tmpAlignInfoNum_Nor2 = peAlignInfo->getAlignInfoVecSize(3);
		int tmpAlignInfoNum_Rcm2 = peAlignInfo->getAlignInfoVecSize(4);
		for(int tmpAlignInfoIndex = 0; tmpAlignInfoIndex < tmpAlignInfoNum_Nor1; 
			tmpAlignInfoIndex ++)
		{
			string tmpAlignInfo_chrNameStr 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(1,
					tmpAlignInfoIndex))->returnAlignChromName();
			int tmpAlignIfno_chrNameInt = indexInfo->convertStringToInt(
				tmpAlignInfo_chrNameStr);
			int tmpAlignInfo_chrMapPos 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(1,
					tmpAlignInfoIndex))->returnAlignChromPos();
			int tmpAlignInfo_genomicRegionIndex
				= this->returnGenomicRegionIndex(tmpAlignIfno_chrNameInt,
					tmpAlignInfo_chrMapPos);
			if(tmpAlignInfo_genomicRegionIndex == fusion_genomicRegionIndex_Nor)
			{
				int tmpAlignInfo_chrMapPos_end 
					= (peAlignInfo->returnAlignInfoInPeAlignInfo(1,
						tmpAlignInfoIndex))->returnEndMatchedPosInChr();
				chrMapPos_end_Nor1.push_back(tmpAlignInfo_chrMapPos_end);;
			}
		}		

		for(int tmpAlignInfoIndex = 0; tmpAlignInfoIndex < tmpAlignInfoNum_Rcm1; 
			tmpAlignInfoIndex ++)
		{
			string tmpAlignInfo_chrNameStr 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(2,
					tmpAlignInfoIndex))->returnAlignChromName();
			int tmpAlignIfno_chrNameInt = indexInfo->convertStringToInt(
				tmpAlignInfo_chrNameStr);
			int tmpAlignInfo_chrMapPos 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(2,
					tmpAlignInfoIndex))->returnAlignChromPos();
			int tmpAlignInfo_genomicRegionIndex
				= this->returnGenomicRegionIndex(tmpAlignIfno_chrNameInt,
					tmpAlignInfo_chrMapPos);
			if(tmpAlignInfo_genomicRegionIndex == fusion_genomicRegionIndex_Rcm)
			{
				chrMapPos_Rcm1.push_back(tmpAlignInfo_chrMapPos);;
			}
		}		

		for(int tmpAlignInfoIndex = 0; tmpAlignInfoIndex < tmpAlignInfoNum_Nor2; 
			tmpAlignInfoIndex ++)
		{
			string tmpAlignInfo_chrNameStr 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(3,
					tmpAlignInfoIndex))->returnAlignChromName();
			int tmpAlignIfno_chrNameInt = indexInfo->convertStringToInt(
				tmpAlignInfo_chrNameStr);
			int tmpAlignInfo_chrMapPos 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(3,
					tmpAlignInfoIndex))->returnAlignChromPos();
			int tmpAlignInfo_genomicRegionIndex
				= this->returnGenomicRegionIndex(tmpAlignIfno_chrNameInt,
					tmpAlignInfo_chrMapPos);
			if(tmpAlignInfo_genomicRegionIndex == fusion_genomicRegionIndex_Nor)
			{
				int tmpAlignInfo_chrMapPos_end 
					= (peAlignInfo->returnAlignInfoInPeAlignInfo(3,
						tmpAlignInfoIndex))->returnEndMatchedPosInChr();
				chrMapPos_end_Nor2.push_back(tmpAlignInfo_chrMapPos_end);;
			}
		}

		for(int tmpAlignInfoIndex = 0; tmpAlignInfoIndex < tmpAlignInfoNum_Rcm2; 
			tmpAlignInfoIndex ++)
		{
			string tmpAlignInfo_chrNameStr 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(4,
					tmpAlignInfoIndex))->returnAlignChromName();
			int tmpAlignIfno_chrNameInt = indexInfo->convertStringToInt(
				tmpAlignInfo_chrNameStr);
			int tmpAlignInfo_chrMapPos 
				= (peAlignInfo->returnAlignInfoInPeAlignInfo(4,
					tmpAlignInfoIndex))->returnAlignChromPos();
			int tmpAlignInfo_genomicRegionIndex
				= this->returnGenomicRegionIndex(tmpAlignIfno_chrNameInt,
					tmpAlignInfo_chrMapPos);
			if(tmpAlignInfo_genomicRegionIndex == fusion_genomicRegionIndex_Rcm)
			{
				chrMapPos_Rcm2.push_back(tmpAlignInfo_chrMapPos);;
			}
		}	

		for(int tmp = 0; tmp < chrMapPos_end_Nor1.size(); tmp++)
		{
			int tmpInnerPos_Nor1 = chrMapPos_end_Nor1[tmp];
			for(int tmp2 = 0; tmp2 < chrMapPos_Rcm2.size(); tmp2++)
			{
				int tmpInnerPos_Rcm2 = chrMapPos_Rcm2[tmp2];
				tmpInnerMapPosPairVec.push_back(pair<int,int>(
					tmpInnerPos_Nor1, tmpInnerPos_Rcm2));
			}
		}

		for(int tmp = 0; tmp < chrMapPos_end_Nor2.size(); tmp++)
		{
			int tmpInnerPos_Nor2 = chrMapPos_end_Nor2[tmp];
			for(int tmp2 = 0; tmp2 < chrMapPos_Rcm1.size(); tmp2++)
			{
				int tmpInnerPos_Rcm1 = chrMapPos_Rcm1[tmp2];
				tmpInnerMapPosPairVec.push_back(pair<int,int>(
					tmpInnerPos_Nor2, tmpInnerPos_Rcm1));				
			}
		}
	}

	void generateAndAddNewInnerMapPosPairVec(PE_Read_Alignment_Info* peAlignInfo,
		fusionRegion_info* fusionRegionInfo, Index_Info* indexInfo)
	{
		vector< pair<int,int> > generatedInnerMapPosPairVec;
		this->generateAndAddNewInnerMapPosPairVec(peAlignInfo,
			fusionRegionInfo, generatedInnerMapPosPairVec ,indexInfo)
		for(int tmp = 0; tmp < generatedInnerMapPosPairVec.size(); tmp++)
		{
			int tmpInnerMapPos_Nor = generatedInnerMapPosPairVec[tmp].first;
			int tmpInnerMapPos_Rcm = generatedInnerMapPosPairVec[tmp].second;
			this->addNewInnerMapPosPair(tmpInnerMapPos_Nor, tmpInnerMapPos_Rcm);
		}
	}
};

#endif