// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXONEENDUNMAPPED_H
#define FIXONEENDUNMAPPED_H

#include <string>
#include <string.h>
//#include "splice_info.h"

using namespace std;

class FixOneEndUnmappedInfo
{
public:
	FixOneEndUnmappedInfo()
	{}

	void fixOneEndUnmapped(PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, 
		Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax,
		bool checkQualSeqForReadSegSeq
		)
	{
		//cout << " ... start of fixOneEndUnmapped ... " << endl;

		int readLength;// = 100; // to debug
		bool End1OrEnd2; 
		bool NorOrRcm;

		vector<int> alignmentInfoVecSize_ori;
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnNorAlignmentInfo_PE_1_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnRcmAlignmentInfo_PE_1_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnNorAlignmentInfo_PE_2_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnRcmAlignmentInfo_PE_2_size());

		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType ++)
		{
			//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;

			End1OrEnd2 = peReadInfo.checkEnd1OrEnd2WithAlignInfoTypeNo(tmpAlignInfoType);
			NorOrRcm = peReadInfo.checkNorOrRcmWithAlignInfoTypeNo(tmpAlignInfoType);
			readLength = peReadInfo.checkReadLengthWithAlignInfoTypeNo(5-tmpAlignInfoType);

			//int peAlignInfoVectorSize = (peAlignInfo->getAlignInfoVecSize(tmpAlignInfoType));
			int peAlignInfoVectorSize = alignmentInfoVecSize_ori[tmpAlignInfoType-1];

			for(int tmpAlignmentNO = 0; 
					tmpAlignmentNO < peAlignInfoVectorSize; 
					tmpAlignmentNO ++)
			{
				//cout << "tmpAlignemntNO: " << tmpAlignmentNO << endl;
				//cout << "start to setUnmapEndInfo" << endl;
				Alignment_Info* tmpAlignInfo 
					= peAlignInfo->getAlignInfo(tmpAlignInfoType, tmpAlignmentNO);
				UnmapEnd_Info* unmapEndInfo = new UnmapEnd_Info();
				unmapEndInfo->setUnmapEndInfo(tmpAlignInfo, End1OrEnd2, NorOrRcm, indexInfo);
				//cout << "start to getIncompleteEndReadSeq" << endl;
				char* unmapEndReadChar =  const_cast<char*>(peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm).c_str());
				Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

				int secondLevelIndexNO = unmapEndInfo->returnSecondLevelIndexNum() - 1;							
				if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO+1))
				{
					delete seg2ndOriInfo;
					delete unmapEndInfo;
					continue;
				}
				//cout << "start to mapMainSecondLevel_compressedIndex" << endl;
				bool unmapEndMapBool;	
				unmapEndMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
					unmapEndReadChar, secondLevelSa[secondLevelIndexNO], 
					secondLevelLcpCompress[secondLevelIndexNO],
					secondLevelChildTab[secondLevelIndexNO],
					secondLevelChrom[secondLevelIndexNO], 
					secondLevelDetChild[secondLevelIndexNO],
					readLength, indexInfo);

				if(!unmapEndMapBool)
				{
					delete seg2ndOriInfo;
					delete unmapEndInfo;
					continue;								
				}
				//cout << "start to do segInfo" << endl;
				Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
					unmapEndInfo->returnMapPosIntervalStart(),
					unmapEndInfo->returnMapPosIntervalEnd(), 
					unmapEndInfo->returnChrPosStartIn2ndLevelIndex(),
					indexInfo, unmapEndInfo->returnChrNameStr());

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
				}
				//cout << "start to do pathInfo" << endl;
				//cout << "SegInfo:\n" << segInfo->segInfoStr(indexInfo) << endl;
				Path_Info* pathInfo = new Path_Info();
				pathInfo->getPossiPathFromSeg_fixOneEndUnmapped(segInfo, spliceJunctionDistanceMax);
				//cout << "PathInfo:\n" << pathInfo->possiPathStr() << endl;
				int pathValidNum = pathInfo->pathValidNumInt();
				if(pathValidNum > 10)
				{
					pathInfo->memoryFree();
					delete pathInfo;
					delete segInfo;
					delete seg2ndOriInfo;
					delete unmapEndInfo;								
					continue;
				}
				//cout << "start to do gapInfo" << endl;
				Gap_Info* gapInfo = new Gap_Info();
				gapInfo->fixGapInPath_phase2(pathInfo, segInfo, 
					indexInfo, peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), 
					readLength, Do_extendHeadTail,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
					spliceJunctionDistanceMax);
				//cout << "start to pushBackPathInfo2PeAlignInfo" << endl;
				peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
				//cout << "end of fixOneEndUnmapped for this alignInfo" << endl;
				delete gapInfo;
				pathInfo->memoryFree();
				delete pathInfo;
				delete segInfo;
				delete seg2ndOriInfo;
				delete unmapEndInfo;
			}
		}
	}

	#ifdef PERSONALIZED_CHR_SEQ
	void fixOneEndUnmapped_includeSNPseqMap(PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, 
		Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax,
		bool checkQualSeqForReadSegSeq,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp, int SNPlocInSyntheticSNPseq	
		)
	{
		//cout << " ... start of fixOneEndUnmapped ... " << endl;

		int readLength;// = 100; // to debug
		bool End1OrEnd2; 
		bool NorOrRcm;

		vector<int> alignmentInfoVecSize_ori;
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnNorAlignmentInfo_PE_1_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnRcmAlignmentInfo_PE_1_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnNorAlignmentInfo_PE_2_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnRcmAlignmentInfo_PE_2_size());

		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType ++)
		{
			//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;

			End1OrEnd2 = peReadInfo.checkEnd1OrEnd2WithAlignInfoTypeNo(tmpAlignInfoType);
			NorOrRcm = peReadInfo.checkNorOrRcmWithAlignInfoTypeNo(tmpAlignInfoType);
			readLength = peReadInfo.checkReadLengthWithAlignInfoTypeNo(5-tmpAlignInfoType);

			//int peAlignInfoVectorSize = (peAlignInfo->getAlignInfoVecSize(tmpAlignInfoType));
			int peAlignInfoVectorSize = alignmentInfoVecSize_ori[tmpAlignInfoType-1];

			for(int tmpAlignmentNO = 0; 
					tmpAlignmentNO < peAlignInfoVectorSize; 
					tmpAlignmentNO ++)
			{
				//cout << "tmpAlignemntNO: " << tmpAlignmentNO << endl;
				//cout << "start to setUnmapEndInfo" << endl;
				Alignment_Info* tmpAlignInfo 
					= peAlignInfo->getAlignInfo(tmpAlignInfoType, tmpAlignmentNO);
				UnmapEnd_Info* unmapEndInfo = new UnmapEnd_Info();
				unmapEndInfo->setUnmapEndInfo(tmpAlignInfo, End1OrEnd2, NorOrRcm, indexInfo);
				//cout << "start to getIncompleteEndReadSeq" << endl;
				char* unmapEndReadChar =  const_cast<char*>(peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm).c_str());
				Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

				int secondLevelIndexNO = unmapEndInfo->returnSecondLevelIndexNum() - 1;							
				if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO+1))
				{
					delete seg2ndOriInfo;
					delete unmapEndInfo;
					continue;
				}
				//cout << "start to mapMainSecondLevel_compressedIndex" << endl;
				bool unmapEndMapBool;	
				unmapEndMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
					unmapEndReadChar, secondLevelSa[secondLevelIndexNO], 
					secondLevelLcpCompress[secondLevelIndexNO],
					secondLevelChildTab[secondLevelIndexNO],
					secondLevelChrom[secondLevelIndexNO], 
					secondLevelDetChild[secondLevelIndexNO],
					readLength, indexInfo);

				if(!unmapEndMapBool)
				{
					delete seg2ndOriInfo;
					delete unmapEndInfo;
					continue;								
				}
				//cout << "start to do segInfo" << endl;
				Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
					unmapEndInfo->returnMapPosIntervalStart(),
					unmapEndInfo->returnMapPosIntervalEnd(), 
					unmapEndInfo->returnChrPosStartIn2ndLevelIndex(),
					indexInfo, unmapEndInfo->returnChrNameStr());

				#ifdef PERSONALIZED_CHR_SEQ
				Seg_Info segInfo_SNPseq;
				int tmpIncompleteSeqLength = (peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm)).length();
				bool tmpMapWithSnpSeqIndexBool = segInfo_SNPseq.greedyMapWithoutPreIndexArray(unmapEndReadChar, 
					sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
					verifyChild_snp, tmpIncompleteSeqLength, indexInfo_snp);
				if(tmpMapWithSnpSeqIndexBool)
					segInfo->update_includeSNPseqMapSegInfo(segInfo_SNPseq, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);
				
				segInfo->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, indexInfo_snp,
					peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), indexInfo, SNPlocInSyntheticSNPseq);
				#endif

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
				}
				//cout << "start to do pathInfo" << endl;
				//cout << "SegInfo:\n" << segInfo->segInfoStr(indexInfo) << endl;
				Path_Info* pathInfo = new Path_Info();
				pathInfo->getPossiPathFromSeg_fixOneEndUnmapped(segInfo, spliceJunctionDistanceMax);
				//cout << "PathInfo:\n" << pathInfo->possiPathStr() << endl;
				int pathValidNum = pathInfo->pathValidNumInt();
				if(pathValidNum > 10)
				{
					pathInfo->memoryFree();
					delete pathInfo;
					delete segInfo;
					delete seg2ndOriInfo;
					delete unmapEndInfo;								
					continue;
				}
				//cout << "start to do gapInfo" << endl;
				Gap_Info* gapInfo = new Gap_Info();
				gapInfo->fixGapInPath_phase2(pathInfo, segInfo, 
					indexInfo, peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), 
					readLength, Do_extendHeadTail,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
					spliceJunctionDistanceMax);
				//cout << "start to pushBackPathInfo2PeAlignInfo" << endl;
				peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
				//cout << "end of fixOneEndUnmapped for this alignInfo" << endl;
				delete gapInfo;
				pathInfo->memoryFree();
				delete pathInfo;
				delete segInfo;
				delete seg2ndOriInfo;
				delete unmapEndInfo;
			}
		}
	}
	#endif
	#ifdef VARY_SNP_MER
	void fixOneEndUnmapped_includeSNPseqMap_varySNPmer(PE_Read_Info& peReadInfo, PE_Read_Alignment_Info* peAlignInfo,
		vector<char*>& secondLevelChrom,
		vector<unsigned int*>& secondLevelSa,
		vector<BYTE*>& secondLevelLcpCompress,
		vector<unsigned int*>& secondLevelChildTab,
		vector<BYTE*>& secondLevelDetChild,
		Index_Info* indexInfo, bool Do_extendHeadTail,
		bool annotation_provided_bool, bool Do_annotation_only_bool, 
		Annotation_Info* annotationInfo,
		int spliceJunctionDistanceMax,
		bool checkQualSeqForReadSegSeq,
		unsigned int* sa_snp, BYTE* lcpCompress_snp, unsigned int* childTab_snp, char* chrom_snp,
		BYTE* verifyChild_snp, Index_Info* indexInfo_snp//, int SNPlocInSyntheticSNPseq	
		)
	{
		//cout << " ... start of fixOneEndUnmapped ... " << endl;

		int readLength;// = 100; // to debug
		bool End1OrEnd2; 
		bool NorOrRcm;

		vector<int> alignmentInfoVecSize_ori;
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnNorAlignmentInfo_PE_1_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnRcmAlignmentInfo_PE_1_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnNorAlignmentInfo_PE_2_size());
		alignmentInfoVecSize_ori.push_back(peAlignInfo->returnRcmAlignmentInfo_PE_2_size());

		for(int tmpAlignInfoType = 1; tmpAlignInfoType <= 4; tmpAlignInfoType ++)
		{
			//cout << "tmpAlignInfoType: " << tmpAlignInfoType << endl;

			End1OrEnd2 = peReadInfo.checkEnd1OrEnd2WithAlignInfoTypeNo(tmpAlignInfoType);
			NorOrRcm = peReadInfo.checkNorOrRcmWithAlignInfoTypeNo(tmpAlignInfoType);
			readLength = peReadInfo.checkReadLengthWithAlignInfoTypeNo(5-tmpAlignInfoType);

			//int peAlignInfoVectorSize = (peAlignInfo->getAlignInfoVecSize(tmpAlignInfoType));
			int peAlignInfoVectorSize = alignmentInfoVecSize_ori[tmpAlignInfoType-1];

			for(int tmpAlignmentNO = 0; 
					tmpAlignmentNO < peAlignInfoVectorSize; 
					tmpAlignmentNO ++)
			{
				//cout << "tmpAlignemntNO: " << tmpAlignmentNO << endl;
				//cout << "start to setUnmapEndInfo" << endl;
				Alignment_Info* tmpAlignInfo 
					= peAlignInfo->getAlignInfo(tmpAlignInfoType, tmpAlignmentNO);
				UnmapEnd_Info* unmapEndInfo = new UnmapEnd_Info();
				unmapEndInfo->setUnmapEndInfo(tmpAlignInfo, End1OrEnd2, NorOrRcm, indexInfo);
				//cout << "start to getIncompleteEndReadSeq" << endl;
				char* unmapEndReadChar =  const_cast<char*>(peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm).c_str());
				Seg2ndOri_Info* seg2ndOriInfo = new Seg2ndOri_Info();

				int secondLevelIndexNO = unmapEndInfo->returnSecondLevelIndexNum() - 1;							
				if(indexInfo->returnInvalidSecondLevelIndexNOset_find_bool(secondLevelIndexNO+1))
				{
					delete seg2ndOriInfo;
					delete unmapEndInfo;
					continue;
				}
				//cout << "start to mapMainSecondLevel_compressedIndex" << endl;
				bool unmapEndMapBool;	
				unmapEndMapBool = seg2ndOriInfo->mapMainSecondLevel_compressedIndex(
					unmapEndReadChar, secondLevelSa[secondLevelIndexNO], 
					secondLevelLcpCompress[secondLevelIndexNO],
					secondLevelChildTab[secondLevelIndexNO],
					secondLevelChrom[secondLevelIndexNO], 
					secondLevelDetChild[secondLevelIndexNO],
					readLength, indexInfo);

				if(!unmapEndMapBool)
				{
					delete seg2ndOriInfo;
					delete unmapEndInfo;
					continue;								
				}
				//cout << "start to do segInfo" << endl;
				Seg_Info* segInfo = new Seg_Info(seg2ndOriInfo, 
					unmapEndInfo->returnMapPosIntervalStart(),
					unmapEndInfo->returnMapPosIntervalEnd(), 
					unmapEndInfo->returnChrPosStartIn2ndLevelIndex(),
					indexInfo, unmapEndInfo->returnChrNameStr());

				//#ifdef PERSONALIZED_CHR_SEQ
				Seg_Info segInfo_SNPseq;
				int tmpIncompleteSeqLength = (peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm)).length();
				bool tmpMapWithSnpSeqIndexBool = segInfo_SNPseq.greedyMapWithoutPreIndexArray(unmapEndReadChar, 
					sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, 
					verifyChild_snp, tmpIncompleteSeqLength, indexInfo_snp);
				
				#ifdef VARY_SNP_MER
				if(tmpMapWithSnpSeqIndexBool)
					segInfo->update_includeSNPseqMapSegInfo_varySNPmer(segInfo_SNPseq, indexInfo_snp, indexInfo);
				segInfo->update_targetMap2SNPseq_varySNPmer(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
					indexInfo_snp, peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), indexInfo);
				#else
				if(tmpMapWithSnpSeqIndexBool)
					segInfo->update_includeSNPseqMapSegInfo(segInfo_SNPseq, SNPlocInSyntheticSNPseq, indexInfo_snp, indexInfo);				
				segInfo->update_targetMap2SNPseq(sa_snp, lcpCompress_snp, childTab_snp, chrom_snp, verifyChild_snp, 
					indexInfo_snp, peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), indexInfo, SNPlocInSyntheticSNPseq);
				#endif

				// fix me: test: checkQualSeqForReadSegSeq
				if(checkQualSeqForReadSegSeq)
				{
					segInfo->filterLowQualitySeg(peReadInfo, tmpAlignInfoType);
				}
				//cout << "start to do pathInfo" << endl;
				//cout << "SegInfo:\n" << segInfo->segInfoStr(indexInfo) << endl;
				Path_Info* pathInfo = new Path_Info();
				pathInfo->getPossiPathFromSeg_fixOneEndUnmapped(segInfo, spliceJunctionDistanceMax);
				//cout << "PathInfo:\n" << pathInfo->possiPathStr() << endl;
				int pathValidNum = pathInfo->pathValidNumInt();
				if(pathValidNum > 10)
				{
					pathInfo->memoryFree();
					delete pathInfo;
					delete segInfo;
					delete seg2ndOriInfo;
					delete unmapEndInfo;								
					continue;
				}
				//cout << "start to do gapInfo" << endl;
				Gap_Info* gapInfo = new Gap_Info();
				gapInfo->fixGapInPath_phase2(pathInfo, segInfo, 
					indexInfo, peReadInfo.getIncompleteEndReadSeq(End1OrEnd2, NorOrRcm), 
					readLength, Do_extendHeadTail,
					annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
					spliceJunctionDistanceMax);
				//cout << "start to pushBackPathInfo2PeAlignInfo" << endl;
				peAlignInfo->pushBackPathInfo2PeAlignInfo(pathInfo, End1OrEnd2, NorOrRcm, indexInfo);
				//cout << "end of fixOneEndUnmapped for this alignInfo" << endl;
				delete gapInfo;
				pathInfo->memoryFree();
				delete pathInfo;
				delete segInfo;
				delete seg2ndOriInfo;
				delete unmapEndInfo;
			}
		}
	}
	#endif
};

#endif