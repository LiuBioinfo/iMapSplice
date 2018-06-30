// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FUSIONALIGN_INFO_H
#define FUSIONALIGN_INFO_H

#include "align_info.h"

using namespace std;

class PE_Fusion_Info
{
public:
	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Nor1_head_nor; // <index_Nor1, alignmentInfo_fusion>
	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Nor1_head_rcm;

	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Nor2_head_nor;
	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Nor2_head_rcm;

	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Rcm1_tail_nor;	
	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Rcm1_tail_rcm;

	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Rcm2_tail_nor;
	vector< pair <int, Alignment_Info*> > fusionAlignInfo_Rcm2_tail_rcm;

	vector<bool> fusionAlignInfo_Nor1_head_nor_bool; // <index_Nor1, alignmentInfo_fusion>
	vector<bool> fusionAlignInfo_Nor2_head_nor_bool;
	vector<bool> fusionAlignInfo_Rcm1_tail_nor_bool;
	vector<bool> fusionAlignInfo_Rcm2_tail_nor_bool;
	vector<bool> fusionAlignInfo_Nor1_head_rcm_bool;
	vector<bool> fusionAlignInfo_Nor2_head_rcm_bool;
	vector<bool> fusionAlignInfo_Rcm1_tail_rcm_bool;
	vector<bool> fusionAlignInfo_Rcm2_tail_rcm_bool;

	vector< pair< int, pair<int, int> > > fusionAlignInfo_Nor1_head_nor_refine; // <index_in_Nor1_head_nor, <firstMatchLength, secondMatchLength> >
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Nor2_head_nor_refine;
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Rcm1_tail_nor_refine;
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Rcm2_tail_nor_refine;
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Nor1_head_rcm_refine;
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Nor2_head_rcm_refine;
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Rcm1_tail_rcm_refine;
	vector< pair< int, pair<int, int> > > fusionAlignInfo_Rcm2_tail_rcm_refine;


	int buffer;
	int LengthOfSeqPerMismatchAllowed;

	PE_Fusion_Info()
	{
		buffer = 3;
		LengthOfSeqPerMismatchAllowed = 10;
	}

	void pushBackFusionTargetRegionAlignInfo(bool Nor1Rcm2OrNor2Rcm1,
			bool unfixedHeadOrTail, bool fusionNorOrRcm,
			int index_Nor1OrNor2OrRcm1OrRcm2, Alignment_Info* alignInfo
			//const string& strandStr
			)
	{
		if(Nor1Rcm2OrNor2Rcm1)
		{
			if(unfixedHeadOrTail)
			{
				if(fusionNorOrRcm)
				{
					fusionAlignInfo_Nor1_head_nor.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
					//fusionAlignInfo_Nor1_head_nor_bool.push_back(false);
					//fusionAlignInfo_Nor1_head_nor_refine.push_back(pair<int,int>(1000,1000))
				}
				else
				{
					fusionAlignInfo_Nor1_head_rcm.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
				}
			}
			else
			{
				if(fusionNorOrRcm)
				{
					fusionAlignInfo_Rcm2_tail_nor.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
				}
				else
				{
					fusionAlignInfo_Rcm2_tail_rcm.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
				}
			}
		}
		else
		{
			if(unfixedHeadOrTail)
			{
				if(fusionNorOrRcm)
				{
					fusionAlignInfo_Nor2_head_nor.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
				}
				else
				{
					fusionAlignInfo_Nor2_head_rcm.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));			
				}
			}
			else
			{
				if(fusionNorOrRcm)
				{
					fusionAlignInfo_Rcm1_tail_nor.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
				}
				else
				{
					fusionAlignInfo_Rcm1_tail_rcm.push_back(pair< int, Alignment_Info* >(
						index_Nor1OrNor2OrRcm1OrRcm2, alignInfo));
				}
			}
		}
	}

	void pushBackFusionAlignInfoRefineResults(bool Nor1Rcm2OrNor2Rcm1,
			bool unfixedHeadOrTail, bool fusionNorOrRcm, 
			int firstMatchLength, int secondMatchLength, int index_num)
	{
		if(Nor1Rcm2OrNor2Rcm1)
		{
			if(unfixedHeadOrTail)
			{
				if(fusionNorOrRcm)
				{
					//fusionAlignInfo_Nor1_head_nor_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Nor1_head_nor_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Nor1_head_nor_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
				}
				else
				{
					//fusionAlignInfo_Nor1_head_rcm_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Nor1_head_rcm_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Nor1_head_rcm_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
				}
			}
			else
			{
				if(fusionNorOrRcm)
				{
					//fusionAlignInfo_Rcm2_tail_nor_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Rcm2_tail_nor_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Rcm2_tail_nor_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
					//fusionAlignInfo_Nor1_head_nor_refine.push_back(pair< int, int >(
					//	firstMatchLength, secondMatchLength));
				}
				else
				{
					//fusionAlignInfo_Rcm2_tail_rcm_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Rcm2_tail_rcm_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Rcm2_tail_rcm_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
				}
			}
		}
		else
		{
			if(unfixedHeadOrTail)
			{
				if(fusionNorOrRcm)
				{
					//fusionAlignInfo_Nor1_head_nor_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Nor1_head_nor_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Nor2_head_nor_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
				}
				else
				{
					//fusionAlignInfo_Nor1_head_rcm_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Nor1_head_rcm_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Nor2_head_rcm_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
				}
			}
			else
			{
				if(fusionNorOrRcm)
				{
					//fusionAlignInfo_Rcm2_tail_nor_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Rcm2_tail_nor_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Rcm1_tail_nor_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
					//fusionAlignInfo_Nor1_head_nor_refine.push_back(pair< int, int >(
					//	firstMatchLength, secondMatchLength));
				}
				else
				{
					//fusionAlignInfo_Rcm2_tail_rcm_refine[index_num].first = firstMatchLength;
					//fusionAlignInfo_Rcm2_tail_rcm_refine[index_num].second = secondMatchLength;
					fusionAlignInfo_Rcm1_tail_rcm_refine.push_back(
						pair< int, pair< int, int > > ( index_num, 
						pair< int, int > (firstMatchLength, secondMatchLength) ));
				}
			}
		}		
	}

	int returnFusionAlignmentSize(bool Nor1Rcm2OrNor2Rcm1,
			bool unfixedHeadOrTail, bool fusionNorOrRcm)
	{
		if(Nor1Rcm2OrNor2Rcm1)
		{
			if(unfixedHeadOrTail)
			{
				if(fusionNorOrRcm)
				{
					return fusionAlignInfo_Nor1_head_nor.size();
				}
				else
				{
					return fusionAlignInfo_Nor1_head_rcm.size();
				}
			}
			else
			{
				if(fusionNorOrRcm)
				{
					return fusionAlignInfo_Nor1_head_rcm.size();
				}
				else
				{
					return fusionAlignInfo_Nor1_head_rcm.size();
				}
			}
		}
		else
		{
			if(unfixedHeadOrTail)
			{
				if(fusionNorOrRcm)
				{
					return fusionAlignInfo_Nor1_head_rcm.size();
				}
				else
				{
					return fusionAlignInfo_Nor1_head_rcm.size();		
				}
			}
			else
			{
				if(fusionNorOrRcm)
				{
					return fusionAlignInfo_Nor1_head_rcm.size();
				}
				else
				{
					return fusionAlignInfo_Nor1_head_rcm.size();
				}
			}
		}
	}

	void refineFusionAlignment_all(PE_Read_Alignment_Info* peAlignInfo, PE_Read_Info* peReadInfo, Index_Info* indexInfo)
	{
		for(int tmp = 0; tmp < fusionAlignInfo_Nor1_head_nor.size(); tmp ++)
		{
			int index_Nor1 = fusionAlignInfo_Nor1_head_nor[tmp].first;
			this->refineFusionAlignemnt(peReadInfo, true, true, true, 
				fusionAlignInfo_Nor1_head_nor[tmp].second, 
				peAlignInfo->returnAlignInfoInPeAlignInfo(true, true, index_Nor1),
				indexInfo, tmp);
		}
		for(int tmp = 0; tmp < fusionAlignInfo_Nor2_head_nor.size(); tmp ++)
		{
			int index_Nor2 = fusionAlignInfo_Nor2_head_nor[tmp].first;
			this->refineFusionAlignemnt(peReadInfo, true, true, true, 
				fusionAlignInfo_Nor2_head_nor[tmp].second,
				peAlignInfo->returnAlignInfoInPeAlignInfo(false, true, index_Nor2),
				 indexInfo, tmp);
		}
		for(int tmp = 0; tmp < fusionAlignInfo_Nor1_head_rcm.size(); tmp ++)
		{

		}
		for(int tmp = 0; tmp < fusionAlignInfo_Nor2_head_rcm.size(); tmp ++)
		{

		}

		for(int tmp = 0; tmp < fusionAlignInfo_Rcm1_tail_nor.size(); tmp ++)
		{

		}
		for(int tmp = 0; tmp < fusionAlignInfo_Rcm2_tail_nor.size(); tmp ++)
		{

		}
		for(int tmp = 0; tmp < fusionAlignInfo_Rcm1_tail_rcm.size(); tmp ++)
		{

		}
		for(int tmp = 0; tmp < fusionAlignInfo_Rcm2_tail_rcm.size(); tmp ++)
		{

		}
	}

	void refineFusionAlignemnt(
		PE_Read_Info* peReadInfo,
		bool Nor1Rcm2OrNor2Rcm1, bool fusionTargetInHeadOrTail, bool NorOrRcmBool,
		Alignment_Info* alignInfo_front, Alignment_Info* alignInfo_rear, Index_Info* indexInfo, int index_num)
	{
		string alignInfo_front_chrom = alignInfo_front->alignChromName;
		int alignInfo_front_chrom_int = indexInfo->convertStringToInt(alignInfo_front_chrom); 
		int alignInfo_front_pos_end = alignInfo_front->getEndMatchedPosInChr(); 

		string alignInfo_rear_chrom = alignInfo_rear->alignChromName;
		int alignInfo_rear_chrom_int = indexInfo->convertStringToInt(alignInfo_rear_chrom);
		int alignInfo_rear_pos_start = alignInfo_rear->alignChromPos;

		int toFixReadLength;

		if(fusionTargetInHeadOrTail)
		{

			toFixReadLength = peReadInfo->returnReadLength(Nor1Rcm2OrNor2Rcm1);

			int targetRegionLength = alignInfo_rear->unfixedHeadLength();

			int firstMatchLength, secondMatchLength;

			if(alignInfo_front->unfixedTailExistsBool())
			{
				int targetRegionUnfixedTailLength = alignInfo_front->unfixedTailLength();
				int tmpBuffer_left = buffer;
				if(targetRegionLength - targetRegionUnfixedTailLength - buffer < 2)
				{
					tmpBuffer_left = 0;
				}
				int tmpBuffer_right = buffer;
				if(targetRegionLength + tmpBuffer_right+2 > toFixReadLength )
				{
					tmpBuffer_right = 0;
				}

				int toFixReadSeqLength = targetRegionUnfixedTailLength + tmpBuffer_left + tmpBuffer_right;
				int toFixLocInRead = targetRegionLength - targetRegionUnfixedTailLength - tmpBuffer_left + 1;

				string readSeq_toFix;
				if(Nor1Rcm2OrNor2Rcm1)
				{
					readSeq_toFix = ((peReadInfo->readInfo_pe1).readSeq).substr(toFixLocInRead-1, toFixReadSeqLength);
				}
				else
				{
					readSeq_toFix = ((peReadInfo->readInfo_pe2).readSeq).substr(toFixLocInRead-1, toFixReadSeqLength);
				}

				int chromSubSeqLengthInProcess = toFixReadSeqLength + 2;
				string left_chrom_seq = indexInfo->chromStr[alignInfo_front_chrom_int].substr(
					alignInfo_front_pos_end + 1 - tmpBuffer_left - 1, chromSubSeqLengthInProcess);
				string right_chrom_seq = indexInfo->chromStr[alignInfo_rear_chrom_int].substr(
					alignInfo_rear_pos_start - 1 + tmpBuffer_right - chromSubSeqLengthInProcess, chromSubSeqLengthInProcess);


				size_t prefix_length = 0;
				size_t max_double_splice_mismatch = toFixReadSeqLength/LengthOfSeqPerMismatchAllowed + 1;
				size_t mismatch_bits = 0;
				size_t comb_bits = 0;
				//bool adjacent_segments = false;
				bool double_anchor_noncanonical = true;//false;//DO_NONCANONICAL; ////debug
				string flank_seq;
				GenomeScan* genome_scan = new GenomeScan;
				bool splice_fixed = (*genome_scan).Double_anchored_score(readSeq_toFix, 
					left_chrom_seq, right_chrom_seq, prefix_length, 
					max_double_splice_mismatch, comb_bits,
		 			//(!adjacent_segments) && 
					double_anchor_noncanonical, flank_seq, mismatch_bits);

				//cout << "prefix_length: " << prefix_length << endl;
				//cout << "splice_fixed: " << splice_fixed << endl;

				if(splice_fixed) 
				{
					firstMatchLength = prefix_length - tmpBuffer_left;
					secondMatchLength = targetRegionUnfixedTailLength - firstMatchLength;

				}
				else
				{

				}

			}
			else
			{
				firstMatchLength = 0;
				secondMatchLength = 0;
			}	
			cout << endl << "......... unfixedHeadLength: " << targetRegionLength << endl << "......... targetRegionUnfixedTailLength: " 
				<< firstMatchLength + secondMatchLength << endl; 	
			cout << "......... firstMatchLength: " << firstMatchLength << " -- secondMatchLength: " << secondMatchLength << endl;
			//if((firstMatchLength != 0)||(secondMatchLength != 0))
			//	exit(0);
			this->pushBackFusionAlignInfoRefineResults(
				Nor1Rcm2OrNor2Rcm1, fusionTargetInHeadOrTail, NorOrRcmBool,
				firstMatchLength, secondMatchLength, index_num);		
		}
		else
		{

		}
	}
};

#endif