// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SEG_INFO_H
#define SEG_INFO_H

#include <stdlib.h>
#include <stdio.h>
#include "local_seg_info.h"

using namespace std;

class Seg_Info_old
{
//public:
private: 
	int repeatRegion_index; // if -1 or < 0, is not mapped to repeatRegion 

	unsigned int segmentNum;
	//unsigned int norSegmentLength[SEGMENTNUM];
	vector<unsigned int> norSegmentLength;
	//unsigned int norSegmentLocInRead[SEGMENTNUM];
	vector<unsigned int> norSegmentLocInRead;
	//unsigned int norSegmentAlignNum[SEGMENTNUM];
	vector<unsigned int> norSegmentAlignNum;
	//unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC];
	vector< vector<unsigned int> > norSegmentAlignLoc;
	int longSegMinLength;

public:
	void assignLongSegMinLength(int newLongSegMinLength)
	{
		longSegMinLength = newLongSegMinLength;
	}
	//void assignSegmentAlignLoc(unsigned int mapPos, int segGroupNO, int segCandiNO)
	//{

	//}
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
	int returnSegmentAlignNum(int segGroupNO)
	{
		return (int)norSegmentAlignNum[segGroupNO];
	}
	unsigned int returnSegmentMapPos(int segGroupNO, int segCandiNO)
	{
		unsigned int mapPos //= *(norSegmentAlignLoc + segGroupNO*CANDALILOC + segCandiNO);
			= (norSegmentAlignLoc[segGroupNO])[segCandiNO];
		return mapPos;
	}
	unsigned int returnSegmentAlignLoc(int segGroupNO, int segCandiNO)
	{
		unsigned int mapPos //= *(norSegmentAlignLoc + segGroupNO*CANDALILOC + segCandiNO);
			= (norSegmentAlignLoc[segGroupNO])[segCandiNO];
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
		//norSegmentLength[segmentNum - 1] = midPartSegLength;
		norSegmentLength.push_back(midPartSegLength);
		//norSegmentLocInRead[segmentNum - 1] = midPartSegLocInRead;
		norSegmentLocInRead.push_back(midPartSegLocInRead);
		//norSegmentAlignNum[segmentNum - 1] = 1;
		norSegmentAlignNum.push_back(1);
		unsigned int tmpMidPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			midPartMapChrInt, midPartMapPosInChr);
		//*( norSegmentAlignLoc + (segmentNum - 1)*CANDALILOC ) = tmpMidPartMapPosInWholeGenome;
		vector<unsigned int> newMapPosVec;
		norSegmentAlignLoc.push_back(newMapPosVec);
		norSegmentAlignLoc[segmentNum - 1].push_back(tmpMidPartMapPosInWholeGenome);
	}

	void addMidPartSeg_incompleteTail(int midPartSegLength, int midPartSegLocInRead, 
		int midPartMapChrInt, int midPartMapPosInChr, Index_Info* indexInfo, Seg_Info* segInfo)
	{
		segmentNum = segInfo->segmentNum + 1;
		//norSegmentLength[0] = midPartSegLength;
		norSegmentLength.push_back(midPartSegLength);
		//norSegmentLocInRead[0] = 1;//midPartSegLocInRead;
		norSegmentLocInRead.push_back(1);
		//norSegmentAlignNum[0] = 1;
		norSegmentAlignNum.push_back(1);
		//int tmpUnfixedTailLocInRead = midPartSegLocInRead + midPartSegLength;

		unsigned int tmpMidPartMapPosInWholeGenome = indexInfo->getWholeGenomeLocation(
			midPartMapChrInt, midPartMapPosInChr);

		//norSegmentAlignLoc[0] = tmpMidPartMapPosInWholeGenome; //midPartMapPosInChr;
		vector<unsigned int> newMapPosVec;
		//newMapPosVec.push_back(tmpMidPartMapPosInWholeGenome);
		norSegmentAlignLoc.push_back(newMapPosVec);
		(norSegmentAlignLoc[0]).push_back(tmpMidPartMapPosInWholeGenome);
	
		for(int tmp = 1; tmp < segmentNum; tmp ++)
		{
			//norSegmentLength[tmp] = segInfo->norSegmentLength[tmp - 1];
			norSegmentLength.push_back(segInfo->norSegmentLength[tmp - 1]);
			//norSegmentLocInRead[tmp] = segInfo->norSegmentLocInRead[tmp - 1] + midPartSegLength;
			norSegmentLocInRead.push_back(segInfo->norSegmentLocInRead[tmp - 1] + midPartSegLength);
			//norSegmentAlignNum[tmp] = segInfo->norSegmentAlignNum[tmp - 1];
			norSegmentAlignNum.push_back(segInfo->norSegmentAlignNum[tmp - 1]);

			vector<unsigned int> newMapPosVec;
			norSegmentAlignLoc.push_back(newMapPosVec);
			for(int tmp2 = 0; tmp2 < norSegmentAlignNum[tmp]; tmp2++)
			{
				//*(norSegmentAlignLoc + (tmp*CANDALILOC) + tmp2) 
				//	= *(segInfo->norSegmentAlignLoc + ((tmp-1) * CANDALILOC) + tmp2); 
				norSegmentAlignLoc[tmp].push_back(segInfo->returnSegmentMapPos(tmp-1, tmp2));
			}
		}
	}

	Seg_Info()
	{
		longSegMinLength = 20;
		segmentNum = 0;
		repeatRegion_index = -1;
	}

	Seg_Info(Seg2ndOri_Info* seg2ndOriInfo, int mapPosIntervalStart, int mapPosIntervalEnd, 
		int chrPosStartIn2ndLevelIndex, Index_Info* indexInfo, const string& chromNameStr)
	{
		repeatRegion_index = -1;

		longSegMinLength = 18;
		segmentNum = seg2ndOriInfo->returnSegmentNum();
		//cout << "segmentNum: " << segmentNum << endl;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			//cout << "tmpSeg: " << tmpSeg << endl;

			//norSegmentLength[tmpSeg] = seg2ndOriInfo->returnSegmentLength(tmpSeg);
			norSegmentLength.push_back(seg2ndOriInfo->returnSegmentLength(tmpSeg));
			//norSegmentLocInRead[tmpSeg] = seg2ndOriInfo->returnSegmentLocInRead(tmpSeg);
			norSegmentLocInRead.push_back(seg2ndOriInfo->returnSegmentLocInRead(tmpSeg));
			
			vector<unsigned int> newMapPosVec;
			norSegmentAlignLoc.push_back(newMapPosVec);

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
						//*(norSegmentAlignLoc + tmpSeg*CANDALILOC + tmpSegCandi) 
						//	= indexInfo->getWholeGenomeLocation(chromNameStr, tmpLoc);
						norSegmentAlignLoc[tmpSeg].push_back(indexInfo->getWholeGenomeLocation(chromNameStr, tmpLoc));
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
			//norSegmentAlignNum[tmpSeg] = tmpSegCandiNum;
			norSegmentAlignNum.push_back(tmpSegCandiNum);
		}
	}		

	bool checkSegLongOrNot(int segGroupNO)
	{
		return (norSegmentLength[segGroupNO] >= longSegMinLength);
	}
	
	int checkSegRelation(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		unsigned int segEndNum1 = segNO_1;
		unsigned int segStartNum2 = segNO_2;

		unsigned int alignLoc1 //= *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			= this->returnSegmentMapPos(segNO_1, segNO_1_candiNO) - norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 //= *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			= this->returnSegmentMapPos(segNO_2, segNO_2_candiNO) - norSegmentLocInRead[segNO_2] + 1;

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
					if (gapInChr > MAX_SPLICE_LENGTH)
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
					if (gapInChr > MAX_SPLICE_LENGTH)
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

		unsigned int alignLoc1 //= *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			= this->returnSegmentMapPos(segNO_1, segNO_1_candiNO)	
				- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 //= *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			= this->returnSegmentMapPos(segNO_2, segNO_2_candiNO)
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
						if(gapInChr <= MAX_SPLICE_LENGTH)
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
					if (gapInChr > MAX_SPLICE_LENGTH)
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
					if (gapInChr > MAX_SPLICE_LENGTH)
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
		int segNO_2, int segNO_2_candiNO)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 //= *(norSegmentAlignLoc + segNO_1 * CANDALILOC + segNO_1_candiNO) 
			= this->returnSegmentMapPos(segNO_1, segNO_1_candiNO)
				- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 //= *(norSegmentAlignLoc + segNO_2 * CANDALILOC + segNO_2_candiNO) 
			= this->returnSegmentMapPos(segNO_2, segNO_2_candiNO)
				- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
			
			if(segmentDistance <= MAX_SPLICE_LENGTH)
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
		int segNO_2, int segNO_2_candiNO)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 = this->returnSegmentMapPos(segNO_1, segNO_1_candiNO)
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = this->returnSegmentMapPos(segNO_2, segNO_2_candiNO)
			- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
			
			if(segmentDistance <= MAX_SPLICE_LENGTH)
				return segmentDistance;
			else
				return segDist;
		}
		else
		{
			segmentDistance = alignLoc1 - alignLoc2;
			if(segmentDistance <= MAX_INSERTION_LENGTH)
				return 0-segmentDistance;
			else
				return segDist;
		}
	}	


	int distanceBetweenSegment_cirRNA(int segNO_1, int segNO_1_candiNO, 
		int segNO_2, int segNO_2_candiNO)
	{
		int segDist = 1000000;
		unsigned int alignLoc1 = this->returnSegmentMapPos(segNO_1, segNO_1_candiNO)
			- norSegmentLocInRead[segNO_1] + 1;
		unsigned int alignLoc2 = this->returnSegmentMapPos(segNO_2, segNO_2_candiNO)
			- norSegmentLocInRead[segNO_2] + 1;

		unsigned int segmentDistance;
		if(alignLoc2 >= alignLoc1)
		{
			segmentDistance = alignLoc2 - alignLoc1;
			
			if(segmentDistance <= MAX_SPLICE_LENGTH)
				return segmentDistance;
			else
				return segDist;
		}
		else
		{
			segmentDistance = alignLoc1 - alignLoc2;
			if(segmentDistance <= MAX_INSERTION_LENGTH)
				return segmentDistance;
			else if(segmentDistance <= MAX_SPLICE_LENGTH)
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
		int segNO_2, int segNO_2_candiNO)
	{
		int segRelation = checkSegRelation(segNO_1, segNO_1_candiNO, segNO_2, segNO_2_candiNO);
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
			if((norSegmentLength[tmpSegNO] >= longSegMinLength)&&(norSegmentAlignNum[tmpSegNO] <= SEGMENTNUM))
			{
				longSegExists = true;
				return tmpSegNO;
			}
		}
		if(!longSegExists)
			return -1;
	}

	bool mapMain_SegInfo_preIndex_repeatRegion(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		unsigned int* valLength, BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray, const string& readStringStr, RepeatRegion_Info* repeatRegionInfo)
	{
		//cout << "mapMain_SegInfo_preIndex function starts ...... "<< endl; 
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		unsigned int norSegmentNum;

		(norSegmentNum) = 0;
		bool mapMain = false;	
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_num = 0;
		//unsigned int segment_length = 0; 
		unsigned int segment_length_max = 0;//used to compare with segment_length for eache segment to get the maximum length segment
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int read_length = readLength; //READ_LENGTH;
		unsigned int interval_begin, interval_end;
		unsigned int n = (indexInfo->returnIndexSize());//size of SA
		unsigned int norAlignLoc;
		unsigned int align_chr_location;
		unsigned int align_chr;
		*valLength = 0;
		char* read_local = read;

		while (stop_loc_overall < read_length) //- 15)
		{
			segment_num++;
			bool queryFound = true;

			/////////////////////////////////////////////////////////////////////////
			/////////////////////////   K-mer search  ///////////////////////////////
			bool firstCheckForHeadOfSegmentBool = true;
			bool KmerSearchFound = false;
			int KmerMappedLength;
			unsigned int KmerSearchIndexIntervalStart = 0;
			unsigned int KmerSearchIndexIntervalEnd = 0;			
			
			if(read_length - stop_loc_overall < INDEX_KMER_LENGTH)
			{
				KmerSearchFound = false;
				firstCheckForHeadOfSegmentBool = false;
			}
			else
			{
				KmerSearchFound = this->getIndexInterval_PreIndex(readStringStr.substr(stop_loc_overall, INDEX_KMER_LENGTH), &KmerMappedLength,
					&KmerSearchIndexIntervalStart, &KmerSearchIndexIntervalEnd, PreIndexMappedLengthArray, 
					PreIndexIntervalStartArray,	PreIndexIntervalEndArray);
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

		   	 		//norSegmentLength[norSegmentNum - 1] = 1;
					norSegmentLength.push_back(1);
					//norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
					norSegmentLocInRead.push_back(stop_loc_overall + 1);
					//norSegmentAlignNum[norSegmentNum - 1] = 0;
					norSegmentAlignNum.push_back(0);
					vector<unsigned int> newMapPosVec;
					norSegmentAlignLoc.push_back(newMapPosVec);

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
		   	 	iterateNum = 0;
	   	 	}
	   	 	else // K-mer found in preIndex base
	   	 	{
	   	 		firstCheckForHeadOfSegmentBool = true;
	   	 		lcp_length = 0;
	   	 		start = 0, end = n-1;
	   	 		c = 0;
	   	 		interval_begin = KmerSearchIndexIntervalStart;
	   	 		interval_end = KmerSearchIndexIntervalEnd;
	    	 	segment_align_SArange[0] = interval_begin;
	   	 		segment_align_SArange[1] = interval_end;
	   	 		segment_align_rangeNum = interval_end - interval_begin + 1;
	   	 		iterateNum = 0;//debug;  	 		
	   	 	}

	   	 	while((c + stop_loc_overall < read_length) && (queryFound == true))
	   	 	{
	   	 		firstCheckForHeadOfSegmentBool = false;
	   	 		iterateNum++;
	   	 		if(iterateNum + stop_loc_overall >read_length)
	   	 		{
	   	 			return false;
	   	 		}
	   	 		unsigned int c_old = c;

				if(interval_begin != interval_end)
				{ 
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);

					Min = min(lcp_length, read_length - stop_loc_overall);
					//cout << "Min: " << Min << endl;
					unsigned int loc_pos = 0;

						int startToCheckPos = 0;
						if(firstCheckForHeadOfSegmentBool)
						{
							 startToCheckPos = KmerMappedLength;
						}
						//cout << "startToCheckPos: " << startToCheckPos << " Min-c_old: " << Min-c_old << endl;
		            	for(loc_pos = startToCheckPos; loc_pos < Min - c_old; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
		            		//cout << "...queryFound: " << queryFound << endl;
		            		if (!queryFound)
		            		{	
		            			break;
		            		}
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

					if (c + stop_loc_overall == read_length)
					{		
						break;			
					}	
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
			    	if(interval_begin > interval_end)
			    	{
			    		queryFound = false;
			    		stop_loc = c-1;
	          			segment_align_SArange[0] = interval_begin_ori;
	            		segment_align_SArange[1] = interval_end_ori;
	            		
	            		segment_align_rangeNum = interval_end_ori - interval_begin_ori + 1;		    			
			    		break;
			    	}
			    	else
			    	{
	          			segment_align_SArange[0] = interval_begin;
	            		segment_align_SArange[1] = interval_end;

	            		segment_align_rangeNum = interval_end - interval_begin + 1;
			    	}
				}//end if
				else 
				{
					unsigned int loc_pos = 0;
		            	for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
		            	{
		            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
		            		if (!queryFound)
		            			break;
		            	}

		    		if(queryFound) 
		    		{
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
		    		}	
	          		segment_align_SArange[0] = interval_begin;
	            	segment_align_SArange[1] = interval_end;
	            	segment_align_rangeNum = interval_end - interval_begin + 1;   	
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

	    		//norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
				norSegmentLength.push_back(tmpSegLength);
	    		//*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		norSegmentLocInRead.push_back(stop_loc_overall + 1);
	    		//*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				norSegmentAlignNum.push_back(segment_align_rangeNum);
				vector<unsigned int> newMapPosVec;
				norSegmentAlignLoc.push_back(newMapPosVec);
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
		    			//*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
		    			(norSegmentAlignLoc[norSegmentNum - 1]).push_back(sa[segment_align_SArange[0] + alignment_num] + 1);
		    		}
	    		}
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
				//norSegmentLength[norSegmentNum - 1] = stop_loc;
				norSegmentLength.push_back(stop_loc);
				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				//norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentLocInRead.push_back(stop_loc_overall + 1);
				//norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;
				norSegmentAlignNum.push_back(segment_align_rangeNum);
				vector<unsigned int> newMapPosVec;
				norSegmentAlignLoc.push_back(newMapPosVec);
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
					for (unsigned int alignment_num = 0; alignment_num < segment_align_rangeNum; alignment_num++)
				    {    			
		    			//*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
		    			norSegmentAlignLoc[norSegmentNum-1].push_back(sa[segment_align_SArange[0] + alignment_num] + 1);
		    		}
		    	}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				//cout << "stop_loc + 1: " << stop_loc +1 << endl;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
			}		

	   	}
		segmentNum = (norSegmentNum);
		mapMain = true;
		return mapMain;
	}

	string segInfoStr(Index_Info* indexInfo)
	{
		string segInfoStr;
		segInfoStr += "\nsegment Info: \n";
		unsigned int align_chr, align_chr_location;
	   	for(unsigned int k1 = 0; k1 < segmentNum; k1++)
	   	{
	  		segInfoStr = segInfoStr + "... segment " + int_to_str(k1+1) + ": " + int_to_str(norSegmentLocInRead[k1]) + "~"   
	  		+ int_to_str(norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) + 
	  		"  Length: " + int_to_str(norSegmentLength[k1]) + " Num: " + int_to_str(norSegmentAlignNum[k1]) + "\n";

	      	if(//(norSegmentLength[k1]>=10)&&
	      		(norSegmentAlignNum[k1] < CANDALILOC))
	      	{
	      		//segInfoStr += "...... Align Location: \n";
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			indexInfo->getChrLocation(
	      				//*(norSegmentAlignLoc + k1*CANDALILOC + k2), 
	      				this->returnSegmentMapPos(k1,k2),
	      				&align_chr, &align_chr_location);
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

class Seg_Info_


#endif
