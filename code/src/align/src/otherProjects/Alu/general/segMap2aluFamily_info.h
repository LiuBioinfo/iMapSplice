// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SEGMAP2ALUFAMILY_INFO_H
#define SEGMAP2ALUFAMILY_INFO_H

#define MAP_TO_ALUFAMILY_SEG_NUM_MAX 10

class SegMap2aluFamily_Info
{
private:
	unsigned int segmentNum;
	unsigned int norSegmentLength[MAP_TO_ALUFAMILY_SEG_NUM_MAX];
	unsigned int norSegmentLocInRead[MAP_TO_ALUFAMILY_SEG_NUM_MAX];
	unsigned int norSegmentAlignNum[MAP_TO_ALUFAMILY_SEG_NUM_MAX];
	bool segMapBool;
public:
	SegMap2aluFamily_Info()
	{
		segmentNum = 0;
		segMapBool = false;
	}

	bool returnSegMapBool(int tmpReadLength)
	{
		if(!segMapBool)
			return false;
		if(segmentNum > MAP_TO_ALUFAMILY_SEG_NUM_MAX)
			return false;
		int tmpSegLenSum = 0;		
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{
			int tmpSegLen = norSegmentLength[tmp];
			if(tmpSegLen < 0)
				return false;
			tmpSegLenSum += tmpSegLen;
		}
		if(tmpSegLenSum > tmpReadLength)
			return false;
		else
			return true;
	}

	void copy2targetSegLenVec(vector<int>& targetSegLenVec)
	{
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{
			targetSegLenVec.push_back(norSegmentLength[tmp]);
		}
	}

	int returnSegMapScore()
	{
		int tmpScoreSum = 0;
		for(int tmp = 0; tmp < segmentNum; tmp++)
		{
			int tmpScore = (norSegmentLength[tmp])*(norSegmentLength[tmp]);
			tmpScoreSum += tmpScore;
		}
		return tmpScoreSum;
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

	void generateSegLen_oneReadOneDir(
		char *read, unsigned int* sa, BYTE* lcpCompress, unsigned int* child, char* chrom, 
		 BYTE* verifyChild, int readLength, Index_Info* indexInfo,
		int* PreIndexMappedLengthArray, unsigned int* PreIndexIntervalStartArray,
		unsigned int* PreIndexIntervalEndArray, const string& readStringStr)
	{
		unsigned int norSegmentNum = 0;
		unsigned int stop_loc = 0; // location in one segment for iterations
		unsigned int stop_loc_overall = 0; //location in the whole read for iterations
		unsigned int segment_num = 0;
		//unsigned int segment_length = 0; 
		unsigned int segment_align_SArange[2] = {0,0};//legal align location in SA
		unsigned int segment_align_rangeNum = 0;
		unsigned int read_length = readLength; //READ_LENGTH;
		unsigned int interval_begin, interval_end;
		unsigned int n = (indexInfo->returnIndexSize());//size of SA
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
		   	 		
		   	 		if(norSegmentNum >= MAP_TO_ALUFAMILY_SEG_NUM_MAX)
		   	 		{
		   	 			segmentNum = MAP_TO_ALUFAMILY_SEG_NUM_MAX;
		   	 			segMapBool = false;
		   	 			return;
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
		   			segMapBool = false;
		   	 		return;
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
		            		break;
		            }
            		//cout << "queryFound: " << queryFound << endl;
	            	if(!queryFound)
	            	{
	            		stop_loc = c_old + loc_pos;
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
						break;			
					unsigned int interval_begin_ori = interval_begin;
					unsigned int interval_end_ori = interval_end;
					//cout << "interval_begin_ori: " << interval_begin_ori << " interval_end_ori: " << interval_end << endl;
			    	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			    		chrom, verifyChild);
			    	if(interval_begin > interval_end)
			    	{
			    		//cout << "interval_begin > interval_end" << endl;
			    		queryFound = false;
			    		stop_loc = c;//c-1;// fixed 11/04/2015
	          			segment_align_SArange[0] = interval_begin_ori;
	            		segment_align_SArange[1] = interval_end_ori;       		
	            		segment_align_rangeNum = interval_end_ori - interval_begin_ori + 1;		    			
			    		break;
			    	}
			    	else
			    	{
			    		//cout << "interval_begin <= interval_end" << endl;
	          			segment_align_SArange[0] = interval_begin;
	            		segment_align_SArange[1] = interval_end;
	            		segment_align_rangeNum = interval_end - interval_begin + 1; // == 1
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
		            		break;
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
	    		if(norSegmentNum > MAP_TO_ALUFAMILY_SEG_NUM_MAX)
	    		{
	    			segmentNum = (MAP_TO_ALUFAMILY_SEG_NUM_MAX);
	    			segMapBool = false;
	    			return;
	    		}
	   	 		if(norSegmentNum > (int)(read_length/5))
	   	 		{
	   	 			segmentNum = (int)(read_length/5);
	   	 			segMapBool = false;
	   	 			return;
	   	 		}
	    		//cout << "segmentNum: " << norSegmentNum << endl;
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;
	    		norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				break;
			}
			else 
			{    
				(norSegmentNum) ++;
				if((norSegmentNum > (int)(read_length/5)) || (norSegmentNum > MAP_TO_ALUFAMILY_SEG_NUM_MAX))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					if((int)(read_length/5) > MAP_TO_ALUFAMILY_SEG_NUM_MAX)
						segmentNum = (MAP_TO_ALUFAMILY_SEG_NUM_MAX);
					else
						segmentNum = (int)(read_length/5);
	   	 			segMapBool = false;
	   	 			return;
				}
				//cout << "stop_loc: " << stop_loc << endl;
				norSegmentLength[norSegmentNum - 1] = stop_loc;
				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				//cout << "stop_loc + 1: " << stop_loc +1 << endl;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;
			}		
	   	}
		segmentNum = (norSegmentNum);
		// Xinan: fixMe: some potential bugs -- segmentLength > read length, to fix
		//int tmpSegLengthSum = 0;
		for(int tmpSeg = 0; tmpSeg < segmentNum; tmpSeg++)
		{
			//tmpSegLengthSum = tmpSegLengthSum + norSegmentLength[tmpSeg]
			if(norSegmentLength[tmpSeg] > readLength)
			{
				segmentNum = 0;
				norSegmentNum = 0;
	   	 		segMapBool = false;
	   	 		return;
	 		}
		}
	   	segMapBool = true;
		return;
	}
};

class SegMap2aluFamily_Info_PeReadBothDir
{
private:
	vector<int> segLenVec_Nor1, segLenVec_Rcm1, segLenVec_Nor2, segLenVec_Rcm2;
	int segMapScore_Nor1, segMapScore_Rcm1, segMapScore_Nor2, segMapScore_Rcm2;
	int segMapScore_max;

	SegMap2aluFamily_Info segMap2aluFamilyInfo_Nor1;
	SegMap2aluFamily_Info segMap2aluFamilyInfo_Rcm1;
	SegMap2aluFamily_Info segMap2aluFamilyInfo_Nor2;
	SegMap2aluFamily_Info segMap2aluFamilyInfo_Rcm2;
public:
	SegMap2aluFamily_Info_PeReadBothDir()
	{}

	string returnSegLenVecStr()
	{
		int tmpSegLenVec_Nor1_size = segLenVec_Nor1.size();
		int tmpSegLenVec_Rcm1_size = segLenVec_Rcm1.size();
		int tmpSegLenVec_Nor2_size = segLenVec_Nor2.size();
		int tmpSegLenVec_Rcm2_size = segLenVec_Rcm2.size();
		string tmpStr_Nor1;
		if(tmpSegLenVec_Nor1_size == 0)
			tmpStr_Nor1 = "NULL";
		else
		{
			for(int tmp = 0; tmp < tmpSegLenVec_Nor1_size; tmp++)
			{	
				tmpStr_Nor1 += (int_to_str(segLenVec_Nor1[tmp]));
				tmpStr_Nor1 += ",";
			}
		}
		string tmpStr_Rcm1;
		if(tmpSegLenVec_Rcm1_size == 0)
			tmpStr_Rcm1 = "NULL";
		else
		{
			for(int tmp = 0; tmp < tmpSegLenVec_Rcm1_size; tmp++)
			{	
				tmpStr_Rcm1 += (int_to_str(segLenVec_Rcm1[tmp]));
				tmpStr_Rcm1 += ",";
			}
		}		
		string tmpStr_Nor2;
		if(tmpSegLenVec_Nor2_size == 0)
			tmpStr_Nor2 = "NULL";
		else
		{
			for(int tmp = 0; tmp < tmpSegLenVec_Nor2_size; tmp++)
			{	
				tmpStr_Nor2 += (int_to_str(segLenVec_Nor2[tmp]));
				tmpStr_Nor2 += ",";
			}
		}
		string tmpStr_Rcm2;
		if(tmpSegLenVec_Rcm2_size == 0)
			tmpStr_Rcm2 = "NULL";
		else
		{
			for(int tmp = 0; tmp < tmpSegLenVec_Rcm2_size; tmp++)
			{	
				tmpStr_Rcm2 += (int_to_str(segLenVec_Rcm2[tmp]));
				tmpStr_Rcm2 += ",";
			}
		}
		string tmpStr = tmpStr_Nor1 + "\n" + tmpStr_Rcm1 + "\n"
			+ tmpStr_Nor2 + "\n" + tmpStr_Rcm2;
		return tmpStr;
	}

	int returnMaxSegMapScore()
	{
		return segMapScore_max;
	}

	void generateMaxSegMapScore()
	{
		int segMapScore_Nor1Rcm2 = segMapScore_Nor1 + segMapScore_Rcm2;
		int segMapScore_Nor2Rcm1 = segMapScore_Nor2 + segMapScore_Rcm1;
		if(segMapScore_Nor1Rcm2 > segMapScore_Nor2Rcm1)
			segMapScore_max = segMapScore_Nor1Rcm2;
		else
			segMapScore_max = segMapScore_Nor2Rcm1;
	}

	void segmentMap(unsigned int* sa, BYTE* lcpCompress, unsigned int* childTab, char* chrom, BYTE* verifyChild, 
		Index_Info* indexInfo, int* preIndexMappedLengthArray, unsigned int* preIndexIntervalStartArray,
		unsigned int* preIndexIntervalEndArray, PE_Read_Info& readInfo)
	{
		//cout << "start to do segmentMap_Nor1" << endl;
		char* read = const_cast<char*>((readInfo.returnReadSeq_1()).c_str());
		segMap2aluFamilyInfo_Nor1.generateSegLen_oneReadOneDir(read, sa, lcpCompress, childTab, chrom, 
			verifyChild, (readInfo.returnReadSeqLength_1()), indexInfo, preIndexMappedLengthArray, 
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnReadSeq_1()));
		//cout << "start to do segmentMap_Rcm1" << endl;
		char* read_RC = const_cast<char*>((readInfo.returnRcmReadSeq_1()).c_str());
		segMap2aluFamilyInfo_Rcm1.generateSegLen_oneReadOneDir(read_RC, sa, lcpCompress, childTab, chrom, 
			verifyChild, (readInfo.returnReadSeqLength_1()), indexInfo, preIndexMappedLengthArray, 
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnRcmReadSeq_1()));
		//cout << "start to do segmentMap_Nor2" << endl;
		char* read_PE = const_cast<char*>((readInfo.returnReadSeq_2()).c_str());
		segMap2aluFamilyInfo_Nor2.generateSegLen_oneReadOneDir(read_PE, sa, lcpCompress, childTab, chrom, 
			verifyChild, (readInfo.returnReadSeqLength_2()), indexInfo, preIndexMappedLengthArray, 
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnReadSeq_2()));		
		//cout << "start to do segmentMap_Rcm2" << endl;
		char* read_RC_PE = const_cast<char*>((readInfo.returnRcmReadSeq_2()).c_str());
		segMap2aluFamilyInfo_Rcm2.generateSegLen_oneReadOneDir(read_RC_PE, sa, lcpCompress, childTab, chrom, 
			verifyChild, (readInfo.returnReadSeqLength_2()), indexInfo, preIndexMappedLengthArray, 
			preIndexIntervalStartArray, preIndexIntervalEndArray, (readInfo.returnRcmReadSeq_2()));
	}

	void generate_segLenVec_segMapScore(PE_Read_Info& readInfo)
	{
		int readLength_1 = readInfo.returnReadSeqLength_1();
		int readLength_2 = readInfo.returnReadSeqLength_2();
		bool segMapBool_Nor1 = segMap2aluFamilyInfo_Nor1.returnSegMapBool(readLength_1);
		if(segMapBool_Nor1)
		{
			segMap2aluFamilyInfo_Nor1.copy2targetSegLenVec(segLenVec_Nor1);
			segMapScore_Nor1 = segMap2aluFamilyInfo_Nor1.returnSegMapScore();
		}
		else
			segMapScore_Nor1 = 0;

		bool segMapBool_Rcm1 = segMap2aluFamilyInfo_Rcm1.returnSegMapBool(readLength_1);
		if(segMapBool_Rcm1)
		{
			segMap2aluFamilyInfo_Rcm1.copy2targetSegLenVec(segLenVec_Rcm1);
			segMapScore_Rcm1 = segMap2aluFamilyInfo_Rcm1.returnSegMapScore();
		}
		else
			segMapScore_Rcm1 = 0;

		bool segMapBool_Nor2 = segMap2aluFamilyInfo_Nor2.returnSegMapBool(readLength_2);
		if(segMapBool_Nor2)
		{
			segMap2aluFamilyInfo_Nor2.copy2targetSegLenVec(segLenVec_Nor2);
			segMapScore_Nor2 = segMap2aluFamilyInfo_Nor2.returnSegMapScore();
		}
		else
			segMapScore_Nor2 = 0;

		bool segMapBool_Rcm2 = segMap2aluFamilyInfo_Rcm2.returnSegMapBool(readLength_2);	
		if(segMapBool_Rcm2)
		{
			segMap2aluFamilyInfo_Rcm2.copy2targetSegLenVec(segLenVec_Rcm2);
			segMapScore_Rcm2 = segMap2aluFamilyInfo_Rcm2.returnSegMapScore();
		}
		else
			segMapScore_Rcm2 = 0;
	}
};
#endif