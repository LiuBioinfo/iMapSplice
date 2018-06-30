// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class Segment_Info
{
public:
	int segmentNum;
	vector<int> segmentLengthVec;
	vector<int> segmentLocInReadVec;
	vector<int> segmentAlignNumVec;
	vector< vector<unsigned int> > segmentAlignLocVec;

	Segment_Info()
	{}

	bool mapMain_segmentInfo(char *read, unsigned int* sa, BYTE* lcpCompress, 
		unsigned int* child, char* chrom, unsigned int* valLength, 
		BYTE* verifyChild, int readLength, Index_Info* indexInfo)
	{
		//input : read, sa, up, down, next, chrom; 
		//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
		
	    unsigned int *norSegmentNum;
	   	*norSegmentNum = 0;  	
		//Denote the length of each segment, assume num of segment <= 20, start from 0
	    unsigned int norSegmentLength[SEGMENTNUM];
	    //Denote the start location of segment in read, start from 1(the first location);
	    unsigned int norSegmentLocInRead[SEGMENTNUM];
	    //Denote hom many locations a segment can align in chromosome, start from 1(first location); 
	    unsigned int norSegmentAlignNum[SEGMENTNUM];
	    //Denote the start location of an alignment for a segment, start from 1;
	    // Note: for convenience, the location is for the whole genome, not for specific chromosomes.
	    unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC]; 	



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
		//unsigned int align_length[102] = {0}; 
		unsigned int n = (indexInfo->indexSize);//size of SA
		unsigned int norAlignLoc;
		unsigned int align_chr_location;
		unsigned int align_chr;
		*valLength = 0;
		//debugln("start mapMain Function!!!");
		char* read_local = read;
		while (stop_loc_overall < read_length) //- 15)
		{
			//debug("stop_loc_overall = "); debugln(stop_loc_overall);
			segment_num++;
			//debug("segmentNum = "); debugln(segment_num);
			bool queryFound = true;

	   	 	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
	   	 	{
	   	 		queryFound = false;
	   	 		//align_length[0] ++;
	   	 		stop_loc = 0;
	   	 		segment_align_SArange[0] = 1;
	   	 		segment_align_SArange[1] = 0;
	   	 		segment_align_rangeNum = 0;
	   	 		queryFound = false;   	 			
	   	 	}
	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	   	 	//debugln("before getFirstInterval");
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
	   	 	//debug("interval_begin = "); debugln(interval_begin);
	   	 	//debug("interval_end = "); debugln(interval_end);
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;
	   	 	//debug("interval_begin = "); debugln(interval_begin);
	   	 	//debug("interval_end = "); debugln(interval_end);

	   	 	unsigned int iterateNum = 0;//debug;
	   	 	while((c < read_length) && (queryFound == true))
	   	 	{
	   	 		iterateNum++;
	   	 		if(iterateNum>read_length)
	   	 		{
	   	 			//debugln("error: interateNum > readLength");
	   	 			return false;
	   	 		}
	   	 		unsigned int c_old = c;
				
				if(interval_begin != interval_end)
				{ 
	           		//Xinan: COMPRESS INDEX
	           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
					//debug("lcp_length = "); debugln(lcp_length);
					Min = min(lcp_length, read_length);
					
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
	            		if (!queryFound)
	            		{	
	            			break;
	            		}
	            	}

	            	if(!queryFound)
	            	{
	            		stop_loc = c_old + loc_pos;
	            		//align_length[stop_loc] ++;	
	            		break;
	            	}
	            	
	            	c = Min;
	            	if(*(read_local+c) == 'N')
	            	{
	            		queryFound = false; 
	            		stop_loc = c;
	            		//align_length[stop_loc] ++; 
	            		break;
	            	}
					start = interval_begin; end = interval_end;
					if (c == read_length)
					{				
						//align_length[101]++;
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
			    		//align_length[stop_loc]++; 
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
	   	 			//debug("interval_begin = "); debugln(interval_begin);
	   	 			//debug("interval_end = "); debugln(interval_end);
				}//end if
				else 
				{
					unsigned int loc_pos = 0;
	            	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
	            	{
	            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
	            		if (!queryFound)
	            			break;
	            	}

		    		if(queryFound) 
		    		{
		    			//align_length[101] ++;
		    		}
		    		else 
		    		{ 
		    			stop_loc = c+loc_pos;
		    			//align_length[stop_loc] ++; 
		    		}	
	          		segment_align_SArange[0] = interval_begin;
	            	segment_align_SArange[1] = interval_end;
	            	segment_align_rangeNum = interval_end - interval_begin + 1;   	

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
	    		(*norSegmentNum) ++;
	    		//debug("segmentNum = "); debugln(*norSegmentNum);
	    		unsigned int tmpSegLength = read_length - stop_loc_overall;

				if(tmpSegLength >= minValSegLength)
				{
					*valLength = *valLength + tmpSegLength;
				}

	    		norSegmentLength[*norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		#ifdef SEGLENGTH
	    			segmentLength1[stop_loc]++;
	    		#endif
	    		*(norSegmentLocInRead + *norSegmentNum - 1) = stop_loc_overall + 1;
	    		*(norSegmentAlignNum + *norSegmentNum - 1) = segment_align_rangeNum;				
				

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
				{    			
	    			*(norSegmentAlignLoc + (*norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				segment_length = read_length-stop_loc_overall;
				
				#ifdef DEBUG
				if(segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall;					
					segment_length_max = segment_length;
				}
				#endif

				break;
			}
			else 
			{    
				(*norSegmentNum) ++;
				//debug("segmentNum = "); debugln(*norSegmentNum);
				if(*norSegmentNum > (int)(read_length/5))
				{
					//debugln("map error, too many segments, there may exist too many Ns");
					return false;
				}
				norSegmentLength[*norSegmentNum - 1] = stop_loc;

				#ifdef SEGLENGTH
					segmentLength1[stop_loc]++;
					segmentLength2[stop_loc]++;
				#endif

				if(stop_loc >= minValSegLength )
				{
					*valLength = *valLength + stop_loc;
				}

				norSegmentLocInRead[*norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[*norSegmentNum - 1] = segment_align_rangeNum;

				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
			    {    			
	    			*(norSegmentAlignLoc + (*norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
				
				#ifdef DEBUG 		
				if (segment_length_max < segment_length)
				{
					norAlignLoc = align_chr_location - stop_loc_overall_ori;
					segment_length_max = segment_length;
				}
				#endif
			}		
	   	}
		
		#ifdef DEBUG
	   	//cout << " " << endl;
	    cout << "# of Segments = " << *norSegmentNum << endl;//; debugln();

	   	for(unsigned int k1 = 0; k1 < *norSegmentNum; k1++)
	   	{
	  		cout << "segment "<< k1+1 << ": " << norSegmentLocInRead[k1] << "~"   
	  		<< (norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) 
	  		<< "  Length: " << norSegmentLength[k1] << " Num: " << norSegmentAlignNum[k1] << endl;
	      	if((norSegmentLength[k1]>=10)&&(norSegmentAlignNum[k1] < 40))
	      	{
	      		cout << "\tAlign Location: " << endl;
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			getChrLocation(*(norSegmentAlignLoc + k1*CANDALILOC + k2), &align_chr, &align_chr_location);
	      			cout << "\t" 
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			<<  " " 
	      			<< chr_name[align_chr] << " " << align_chr_location << endl;
	      		}
	      	}
	   	}
		cout << ("segment_length_max = ") << segment_length_max << endl;
		#endif 
		
		segmentNum = *norSegmentNum;
		for(int segmentNumTmp = 0; segmentNumTmp < segmentNum; segmentNumTmp++)
		{
			vector<unsigned int> tmpAlignLocVec;

			segmentLengthVec.push_back(norSegmentLength[segmentNumTmp]);
			segmentLocInReadVec.push_back(norSegmentLocInRead[segmentNumTmp]);
			segmentAlignNumVec.push_back(norSegmentAlignNum[segmentNumTmp]);
			for(int alignLocTmp = 0; alignLocTmp < norSegmentAlignNum[segmentNumTmp]; alignLocTmp++)
			{
				tmpAlignLocVec.push_back(*(norSegmentAlignLoc + segmentNumTmp*CANDALILOC + alignLocTmp));
			}

			segmentAlignLocVec.push_back(tmpAlignLocVec);
		}

		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}

};

class Fragment_Info
{
public:
	int mapLabel;
	vector<int> segMapRangeStartVec;
	vector<int> segMapRangeEndVec;
	vector<unsigned int> segMapLocVec;

	Fragment_Info()
	{}

	bool detAliLoc_fragmentInfo(/*unsigned int norSegmentNum, unsigned int* norSegmentLength, 
		unsigned int* norSegmentLocInRead,
		unsigned int* norSegmentAlignNum, unsigned int* norSegmentAlignLoc,*/ 
		Segment_Info* segmentInfo,
		char* read, char* chrom, 
		/*int* finalMapLabel, unsigned int* segMapRangeStart, 
		unsigned int* segMapRangeEnd, unsigned int* segMapLoc,*/
		Index_Info* indexInfo)
	{
		#ifdef DEBUG
		cout << "debug -- start detAliLoc function..." << endl;
		#endif
		bool detAliLoc = false;
		unsigned int segmentNum_max = SEGMENTNUM;
		unsigned int segmentAliNum_max = CANDALILOC;
		
		unsigned int norSegMapRangeStart[SEGMENTNUM * CANDALILOC] = {0};
		unsigned int norSegMapRangeEnd[SEGMENTNUM * CANDALILOC] = {0};
		
		unsigned int norSegMapLoc[SEGMENTNUM * CANDALILOC] = {0};

		unsigned int norValSegRange[2] = {50, 0};

		map <unsigned int, unsigned int> candAliLocMap;
		map <unsigned int, unsigned int> ::iterator candAliLoc_it;
		unsigned int mapLabel = 0; // normal
		unsigned int candAliLocWeight[segmentNum_max * segmentAliNum_max]; //normal

		#ifdef DEBUG
		cout << "debug -- start detAliLoc for alignment..." << endl;
		#endif
		for (unsigned int tmpSegNum = 1; tmpSegNum <= segmentInfo->norSegmentNum; tmpSegNum++)
		{
			//if ((*(norSegmentLength + tmpSegNum - 1) < minValSegLength)||
			//	(*(norSegmentAlignNum + tmpSegNum - 1) > POSSIBLE_MAP_CASES_MAX))
			if( ( segmentInfo->segmentLengthVec[tmpSegNum - 1] < minValSegLength ) ||
				( segmentInfo->segmentAlignNumVec[tmpSegNum - 1] > POSSIBLE_MAP_CASES_MAX ) )
			{
				continue;
			}			

			for(unsigned int tmpSegAliLoc = 1; 
				tmpSegAliLoc <= segmentInfo->segmentAlignNumVec[tmpSegNum-1]; tmpSegAliLoc++) 
			{				
				/*unsigned int tmpAlignLoc = *(norSegmentAlignLoc + (tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1) 
				- *(norSegmentLocInRead + tmpSegNum - 1) + 1;*/

				unsigned int tmpAlignLoc 
					= (segmentInfo->segmentAlignLocVec[tmpSegNum-1])[tmpSegAliLoc-1]
						- segmentInfo->segmentLocInReadVec[tmpSegNum-1] + 1;

				if(tmpAlignLoc > (indexInfo->indexSize))
				{
					continue;
				}
				candAliLoc_it = candAliLocMap.find(tmpAlignLoc);

				//if ((*(norSegmentLength + tmpSegNum - 1) >= minValSegLength) 
				//	&& (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
				if( (segmentInfo->segmentLengthVec[tmpSegNum-1] >= minValSegLength)
					&& (candAliLoc_it == candAliLocMap.end()) )  //cannot find, insert
				{				
					candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
					//candAliLocWeight[mapLabel] = *(norSegmentLength + tmpSegNum - 1);
					candAliLocWeight[mapLabel] = segmentInfo->segmentLengthVec[tmpSegNum-1];

					norSegMapRangeStart[mapLabel] = tmpSegNum;
					norSegMapRangeEnd[mapLabel] = tmpSegNum;
					norSegMapLoc[mapLabel] = tmpAlignLoc;
					mapLabel ++;
				}
				//else if (*(norSegmentLength + tmpSegNum - 1) >= minValSegLength)// find, add corresponding weight 
				else if( segmentInfo->segmentLengthVec[tmpSegNum-1] >= minValSegLength ) // find, add corresponding weight 
				{
					candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + segmentInfo->segmentLengthVec[tmpSegNum-1]; //*(norSegmentLength + tmpSegNum - 1);
					norSegMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
				}
				else
				{
				}
			}
		}
		candAliLocMap.clear();

		if((mapLabel <= 0) || (mapLabel > (POSSIBLE_MAP_CASES_MAX/* * (READ_LENGTH/minValSegLength)*/)))
		{
			#ifdef DEBUG
			cout << "debug -- mapLabel too small or large : " << mapLabel << endl;
			#endif
			(*finalMapLabel) = 0;
			return false;
		}

		#ifdef DEBUG
		cout << "# of Label is : " << mapLabel << endl;
		for (unsigned int tmp_label = 0; tmp_label < mapLabel; tmp_label++)
		{
			unsigned int tmpSegMapChr, tmpSegMapChrPos;
			getChrLocation(norSegMapLoc[tmp_label], &tmpSegMapChr, &tmpSegMapChrPos);

			cout << "label = " << tmp_label << ", loc = " << //norSegMapLoc[tmp_label] 
			chrNameStr[tmpSegMapChr] << " " << tmpSegMapChrPos << ", map_range = " 
				<< norSegMapRangeStart[tmp_label] << "~" << norSegMapRangeEnd[tmp_label] << endl;
		}
		#endif

		//////////////try second Level detAlignLoc

		(*finalMapLabel) = detAliLocShortSegIncluded_fragmentInfo(
			/*norSegmentNum, norSegmentLength, norSegmentLocInRead, 
			norSegmentAlignNum, 
			norSegmentAlignLoc, */
			segmentInfo,
			mapLabel, norSegMapLoc, 
			segMapRangeStart, segMapRangeEnd, segMapLoc, indexInfo);	

		if(((*finalMapLabel) <= 0) || ((*finalMapLabel) > POSSIBLE_MAP_CASES_MAX/* * (READ_LENGTH/minValSegLength)*/))
		{
			#ifdef DEBUG
			cout << "debug -- mapLabel too small or large : " << (*finalMapLabel) << endl;
			#endif
			(*finalMapLabel) = 0;
			return false;
		}	
		
		#ifdef DEBUG
		cout << "# of Label is : " << (*finalMapLabel) << endl;
		for (unsigned int tmp_label = 0; tmp_label < (*finalMapLabel); tmp_label++)
		{
			unsigned int tmpSegMapChr, tmpSegMapChrPos;
			getChrLocation(segMapLoc[tmp_label], &tmpSegMapChr, &tmpSegMapChrPos);

			cout << "label = " << tmp_label << ", loc = " << //segMapLoc[tmp_label] 
			chrNameStr[tmpSegMapChr] << " " << tmpSegMapChrPos << ", map_range = " 
				<< segMapRangeStart[tmp_label] << "~" << segMapRangeEnd[tmp_label] << endl;
		}
		#endif



		mapLabel = *finalMapLabel;
		for(int mapLabelTmp = 0; mapLabelTmp < mapLabel; mapLabelTmp++)
		{
			segMapRangeStartVec.push_back();
			segMapRangeEndVec.push_back();
			segMapLocVec.push_back();			
		}


		detAliLoc = true;
		return detAliLoc;
	}

	int detAliLocShortSegIncluded_fragmentInfo(/*unsigned int segmentNum, 
		unsigned int* segmentLength, unsigned int* segmentLocInRead,
		unsigned int* segmentAlignNum, unsigned int* segmentAlignLoc, */
		Segment_Info* segmentInfo,

		unsigned int oldMapLabel, 
		//unsigned int* oldSegMapRangeStart, unsigned int* oldSegMapRangeEnd, 
		unsigned int* oldSegMapPos, 
		//unsigned int* mapLabel, 
		unsigned int* segMapRangeStart, unsigned int* segMapRangeEnd, 
		unsigned int* segMapLoc,
		Index_Info* indexInfo
		)
	{
		bool detAliLocShortSegIncluded = false;
		unsigned int segmentNum_max = SEGMENTNUM;
		unsigned int segmentAliNum_max = CANDALILOC;
		map <unsigned int, unsigned int> candAliLocMap;
		map <unsigned int, unsigned int> ::iterator candAliLoc_it;
		unsigned int mapLabel = 0; // normal
		unsigned int candAliLocWeight[segmentNum_max * segmentAliNum_max]; //normal


		#ifdef DEBUG
		cout << "debug -- start second level detAliLoc for alignment..." << endl;
		#endif
		for (unsigned int tmpSegNum = 1; tmpSegNum <= segmentNum; tmpSegNum++)
		{
			if (//(*(norSegmentLength + tmpSegNum - 1) < minValSegLength)||
				(*(segmentAlignNum + tmpSegNum - 1) > POSSIBLE_MAP_CASES_MAX))
			{
				continue;
			}		

			for(unsigned int tmpSegAliLoc = 1; tmpSegAliLoc <= segmentAlignNum[tmpSegNum-1]; tmpSegAliLoc++) 
			{				
				unsigned int tmpAlignLoc = *(segmentAlignLoc + (tmpSegNum-1)*segmentAliNum_max + tmpSegAliLoc - 1) 
				- *(segmentLocInRead + tmpSegNum - 1) + 1;
				if(tmpAlignLoc > indexInfo->indexSize)
				{
					continue;
				}
				candAliLoc_it = candAliLocMap.find(tmpAlignLoc);

				if ((*(segmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
				{				
					candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
					candAliLocWeight[mapLabel] = *(segmentLength + tmpSegNum - 1);

					segMapRangeStart[mapLabel] = tmpSegNum;
					segMapRangeEnd[mapLabel] = tmpSegNum;
					segMapLoc[mapLabel] = tmpAlignLoc;
					mapLabel ++;
				}
				else if ((*(segmentLength + tmpSegNum - 1) < minValSegLength) && (candAliLoc_it == candAliLocMap.end()))  //cannot find, insert
				{	
					if(checkShortSegWithCandLoc(tmpAlignLoc, oldMapLabel, oldSegMapPos, indexInfo))
					{			
						candAliLocMap.insert(pair<unsigned int, unsigned int>(tmpAlignLoc, mapLabel));
						candAliLocWeight[mapLabel] = *(segmentLength + tmpSegNum - 1);

						segMapRangeStart[mapLabel] = tmpSegNum;
						segMapRangeEnd[mapLabel] = tmpSegNum;
						segMapLoc[mapLabel] = tmpAlignLoc;
						mapLabel ++;
					}
				}
				else if ((*(segmentLength + tmpSegNum - 1) >= minValSegLength) && (candAliLoc_it != candAliLocMap.end()))// find, add corresponding weight 
				{
					candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + *(segmentLength + tmpSegNum - 1);
					segMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
				}
				else //((*(segmentLength + tmpSegNum - 1) < minValSegLength) && (candAliLoc_it != candAliLocMap.end()))// find, add corresponding weight 
				{
					candAliLocWeight[candAliLoc_it->second] = candAliLocWeight[candAliLoc_it->second] + *(segmentLength + tmpSegNum - 1);
					segMapRangeEnd[candAliLoc_it->second] = tmpSegNum;
				}
			}
		}
		candAliLocMap.clear();	

		return mapLabel;
	}

};