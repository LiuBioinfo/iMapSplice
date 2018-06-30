// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef LOCAL_SEG_INFO_H
#define LOCAL_SEG_INFO_H

#include <stdlib.h>
#include <stdio.h>

using namespace std;

class Seg2ndOri_Info
{
//public:
private: 
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
	Seg2ndOri_Info()
	{}
	void assignSegmentMapPos(unsigned int tmpMapPos, int tmpSegGroup, int tmpSegCandi)
	{
		//*(norSegmentAlignLoc + tmpSegGroup*CANDALILOC + tmpSegCandi);
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
		//unsigned int mapPos = *(norSegmentAlignLoc + segGroupNO*CANDALILOC + segCandiNO);
		unsigned int mapPos = (norSegmentAlignLoc[segGroupNO])[segCandiNO];
		return mapPos;
	}
	unsigned int returnSegmentAlignLoc(int segGroupNO, int segCandiNO)
	{
		unsigned int mapPos = (norSegmentAlignLoc[segGroupNO])[segCandiNO];
		return mapPos;		
	}
	int returnLongSegMinLength()
	{
		return longSegMinLength;
	}
	bool mapMainSecondLevel_compressedIndex(
		char *read, unsigned int* sa, 
		BYTE* lcpCompress, 
		unsigned int* child, char* chrom, 
		BYTE* verifyChild, int readLength, 
		Index_Info* indexInfo)
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
	   	 	unsigned int lcp_length = 0;
	   	 	unsigned int start = 0, end = n-1;
	   	 	unsigned int Min;
	   	 	unsigned int c = 0;
	   	 	//cout << "char: " << (*read_local) << endl;
	   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
	   	 	//cout << "new Segment: " << segment_num << endl;
	   	 	//cout << "char: " << (*read_local) << " firstInterval: " << interval_begin << " ~ " << interval_end << endl;
	   	 	segment_align_SArange[0] = interval_begin;
	   	 	segment_align_SArange[1] = interval_end;
	   	 	segment_align_rangeNum = interval_end - interval_begin + 1;

	   	 	unsigned int iterateNum = 0;//debug;
	   	 	while((c + stop_loc_overall< read_length) && (queryFound == true))
	   	 	{
	   	 		iterateNum++;
	   	 		if(iterateNum>read_length)
	   	 		{
	   	 			return false;
	   	 		}
	   	 		unsigned int c_old = c;
				
				if(interval_begin != interval_end)
				{ 
	           		//Xinan: COMPRESS INDEX
	           		//lcp_length = getlcp(interval_begin, interval_end, lcp, child_up, child_down);
	 				lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild//child_up, child_down
	 					);
	 				//cout << "lcp_length: " << lcp_length << endl;
					//Min = min(lcp_length, read_length);
					Min = min(lcp_length, read_length - stop_loc_overall);
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
					{				
						break;			
					}	
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
	            	for(loc_pos = 0; loc_pos < read_length - c- stop_loc_overall; loc_pos++)
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

	    		//norSegmentLength[norSegmentNum - 1] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
	    		norSegmentLength.push_back(tmpSegLength);
	    		//*(norSegmentLocInRead + norSegmentNum - 1) = stop_loc_overall + 1;
	    		norSegmentLocInRead.push_back(stop_loc_overall + 1);
	    		//*(norSegmentAlignNum + norSegmentNum - 1) = segment_align_rangeNum;				
				norSegmentAlignNum.push_back(segment_align_rangeNum);
				vector<unsigned int> newMapPosVec;
				norSegmentAlignLoc.push_back(newMapPosVec);
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
				{    			
	    			//*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    			norSegmentAlignLoc[norSegmentNum-1].push_back(sa[segment_align_SArange[0] + alignment_num] + 1);
	    		}

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
				//norSegmentLength[norSegmentNum - 1] = stop_loc;
				norSegmentLength.push_back(stop_loc);
				//norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentLocInRead.push_back(stop_loc_overall + 1);
				//norSegmentAlignNum[norSegmentNum - 1] = segment_align_rangeNum;
				norSegmentAlignNum.push_back(segment_align_rangeNum);
				vector<unsigned int> newMapPosVec;
				norSegmentAlignLoc.push_back(newMapPosVec);
				for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,CANDALILOC); alignment_num++)
			    {    			
	    			//*(norSegmentAlignLoc + (norSegmentNum - 1) * CANDALILOC + alignment_num) = sa[segment_align_SArange[0] + alignment_num] + 1;
	    			(norSegmentAlignLoc[norSegmentNum-1]).push_back(sa[segment_align_SArange[0] + alignment_num] + 1);
	    		}

				unsigned int stop_loc_overall_ori = stop_loc_overall;
				read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
				stop_loc_overall = stop_loc_overall + stop_loc + 1;

	    		segment_length = stop_loc;
 
			}		
	   	}

		
		segmentNum = (norSegmentNum);
		mapMain = true;
		//debugln("mapMain ended!!!");
		return mapMain;
	}

	string segInfoStr(Index_Info* indexInfo, int chrPosStartIn2ndLevelIndex, string chrNameStr)
	{
		string segInfoStr;
		segInfoStr += "\noriginal 2nd level segment Info: \n";
		segInfoStr = segInfoStr + "SegmentNum: " + int_to_str(segmentNum) + "\n";
		unsigned int align_chr, align_chr_location;
	   	for(unsigned int k1 = 0; k1 < segmentNum; k1++)
	   	{
	  		segInfoStr = segInfoStr + "... segment " + int_to_str(k1+1) + ": " + int_to_str(norSegmentLocInRead[k1]) + "~"   
	  		+ int_to_str(norSegmentLocInRead[k1] + norSegmentLength[k1] - 1) + 
	  		"  Length: " + int_to_str(norSegmentLength[k1]) + " Num: " + int_to_str(norSegmentAlignNum[k1]) + "\n";

	      	if(//(norSegmentLength[k1]>=10)&&
	      		(norSegmentAlignNum[k1] <= CANDALILOC))
	      	{
	      		//segInfoStr += "...... Align Location: \n";
	      		for(unsigned int k2 = 0; k2 < norSegmentAlignNum[k1]; k2++)
	      		{
	      			segInfoStr = segInfoStr + "\t" + int_to_str(k2+1) +
	      			//<< *(norSegmentAlignLoc + k1*CANDALILOC + k2) 
	      			+ ". " 
	      			+ chrNameStr//(indexInfo->chrNameStr)[align_chr] 
	      			+ " " + int_to_str((int)
	      				((
	      					//*(norSegmentAlignLoc + k1*CANDALILOC + k2)
	      					(norSegmentAlignLoc[k1])[k2]
	      					) + chrPosStartIn2ndLevelIndex - 1)
	      				) 
	      			+ "\n";
	      		}
	      	}
	   	}
		return segInfoStr;//+"\n";
	}	
};

#endif