// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
/*    
 *    buildWholeGenomeIndex.cpp
 *	  MapSplice3
 *
 *    Copyright (C) 2016 University of Kentucky and
 *                       Xinan Liu
 *
 *    Authors: Xinan Liu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <sys/types.h>    
#include <dirent.h>    
#include <stdio.h>    
#include <errno.h>
#include <set>
#include "buildIndexParameter.h"
//#include "buildIndex_info.h"

#define Range 6
#define LONGLCP 255
#define NULL_NUM 4294967290 
#define PREINDEX_STRINGLENGTH 14
#define CANDALILOC 100
#define SEGMENTNUM 20
#define minValSegLength 20
#define PreIndexSize 268435456

using namespace std; 

typedef unsigned char BYTE;

int read_num = 0; // to calculate the read num;

//#include "seg_info.h"

int baseChar2intArray[26] = {0, 100, 1, 100, 100, 100, 2,
			100, 100, 100, 100, 100, 100, 100,
			100, 100, 100, 100, 100, 3, 
			100, 100, 100, 100, 100, 100};


inline string int_to_str(int numerical)
{
		char c[100];
		sprintf(c,"%d",numerical);
		string str(c);
		return str;
}


unsigned int getPreIndexNO(const string& readPreStr)
{
	int preIndexStrSize = readPreStr.length();

	unsigned int preIndexNO = 0;

	int baseForCount = 1;

	for(int tmp = preIndexStrSize - 1; tmp >= 0; tmp--)
	{
		//char tmpChar = readPreStr.at(tmp);
		//unsigned int tmpArrayIndex = tmpChar - 'A';
		//baseChar2int(tmpChar);
		//preIndexNO = preIndexNO + baseChar2int(tmpChar) * baseForCount;
		preIndexNO = preIndexNO + baseChar2intArray[(readPreStr.at(tmp) - 'A')] * baseForCount;
		baseForCount = baseForCount * 4;
	}		
	return preIndexNO;
}

unsigned int getChildUpValue(unsigned int *child, BYTE *verifyChild, unsigned int index)
{
	return ((verifyChild[index-1]==1)*child[index-1]);
}

unsigned int getChildDownValue(unsigned int *child, BYTE *verifyChild, unsigned int index)
{
	if(verifyChild[index] == 4)
		return child[child[index]-1];
	else if (verifyChild[index] == 2)
		return child[index];
	else 
		return 0;
}

unsigned int getChildNextValue(unsigned int *child, BYTE *verifyChild, unsigned int index)
{
	return ((verifyChild[index]>2)*child[index]);
}


void getFirstInterval(char ch, unsigned int *interval_begin, unsigned int *interval_end, 
	unsigned int* child, BYTE* verifyChild)
{
	if (ch == 'C')
	{
		*interval_begin = ((verifyChild[0]>2)*child[0]);
		*interval_end = ((verifyChild[*interval_begin]>2)*child[*interval_begin])-1;//getChildNextValue(child, verifyChild, *interval_begin) - 1;
	}
	else if(ch == 'A') //ch == 'A'
	{
		*interval_begin = 0;
		*interval_end = ((verifyChild[0]>2)*child[0]) - 1;
	}
	else //ch == 'G' or ch == 'T'
	{
		if (ch == 'G')
		{
			unsigned int child_next_tmp = ((verifyChild[0]>2)*child[0]);
			*interval_begin = ((verifyChild[child_next_tmp]>2)*child[child_next_tmp]);
			*interval_end = ((verifyChild[*interval_begin]>2)*child[*interval_begin]) - 1;
		}
		else
		{
			//cout << " Yes ! T " << endl;
			unsigned int child_next_tmp = ((verifyChild[0]>2)*child[0]);
			unsigned int child_next_tmp2 = ((verifyChild[child_next_tmp]>2)*child[child_next_tmp]);
			*interval_begin = ((verifyChild[child_next_tmp2]>2)*child[child_next_tmp2]);
			*interval_end = ((verifyChild[*interval_begin]>2)*child[*interval_begin]) - 1;
		}
	}
}

void getInterval_old(unsigned int start, unsigned int end, 
	unsigned int position, char ch, 
	unsigned int *interval_begin, unsigned int *interval_end, 
	unsigned int* sa, 
	unsigned int* child,
	char* chrom,
	BYTE* verifyChild)
{   
	unsigned int index_begin;
	unsigned int pos;
	*interval_end = 0;
	*interval_begin = 1;
	if(ch == 'C')
	{
		unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
		if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
			index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
		else index_begin = getChildDownValue(child, verifyChild, start);

		if(chrom[sa[index_begin]+position] == 'C')
		{
			*interval_begin = index_begin;
			*interval_end = ((verifyChild[index_begin]>2)*child[index_begin]) - 1; 
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}	
		else if(chrom[sa[start]+position] == 'C')
		{
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}		
	}
	else if(ch == 'A') // ch == 'A'
	{
		unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
		if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
			index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
		else index_begin = getChildDownValue(child, verifyChild, start);

		if(chrom[sa[start]+position] == 'A')
		{
			*interval_begin = start;
			*interval_end = index_begin - 1;
			if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
		}	
	}
	else // ch == 'G' or ch == 'T'
	{	
		if(ch == 'G')
		{
			unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(child, verifyChild, start);
			
			pos = ((verifyChild[index_begin]>2)*child[index_begin]);
			if(chrom[sa[pos]+position] == 'G')
			{
				*interval_begin = pos;
				*interval_end = ((verifyChild[pos]>2)*child[pos]) - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}
			else if(chrom[sa[start]+position] == 'G')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[index_begin]+position] == 'G')
			{
				*interval_begin = index_begin;
				*interval_end = ((verifyChild[index_begin]>2)*child[index_begin]) - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}				
		}
		else //(ch == 'T')
		{

			unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(child, verifyChild, start);

			if(chrom[sa[start]+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[index_begin]+position] == 'T')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((verifyChild[index_begin]>2)*child[index_begin]) - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else 
			{
				//pos = child_next[index_begin];
				pos = ((verifyChild[index_begin]>2)*child[index_begin]);
				if(chrom[sa[pos]+position] == 'T')
				{
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((verifyChild[pos]>2)*child[pos]) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					//pos = child_next[pos];
					pos = ((verifyChild[pos]>2)*child[pos]);
					if(chrom[sa[pos]+position] == 'T')
					{
						*interval_begin = pos;
						//*interval_end = child_next[pos] - 1;
						*interval_end = ((verifyChild[pos]>2)*child[pos]) - 1;
						if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
					}
				}						
			}
		}
	}			
}

void getInterval(unsigned int start, unsigned int end, 
	unsigned int position, char ch, 
	unsigned int *interval_begin, unsigned int *interval_end, 
	unsigned int* sa, 
	unsigned int* child,
	char* chrom,
	BYTE* verifyChild)
{   
	//#ifdef CAL_TIME
	//getInterval_begin = clock();
	//#endif
	unsigned int index_begin;
	unsigned int pos;
	*interval_end = 0;
	*interval_begin = 1;
	//cout << "......getInterval: " << endl;
	//cout << "start: " << start << " end: " << end << endl;
	//cout << "ch: " << ch << endl;
	switch(ch)
	{
		case 'C':
		{
			//if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
	       	//	index_begin = child_up[end+1];
			//else index_begin = child_down[start];
			unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(child, verifyChild, start);

	
			if(chrom[sa[start]+position] == 'C')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else if(chrom[sa[index_begin]+position] == 'C')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((verifyChild[index_begin]>2)*child[index_begin]) - 1; 
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}
			else
			{
				//cout << "no C found !" << endl;
			}
		}
			break;
		case 'A':
		{
			//if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
	       	//	index_begin = child_up[end+1];
			//else index_begin = child_down[start];
			unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(child, verifyChild, start);

			if(chrom[sa[start]+position] == 'A')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}
			else
			{
				//cout << "no A found !" << endl;
			}	
		}	
			break;
		case 'G':
			{
				//if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
	       		//	index_begin = child_up[end+1];
				//else index_begin = child_down[start];
				unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
				if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
					index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
				else index_begin = getChildDownValue(child, verifyChild, start);

				//pos = child_next[index_begin];
				pos = ((verifyChild[index_begin]>2)*child[index_begin]);
				
				//cout << "case 'G' !" << endl;
				//cout << "pos = " << pos << endl;
				//cout << "start = " << start << endl;
				//cout << "index_begin = " << index_begin << endl;


				if(chrom[sa[start]+position] == 'G')
				{
					//cout << " (chrom[sa[start]+position] == 'G') " << endl;
					*interval_begin = start;
					*interval_end = index_begin - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}	
				else if(chrom[sa[index_begin]+position] == 'G')
				{
					//cout << " (chrom[sa[index_begin]+position] == 'G') " << endl;
					*interval_begin = index_begin;
					//*interval_end = child_next[index_begin]-1;
					*interval_end = ((verifyChild[index_begin]>2)*child[index_begin]) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else if(chrom[sa[pos]+position] == 'G')
				{
					//cout << " (chrom[sa[pos]+position] == 'G') " << endl;
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((verifyChild[pos]>2)*child[pos]) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					//cout << "no G found !" << endl;
				}				
			}
			break;
		case 'T':
		default:
		{
			unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(child, verifyChild, start);

			if(chrom[sa[start]+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[sa[index_begin]+position] == 'T')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((verifyChild[index_begin]>2)*child[index_begin]) - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else 
			{
				//pos = child_next[index_begin];
				pos = ((verifyChild[index_begin]>2)*child[index_begin]);
				if(chrom[sa[pos]+position] == 'T')
				{
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((verifyChild[pos]>2)*child[pos]) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					//pos = child_next[pos];
					pos = ((verifyChild[pos]>2)*child[pos]);
					if(chrom[sa[pos]+position] == 'T')
					{
						*interval_begin = pos;
						//*interval_end = child_next[pos] - 1;
						*interval_end = ((verifyChild[pos]>2)*child[pos]) - 1;
						if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
					}
				}						
			}
		}	
			break;

	}		
	//#ifdef CAL_TIME
	//getInterval_end = clock();
	//getInterval_cost = getInterval_cost + getInterval_end - getInterval_begin;
	//#endif
}



unsigned int getlcp(unsigned int start, unsigned int end, BYTE* lcpCompress, unsigned int* child, BYTE* verifyChild)
{
	unsigned int tmpIndex;
	
	unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
	if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
		tmpIndex = child_up_value_tmp;
	else
		tmpIndex = getChildDownValue(child, verifyChild, start);


	return lcpCompress[tmpIndex];
}

unsigned int min(unsigned int a, unsigned int b)
{
   unsigned int min;
   if(a >= b)
      min = b;
   else
	  min = a;
   return min;
}


bool mapMainForIndexStringHash(
	char *read, unsigned int* sa, BYTE* lcpCompress, unsigned int* child, char* chrom, int* mappedLength, 
	unsigned int* indexIntervalStart, unsigned int* indexIntervalEnd, BYTE* verifyChild, int readLength, int MAX)
{
	//input : read, sa, up, down, next, chrom; 
	//output: segmentNum, segmentLength, segmentLocInread, segmentAlignNum, segmentAlignLoc 
	bool mapMain = false;	
	unsigned int stop_loc = 0; // location in one segment for iterations
	unsigned int read_length = readLength; //READ_LENGTH;
	unsigned int interval_begin, interval_end;
	//unsigned int align_length[102] = {0}; 
	unsigned int n = MAX;//size of SA
	char* read_local = read;
	bool queryFound = true;

   	if((*read_local != 'A')&&(*read_local != 'C')&&(*read_local != 'G')&&(*read_local != 'T')) 
   	{
		(*mappedLength) = 0;
		return true;		 			
   	}
   	unsigned int lcp_length = 0;
   	unsigned int start = 0, end = n-1;
   	unsigned int Min;
   	unsigned int c = 0;
   	//cout << "read " << (*read_local) << endl;
   	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);
   	//cout << "firstInterval: " << interval_begin << " ~ " << interval_end << endl;

   	unsigned int iterateNum = 0;//debug;
    //cout << "interval_begin = " << interval_begin << endl << "interval_end = " << interval_end << endl; 
   	while((c < read_length) && (queryFound == true))
    {
   	 	iterateNum++;
   	 	//cout << "iterateNum: " << iterateNum << endl;
   	 	if(iterateNum>read_length)
   	 	{
   	 			//debugln("error: interateNum > readLength");
   	 			return false;
   	 	}
   	 	unsigned int c_old = c;
			
		if(interval_begin != interval_end)
		{ 
 			lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild);
			Min = min(lcp_length, read_length);
			//cout << "lcp: " << lcp_length << endl;


			unsigned int loc_pos = 0;
            for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
            {
            	queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
            	if (!queryFound)
            	{	
            		break;
            	}
            }
            //cout << "queryFound: " << queryFound << endl;
            if(!queryFound)
            {
            	stop_loc = c_old + loc_pos;
            	(*mappedLength) = stop_loc;
            	(*indexIntervalStart) = interval_begin;
            	(*indexIntervalEnd) = interval_end;
            	return true;
            	//break;
            }
            	
            c = Min;
            if(*(read_local+c) == 'N')
            {
            	//queryFound = false;             	
            	stop_loc = c;
            	(*mappedLength) = stop_loc;
            	(*indexIntervalStart) = interval_begin;
            	(*indexIntervalEnd) = interval_end;
            	return true;
            }
			start = interval_begin; end = interval_end;
			if (c == read_length)
			{	
				(*mappedLength) = read_length;		
            	(*indexIntervalStart) = interval_begin;
            	(*indexIntervalEnd) = interval_end;
            	return true;	
				//break;			
			}	
			unsigned int interval_begin_ori = interval_begin;
			unsigned int interval_end_ori = interval_end;
		    getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
		    		chrom, verifyChild);

		    //cout << "$1 Interval: " << interval_begin << " ~ " << interval_end << endl;
		    if(interval_begin > interval_end)
		    {
		    	queryFound = false;
		    	stop_loc = c;
          		//cout << "$2 interval_begin > interval_end" << endl; cout << "stop_loc = " << stop_loc << endl; 
				(*mappedLength) = stop_loc;		
            	(*indexIntervalStart) = interval_begin_ori;
            	(*indexIntervalEnd) = interval_end_ori;
            	return true;	          		 			
		    }
		    else
		    {
		    
		    }
 		}//end if
		else 
		{
			//cout << "$3 interval_begin == interval_end" << endl;
			unsigned int loc_pos = 0;
           	for(loc_pos = 0; loc_pos < read_length - c; loc_pos++)
           	{
           		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
           		if (!queryFound)
           			break;
           	}
	    	if(queryFound) 
	    	{
	    		(*mappedLength) = read_length;
	    		(*indexIntervalStart) = interval_begin;
	    		(*indexIntervalEnd) = interval_end;
	    		return true;
	    		//align_length[101] ++;
	    	}
	    	else 
	    	{ 
	    		stop_loc = c+loc_pos;
	    		(*mappedLength) = stop_loc;
	    		(*indexIntervalStart) = interval_begin;
	    		(*indexIntervalEnd) = interval_end;
	    		return true;
	    	}	
	    }
	} //end while

	return true;
}

void build_up_down(unsigned int *lcptab, unsigned int *up, unsigned int *down, unsigned int n)
{
    unsigned int lastIndex = NULL_NUM;	
	stack<unsigned int> up_down;
	up_down.push(0);
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		while (lcptab[i] < lcptab[up_down.top()])
		{	lastIndex = up_down.top();
			up_down.pop();
			if((lcptab[i] <= lcptab[up_down.top()]) && (lcptab[up_down.top()] != lcptab[lastIndex]))
				down[up_down.top()] = lastIndex;
		}
		// now lcptab[i] >= lcptab[up_down.top()] holds
		if(lastIndex != NULL_NUM)
		{
			up[i] = lastIndex;
			lastIndex = NULL_NUM;
		}
		up_down.push(i);	
	}
	return;
}

void build_next(unsigned int *lcptab, unsigned int *next, unsigned int n)
{
	stack<unsigned int> nextIndex;
	unsigned int lastIndex;
	unsigned int j;
	nextIndex.push(0);
	for(j = 0; j < n; j++)
	{
		while(lcptab[j] < lcptab[nextIndex.top()])
			nextIndex.pop();
		if(lcptab[j] == lcptab[nextIndex.top()])
		{
			lastIndex = nextIndex.top();
			nextIndex.pop();
			next[lastIndex] = j;
		}
		nextIndex.push(j);
	}	
	return;
}

void build_lcp(unsigned int *r, unsigned int *sa, unsigned int *lcp, unsigned int *rank, unsigned int n)
{

	unsigned int i, j; 
	int k=0;
	for (i = 0; i < n; i++) rank[sa[i]] = i;
		cout << "build lcp - step 1" << endl;
	for (i = 0; i < n; lcp[rank[i++]] = k) 
	for (k?k--:0, (rank[i] == 0)?(j=0):(j=sa[rank[i]-1]); r[i+k] == r[j+k]; k++);
	lcp[0] = 0;
		cout << "build lcp - step 2" << endl; 
	return;
}

unsigned int cmp(unsigned int *r,unsigned int a,unsigned int b,unsigned int l)
{return r[a]==r[b]&&r[a+l]==r[b+l];}  

const unsigned int compare(unsigned int *r, unsigned int a, unsigned int b,
		unsigned int l) {
	return r[a] == r[b] && r[a + l] == r[b + l];
}

inline bool leq(int a1, int a2, int b1, int b2) // lexicographic order
{
	return (a1 < b1 || (a1 == b1 && a2 <= b2));
} // for pairs

inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3)
{
	return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
} // and triples

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(unsigned int* a, unsigned int* b, unsigned int* r, unsigned int n, int K)
{
	// count occurrences
	unsigned int* c = new unsigned int[K + 1]; // counter array

	for (int i = 0; i <= K; i++)
		c[i] = 0; // reset counters

	for (unsigned int i = 0; i < n; i++)
		c[r[a[i]]]++; // count occurrences

	for (int i = 0, sum = 0; i <= K; i++) // exclusive prefix sums
	{
		int t = c[i];
		c[i] = sum;
		sum += t;
	}

	for (unsigned int i = 0; i < n; i++)
		b[c[r[a[i]]]++] = a[i]; // sort

	delete[] c;
}

// find the suffix array SA of T[0..n-1] in {1..K}^n
// require T[n]=T[n+1]=T[n+2]=0, n>=2
void suffixArray(unsigned int* T, unsigned int* SA, unsigned int n, int K) {
	unsigned int n0 = (n + 2) / 3,
		n1 = (n + 1) / 3,
		n2 = n / 3,
		n02 = n0 + n2;

	unsigned int* R = new unsigned int[n02 + 3];
	R[n02] = R[n02 + 1] = R[n02 + 2] = 0;
	unsigned int* SA12 = new unsigned int[n02 + 3];
	SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
	unsigned int* R0 = new unsigned int[n0];
	unsigned int* SA0 = new unsigned int[n0];

	//******* Step 0: Construct sample ********
	// generate positions of mod 1 and mod 2 suffixes
	// the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
	for (unsigned int i = 0, j = 0; i < n + (n0 - n1); i++)
		if (i % 3 != 0)
			R[j++] = i;

	//******* Step 1: Sort sample suffixes ********
	// lsb radix sort the mod 1 and mod 2 triples
	radixPass(R, SA12, T + 2, n02, K);
	radixPass(SA12, R, T + 1, n02, K);
	radixPass(R, SA12, T, n02, K);

	// find lexicographic names of triples and
	// write them to correct places in R
	unsigned int name = 0, c0 = 0-1, c1 = 0-1, c2 = 0-1;
	for (unsigned int i = 0; i < n02; i++)
	{
		if (T[SA12[i]] != c0 || T[SA12[i] + 1] != c1 || T[SA12[i] + 2] != c2)
		{
			name++;
			c0 = T[SA12[i]];
			c1 = T[SA12[i] + 1];
			c2 = T[SA12[i] + 2];
		}

		if (SA12[i] % 3 == 1)
		{
			// write to R1
			R[SA12[i] / 3] = name;
		}
		else
		{
			// write to R2
			R[SA12[i] / 3 + n0] = name;
		}
	}

	// recurse if names are not yet unique
	if (name < n02)
	{
		suffixArray(R, SA12, n02, name);

		// store unique names in R using the suffix array
		for (unsigned int i = 0; i < n02; i++)
			R[SA12[i]] = i + 1;
	}
	else // generate the suffix array of R directly
		for (unsigned int i = 0; i < n02; i++)
			SA12[R[i] - 1] = i;

	//******* Step 2: Sort nonsample suffixes ********
	// stably sort the mod 0 suffixes from SA12 by their first character
	for (unsigned int i = 0, j = 0; i < n02; i++)
		if (SA12[i] < n0)
			R0[j++] = 3 * SA12[i];
	radixPass(R0, SA0, T, n0, K);

	//******* Step 3: Merge ********
	// merge sorted SA0 suffixes and sorted SA12 suffixes
	for (unsigned int p = 0, t = n0 - n1, k = 0; k < n; k++)
	{
		#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
		unsigned int i = GetI(); // pos of current offset 12 suffix
		unsigned int j = SA0[p]; // pos of current offset 0 suffix
		if (SA12[t] < n0 // different compares for mod 1 and mod 2 suffixes
				? leq(T[i], R[SA12[t] + n0], T[j], R[j / 3])
				: leq(T[i], T[i + 1], R[SA12[t] - n0 + 1], T[j], T[j + 1], R[j / 3 + n0]))
		{ 	// suffix from SA12 is smaller

			SA[k] = i;
			t++;
			if (t == n02) // done --- only SA0 suffixes left
				for (k++; p < n0; p++, k++)
					SA[k] = SA0[p];
		}
		else // suffix from SA0 is smaller
		{
			SA[k] = j;
			p++;
			if (p == n0) // done --- only SA12 suffixes left
				for (k++; t < n02; t++, k++)
					SA[k] = GetI();
		}
	}
	delete[] R;
	delete[] SA12;
	delete[] SA0;
	delete[] R0;
}
void da(unsigned int *r,unsigned int *sa,unsigned int n,unsigned int m)
{   
	unsigned int maxn = n;

    unsigned int *wa = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *wb = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *wv = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *ws = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int i,j,p,*x=wa,*y=wb,*t;
   
    for(i=0;i<m;i++) ws[i]=0;

    for(i=0;i<n;i++) ws[x[i]=r[i]]++; 

	for(i=1;i<m;i++) ws[i]+=ws[i-1];

    for(i=n-1;i>=0;i--) 
    {
    	sa[--ws[x[i]]]=i; 
		if(i == 0)
			break;
	}

    for(j=1,p=1;p<n;j*=2,m=p)
    {
        for(p=0,i=n-j;i<n;i++) y[p++]=i;  
        for(i=0;i<n;i++) if(sa[i]>=j) y[p++]=sa[i]-j; 
        for(i=0;i<n;i++) wv[i]=x[y[i]];  
        for(i=0;i<m;i++) ws[i]=0;
        for(i=0;i<n;i++) ws[wv[i]]++;
        for(i=1;i<m;i++) ws[i]+=ws[i-1];
        for(i=n-1;i>=0;i--) 
        {
        	sa[--ws[wv[i]]]=y[i];  
			if(i == 0)
				break;
		}
		for(t=x,x=y,y=t,p=1,x[sa[0]]=0,i=1;i<n;i++)
        x[sa[i]]=cmp(y,sa[i-1],sa[i],j)?p-1:p++; 
    }
    return;
}

//set<string> chrNameSet;

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		cout << "Executable <InputChromosomesFolder> <outputIndexFolder>" << endl;
		exit(0);
	}	
	
	BuildIndexParameter_info* buildIndexInfo = new BuildIndexParameter_info();

	string InputChrFolder = argv[1];
	string OutputIndexFolder = argv[2];

	///////////// initiate index files //////////////
	//string outputIndexFileStr = OutputIndexFolder + "/WholeGenomeIndex";
	buildIndexInfo->InputChromFolderStr = InputChrFolder;
	buildIndexInfo->OutputIndexFolderStr = OutputIndexFolder + "/";
	
	cout << "mkdir " << endl;
	string mkdir_cmd_outputIndex = "mkdir -p " + OutputIndexFolder;
	system(mkdir_cmd_outputIndex.c_str());

	string parameter_file = buildIndexInfo->OutputIndexFolderStr + "/_parameter";
	ofstream parameter_file_ofs(parameter_file.c_str());
	string chrom_file = buildIndexInfo->OutputIndexFolderStr + "/_chromMerge";
	ofstream chrom_file_ofs(chrom_file.c_str(),ios::binary);	

	string log_file = buildIndexInfo->OutputIndexFolderStr + "/_progress.log";
	ofstream log_ofs(log_file.c_str());
	//buildIndexInfo->initiateAllOutputIndexFile();
	/////////////////////////////////////////////////////////////////////////////////
	////////////// scan all chromosomes, generate index_parameter, index_chrom //////
	/////////////////////////////////////////////////////////////////////////////////
	log_ofs << endl << "start to look for chromosomes ..." << endl << endl;
	DIR *dp;    
    struct dirent *dirp;    
    if((dp=opendir(argv[1])) == NULL)    
    {
    	log_ofs << "can't open " << argv[1] << endl;
        //printf("can't open %s",argv[1]);   
    }
    log_ofs << endl << " ******************  start to look for chromosome files ...  ***************************" << endl << endl;
    while ((dirp=readdir(dp))!=NULL)    
    {     
    	//log_ofs << "read dir ..." << " dirp->d_type: " << dirp->d_type << " dirp->d_name: " << dirp->d_name << endl;
    	//if (dirp->d_type==8)
	    //{
    	string fileName = dirp->d_name;
	    if(fileName.length() < 3)
	    	continue;
	    if(fileName.substr(fileName.length()-3) == ".fa")
	    {  
	    	(buildIndexInfo->chrNameStrSet).insert((dirp->d_name));
	    	//log_ofs << "chromsome name: " << endl << dirp->d_name << endl;
    	}
    	//}
    }    
    closedir(dp);    

    log_ofs << endl << "chromFileSetNum: " << (buildIndexInfo->chrNameStrSet).size() << endl << endl;
    
    for(set<string>::iterator tmpIter = (buildIndexInfo->chrNameStrSet).begin(); 
   		 	tmpIter != (buildIndexInfo->chrNameStrSet).end(); tmpIter++)
    {   
		//log_ofs << "chromNameStr: " << (*tmpIter) << endl;
		(buildIndexInfo->chrNameStrVec).push_back((*tmpIter));
	}
	log_ofs << endl << "********************   start to read genome files ***********************" << endl << endl;
	
	for(int tmpChromNum = 0; tmpChromNum < (buildIndexInfo->chrNameStrVec).size();
		tmpChromNum ++)
	{
		string tmpChromFileName = buildIndexInfo->InputChromFolderStr 
			+ "/" + (buildIndexInfo->chrNameStrVec)[tmpChromNum];
		ifstream tmpChrInput_ifs(tmpChromFileName.c_str());
		//FILE *fp_tmpChr = fopen(tmpChromFileName.c_str(), "r");
		//log_ofs << "readFileName " << tmpChromNum + 1 << ": " << tmpChromFileName << endl;		
		string tmpChrNameLineStr;
		getline(tmpChrInput_ifs, tmpChrNameLineStr);
		string tmpInputLineStr;
		int tmpChrSeqLength = 0;
		int tmpLineLength = 0;
		while(getline(tmpChrInput_ifs, tmpInputLineStr))
		{	
			int tmpLineLength = tmpInputLineStr.length();
			//cout << "tmpLineLength: " << tmpLineLength << endl;
			tmpChrSeqLength += tmpLineLength;
		}
		(buildIndexInfo->chromLengthVec).push_back(tmpChrSeqLength);
	}

	log_ofs << endl << " ****************   getChrEndPosVec() starts ... ***********************" << endl << endl;
	buildIndexInfo->getChrEndPosVec(log_ofs);
	log_ofs << endl << " ****************   outputChromSeq() starts ...   **************************" << endl << endl; 
	buildIndexInfo->outputChromSeq(chrom_file_ofs, log_ofs);
	log_ofs << endl << " ****************   outputIndexParameter() starts ... **********************" << endl << endl;	
	buildIndexInfo->outputIndexParameter(parameter_file_ofs);

	/////////////////////////////////////////////////
	//////generate original size index
	/////////////////////////////////////////////////
	log_ofs << endl << " *****************   start to generate original size index   **********************" << endl << endl;

	unsigned int MAX = (buildIndexInfo->chrEndPosInGenomeVec)[(buildIndexInfo->chrNameStrVec).size()-1] + 1 + 1;
	
	log_ofs << "MAX: " << MAX << endl;

	string chrom_file_str = chrom_file;

	string SA_file = buildIndexInfo->OutputIndexFolderStr; SA_file.append("/_SA"); ofstream SA_file_ofs(SA_file.c_str(),ios::binary); 
	string lcp_file = buildIndexInfo->OutputIndexFolderStr; lcp_file.append("/_lcp"); ofstream lcp_file_ofs(lcp_file.c_str(),ios::binary);
	string up_file = buildIndexInfo->OutputIndexFolderStr; up_file.append("/_up"); ofstream up_file_ofs(up_file.c_str(),ios::binary);
	string down_file = buildIndexInfo->OutputIndexFolderStr; down_file.append("/_down"); ofstream down_file_ofs(down_file.c_str(),ios::binary);
	string next_file = buildIndexInfo->OutputIndexFolderStr; next_file.append("/_next"); ofstream next_file_ofs(next_file.c_str(),ios::binary);
	string chrom_bit_file = buildIndexInfo->OutputIndexFolderStr; chrom_bit_file.append("/_chrom"); ofstream chrom_bit_file_ofs(chrom_bit_file.c_str(),ios::binary);
	
  	FILE *fp = fopen(chrom_file_str.c_str(), "r");
    unsigned int *r = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	char ch;
	char head[100];
	char base[Range] = {'X','A','C','G','T','N'};
	char *chrom = (char*)malloc(MAX * sizeof(char));	

	unsigned int chrom_base_num = 0;
	while((ch = fgetc(fp)) != EOF)
	{   
		//printf("ch = %c\n",ch); 
		if((ch == 'A')||(ch == 'a')) {chrom[chrom_base_num] = 'A'; r[chrom_base_num] = 1; chrom_base_num++;}
		else if((ch == 'C')||(ch == 'c')) {chrom[chrom_base_num] = 'C'; r[chrom_base_num] = 2; chrom_base_num++;}
		else if((ch == 'G')||(ch == 'g')) {chrom[chrom_base_num] = 'G'; r[chrom_base_num] = 3; chrom_base_num++;}
		else if((ch == 'T')||(ch == 't')) {chrom[chrom_base_num] = 'T'; r[chrom_base_num] = 4; chrom_base_num++;}
		else if((ch == 'N')||(ch == 'n')) 
		{
			chrom[chrom_base_num] = 'N'; r[chrom_base_num] = 5; chrom_base_num++;
		}
		else if((ch == 'X')) {chrom[chrom_base_num] = 'X'; r[chrom_base_num] = 6; chrom_base_num++;}
		else if((ch == '\t')||(ch == '\n')) {continue;}
		else 
		{
			printf("\n found illegal input is '%c', For now, replace it with N",ch); 
			log_ofs << endl << "illegal input in provided reference Fasta file: " << ch 
				<< endl << "For now, replace it with N" << endl;
			//break;
			chrom[chrom_base_num] = 'N'; r[chrom_base_num] = 5; chrom_base_num++;
		}
	}
	log_ofs << "the number of bases in Chromo is "<< chrom_base_num << endl;
	log_ofs << "chrom is ready" << endl;
	r[MAX-1] = Range;
	chrom[MAX-1] = 'X';	

	log_ofs << "chrom[MAX-3]: " << chrom[MAX-3] << endl;	
	log_ofs << "chrom[MAX-2]: " << chrom[MAX-2] << endl;
	log_ofs << "chrom[MAX-1]: " << chrom[MAX-1] << endl;

	chrom_bit_file_ofs.write((const char*) chrom, MAX * sizeof(char));
	//free(chrom);



	
    unsigned int *sa = (unsigned int*)malloc(MAX * sizeof(unsigned int));
   
	log_ofs << "start to build SA array" << endl;

	//Xinan's SA generation algorithm
	//da(r,sa,MAX,Range+1); 

	//Kyle's SA generation algorithm
	suffixArray(r, sa, MAX, Range+1);

	log_ofs << "SA is ready" << endl;	
    //////////////////////////////////////////////////
    //unsigned int rank[MAX]={0}, lcp[MAX]={0}, up[MAX] = {0}, down[MAX] = {0}, next[MAX] = {0};  
    unsigned int *rank = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *lcp = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *up = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *down = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *next = (unsigned int*)malloc(MAX * sizeof(unsigned int));

    log_ofs << "start to build LCP array"<< endl;
	build_lcp(r, sa, lcp, rank, MAX); //build lcp array
	free(r); free(rank);


	log_ofs << "lcp is ready" << endl;	
	log_ofs << "start to build up & down array" << endl;
	build_up_down(lcp, up, down, MAX); // build up and down array
	log_ofs << "up & down are ready" << endl;
	log_ofs << "start to build next array" << endl;		
	build_next(lcp, next, MAX); //build next array
	log_ofs << "next is ready" << endl;	

	log_ofs << "start to output Index to file "<< endl;
	SA_file_ofs.write((const char*) sa, MAX * sizeof(unsigned int));
	lcp_file_ofs.write((const char*) lcp, MAX * sizeof(unsigned int));
	up_file_ofs.write((const char*) up, MAX * sizeof(unsigned int));
	down_file_ofs.write((const char*) down, MAX * sizeof(unsigned int));
	next_file_ofs.write((const char*) next, MAX * sizeof(unsigned int));
	//free(sa);
	log_ofs << "finish writing index to files " << endl;
	/////////////////////////////////////////////////////////////////////
	/////  compress index for whole genome original size index
	/////////////////////////////////////////////////////////////////////

	unsigned int indexSize = MAX;

	string childTab_file = buildIndexInfo->OutputIndexFolderStr; 
	childTab_file.append("_childTab"); ofstream childTab_file_ofs(childTab_file.c_str(),ios::binary);
	string detChild_file = buildIndexInfo->OutputIndexFolderStr; 
	detChild_file.append("_detChild"); ofstream detChild_file_ofs(detChild_file.c_str(), ios::binary);

    unsigned int *childTab;
    childTab = (unsigned int*)malloc(indexSize * sizeof(unsigned int));

   	BYTE *verifyChild;
	verifyChild = (BYTE*)malloc(indexSize * sizeof(BYTE));

	log_ofs << "start to compress original size index " << endl;
	BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();

	tmpBuildIndexInfo->compressUpDownNext2ChildtabVerifyChild(
		up, down, next, childTab, verifyChild, indexSize);
	//log_ofs << "finish compressing index" << endl;
	childTab_file_ofs.write((const char*) childTab, indexSize * sizeof(unsigned int));
	detChild_file_ofs.write((const char*) verifyChild, indexSize * sizeof(BYTE));	

	childTab_file_ofs.close();
	detChild_file_ofs.close();

	log_ofs << "finish compressing index" << endl; 
	free(up); free(down); free(next); //free(childTab); free(verifyChild);

	string lcpCompress_file = buildIndexInfo->OutputIndexFolderStr; 
	lcpCompress_file.append("_lcpCompress"); ofstream lcpCompress_file_ofs(lcpCompress_file.c_str(),ios::binary);

	//BuildIndex_Info* tmpBuildIndexInfo = new BuildIndex_Info();
	log_ofs << "start to compress Lcp array" << endl;
	BYTE *lcpCompress = (BYTE*)malloc(indexSize * sizeof(BYTE));
	tmpBuildIndexInfo->compressLcp2Lcpcompress(lcp, 
		lcpCompress, indexSize);
	log_ofs << "finish compressing Lcp array " << endl;
	lcpCompress_file_ofs.write((const char*) lcpCompress, indexSize * sizeof(BYTE));	

	free(lcp);
	//free(lcpCompress);
	lcpCompress_file_ofs.close();

	////////////////////////////////////////////////////////
	/////  start to generate preIndex String 
	////////////////////////////////////////////////////////
	
	int preIndexStringLength = PREINDEX_STRINGLENGTH;
	string generateStringFilePrefix = buildIndexInfo->OutputIndexFolderStr + "/_preIndexString";
	int stringLengthInHash = preIndexStringLength;
	
	
	log_ofs << "start to write preIndex string files" << endl;

	string baseStr[4] = {"A", "C", "G", "T"};
	string generateStringFile[stringLengthInHash];// = argv[1];
	//string generateStringFilePrefix = argv[1];
	generateStringFile[0] = generateStringFilePrefix + ".1";

	ofstream GenerateStringFile_ofs(generateStringFile[0].c_str());
	GenerateStringFile_ofs << "A" << endl << "C" << endl
		<< "G" << endl << "T"; // need to debug ...
	GenerateStringFile_ofs.close();
	string readString;
	char readChar[20];	
	for(int tmp = 1; tmp < stringLengthInHash; tmp++)
	{
		char tmpChar[2];
		sprintf(tmpChar, "%d", tmp+1);
		string tmpString = tmpChar;
		generateStringFile[tmp] = generateStringFilePrefix + "." + tmpString;

		FILE *fp_in = fopen(generateStringFile[tmp-1].c_str(), "r");
		ofstream GenerateStringFile_ofs(generateStringFile[tmp].c_str());

    	while(!feof(fp_in))
    	{
    		fgets(readChar, sizeof(readChar), fp_in);
    		readString = readChar;
    		readString = readString.substr(0, tmp);
    		for(int tmpBase = 0; tmpBase < 4; tmpBase ++)
    		{
    			string resultString = readString + baseStr[tmpBase];
    			GenerateStringFile_ofs << resultString << endl;
    		}
    	}  fclose(fp_in);
    	GenerateStringFile_ofs.close();
    	//remove(generateStringFile[tmp-1].c_str());
	}
	log_ofs << "finish writing preIndex string files" << endl;

	log_ofs << "start to generate preIndex string mapping record" << endl;
	string preIndexStringFileStr = generateStringFilePrefix + "." + int_to_str(preIndexStringLength);
	FILE *fp_in_2 = fopen(preIndexStringFileStr.c_str(), "r"); 

	string stringHashRecordFile = buildIndexInfo->OutputIndexFolderStr;
	stringHashRecordFile += "_preIndexRecord";
	ofstream StringHashRecordFile_ofs(stringHashRecordFile.c_str());

	char read[20];
	int readLength = stringLengthInHash;
	int preIndexStringNum = 0;
	while(!feof(fp_in_2))
	{
		preIndexStringNum++;
		//log_ofs << "preIndexStringNum: " << preIndexStringNum << endl;
		fgets(read, sizeof(read), fp_in_2);

		int mappedLength;
		unsigned int indexIntervalStart, indexIntervalEnd;
		bool mapMain = mapMainForIndexStringHash(read, sa, lcpCompress, childTab, chrom, &mappedLength, 
		&indexIntervalStart, &indexIntervalEnd, verifyChild, readLength, MAX);
		string readString = read;
		readString = readString.substr(0, readLength);
		if(mapMain)
		{
			StringHashRecordFile_ofs << readString << "\t" << mappedLength << "\t" << indexIntervalStart << "\t" << indexIntervalEnd << endl;
		}
		else
		{
			log_ofs << "mapMain error " << endl;
			exit(1);
		}
		//exit(1);
	}

	fclose(fp_in_2);
	StringHashRecordFile_ofs.close();
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);free(chrom);free(verifyChild);	
	log_ofs << "finish writing preIndex string mapping record !" << endl;

	log_ofs << "start to generate preIndex array" << endl;
	string preIndexFileStr = buildIndexInfo->OutputIndexFolderStr;//"/data/homes/lxauky/adSA_result/chrAll/result_0222/preIndexRecord";
	preIndexFileStr.append("_preIndexRecord");
	FILE* fp_in_preIndex = fopen(preIndexFileStr.c_str(), "r"); 

	log_ofs << "start to creat mapLength, intervalBegin, intervalEnd files" << endl;
	string preIndexArrayPreStr = buildIndexInfo->OutputIndexFolderStr;

	string preIndexMapLengthArrayStr = preIndexArrayPreStr;
	preIndexMapLengthArrayStr.append("_MapLength"); 
	ofstream preIndexMapLengthArray_ofs(preIndexMapLengthArrayStr.c_str(), ios::binary);

	string preIndexIntervalStartArrayStr = preIndexArrayPreStr;
	preIndexIntervalStartArrayStr.append("_IntervalStart");
	ofstream preIndexIntervalStartArray_ofs(preIndexIntervalStartArrayStr.c_str(), ios::binary);

	string preIndexIntervalEndArrayStr = preIndexArrayPreStr;
	preIndexIntervalEndArrayStr.append("_IntervalEnd");
	ofstream preIndexIntervalEndArray_ofs(preIndexIntervalEndArrayStr.c_str(), ios::binary);

	log_ofs << "finish creating mapLength, intervalBegin, intervalEnd files" << endl;

	int* preIndexMapLengthArray;
	preIndexMapLengthArray = (int*)malloc(PreIndexSize * sizeof(int));

	unsigned int *preIndexIntervalStartArray;
    preIndexIntervalStartArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));

	unsigned int *preIndexIntervalEndArray;
    preIndexIntervalEndArray = (unsigned int*)malloc(PreIndexSize * sizeof(unsigned int));

	char preIndexRecordChar[100];
	char preIndexStrChar[20];
	char preIndexMapLengthChar[10];
	char preIndexIntervalStartChar[20];
	char preIndexIntervalEndChar[20];

	int preIndexMapLengthInt;
	unsigned int preIndexIntervalStartInt;
	unsigned int preIndexIntervalEndInt;

	unsigned int preIndexNO = 0;

	log_ofs << "start to extract record from record file" << endl;
	while(!feof(fp_in_preIndex))
	{
		//log_ofs << "preIndexNO: " << preIndexNO;// << endl;
		fgets(preIndexRecordChar, sizeof(preIndexRecordChar), fp_in_preIndex);
		//log_ofs << "record: " << preIndexRecordChar << endl;

		sscanf(preIndexRecordChar, "%s\t%s\t%s\t%s", preIndexStrChar, preIndexMapLengthChar,
			preIndexIntervalStartChar, preIndexIntervalEndChar);
		preIndexMapLengthInt = atoi(preIndexMapLengthChar);
		preIndexIntervalStartInt = strtoul(preIndexIntervalStartChar, NULL, 10);
		preIndexIntervalEndInt = strtoul(preIndexIntervalEndChar, NULL, 10);
		
		string tmpPreIndexStr = preIndexStrChar;
		preIndexNO = getPreIndexNO(tmpPreIndexStr);

		preIndexMapLengthArray[preIndexNO] = preIndexMapLengthInt;
		preIndexIntervalStartArray[preIndexNO] = preIndexIntervalStartInt;
		preIndexIntervalEndArray[preIndexNO] = preIndexIntervalEndInt;

 	}
 	log_ofs << "finish extracting record from record file" << endl;
 	//delete(preIndexSegInfo);
 	fclose(fp_in_preIndex);

 	preIndexMapLengthArray_ofs.write((const char*) preIndexMapLengthArray, PreIndexSize * sizeof(int));
 	preIndexIntervalStartArray_ofs.write((const char*) preIndexIntervalStartArray, PreIndexSize * sizeof(unsigned int));
 	preIndexIntervalEndArray_ofs.write((const char*) preIndexIntervalEndArray, PreIndexSize * sizeof(unsigned int));

 	//cout << "mapLengthArray[0] = " << preIndexMapLengthArray[0] << endl;
 	//cout << "mapLengthArray[PreIndexSize-2] = " << preIndexMapLengthArray[PreIndexSize-2] << endl;
 	//cout << "mapLengthArray[PreIndexSize-1] = " << preIndexMapLengthArray[PreIndexSize-1] << endl;

 	//cout << "preIndexIntervalStartArray[0] = " << preIndexIntervalStartArray[0] << endl;
 	//cout << "preIndexIntervalStartArray[PreIndexSize-2] = " << preIndexIntervalStartArray[PreIndexSize-2] << endl;
 	//cout << "preIndexIntervalStartArray[PreIndexSize-1] = " << preIndexIntervalStartArray[PreIndexSize-1] << endl;

 	//cout << "preIndexIntervalEndArray[0] = " << preIndexIntervalEndArray[0] << endl;
 	//cout << "preIndexIntervalEndArray[PreIndexSize-2] = " << preIndexIntervalEndArray[PreIndexSize-2] << endl;
 	//cout << "preIndexIntervalEndArray[PreIndexSize-1] = " << preIndexIntervalEndArray[PreIndexSize-1] << endl;

 	free(preIndexMapLengthArray);
 	free(preIndexIntervalStartArray);
 	free(preIndexIntervalEndArray);

 	log_ofs << "finish generating preIndex array" << endl;
 	log_ofs << endl << "all index building jobs done" << endl;
	return 0;
}