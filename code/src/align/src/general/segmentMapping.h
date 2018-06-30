// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SEGMENTMAPPING_H
#define SEGMENTMAPPING_H

#include "enhanced_suffix_array_info.h"

using namespace std;


unsigned int getChildNextValue(unsigned int *child, BYTE *verifyChild, unsigned int index)
{
	return ((verifyChild[index]>2)*child[index]);
}
unsigned int getChildUpValue(unsigned int *child, BYTE *verifyChild, unsigned int index)
{
	return ((verifyChild[index-1]==1) * child[index-1]);
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

unsigned int getChildNextValue(Enhanced_Suffix_Array_Info* esaInfo, BYTE *verifyChild, unsigned int index)
{
	return ((verifyChild[index]>2) * esaInfo->returnChild(index));
}
unsigned int getChildUpValue(Enhanced_Suffix_Array_Info* esaInfo, BYTE *verifyChild, unsigned int index)
{
	return ((verifyChild[index-1]==1) * esaInfo->returnChild(index-1));
}
unsigned int getChildDownValue(Enhanced_Suffix_Array_Info* esaInfo, BYTE *verifyChild, unsigned int index)
{
	if(verifyChild[index] == 4)
		return esaInfo->returnChild(esaInfo->returnChild(index)-1);
	else if (verifyChild[index] == 2)
		return esaInfo->returnChild(index);
	else 
		return 0;
}

unsigned int getChildNextValue(Enhanced_Suffix_Array_Info* esaInfo, unsigned int index)
{
	return ((esaInfo->returnVerifyChild(index)>2) * esaInfo->returnChild(index));
}
unsigned int getChildUpValue(Enhanced_Suffix_Array_Info* esaInfo, unsigned int index)
{
	return ((esaInfo->returnVerifyChild(index-1)==1) * esaInfo->returnChild(index-1));
}
unsigned int getChildDownValue(Enhanced_Suffix_Array_Info* esaInfo, unsigned int index)
{
	if(esaInfo->returnVerifyChild(index) == 4)
		return esaInfo->returnChild(esaInfo->returnChild(index)-1);
	else if (esaInfo->returnVerifyChild(index) == 2)
		return esaInfo->returnChild(index);
	else 
		return 0;
}


void getFirstInterval(char ch, unsigned int *interval_begin, unsigned int *interval_end, 
	unsigned int* child, BYTE* verifyChild)
{
	switch(ch)
	{
		case 'C':
		{
			*interval_begin = ((verifyChild[0]>2)*child[0]);//getChildNextValue(child, verifyChild, 0);
			//*interval_end = child_next[*interval_begin] - 1; 
			*interval_end = ((verifyChild[*interval_begin]>2)*child[*interval_begin])-1;//getChildNextValue(child, verifyChild, *interval_begin) - 1;

		}
			break;
		case 'A':
		{
			*interval_begin = 0;
			//*interval_end = child_next[0] - 1;
			*interval_end = ((verifyChild[0]>2)*child[0]) - 1;
			
		}
			break;
		case 'G':
		{
			unsigned int child_next_tmp = ((verifyChild[0]>2)*child[0]);
			//cout << "child_next_tmp: " << child_next_tmp << endl;
			*interval_begin = ((verifyChild[child_next_tmp]>2)*child[child_next_tmp]);
			*interval_end = ((verifyChild[*interval_begin]>2)*child[*interval_begin]) - 1;
			
		}
			break;
		case 'T':
		default:
			unsigned int child_next_tmp = ((verifyChild[0]>2)*child[0]);
			unsigned int child_next_tmp2 = ((verifyChild[child_next_tmp]>2)*child[child_next_tmp]);
			*interval_begin = ((verifyChild[child_next_tmp2]>2)*child[child_next_tmp2]);
			*interval_end = ((verifyChild[*interval_begin]>2)*child[*interval_begin]) - 1;
			break;		
	}
}

void getFirstInterval(Enhanced_Suffix_Array_Info* esaInfo,
	char ch, unsigned int *interval_begin, unsigned int *interval_end, 
	//unsigned int* child, 
	BYTE* verifyChild)
{
	switch(ch)
	{
		case 'C':
		{
			*interval_begin = ((verifyChild[0]>2)*(esaInfo->returnChild(0)));//getChildNextValue(child, verifyChild, 0);
			//*interval_end = child_next[*interval_begin] - 1; 
			*interval_end = ((verifyChild[*interval_begin]>2)* (esaInfo->returnChild(*interval_begin)) )-1;//getChildNextValue(child, verifyChild, *interval_begin) - 1;

		}
			break;
		case 'A':
		{
			*interval_begin = 0;
			//*interval_end = child_next[0] - 1;
			*interval_end = ((verifyChild[0]>2) * (esaInfo->returnChild(0)) ) - 1;
			
		}
			break;
		case 'G':
		{
			unsigned int child_next_tmp = ((verifyChild[0]>2) * (esaInfo->returnChild(0)));
			//cout << "child_next_tmp: " << child_next_tmp << endl;
			*interval_begin = ((verifyChild[child_next_tmp]>2)* (esaInfo->returnChild(child_next_tmp)) );
			*interval_end = ((verifyChild[*interval_begin]>2)* (esaInfo->returnChild(*interval_begin)) ) - 1;
			
		}
			break;
		case 'T':
		default:
			unsigned int child_next_tmp = ((verifyChild[0]>2) * (esaInfo->returnChild(0)) );
			unsigned int child_next_tmp2 = ((verifyChild[child_next_tmp]>2) * esaInfo->returnChild(child_next_tmp) );
			*interval_begin = ((verifyChild[child_next_tmp2]>2) * (esaInfo->returnChild(child_next_tmp2)) );
			*interval_end = ((verifyChild[*interval_begin]>2) * (esaInfo->returnChild(*interval_begin))) - 1;
			break;		
	}
}

void getFirstInterval(Enhanced_Suffix_Array_Info* esaInfo,
	char ch, unsigned int *interval_begin, unsigned int *interval_end)
{
	switch(ch)
	{
		case 'C':
		{
			*interval_begin = ((esaInfo->returnVerifyChild(0)>2)*(esaInfo->returnChild(0)));//getChildNextValue(child, verifyChild, 0);
			//*interval_end = child_next[*interval_begin] - 1; 
			*interval_end = ((esaInfo->returnVerifyChild(*interval_begin)>2)* (esaInfo->returnChild(*interval_begin)) )-1;//getChildNextValue(child, verifyChild, *interval_begin) - 1;

		}
			break;
		case 'A':
		{
			*interval_begin = 0;
			//*interval_end = child_next[0] - 1;
			*interval_end = ((esaInfo->returnVerifyChild(0)>2) * (esaInfo->returnChild(0)) ) - 1;
			
		}
			break;
		case 'G':
		{
			unsigned int child_next_tmp = ((esaInfo->returnVerifyChild(0)>2) * (esaInfo->returnChild(0)));
			//cout << "child_next_tmp: " << child_next_tmp << endl;
			*interval_begin = ((esaInfo->returnVerifyChild(child_next_tmp)>2)* (esaInfo->returnChild(child_next_tmp)) );
			*interval_end = ((esaInfo->returnVerifyChild(*interval_begin)>2)* (esaInfo->returnChild(*interval_begin)) ) - 1;
			
		}
			break;
		case 'T':
		default:
			unsigned int child_next_tmp = ((esaInfo->returnVerifyChild(0)>2) * (esaInfo->returnChild(0)) );
			unsigned int child_next_tmp2 = ((esaInfo->returnVerifyChild(child_next_tmp)>2) * esaInfo->returnChild(child_next_tmp) );
			*interval_begin = ((esaInfo->returnVerifyChild(child_next_tmp2)>2) * (esaInfo->returnChild(child_next_tmp2)) );
			*interval_end = ((esaInfo->returnVerifyChild(*interval_begin)>2) * (esaInfo->returnChild(*interval_begin))) - 1;
			break;		
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

void getInterval(Enhanced_Suffix_Array_Info* esaInfo,
	unsigned int start, unsigned int end, 
	unsigned int position, char ch, 
	unsigned int *interval_begin, unsigned int *interval_end, 
	//unsigned int* sa, 
	//unsigned int* child,
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
			unsigned int child_up_value_tmp = ((verifyChild[end]==1) * (esaInfo->returnChild(end)) );
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(esaInfo, verifyChild, start);

	
			if(chrom[(esaInfo->returnSA(start))+position] == 'C')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else if(chrom[ (esaInfo->returnSA(index_begin)) + position] == 'C')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((verifyChild[index_begin]>2) * (esaInfo->returnChild(index_begin)) ) - 1; 
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
			unsigned int child_up_value_tmp = ((verifyChild[end]==1) * (esaInfo->returnChild(end)) );
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(esaInfo, verifyChild, start);

			if(chrom[(esaInfo->returnSA(start))+position] == 'A')
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
				unsigned int child_up_value_tmp = ((verifyChild[end]==1) * (esaInfo->returnChild(end)));
				if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
					index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
				else index_begin = getChildDownValue(esaInfo, verifyChild, start);

				//pos = child_next[index_begin];
				pos = ((verifyChild[index_begin]>2) * (esaInfo->returnChild(index_begin)));
				
				//cout << "case 'G' !" << endl;
				//cout << "pos = " << pos << endl;
				//cout << "start = " << start << endl;
				//cout << "index_begin = " << index_begin << endl;


				if(chrom[esaInfo->returnSA(start)+position] == 'G')
				{
					//cout << " (chrom[sa[start]+position] == 'G') " << endl;
					*interval_begin = start;
					*interval_end = index_begin - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}	
				else if(chrom[esaInfo->returnSA(index_begin)+position] == 'G')
				{
					//cout << " (chrom[sa[index_begin]+position] == 'G') " << endl;
					*interval_begin = index_begin;
					//*interval_end = child_next[index_begin]-1;
					*interval_end = ((verifyChild[index_begin]>2) * (esaInfo->returnChild(index_begin))) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else if(chrom[esaInfo->returnSA(pos)+position] == 'G')
				{
					//cout << " (chrom[sa[pos]+position] == 'G') " << endl;
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((verifyChild[pos]>2) * (esaInfo->returnChild(pos))) - 1;
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
			unsigned int child_up_value_tmp = ((verifyChild[end]==1) * (esaInfo->returnChild(end)));
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(esaInfo, verifyChild, start);

			if(chrom[esaInfo->returnSA(start)+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[esaInfo->returnSA(index_begin)+position] == 'T')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((verifyChild[index_begin]>2) * (esaInfo->returnChild(index_begin))) - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else 
			{
				//pos = child_next[index_begin];
				pos = ((verifyChild[index_begin]>2) * esaInfo->returnChild(index_begin));
				if(chrom[esaInfo->returnSA(pos)+position] == 'T')
				{
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((verifyChild[pos]>2) * esaInfo->returnChild(pos)) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					//pos = child_next[pos];
					pos = ((verifyChild[pos]>2) * esaInfo->returnChild(pos));
					if(chrom[esaInfo->returnSA(pos)+position] == 'T')
					{
						*interval_begin = pos;
						//*interval_end = child_next[pos] - 1;
						*interval_end = ((verifyChild[pos]>2) * esaInfo->returnChild(pos)) - 1;
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

void getInterval(Enhanced_Suffix_Array_Info* esaInfo,
	unsigned int start, unsigned int end, 
	unsigned int position, char ch, 
	unsigned int *interval_begin, unsigned int *interval_end, 
	char* chrom)
{   
	/*#ifdef CAL_TIME
	getInterval_begin = clock();
	#endif*/
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
			unsigned int child_up_value_tmp = ((esaInfo->returnVerifyChild(end)==1) * (esaInfo->returnChild(end)) );
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(esaInfo, start);

	
			if(chrom[(esaInfo->returnSA(start))+position] == 'C')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else if(chrom[ (esaInfo->returnSA(index_begin)) + position] == 'C')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((esaInfo->returnVerifyChild(index_begin)>2) * (esaInfo->returnChild(index_begin)) ) - 1; 
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
			unsigned int child_up_value_tmp = ((esaInfo->returnVerifyChild(end)==1) * (esaInfo->returnChild(end)) );
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(esaInfo, start);

			if(chrom[(esaInfo->returnSA(start))+position] == 'A')
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
				unsigned int child_up_value_tmp = ((esaInfo->returnVerifyChild(end)==1) * (esaInfo->returnChild(end)));
				if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
					index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
				else index_begin = getChildDownValue(esaInfo, start);

				//pos = child_next[index_begin];
				pos = ((esaInfo->returnVerifyChild(index_begin)>2) * (esaInfo->returnChild(index_begin)));

				if(chrom[esaInfo->returnSA(start)+position] == 'G')
				{
					//cout << " (chrom[sa[start]+position] == 'G') " << endl;
					*interval_begin = start;
					*interval_end = index_begin - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}	
				else if(chrom[esaInfo->returnSA(index_begin)+position] == 'G')
				{
					//cout << " (chrom[sa[index_begin]+position] == 'G') " << endl;
					*interval_begin = index_begin;
					//*interval_end = child_next[index_begin]-1;
					*interval_end = ((esaInfo->returnVerifyChild(index_begin)>2) * (esaInfo->returnChild(index_begin))) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else if(chrom[esaInfo->returnSA(pos)+position] == 'G')
				{
					//cout << " (chrom[sa[pos]+position] == 'G') " << endl;
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((esaInfo->returnVerifyChild(pos)>2) * (esaInfo->returnChild(pos))) - 1;
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
			unsigned int child_up_value_tmp = ((esaInfo->returnVerifyChild(end)==1) * (esaInfo->returnChild(end)));
			if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
				index_begin = child_up_value_tmp;//getChildUpValue(child, verifyChild, end+1);
			else index_begin = getChildDownValue(esaInfo, start);

			if(chrom[esaInfo->returnSA(start)+position] == 'T')
			{
				*interval_begin = start;
				*interval_end = index_begin - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}	
			else if(chrom[esaInfo->returnSA(index_begin)+position] == 'T')
			{
				*interval_begin = index_begin;
				//*interval_end = child_next[index_begin]-1;
				*interval_end = ((esaInfo->returnVerifyChild(index_begin)>2) * (esaInfo->returnChild(index_begin))) - 1;
				if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
			}		
			else 
			{
				//pos = child_next[index_begin];
				pos = ((esaInfo->returnVerifyChild(index_begin)>2) * esaInfo->returnChild(index_begin));
				if(chrom[esaInfo->returnSA(pos)+position] == 'T')
				{
					*interval_begin = pos;
					//*interval_end = child_next[pos] - 1;
					*interval_end = ((esaInfo->returnVerifyChild(pos)>2) * esaInfo->returnChild(pos)) - 1;
					if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
				}
				else
				{
					//pos = child_next[pos];
					pos = ((esaInfo->returnVerifyChild(pos)>2) * esaInfo->returnChild(pos));
					if(chrom[esaInfo->returnSA(pos)+position] == 'T')
					{
						*interval_begin = pos;
						//*interval_end = child_next[pos] - 1;
						*interval_end = ((esaInfo->returnVerifyChild(pos)>2) * esaInfo->returnChild(pos)) - 1;
						if ((*interval_end < start) || (*interval_end > end)) *interval_end = end;
					}
				}						
			}
		}	
			break;

	}		
	/*#ifdef CAL_TIME
	getInterval_end = clock();
	getInterval_cost = getInterval_cost + getInterval_end - getInterval_begin;
	#endif*/
}


unsigned int getlcp(unsigned int start, unsigned int end, BYTE* lcpCompress, unsigned int* child, BYTE* verifyChild)
{
	//getLcpBegin = clock();
	//unsigned int lcp_length;
	unsigned int tmpIndex;
	
	//if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
    //    tmpIndex = child_up[end+1];
	//else 
	//	tmpIndex = child_down[start];
	unsigned int child_up_value_tmp = ((verifyChild[end]==1)*child[end]);
	if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
		tmpIndex = child_up_value_tmp;
	else
		tmpIndex = getChildDownValue(child, verifyChild, start);

	//lcp_length = (unsigned int)lcpCompress[tmpIndex];
	//lcp_length = lcp[tmpIndex];

	/*if(lcp_length == LONG_LCP)
	{
		debugln("long lcp!");
		longLcpMapIter = longLcpMap.find(tmpIndex);
		if(longLcpMapIter == longLcpMap.end())
			lcp_length = LONG_LCP;
		else
		{
			debugln("longLcp>255");
			lcp_length = longLcpMapIter -> second; 
		}
	}*/
	//debugln("lcp_length = "); debugln(lcp_length);
	//getLcpEnd = clock();
	//getLcpCost = getLcpCost + getLcpEnd - getLcpBegin;
	return lcpCompress[tmpIndex];
}

unsigned int getlcp(Enhanced_Suffix_Array_Info* esaInfo, unsigned int start, unsigned int end, BYTE* lcpCompress, BYTE* verifyChild)
{
	//getLcpBegin = clock();
	//unsigned int lcp_length;
	unsigned int tmpIndex;
	
	//if ((start < child_up[end + 1]) && (end >= child_up[end +1]))
    //    tmpIndex = child_up[end+1];
	//else 
	//	tmpIndex = child_down[start];
	unsigned int child_up_value_tmp = ((verifyChild[end]==1) * esaInfo->returnChild(end));
	if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
		tmpIndex = child_up_value_tmp;
	else
		tmpIndex = getChildDownValue(esaInfo, verifyChild, start);

	return lcpCompress[tmpIndex];
}

unsigned int getlcp(Enhanced_Suffix_Array_Info* esaInfo, unsigned int start, unsigned int end)
{
	/*#ifdef CAL_TIME
	getLcp_begin_1 = clock();
	#endif*/
	unsigned int tmpIndex;
	
	unsigned int child_up_value_tmp = ((esaInfo->returnVerifyChild(end)==1) * esaInfo->returnChild(end));

	if((start < child_up_value_tmp) && (end >= child_up_value_tmp))
		tmpIndex = child_up_value_tmp;
	else
		tmpIndex = getChildDownValue(esaInfo, start);

	int tmpLcp = esaInfo->returnLcpCompress(tmpIndex);
	/*#ifdef CAL_TIME
	getLcp_end_1 = clock();
	getLcp_cost_1 = getLcp_cost_1 + getLcp_end_1 - getLcp_begin_1;
	#endif*/	
	return tmpLcp;
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


bool mapMainSecondLevelForTargetMapping_compressedIndex(char *read, 
	unsigned int* sa, 
	BYTE* lcpCompress, unsigned int* child, BYTE* verifyChild,
	char* chrom, int readLength, unsigned int indexSize, 
	unsigned int midPartMapPosForLongHead, 
	unsigned int midPartMapPosForLongHeadInSecondLevelIndex,
	int* targetMappingAlignNum, unsigned int* targetMappingAlignLoc, 
	Index_Info* indexInfo)
{
	//return false;
	//cout << "start target mapping " << endl;
	//cout << "read: " << endl;
	//cout << read << endl;
    unsigned int norSegmentNum = 0;  	
    //cout << "*norSegmentNum = " << endl; cout << (*norSegmentNum) << endl;
   	unsigned int norSegmentLength[SEGMENTNUM];// = (unsigned int*)malloc(SEGMENTNUM * sizeof(unsigned int));
   	unsigned int norSegmentLocInRead[SEGMENTNUM];// = (unsigned int*)malloc(SEGMENTNUM * sizeof(unsigned int));
   	unsigned int norSegmentAlignNum[SEGMENTNUM];// = (unsigned int*)malloc(SEGMENTNUM * sizeof(unsigned int));
   	//Denote the start location of an alignment for a segment, start from 1;
   	unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC];// = (unsigned int*)malloc(SEGMENTNUM * CANDALILOC * sizeof(unsigned int)); 
   	//cout << "finish to set space " << endl;
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
	unsigned int n = indexInfo->returnIndexSize();//size of SA
	unsigned int norAlignLoc;
	unsigned int align_chr_location;
	unsigned int align_chr;
	char* read_local = read;
	//return false;
	//cout << "start to check every base......" << endl;
	while (stop_loc_overall < read_length) //- 15)
	{
		//cout << "stop_loc_overall = " << stop_loc_overall << endl;
		segment_num++;
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
   	 	
   	 	getFirstInterval(*read_local, &interval_begin, &interval_end, child, verifyChild);   	 	

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

 				//lcp_length = getlcp(interval_begin, interval_end, lcp, up, down);
	 			lcp_length = getlcp(interval_begin, interval_end, lcpCompress, child, verifyChild);

				Min = min(lcp_length, read_length - stop_loc_overall);
				
				unsigned int loc_pos = 0;
            	for(loc_pos = 0; loc_pos < Min - c_old; loc_pos++)
            	{
            		queryFound = (*(read_local+c_old+loc_pos) == *(chrom+sa[interval_begin]+c_old+loc_pos));
            		if (!queryFound)
            		{	
            			return false;
            			break;
            		}
            	}
            	
            	c = Min;
            	if(*(read_local+c) == 'N')
            	{
            		return false;
            	}
				start = interval_begin; end = interval_end;
				if (c + stop_loc_overall == read_length)
				{				
					break;			
				}	
				//cout << "to get interval" << endl;
				unsigned int interval_begin_ori = interval_begin;
				unsigned int interval_end_ori = interval_end;
	    	
			   	getInterval(start, end, c, *(read_local+c), &interval_begin, &interval_end, sa, child,//child_up, child_down, child_next, 
			   		chrom, verifyChild);		    	

		    	if(interval_begin > interval_end)
		    	{
		    		return false;
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
				//cout << "interval_begin == interval_end " << endl;
				unsigned int loc_pos = 0;
            	for(loc_pos = 0; loc_pos < read_length - c - stop_loc_overall; loc_pos++)
            	{
            		queryFound = (*(read_local+c+loc_pos) == *(chrom+sa[interval_begin]+c+loc_pos));
            		if (!queryFound)
            			return false;
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
    		unsigned int tmpSegLength = read_length - stop_loc_overall;
    		if((norSegmentNum > 1)||(tmpSegLength < read_length))
    		{	
    			return false;
    		}
    		norSegmentLength[0] = tmpSegLength;//READ_LENGTH - stop_loc_overall;
    		norSegmentLocInRead[0] = stop_loc_overall + 1;
    		norSegmentAlignNum[0] = segment_align_rangeNum;				
			
			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
			{    			
    			norSegmentAlignLoc[alignment_num] 
    				= sa[segment_align_SArange[0] + alignment_num] + 1;
    				//- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead;
    		}

			segment_length = read_length-stop_loc_overall;
			break;
		}
		else 
		{    
			(norSegmentNum) ++;
    		if((norSegmentNum > 1)||(stop_loc < read_length))
    		{	
    			return false;
    		}			
			norSegmentLength[0] = stop_loc;
			norSegmentLocInRead[0] = stop_loc_overall + 1;
			norSegmentAlignNum[0] = segment_align_rangeNum;

			for (unsigned int alignment_num = 0; alignment_num < min(segment_align_rangeNum,100); alignment_num++)
		    {    			
    			norSegmentAlignLoc[alignment_num]
    				= sa[segment_align_SArange[0] + alignment_num] + 1; 
    				//- midPartMapPosForLongHeadInSecondLevelIndex + midPartMapPosForLongHead;
    		}

			unsigned int stop_loc_overall_ori = stop_loc_overall;
			read_local = read_local + stop_loc + 1; // if read-align failed at some location, then restart from that location //align_length[c] ++;
			stop_loc_overall = stop_loc_overall + stop_loc + 1;
    		segment_length = stop_loc;	
		}		
   	}
	//return false;
	//cout << "end of searching for segments" << endl;
	mapMain = (((norSegmentNum) == 1) && (norSegmentLength[0]==readLength) && ((norSegmentAlignNum[0])<40));
	//cout << "mapMain: " << mapMain << endl;
	//return false;

	if(mapMain)
	{
		*targetMappingAlignNum = (int)(norSegmentAlignNum[0]);
		//cout << "targetMappingAlignNum: " << (*targetMappingAlignNum) << endl;
		for(int tmp = 0; tmp < (*targetMappingAlignNum); tmp++)
		{
			targetMappingAlignLoc[tmp] = norSegmentAlignLoc[tmp] + midPartMapPosForLongHead - midPartMapPosForLongHeadInSecondLevelIndex;
		}
	}
	//cout << "start to free memory" << endl;
   	//cout << "finish all memory free" << endl;
	return mapMain;
}


	bool targetMapWithoutPreIndexArray(
		char *read, unsigned int* sa, BYTE* lcpCompress, unsigned int* child, char* chrom, 
		BYTE* verifyChild, int readLength, Index_Info* indexInfo, int& targetMapAlignNum, 
		vector<unsigned int>& targetMapAlignLocVec)
	{
		unsigned int segmentNum;
		unsigned int norSegmentLength[SEGMENTNUM];
		unsigned int norSegmentLocInRead[SEGMENTNUM];
		unsigned int norSegmentAlignNum[SEGMENTNUM];
		unsigned int norSegmentAlignLoc[SEGMENTNUM * CANDALILOC];

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
			//// target mapping, segmentMap === 1. 
			if(segment_num > 1)
				return false;
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
				//// target mapping, segmentMap === 1. 
				if(norSegmentNum > 1)
					return false;
	   	 		norSegmentLength[norSegmentNum - 1] = 1;
				norSegmentLocInRead[norSegmentNum - 1] = stop_loc_overall + 1;
				norSegmentAlignNum[norSegmentNum - 1] = 0;				
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
			///////////////////////////////////////////////////////////////////////////////////////////////   	 
			/////////////////////////////////////////SEGMENT MAP RESULT////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////

	    	if (queryFound && (interval_end >= interval_begin)) 
	    	{
	    		(norSegmentNum) ++;
				//// target mapping, segmentMap === 1. 
				if(norSegmentNum > 1)
					return false;	    		
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
				//// target mapping, segmentMap === 1. 
				if(norSegmentNum > 1)
					return false;
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
		mapMain = (((norSegmentNum) == 1) && (norSegmentLength[0]==readLength) && ((norSegmentAlignNum[0]) <= CANDALILOC));
		if(mapMain)
		{
			targetMapAlignNum = norSegmentAlignNum[0];
			for(int tmp = 0; tmp < targetMapAlignNum; tmp++)
				targetMapAlignLocVec.push_back(norSegmentAlignLoc[tmp]);
			return true;
		}
	}	

#endif
