// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

//#define stringLengthInHash 5
//#define MAX 2654911539 
#define CANDALILOC 100
#define SEGMENTNUM 20
#define minValSegLength 20

using namespace std;

typedef unsigned char BYTE;

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

void getInterval(unsigned int start, unsigned int end, unsigned int position, char ch, unsigned int *interval_begin, unsigned int *interval_end, 
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

		    //cout << "firstInterval: " << interval_begin << " ~ " << interval_end << endl;
		    if(interval_begin > interval_end)
		    {
		    	queryFound = false;
		    	stop_loc = c-1;
          		//cout << "interval_begin > interval_end" << endl; cout << "stop_loc = " << stop_loc << endl; 
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

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "Executable <preIndexStringPrefix> <preIndexRecordPrefix> <indexPrefix> <preIndexStringLength> <IndexSize>" << endl;
		exit(0);
	}
	string indexHashStringLengthStr = argv[4];
	int stringLengthInHash = atoi(indexHashStringLengthStr.c_str());

	string indexSizeStr = argv[5];
	unsigned int indexSize = strtoul(indexSizeStr.c_str(), NULL, 10);
	unsigned int MAX = indexSize;

	string preIndexStringPrefix = argv[1];
	string preIndexStringFileStr = preIndexStringPrefix + "_preIndexString." + indexHashStringLengthStr; 

	string SA_file = argv[3]; SA_file.append("_SA"); ifstream SA_file_ifs(SA_file.c_str(),ios::binary); 
	//CompressIndex
	string lcpCompress_file = argv[3]; lcpCompress_file.append("_lcpCompress"); ifstream lcpCompress_file_ifs(lcpCompress_file.c_str(),ios::binary);
	string childTab_file = argv[3]; childTab_file.append("_childTab"); ifstream childTab_file_ifs(childTab_file.c_str(),ios::binary);
	string verifyChild_file = argv[3]; verifyChild_file.append("_detChild"); ifstream verifyChild_file_ifs(verifyChild_file.c_str(),ios::binary);
	string chrom_bit_file = argv[3]; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	cout << "start to load whole genome" << endl;
	char *chrom;
	chrom = (char*)malloc(MAX * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, MAX * sizeof(char)); 

	cout << "start to load SA" << endl;
    unsigned int *sa;
    sa = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    SA_file_ifs.read((char*)sa, MAX * sizeof(unsigned int));

	cout << "start to load lcpCompress" << endl;
	BYTE *lcpCompress;
	lcpCompress = (BYTE*)malloc(MAX * sizeof(BYTE));
	lcpCompress_file_ifs.read((char*)lcpCompress, MAX * sizeof(BYTE));	

	cout << "start to load childTab " << endl;
	unsigned int *childTab;
	childTab = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	childTab_file_ifs.read((char*)childTab, MAX * sizeof(unsigned int));

	cout << "start to load detChild" << endl;
	BYTE *verifyChild;
	verifyChild = (BYTE*)malloc(MAX * sizeof(BYTE));
	verifyChild_file_ifs.read((char*)verifyChild, MAX * sizeof(BYTE));
	
	cout << "All index files loaded" << endl;

	cout << "start to generate indexStringHashRecord ... " << endl;
	FILE *fp_in_2 = fopen(preIndexStringFileStr.c_str(), "r"); 

	string stringHashRecordFile = argv[2];
	stringHashRecordFile += "_preIndexRecord";
	ofstream StringHashRecordFile_ofs(stringHashRecordFile.c_str());

	char read[20];
	int readLength = stringLengthInHash;
	int preIndexStringNum = 0;
	while(!feof(fp_in_2))
	{
		preIndexStringNum++;
		//cout << "preIndexStringNum: " << preIndexStringNum << endl;
		/*if(preIndexStringNum > 2)
		{
			break;
		}*/
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
			cout << "mapMain error " << endl;
		}
	}
	fclose(fp_in_2);
	StringHashRecordFile_ofs.close();
	free(sa);free(lcpCompress);//free(child_up);free(child_down);free(child_next);
	free(childTab);free(chrom);	

	cout << "all done !" << endl;
	return 0;
	//return 0;
}