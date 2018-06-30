// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
//used to check mappedLength's distribution of primary alignments
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>

using namespace std;
/*
int getFlagFromSamStr(const string& tmpSamStr)
{
	int tmpFlag;

	return tmpFlag;
}*/

bool mappedOrNot(int tmpFlag)
{
	if(tmpFlag & 0x4)
		return false;
	else
		return true;
}


void unmappedAlignemnt2ReadFile(const string& SAMfile)
{
	string outputFileStr_mappedAlignment = SAMfile + ".mapped.alignment.sam";
	string outputFileStr = SAMfile + ".unmapped.read";
	string outputFileStr_read_end1 = outputFileStr + ".1.fq";
	string outputFileStr_read_end2 = outputFileStr + ".2.fq"; 
	string outputFileStr_log = outputFileStr + ".log";	

	FILE *fp;
	fp = fopen(SAMfile.c_str(), "r");
	ofstream output_mappedAlignment_ofs(outputFileStr_mappedAlignment.c_str());
	ofstream output_read_end1_ofs(outputFileStr_read_end1.c_str());
	ofstream output_read_end2_ofs(outputFileStr_read_end2.c_str());
	ofstream output_log_ofs(outputFileStr_log.c_str());

	int oriAlignmentNum = 0;
	int unmappedAlignmentNum = 0;

	char tmpSamStrChar[2000];
	char read_name[100], flag[3], rname[50],pos[20],mapq[3],cigar[100],rnext[50];
	char pnext[20], tlen[10], seq[700], qual[700], others[200];

	bool end1_or_end2_bool = false;
	while(!feof(fp))
	{
		fgets(tmpSamStrChar, sizeof(tmpSamStrChar),fp);
		if(feof(fp))
			break;	
		if(tmpSamStrChar[0] == '@')
			continue;
		oriAlignmentNum ++;
		sscanf(tmpSamStrChar,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,others);  
		end1_or_end2_bool = (!end1_or_end2_bool);
		int tmpFlag = atoi(flag);
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		if(mappedOrNot_bool)
		{
			string tmpSamString = tmpSamStrChar;
			output_mappedAlignment_ofs << tmpSamString << endl;
		}
		else
		{
			unmappedAlignmentNum ++;
			string tmpReadName = read_name;
			string tmpReadSeq = seq;
			string tmpQualSeq = qual;
			if(tmpQualSeq.length() != tmpReadSeq.length())
			{
				tmpQualSeq = "I";
				for(int tmp = 0;  tmp < tmpReadSeq.length() - 1; tmp++)
				{
					tmpQualSeq += "I";
				}
			}
			if(end1_or_end2_bool)
			{
				output_read_end1_ofs << "@" << tmpReadName << endl
					<< tmpReadSeq << endl << "+" << endl << tmpQualSeq << endl; 
			}
			else
			{
				output_read_end2_ofs << "@" << tmpReadName << endl
					<< tmpReadSeq << endl << "+" << endl << tmpQualSeq << endl; 
			}
		}
	}

	output_log_ofs << "oriAlignmentNum: " << oriAlignmentNum << endl;
	output_log_ofs << "unmappedAlignmentNum: " << unmappedAlignmentNum 
		<< "\tperc: " << ((double)unmappedAlignmentNum/(double)oriAlignmentNum)*100 << endl;
	output_log_ofs << "generated read pair number: " << unmappedAlignmentNum/2 << endl;
	fclose(fp);
	output_log_ofs.close();
	output_read_end1_ofs.close();
	output_read_end2_ofs.close();
	output_mappedAlignment_ofs.close();
	return;	
}


int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout << "Exetuable inputSamFile !" << endl;
		exit(1);
	}

	string inputFileStr = argv[1];

	unmappedAlignemnt2ReadFile(inputFileStr);
	/*string outputFileStr_read_end1 = outputFileStr + ".1.fq";
	string outputFileStr_read_end2 = outputFileStr + ".2.fq"; 
	string outputFileStr_log = outputFileStr + ".log";

	FILE *fp;
	fp = fopen(inputFileStr.c_str(), "r");
	ofstream output_read_end1_ofs(outputFileStr_read_end1.c_str());
	ofstream output_read_end2_ofs(outputFileStr_read_end2.c_str());
	ofstream output_log_ofs(outputFileStr_log.c_str());

	int oriAlignmentNum = 0;
	int unmappedAlignmentNum = 0;

	char tmpSamStrChar[2000];
	char read_name[100], flag[3], rname[50],pos[20],mapq[3],cigar[100],rnext[50];
	char pnext[20], tlen[10], seq[700], qual[700], others[200];

	bool end1_or_end2_bool = false;
	while(!feof(fp))
	{
		fgets(tmpSamStrChar, sizeof(tmpSamStrChar),fp);
		if(feof(fp))
			break;	
		if(tmpSamStrChar[0] == '@')
			continue;
		oriAlignmentNum ++;
		sscanf(tmpSamStrChar,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,others);  
		end1_or_end2_bool = (!end1_or_end2_bool);
		int tmpFlag = atoi(flag);
		bool mappedOrNot_bool = mappedOrNot(tmpFlag);
		if(!mappedOrNot_bool)
		{
			unmappedAlignmentNum ++;
			string tmpReadName = read_name;
			string tmpReadSeq = seq;
			string tmpQualSeq = qual;
			if(tmpQualSeq.length() != tmpReadSeq.length())
			{
				tmpQualSeq = "I";
				for(int tmp = 0;  tmp < tmpReadSeq.length() - 1; tmp++)
				{
					tmpQualSeq += "I";
				}
			}
			if(end1_or_end2_bool)
			{
				output_read_end1_ofs << "@" << tmpReadName << endl
					<< tmpReadSeq << endl << "+" << endl << tmpQualSeq << endl; 
			}
			else
			{
				output_read_end2_ofs << "@" << tmpReadName << endl
					<< tmpReadSeq << endl << "+" << endl << tmpQualSeq << endl; 
			}
		}
	}

	output_log_ofs << "oriAlignmentNum: " << oriAlignmentNum << endl;
	output_log_ofs << "unmappedAlignmentNum: " << unmappedAlignmentNum 
		<< "\tperc: " << ((double)unmappedAlignmentNum/(double)oriAlignmentNum)*100 << endl;
	output_log_ofs << "generated read pair number: " << unmappedAlignmentNum/2 << endl;
	fclose(fp);
	output_log_ofs.close();
	output_read_end1_ofs.close();
	output_read_end2_ofs.close();*/


	return 0;
}