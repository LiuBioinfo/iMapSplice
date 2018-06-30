// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef NORMAL_TRANSCRIPT_INFO_H
#define NORMAL_TRANSCRIPT_INFO_H

#include <string>
#include <string.h>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <stdlib.h>    
#include <time.h> 
#include <stdio.h>
#include <math.h>
#include <algorithm> 
//#include "index_info.h"
//#include "SNP_set.h"

using namespace std;

bool sortExonPosOrderBool(int i, int j) { return (i < j);}

class Normal_Transcript_Info
{
private:
	string transcript_name;
	string chrom_name;
	vector<int> exon_start_pos;
	vector<int> exon_end_pos;
	int exon_num;
	string strand;

	vector< pair<int,int> > spliceSiteInTranscript; //50M1000N50M : splice site is 51

public:
	int returnMapPosInChrWithMapPosInTranscript(int mapPosInTranscript)
	{
		if(exon_num == 1)
		{
			return mapPosInTranscript + exon_start_pos[0] - 1;
		}
		for(int tmp = 0; tmp < spliceSiteInTranscript.size(); tmp++)
		{
			int tmp_splice_site = spliceSiteInTranscript[tmp].first;
			if(mapPosInTranscript < tmp_splice_site)
			{
				int tmp_dist2exonEnd = tmp_splice_site - 1 - mapPosInTranscript;
				int tmp_exon_end_pos = exon_end_pos[tmp];
				return tmp_exon_end_pos - tmp_dist2exonEnd;
			}
		}
		int last_splice_site = spliceSiteInTranscript[exon_num - 2].first;
		int last_exon_start_pos_chrom = exon_start_pos[exon_num - 1];
		return (last_exon_start_pos_chrom + mapPosInTranscript - last_splice_site);
	}

	void generateSJsiteAndSizeInAlignment(int alignmentStartPosInTranscript, int alignmentEndPosInTranscript,
		vector< pair<int,int> >& SJsiteAndSize)
	{
		for(int tmp = 0; tmp < spliceSiteInTranscript.size(); tmp++)
		{
			int tmp_splice_site = spliceSiteInTranscript[tmp].first;
			int tmp_splice_size = spliceSiteInTranscript[tmp].second;

			if((tmp_splice_site > alignmentStartPosInTranscript)&&(tmp_splice_site <= alignmentEndPosInTranscript))
			{    
				SJsiteAndSize.push_back(pair<int,int>(tmp_splice_site, tmp_splice_size));
			}
		}
	}

	void outputRefSeq(const string& outputRefFolderStr, Index_Info* indexInfo, int baseNumInEachLine)
	{
		string transcript_seq = this->getTranscriptSeq_noDir(indexInfo);
		int transcript_length = transcript_seq.length();
		string refSeq_file_name = outputRefFolderStr + "/" + transcript_name + ".fa";
		ofstream refSeq_ofs(refSeq_file_name.c_str());
		refSeq_ofs << ">" << transcript_name << endl;
		int tmp_start_pos = 0;
		int tmp_left_length = transcript_length;
		while(tmp_left_length > baseNumInEachLine)
		{
			refSeq_ofs << transcript_seq.substr(tmp_start_pos, baseNumInEachLine) << endl;
			tmp_start_pos += baseNumInEachLine;
			tmp_left_length = tmp_left_length - baseNumInEachLine;
		}
		if(tmp_left_length > 0)
			refSeq_ofs << transcript_seq.substr(tmp_start_pos) << endl;
	}

	void outputTranscriptInfo(ofstream& transcriptInfo_ofs)
	{
		transcriptInfo_ofs << transcript_name << "\t" << chrom_name << "\t" 
		<< strand << "\t" << exon_num << "\t";
		for(int tmp = 0; tmp < exon_num; tmp++)
		{
			transcriptInfo_ofs << exon_start_pos[tmp] << ",";
		}
		transcriptInfo_ofs << "\t";
		for(int tmp = 0; tmp < exon_num; tmp++)
		{
			transcriptInfo_ofs << exon_end_pos[tmp] << ",";
		}		
		transcriptInfo_ofs << endl;
	}

	string returnNormalTranscriptInfo()
	{
		string tmpStr;
		tmpStr = "transcript name: " + transcript_name + "\n"
			+ "strand: " + strand + "\n";
		tmpStr = tmpStr + "spliceSiteInTranscript: " + "\n";
		for(int tmp = 0; tmp < exon_num - 1; tmp++)
		{
			tmpStr = tmpStr + "splice: " + int_to_str(spliceSiteInTranscript[tmp].first)
				+ " " + int_to_str(spliceSiteInTranscript[tmp].second) + "\n";
		}
		return tmpStr;
	}

	void getSpliceSiteInTranscript()
	{
		vector<int> exon_length;		
		for(int tmp = 0; tmp < exon_num; tmp++)
		{
			int tmp_exon_start_pos = exon_start_pos[tmp];
			int tmp_exon_end_pos = exon_end_pos[tmp];
			int tmp_exon_length = tmp_exon_end_pos - tmp_exon_start_pos + 1;
			exon_length.push_back(tmp_exon_length);
		}
		
		int currentLocInRead = 0;
		for(int tmp = 0; tmp < exon_num-1; tmp++)
		{
			currentLocInRead += exon_length[tmp];
			int tmp_splice_size = exon_start_pos[tmp+1] - exon_end_pos[tmp] - 1;
			spliceSiteInTranscript.push_back( pair<int, int> (currentLocInRead+1, tmp_splice_size)); //50M1000N50M : splice site is 51
		}
	}

	Normal_Transcript_Info()
	{}

	bool Normal_Transcript_Info_parseBEER(const string& geneNameStr, const string& exonStartStr,
		const string& exonEndStr, const string& strandStr, const string& chrStr, Index_Info* indexInfo)
	{
		/// transcript name ////
		int equalSignLoc = geneNameStr.find("=");
		transcript_name = geneNameStr.substr(equalSignLoc + 2);

		// strand ///
		equalSignLoc = strandStr.find("=");
		strand = strandStr.substr(equalSignLoc + 2);

		// chr //
		equalSignLoc = chrStr.find("=");
		chrom_name = chrStr.substr(equalSignLoc + 2);

		// exonStart //
		equalSignLoc = exonStartStr.find("=");
		string tmp_exonStartStr = exonStartStr.substr(equalSignLoc + 2);
		int searchStartLoc_exonStart = 0;
		int commaLoc_exonStart = 0;
		while(1)
		{
			commaLoc_exonStart = tmp_exonStartStr.find(',', searchStartLoc_exonStart);
			if(commaLoc_exonStart == string::npos)
				break;
			string tmpExonStartPosStr = tmp_exonStartStr.substr(
				searchStartLoc_exonStart, commaLoc_exonStart - 1 - searchStartLoc_exonStart + 1);
			//cout << "tmpExonStartPosStr: " << endl << tmpExonStartPosStr << endl;
			int tmpExonStartPosInt = atoi(tmpExonStartPosStr.c_str())+1;
			exon_start_pos.push_back(tmpExonStartPosInt);   // initiate exon_start_pos

			searchStartLoc_exonStart = commaLoc_exonStart + 1; 
		}
		if(searchStartLoc_exonStart < tmp_exonStartStr.length())
		{
			string tmpExonStartPosStr_last = tmp_exonStartStr.substr(searchStartLoc_exonStart);
			int tmpExonStartPosInt_last = atoi(tmpExonStartPosStr_last.c_str()) + 1; // FIX ME: here
			exon_start_pos.push_back(tmpExonStartPosInt_last);
		}

		// exonEnd //
		equalSignLoc = exonEndStr.find("=");
		string tmp_exonEndStr = exonEndStr.substr(equalSignLoc + 2);
		int searchStartLoc_exonEnd  = 0;
		int commaLoc_exonEnd = 0;
		while(1)
		{
			commaLoc_exonEnd = tmp_exonEndStr.find(',', searchStartLoc_exonEnd);
			if(commaLoc_exonEnd == string::npos)
				break;
			string tmpExonEndPosStr = tmp_exonEndStr.substr(
				searchStartLoc_exonEnd, commaLoc_exonEnd - 1 - searchStartLoc_exonEnd + 1);
			//cout << "tmpExonEndPosStr: " << endl << tmpExonEndPosStr << endl;
			int tmpExonEndPosInt = atoi(tmpExonEndPosStr.c_str());
			exon_end_pos.push_back(tmpExonEndPosInt);    // initiate exon_end_pos

			searchStartLoc_exonEnd = commaLoc_exonEnd + 1;
		}	
		if(searchStartLoc_exonEnd < tmp_exonEndStr.length())
		{		
			string tmpExonEndPosStr_last = tmp_exonEndStr.substr(searchStartLoc_exonEnd);
			int tmpExonEndPosInt_last = atoi(tmpExonEndPosStr_last.c_str());
			exon_end_pos.push_back(tmpExonEndPosInt_last);
		}

		// exon_num //
		exon_num = exon_start_pos.size();

		int chrom_int = indexInfo->convertStringToInt(chrom_name);

		if((chrom_int < 0) || (chrom_int >= indexInfo->returnChromNum()))
		{
			return false;
		}
		//this->getTranscriptSeq(indexInfo);
		//transcript_length = transcript_seq.length();; // initiate transcript_length
		this->getSpliceSiteInTranscript();
		//cout << "transcript_name: " << transcript_name << endl;
		return true;		
	}

	bool Normal_Transcript_Info_parseGAF(const string& GAFentry, Index_Info* indexInfo)
	{
		//cout << "GAFentry: " << endl << GAFentry << endl;
		string chromStr;
		string exonCountStr;
		string exonStartsStr, exonEndsStr;

		vector<string> fieldStrVec;
		int searchStartLoc = 0;
		int tabLoc = 0;
		for(int tmp = 0; tmp < 10; tmp++)
		{
			tabLoc = GAFentry.find('\t', searchStartLoc);
			//cout << "tabLoc: " << tabLoc << endl;
			if(tabLoc == string::npos)
				return false;
			string tmpFieldStr = GAFentry.substr(searchStartLoc, tabLoc - 1 - searchStartLoc + 1);
			//cout << "tmpFieldStr: " << tmpFieldStr << endl;
			fieldStrVec.push_back(tmpFieldStr);
			searchStartLoc = tabLoc + 1;
		}

		transcript_name = fieldStrVec[0]; // initiate transcript_name
		chrom_name = fieldStrVec[1];   // initiate chrom_name
		strand = fieldStrVec[2]; // initiate strand
		exonCountStr = fieldStrVec[7];
		exonStartsStr = fieldStrVec[8];
		exonEndsStr = fieldStrVec[9];
		
		/*cout << "transcript_name: " << endl << transcript_name << endl 
			<< "chrom_name: " << endl << chrom_name << endl
			<< "strand: " << endl << strand << endl
			<< "exonCountStr: " << endl << exonCountStr << endl
			<< "exonStartsStr: " << endl << exonStartsStr << endl
			<< "exonEndsStr: " << endl << exonEndsStr << endl << endl;*/

		exon_num = atoi(exonCountStr.c_str());  // initiate exon_num

		int searchStartLoc_exonStart = 0;
		int commaLoc_exonStart = 0;
		for(int tmp = 0; tmp < exon_num; tmp++)
		{
			commaLoc_exonStart = exonStartsStr.find(',', searchStartLoc_exonStart);
			
			string tmpExonStartPosStr = exonStartsStr.substr(
				searchStartLoc_exonStart, commaLoc_exonStart - 1 - searchStartLoc_exonStart + 1);
			//cout << "tmpExonStartPosStr: " << endl << tmpExonStartPosStr << endl;
			int tmpExonStartPosInt = atoi(tmpExonStartPosStr.c_str())+1;
			exon_start_pos.push_back(tmpExonStartPosInt);   // initiate exon_start_pos

			searchStartLoc_exonStart = commaLoc_exonStart + 1; 
		}

		int searchStartLoc_exonEnd  = 0;
		int commaLoc_exonEnd = 0;
		for(int tmp = 0; tmp < exon_num; tmp++)
		{
			commaLoc_exonEnd = exonEndsStr.find(',', searchStartLoc_exonEnd);

			string tmpExonEndPosStr = exonEndsStr.substr(
				searchStartLoc_exonEnd, commaLoc_exonEnd - 1 - searchStartLoc_exonEnd + 1);
			//cout << "tmpExonEndPosStr: " << endl << tmpExonEndPosStr << endl;
			int tmpExonEndPosInt = atoi(tmpExonEndPosStr.c_str());
			exon_end_pos.push_back(tmpExonEndPosInt);    // initiate exon_end_pos

			searchStartLoc_exonEnd = commaLoc_exonEnd + 1;
		}		

		//transcript_start_pos = exon_start_pos[0];  // initiate transcript_start_pos
		//transcript_end_pos = exon_end_pos[exon_num - 1]; // initiate transcript_end_pos
		/// ***************************************************************** //////
		/// ***************      sort exon positions     ******************** //////
		sort(exon_start_pos.begin(), exon_start_pos.end(), sortExonPosOrderBool);
		sort(exon_end_pos.begin(), exon_end_pos.end(), sortExonPosOrderBool);		
		/// ***************************************************************** //////
		/// ***************************************************************** //////
		
		int chrom_int = indexInfo->convertStringToInt(chrom_name);
		//cout << "chrom_int: " << endl;
		if((chrom_int < 0) || (chrom_int >= indexInfo->returnChromNum()))
		{
			return false;
		}
		//this->getTranscriptSeq(indexInfo);
		//transcript_length = transcript_seq.length();; // initiate transcript_length
		this->getSpliceSiteInTranscript();
		return true;
	}

	string getTranscriptSeq_noDir(Index_Info* indexInfo)
	{
		string tmp_transcript_seq;
		//int tmp_chrom_int = indexInfo->convertStringToInt(chrom_name);
		for(int tmpExon = 0; tmpExon < exon_num; tmpExon ++)
		{
			int tmp_exon_start_pos = exon_start_pos[tmpExon];
			int tmp_exon_end_pos = exon_end_pos[tmpExon];
			tmp_transcript_seq = tmp_transcript_seq + indexInfo->getReferenceGenomeSubstr(
				chrom_name, tmp_exon_start_pos, tmp_exon_end_pos);
			//cout << "tmp_transcript_seq: " << endl << tmp_transcript_seq << endl;
		}
		return tmp_transcript_seq;
	}

	/*
	string getTranscriptSeq_withDir(Index_Info* indexInfo)
	{
		string tmp_transcript_seq;
		//int tmp_chrom_int = indexInfo->convertStringToInt(chrom_name);
		for(int tmpExon = 0; tmpExon < exon_num; tmpExon ++)
		{
			int tmp_exon_start_pos = exon_start_pos[tmpExon];
			int tmp_exon_end_pos = exon_end_pos[tmpExon];
			tmp_transcript_seq = tmp_transcript_seq + indexInfo->getReferenceGenomeSubstr(
				chrom_name, tmp_exon_start_pos, tmp_exon_end_pos);
			//cout << "tmp_transcript_seq: " << endl << tmp_transcript_seq << endl;
		}
		if(strand == "+")
		{
			return tmp_transcript_seq;
		}
		else
		{
			return getRcmSeq(tmp_transcript_seq);
		}
	}*/

	//string returnTranscriptSeq()
	//{
	//	return transcript_seq;
	//}
	//int returnTranscriptLength()
	//{
	//	return transcript_length;
	//}
	string returnTranscriptName()
	{
		return transcript_name;
	}
	string returnChromName()
	{
		return chrom_name;
	}
	int returnExonNum()
	{
		return exon_num;
	}
	/*int returnTranscriptStartPos()
	{
		return transcript_start_pos;
	}
	int returnTranscriptEndPos()
	{
		return transcript_end_pos;
	}*/
	int returnTranscriptStartPosInChr()
	{
		return exon_start_pos[0];
	}
	int returnTranscriptEndPosInChr()
	{
		return exon_end_pos[exon_num-1];
	}
	int returnExonStartPos(int tmpExonNum) // exonNum >= 1
	{
		if(( tmpExonNum < 0 )||(tmpExonNum > exon_num))
		{
			cout << "error in returnExonStartPos !";
		}
		return exon_start_pos[tmpExonNum-1];
	}
	int returnExonEndPos(int tmpExonNum) // exonNum >= 1
	{
		if(( tmpExonNum < 0 )||(tmpExonNum > exon_num))
		{
			cout << "error in returnExonEndPos !";
		}
		return exon_end_pos[tmpExonNum-1];
	}

};

#endif