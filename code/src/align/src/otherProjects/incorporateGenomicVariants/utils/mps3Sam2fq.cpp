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
#include <sstream>

using namespace std;

string covertCharToReverseComplement(const string& Ori_Char)
{
	if(Ori_Char == "A")
	{
		return "T";
	}
	else if(Ori_Char == "T")
	{
		return "A";
	}
	else if(Ori_Char == "G")
	{
		return "C";
	}
	else if(Ori_Char == "C")
	{
		return "G";
	}
	else if(Ori_Char == "N")
	{
		return "N";
	}
	else
	{
		cout << "incorrect Ori_Char in covertCharToReverseComplement" << endl;
		exit(1);
		return "X";
	}
}

string convertStringToReverseComplement(const string& originalString)
{
	//cout << "convertStringToReverseComplement starts ..." << endl;
	int stringLength = originalString.size();
	//cout << "stringLength: " << stringLength << endl;
	string resultString = covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + covertCharToReverseComplement(
			originalString.substr(stringLength-1-tmp, 1));
	}
	//cout << "resultString: " << resultString << endl;
	return resultString;
}

string convertQualityScoreString2Reverse(const string& originalQualityScoreString)
{
	int stringLength = originalQualityScoreString.size();
	//cout << "qualSeqLength: " << stringLength << endl;
	string resultString = originalQualityScoreString.substr(stringLength-1, 1);//covertCharToReverseComplement(originalString.substr(stringLength-1, 1));
	for (int tmp = 1; tmp < stringLength; tmp++)
	{
		resultString = resultString + originalQualityScoreString.substr(stringLength-1-tmp, 1);
			//covertCharToReverseComplement(originalString.substr(stringLength-1-tmp, 1));
	}
	//cout << "qualSeq: " << resultString << endl;
	return resultString;
}

void getFlagSeqQualFromSam(string& sam_1, string& sam_2,
	int& sam_flagInt_1, int& sam_flagInt_2, string& sam_id_1, string& sam_id_2, 
	string& sam_seq_1, string& sam_seq_2, string& sam_qual_1, string& sam_qual_2)
{
		vector<string> samFieldVec_1;
		vector<string> samFieldVec_2;
		int startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = sam_1.find("\t", startLoc);
			string tmpSamField = sam_1.substr(startLoc, tabLoc-startLoc);
			samFieldVec_1.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}	
		startLoc = 0;
		for(int tmp = 0; tmp < 11; tmp++)
		{
			int tabLoc = sam_2.find("\t", startLoc);
			string tmpSamField = sam_2.substr(startLoc, tabLoc-startLoc);
			samFieldVec_2.push_back(tmpSamField);
			startLoc = tabLoc + 1;
		}
		sam_id_1 = samFieldVec_1[0];
		sam_id_2 = samFieldVec_2[0];
		string sam_flagInt_1_str = samFieldVec_1[1];
		string sam_flagInt_2_str = samFieldVec_2[1];
		sam_flagInt_1 = atoi(sam_flagInt_1_str.c_str());
		sam_flagInt_2 = atoi(sam_flagInt_2_str.c_str());
		sam_seq_1 = samFieldVec_1[9];
		sam_seq_2 = samFieldVec_2[9];
		sam_qual_1 = samFieldVec_1[10];
		sam_qual_2 = samFieldVec_2[10];
}

bool getFastqStr(
	int& sam_flagInt_1, int& sam_flagInt_2, string& sam_id_1, string& sam_id_2, 
	string& sam_seq_1, string& sam_seq_2, string& sam_qual_1, string& sam_qual_2,

	string& fq_id_end1, string& fq_seq_end1, string& fq_com_end1, string& fq_qual_end1,
	string& fq_id_end2, string& fq_seq_end2, string& fq_com_end2, string& fq_qual_end2)
{
	//cout << "getFastqStr: " << getFastqStr << endl;
	fq_com_end1 = "+";
	fq_com_end2 = "+";
	bool read1_end1_or_end2_bool = (sam_flagInt_1&0x40);
	bool read1_end2_or_end1_bool = (sam_flagInt_1&0x80);
	bool read2_end1_or_end2_bool = (sam_flagInt_2&0x40);
	bool read2_end2_or_end1_bool = (sam_flagInt_2&0x80);	

	if((read1_end1_or_end2_bool && read1_end2_or_end1_bool)
		||((!read1_end1_or_end2_bool)&&(!read1_end2_or_end1_bool))
		||(read1_end1_or_end2_bool && read2_end1_or_end2_bool)
		||((!read1_end1_or_end2_bool) && (!read2_end1_or_end2_bool)))
	{
		cout << "conflict between read1_end1_or_end2_bool and read1_end2_or_end1_bool!" << endl;
		cout << "read1_end1_or_end2_bool: " << read1_end1_or_end2_bool << endl;
		cout << "read1_end2_or_end1_bool: " << read1_end2_or_end1_bool << endl;
		cout << "read2_end1_or_end2_bool: " << read2_end1_or_end2_bool << endl;
		cout << "read2_end2_or_end1_bool: " << read2_end2_or_end1_bool << endl;		
		return false; 
	}
	//cout << "read1_end1_or_end2_bool: " << read1_end1_or_end2_bool << endl;
	if(read1_end1_or_end2_bool) // read1 -- end1, read2 -- end2
	{
		// read1 -- end1
		fq_id_end1 = sam_id_1;
		if(sam_flagInt_1 & 0x4)// unmapped
		{
			fq_seq_end1 = sam_seq_1;
			fq_qual_end1 = sam_qual_1;
		}
		else // mapped
		{
			if(!(sam_flagInt_1 & 0x10)) // forward
			{
				//cout << "forward!" << endl; 
				fq_seq_end1 = sam_seq_1;
				fq_qual_end1 = sam_qual_1;				
			}
			else // reverse
			{
				//cout << "reverse!" << endl;
				fq_seq_end1 = convertStringToReverseComplement(sam_seq_1);
				fq_qual_end1 = convertQualityScoreString2Reverse(sam_qual_1);
			}
		}
	
		// read2 -- end2
		fq_id_end2 = sam_id_2;
		if(sam_flagInt_2 & 0x4)
		{
			fq_seq_end2 = sam_seq_2;
			fq_qual_end2 = sam_qual_2;
		}
		else // mapped
		{
			if(!(sam_flagInt_2 & 0x10)) // forward
			{
				//cout << "forward!" << endl;
				fq_seq_end2 = sam_seq_2;
				fq_qual_end2 = sam_qual_2;				
			}
			else // reverse
			{
				//cout << "reverse!" << endl;
				fq_seq_end2 = convertStringToReverseComplement(sam_seq_2);
				fq_qual_end2 = convertQualityScoreString2Reverse(sam_qual_2);
			}
		}
		// cout << "fq_seq_end1: " << fq_seq_end1 << endl;
		// cout << "fq_seq_end2: " << fq_seq_end2 << endl;
		// cout << "fq_qual_end1: " << fq_qual_end1 << endl;
		// cout << "fq_qual_end2: " << fq_qual_end2 << endl;		
	}
	else // read2 -- end1, read1 -- end2
	{
		// read2 -- end1
		fq_id_end1 = sam_id_2;
		if(sam_flagInt_2 & 0x4)// unmapped
		{
			fq_seq_end1 = sam_seq_2;
			fq_qual_end1 = sam_qual_2;
		}
		else // mapped
		{
			if(!(sam_flagInt_2 & 0x10)) // forward
			{
				fq_seq_end1 = sam_seq_2;
				fq_qual_end1 = sam_qual_2;				
			}
			else // reverse
			{
				fq_seq_end1 = convertStringToReverseComplement(sam_seq_2);
				fq_qual_end1 = convertQualityScoreString2Reverse(sam_qual_2);
			}
		}
	
		// read1 -- end2
		fq_id_end2 = sam_id_1;
		if(sam_flagInt_1 & 0x4)
		{
			fq_seq_end2 = sam_seq_1;
			fq_qual_end2 = sam_qual_1;
		}
		else // mapped
		{
			if(!(sam_flagInt_1 & 0x10)) // forward
			{
				fq_seq_end2 = sam_seq_1;
				fq_qual_end2 = sam_qual_1;				
			}
			else // reverse
			{
				fq_seq_end2 = convertStringToReverseComplement(sam_seq_1);
				fq_qual_end2 = convertQualityScoreString2Reverse(sam_qual_1);
			}
		}	
	}
	return true;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "#0 Executable" << endl;
		cout << "#1 inputSam" << endl;
		cout << "#2 outputReadFilePrefix" << endl;
		exit(1);
	}
	string inputSam = argv[1];
	string outputReadFilePrefix = argv[2];
	string outputReadFile_end1 = outputReadFilePrefix + ".1.fq";
	string outputReadFile_end2 = outputReadFilePrefix + ".2.fq";
	string outputReadFile_log = outputReadFilePrefix + ".log";
	ifstream sam_ifs(inputSam.c_str());
	ofstream fq1_ofs(outputReadFile_end1.c_str());
	ofstream fq2_ofs(outputReadFile_end2.c_str());
	ofstream log_ofs(outputReadFile_log.c_str());
	while(!sam_ifs.eof())
	{
			string tmpAlign_1, tmpAlign_2;
			getline(sam_ifs, tmpAlign_1);
			//cout << "tmpAlign_1: " << tmpAlign_1 << endl;
			if(tmpAlign_1.at(0) == '@')
				continue;
			if(tmpAlign_1 == "")
				break;
			getline(sam_ifs, tmpAlign_2);
			//cout << "tmpAlign_2: " << tmpAlign_2 << endl;
			//cout << "tmpAlign_1: " << endl << "tmpAlign_2: " << tmpAlign_2 << endl;
			int sam_flagInt_1, sam_flagInt_2; 
			string sam_id_1, sam_id_2, sam_seq_1, sam_seq_2, sam_qual_1, sam_qual_2;
			getFlagSeqQualFromSam(tmpAlign_1, tmpAlign_2, sam_flagInt_1, sam_flagInt_2, 
				sam_id_1, sam_id_2, sam_seq_1, sam_seq_2, sam_qual_1, sam_qual_2);
			//cout << "sam_flagInt_1: " << sam_flagInt_1 << endl << "sam_flagInt_2: " << sam_flagInt_2 << endl;
			if((sam_flagInt_1&0x100)||(sam_flagInt_2&0x100)) // mapped, non-primary alignment
				continue;
			string fq_id_end1, fq_id_end2, fq_seq_end1, fq_seq_end2, fq_qual_end1, fq_qual_end2, fq_com_end1, fq_com_end2;
			bool getFastqStrBool = getFastqStr(sam_flagInt_1, sam_flagInt_2, sam_id_1, sam_id_2, sam_seq_1, sam_seq_2, 
				sam_qual_1, sam_qual_2, fq_id_end1, fq_seq_end1, fq_com_end1, fq_qual_end1, 
				fq_id_end2, fq_seq_end2, fq_com_end2, fq_qual_end2);
			if(!getFastqStrBool)
				continue;
			//cout << "getFastqStr done!" << endl;
			fq1_ofs << "@" << fq_id_end1 << endl << fq_seq_end1 << endl << fq_com_end1 << endl << fq_qual_end1 << endl;
			fq2_ofs << "@" << fq_id_end2 << endl << fq_seq_end2 << endl << fq_com_end2 << endl << fq_qual_end2 << endl;
	}
	sam_ifs.close();
	log_ofs.close();
	fq1_ofs.close();
	fq2_ofs.close();
	return 0;
}