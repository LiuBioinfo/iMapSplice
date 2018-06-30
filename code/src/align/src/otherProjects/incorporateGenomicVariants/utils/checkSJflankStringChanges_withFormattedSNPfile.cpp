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

#include "../../../general/read_block_test.h"
#include "../../../general/index_info.h"

using namespace std;

void extractSJchrNamePosSupNumFromStr(string& tmpSJstr, string& tmpSJ_chrName, 
	int& tmpSJ_startPos, int& tmpSJ_endPos, int& tmpSJ_supNum, string& tmpOtherStr)
{
	int tabLoc_1 = tmpSJstr.find("\t");
	int tabLoc_2 = tmpSJstr.find("\t", tabLoc_1 + 1);
	int tabLoc_3 = tmpSJstr.find("\t", tabLoc_2 + 1);
	int tabLoc_4 = tmpSJstr.find("\t", tabLoc_3 + 1);
	if(tabLoc_4 == string::npos)
	{
		cout << "incorrect SJ format" << endl;
		exit(1);
	}
	int tabLoc_5 = tmpSJstr.find("\t", tabLoc_4 + 1);
	string tmpSupNumStr;
	if(tabLoc_5 == string::npos)
	{	
		tmpSupNumStr = tmpSJstr.substr(tabLoc_4 + 1);
		tmpOtherStr = "";
	}
	else
	{	
		tmpSupNumStr = tmpSJstr.substr(tabLoc_4 + 1, tabLoc_5 - tabLoc_4 - 1);	
		tmpOtherStr = tmpSJstr.substr(tabLoc_5 + 1);
	}
	tmpSJ_chrName = tmpSJstr.substr(0,tabLoc_1);
	string tmpStartPosStr = tmpSJstr.substr(tabLoc_1 + 1, tabLoc_2 - tabLoc_1 - 1);
	string tmpEndPosStr = tmpSJstr.substr(tabLoc_2 + 1, tabLoc_3 - tabLoc_2 - 1);	
	tmpSJ_startPos = atoi(tmpStartPosStr.c_str());
	tmpSJ_endPos = atoi(tmpEndPosStr.c_str());
	tmpSJ_supNum = atoi(tmpSupNumStr.c_str());
}

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cout << "Executable inputSdIndex inputPsSNPfile inputSJ outputFlankStringChangesFolder SJ_size_min SJ_size_max" << endl;
		exit(1);
	}
	string SJ_size_min_str = argv[5];
	string SJ_size_max_str = argv[6];
	int SJ_size_min = atoi(SJ_size_min_str.c_str());
	int SJ_size_max = atoi(SJ_size_max_str.c_str());

	cout << "creating results folder ...." << endl;
	string outputFolderStr = argv[4];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_path = outputFolderStr + "log.txt";
	ofstream log_ofs(log_path.c_str());	

	cout << "start to initiate indexInfo for both sd and ps genome" << endl;
	log_ofs << "start to initiate indexInfo for both sd and ps genome" << endl;
	string indexFolderPath_sd = argv[1];
	string indexFolderPath_ps = argv[1];
	indexFolderPath_sd += "/";
	indexFolderPath_ps += "/";
	string chrom_bit_file_sd = indexFolderPath_sd; chrom_bit_file_sd.append("_chrom");
	string chrom_bit_file_ps = indexFolderPath_ps; chrom_bit_file_ps.append("_chrom"); 
	ifstream chrom_bit_file_ifs_sd(chrom_bit_file_sd.c_str(),ios::binary);
	ifstream chrom_bit_file_ifs_ps(chrom_bit_file_ps.c_str(),ios::binary);	
	string indexParameterFileStr_sd = indexFolderPath_sd + "/_parameter";
	string indexParameterFileStr_ps = indexFolderPath_ps + "/_parameter";
	ifstream parameter_ifs_sd(indexParameterFileStr_sd.c_str());
	ifstream parameter_ifs_ps(indexParameterFileStr_ps.c_str());
	Index_Info* indexInfo_sd = new Index_Info(parameter_ifs_sd);
	Index_Info* indexInfo_ps = new Index_Info(parameter_ifs_ps);
	int chromNum_sd = indexInfo_sd->returnChromNum();
	int chromNum_ps = indexInfo_ps->returnChromNum();
	if(chromNum_sd != chromNum_ps)
	{
		cout << "different chrom # for sd and ps genome" << endl;
		exit(1);
	}
	char *chrom_sd; chrom_sd = (char*)malloc((indexInfo_sd->returnIndexSize()) * sizeof(char));
	char *chrom_ps; chrom_ps = (char*)malloc((indexInfo_ps->returnIndexSize()) * sizeof(char)); 
	chrom_bit_file_ifs_sd.read((char*)chrom_sd, (indexInfo_sd->returnIndexSize()) * sizeof(char));
	chrom_bit_file_ifs_ps.read((char*)chrom_ps, (indexInfo_ps->returnIndexSize()) * sizeof(char));
	indexInfo_sd->readGenome(chrom_sd);
	indexInfo_ps->readGenome(chrom_ps);
	indexInfo_sd->initiate();
	indexInfo_ps->initiate();
	indexInfo_sd->initiateChrNameIndexArray(1000);
	indexInfo_ps->initiateChrNameIndexArray(1000);
	cout << "insert SNPs into indexInfo_ps" << endl;
	log_ofs << "insert SNPs into indexInfo_ps" << endl;
	string formattedSNPfile = argv[2];
	indexInfo_ps->insertSNP2chromStr(formattedSNPfile, log_ofs);
	cout << "end of initiating indexInfo" << endl;
	log_ofs << "end of initiating indexInfo" << endl;

	string output_flankStringChanged = outputFolderStr + "SJ_flankString_changed.txt";
	string output_flankStringChanged_2canonicalSJ = outputFolderStr + "SJ_flankString_changed_2canonicalSJ.txt";
	string output_flankStringChanged_others = outputFolderStr + "SJ_flankString_changed_others.txt";
	string output_flankStringKept = outputFolderStr + "SJ_flankString_kept.txt";
	string output_compare = outputFolderStr + "compare.txt";
	ofstream changed_ofs(output_flankStringChanged.c_str());
	ofstream changed_2canonicalSJ_ofs(output_flankStringChanged_2canonicalSJ.c_str());
	ofstream changed_others_ofs(output_flankStringChanged_others.c_str());
	ofstream kept_ofs(output_flankStringKept.c_str());
	ofstream compare_ofs(output_compare.c_str());
	string inputSJpath = argv[3];
	ifstream SJ_ifs(inputSJpath.c_str());
	int SJnum_flankStringChanged = 0;
	int SJnum_flankStringChanged_2canonicalSJ = 0;
	int SJnum_flankStringChanged_others = 0;
	int SJnum_flankStringkept = 0;
	while(!SJ_ifs.eof())
	{
		string tmpSJstr;
		getline(SJ_ifs, tmpSJstr);
		if(tmpSJstr == "")
			break;
		string tmpSJ_chrName;
		int tmpSJ_startPos;
		int tmpSJ_endPos;
		int tmpSJ_supNum;
		string tmpOtherStr;
		extractSJchrNamePosSupNumFromStr(tmpSJstr, tmpSJ_chrName, tmpSJ_startPos, 
			tmpSJ_endPos, tmpSJ_supNum, tmpOtherStr);
		int tmpSJ_size = tmpSJ_endPos - tmpSJ_startPos - 1;
		if((tmpSJ_size < SJ_size_min)||(tmpSJ_size > SJ_size_max))
			continue;
		int tmpSJ_chrNameInt_sd = indexInfo_sd->convertStringToInt(tmpSJ_chrName);
		int tmpSJ_chrNameInt_ps = indexInfo_ps->convertStringToInt(tmpSJ_chrName);
		if(tmpSJ_chrNameInt_sd != tmpSJ_chrNameInt_ps)
		{
			cout << "inconsistent SJ chrName !" << endl;
			exit(1);
		}
		if(tmpSJ_chrNameInt_sd < 0)
		{
			cout << "SJ chrNameInt < 0 !" << endl;
			exit(1);
		}
		int tmpSJ_chrNameInt = tmpSJ_chrNameInt_sd;
		string tmpSJ_flankString_sd = indexInfo_sd->returnFlankString(
			tmpSJ_chrNameInt, tmpSJ_startPos, tmpSJ_endPos);
		string tmpSJ_flankString_ps = indexInfo_ps->returnFlankString(
			tmpSJ_chrNameInt, tmpSJ_startPos, tmpSJ_endPos);
		if(tmpSJ_flankString_sd == tmpSJ_flankString_ps)
		{
			SJnum_flankStringkept ++;
			kept_ofs << tmpSJ_chrName << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos
				<< "\t" << tmpSJ_flankString_sd << "\t" << tmpSJ_supNum << "\t" << tmpOtherStr << endl;
		}
		else
		{
			SJnum_flankStringChanged ++;
			changed_ofs << tmpSJ_chrName << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos
				<< "\t" << tmpSJ_flankString_sd << "->" << tmpSJ_flankString_ps << "\t" << tmpSJ_supNum 
				<< "\t" << tmpOtherStr  << endl;
			if((tmpSJ_flankString_ps == "GTAG")||(tmpSJ_flankString_ps == "CTAC"))
			{
				SJnum_flankStringChanged_2canonicalSJ ++;
				changed_2canonicalSJ_ofs << tmpSJ_chrName << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos
					<< "\t" << tmpSJ_flankString_sd << "->" << tmpSJ_flankString_ps << "\t" << tmpSJ_supNum 
					<< "\t" << tmpOtherStr  << endl;
			}
			else
			{
				SJnum_flankStringChanged_others ++;
				changed_others_ofs << tmpSJ_chrName << "\t" << tmpSJ_startPos << "\t" << tmpSJ_endPos
					<< "\t" << tmpSJ_flankString_sd << "->" << tmpSJ_flankString_ps << "\t" << tmpSJ_supNum 
					<< "\t" << tmpOtherStr  << endl;
			}
		}
	}
	int totalSJnum = SJnum_flankStringkept + SJnum_flankStringChanged;
	double SJperc_flankStringChanged = ((double)SJnum_flankStringChanged / (double)totalSJnum) * 100;
	double SJperc_flankStringChanged_2canonicalSJ = ((double)SJnum_flankStringChanged_2canonicalSJ / (double)totalSJnum) * 100;
	double SJperc_flankStringChanged_others = ((double)SJnum_flankStringChanged_others / (double)totalSJnum) * 100;
	double SJperc_flankStringKept = ((double)SJnum_flankStringkept / (double)totalSJnum) * 100;
	double SJperc_2canonicalSJ_amongChangedSJ = ((double)SJnum_flankStringChanged_2canonicalSJ / (double)SJnum_flankStringChanged) * 100;
	double SJperc_others_amongChangedSJ = ((double)SJnum_flankStringChanged_others / (double)SJnum_flankStringChanged) * 100;

	cout << endl << "        totalSJnum: " << totalSJnum << endl << endl;
	cout << "   flankStringkept: " << SJnum_flankStringkept << " -- " << SJperc_flankStringKept << "%" << endl << endl;	
	cout << "flankStringChanged: " << SJnum_flankStringChanged << " -- " << SJperc_flankStringChanged << "%" << endl;
	cout << "      2canonicalSJ: " << SJnum_flankStringChanged_2canonicalSJ << " -- " 
		<< SJperc_flankStringChanged_2canonicalSJ << "% -- " << SJperc_2canonicalSJ_amongChangedSJ << "%" << endl;
	cout << "            others: " << SJnum_flankStringChanged_others << " -- " << SJperc_flankStringChanged_others 
		<< "% -- " << SJperc_others_amongChangedSJ << "%" << endl;

	compare_ofs << "        totalSJnum: " << totalSJnum << endl << endl;
	compare_ofs << "   flankStringkept: " << SJnum_flankStringkept << " -- " << SJperc_flankStringKept << "%" << endl << endl;	
	compare_ofs << "flankStringChanged: " << SJnum_flankStringChanged << " -- " << SJperc_flankStringChanged << "%" << endl;
	compare_ofs << "      2canonicalSJ: " << SJnum_flankStringChanged_2canonicalSJ << " -- " 
		<< SJperc_flankStringChanged_2canonicalSJ << "% -- " << SJperc_2canonicalSJ_amongChangedSJ << "%" << endl;
	compare_ofs << "            others: " << SJnum_flankStringChanged_others << " -- " << SJperc_flankStringChanged_others 
		<< "% -- " << SJperc_others_amongChangedSJ << "%" << endl;

	SJ_ifs.close();
	compare_ofs.close();
	kept_ofs.close();
	changed_others_ofs.close();
	changed_2canonicalSJ_ofs.close();
	changed_ofs.close();
	delete indexInfo_ps;
	delete indexInfo_sd;
	log_ofs.close();
	parameter_ifs_ps.close();
	parameter_ifs_sd.close();
	return 0;
}