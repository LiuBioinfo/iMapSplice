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

int baseStr2int(string& base)
{
	if(base == "A")
		return 0;
	else if(base == "C")
		return 1;
	else if(base == "G")
		return 2;
	else if(base == "T")
		return 3;
	else
	{
		cout << "base is not A, T, C or G" << endl;
		exit(1);
	}
}

int flankString2int(string& tmpFlkStr)
{
	string tmpFlkStr_1stBase = tmpFlkStr.substr(0,1);
	string tmpFlkStr_2ndBase = tmpFlkStr.substr(1,1);
	string tmpFlkStr_3rdBase = tmpFlkStr.substr(2,1);
	string tmpFlkStr_4thBase = tmpFlkStr.substr(3,1);
	int int_1stBase = baseStr2int(tmpFlkStr_1stBase);
	int int_2ndBase = baseStr2int(tmpFlkStr_2ndBase);
	int int_3rdBase = baseStr2int(tmpFlkStr_3rdBase);
	int int_4thBase = baseStr2int(tmpFlkStr_4thBase);
	int tmpScore = int_1stBase * 64 + int_2ndBase * 16 + int_3rdBase * 4 + int_4thBase;
	return tmpScore;
}

string getRcmBase(string& base)
{
	if(base == "A")
		return "T";
	else if(base == "C")
		return "G";
	else if(base == "G")
		return "C";
	else if(base == "T")
		return "A";
	else
	{
		cout << "base is not A, T, C or G" << endl;
		exit(1);
	}
}

string getRcmFlkStr(string& tmpFlkStr)
{
	string tmpFlkStr_1stBase = tmpFlkStr.substr(0,1);
	string tmpFlkStr_2ndBase = tmpFlkStr.substr(1,1);
	string tmpFlkStr_3rdBase = tmpFlkStr.substr(2,1);
	string tmpFlkStr_4thBase = tmpFlkStr.substr(3,1);
	string tmpFlkStr_1stBase_rcm = getRcmBase(tmpFlkStr_1stBase);
	string tmpFlkStr_2ndBase_rcm = getRcmBase(tmpFlkStr_2ndBase);
	string tmpFlkStr_3rdBase_rcm = getRcmBase(tmpFlkStr_3rdBase);
	string tmpFlkStr_4thBase_rcm = getRcmBase(tmpFlkStr_4thBase);
	string tmpRcmFlkStr = tmpFlkStr_4thBase_rcm + tmpFlkStr_3rdBase_rcm + tmpFlkStr_2ndBase_rcm + tmpFlkStr_1stBase_rcm;
	return tmpRcmFlkStr;
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputSJfile flankStringFieldNO outputFile" << endl;
		exit(1);
	}
	string inputSJfile = argv[1];
	string flankStringFieldNOstr = argv[2];
	int flankStringFieldNO = atoi(flankStringFieldNOstr.c_str());
	string outputFile = argv[3];
	ifstream SJ_ifs(inputSJfile.c_str());

	vector<string> baseVec;
	baseVec.push_back("A");
	baseVec.push_back("C");
	baseVec.push_back("G");
	baseVec.push_back("T");

	vector<int> flkFreqVec;
	for(int tmp = 0; tmp < 256; tmp++)
		flkFreqVec.push_back(0);

	vector<string> flkStrVec;
	for(int tmp1 = 0; tmp1 < 4; tmp1++)
	{
		string tmpBase_1st = baseVec[tmp1];
		for(int tmp2 = 0; tmp2 < 4; tmp2++)
		{
			string tmpBase_2nd = baseVec[tmp2];
			for(int tmp3 = 0; tmp3 < 4; tmp3++)
			{
				string tmpBase_3rd = baseVec[tmp3];
				for(int tmp4 = 0; tmp4 < 4; tmp4++)
				{
					string tmpBase_4th = baseVec[tmp4];
					string tmpFlkStr = tmpBase_1st + tmpBase_2nd + tmpBase_3rd + tmpBase_4th;
					flkStrVec.push_back(tmpFlkStr);
				}
			}	
		}		
	}

	while(!SJ_ifs.eof())
	{
		string tmpStr;
		getline(SJ_ifs, tmpStr);
		if(tmpStr == "")
			break;
		vector<string> strVec;
		int startLoc = 0;
		for(int tmp = 0; tmp < flankStringFieldNO; tmp++)
		{
			int tabLoc = tmpStr.find("\t", startLoc);
			string tmpFieldStr = tmpStr.substr(startLoc, tabLoc - startLoc);
			strVec.push_back(tmpFieldStr);
			startLoc = tabLoc + 1;
		}
		string tmpFlkStr = strVec[flankStringFieldNO - 1];
		if(tmpFlkStr.length() != 4)
		{
			cout << "tmpFlkStr.leng() != 4 " << endl;
			exit(1);
		}
		int tmpFlkStr_score = flankString2int(tmpFlkStr);
		flkFreqVec[tmpFlkStr_score] ++;
	}

	ofstream flkFreq_ofs(outputFile.c_str());
	// for(int tmp = 0; tmp < 256; tmp++)
	// {
	// 	flkFreq_ofs << flkStrVec[tmp] << "\t" << flkFreqVec[tmp] << endl;
	// }
	for(int tmp = 0; tmp < 256; tmp++)
	{
		string tmpFlkStr = flkStrVec[tmp];
		string tmpFlkStr_rcm = getRcmFlkStr(tmpFlkStr);
		int tmpFlkStr_score = flankString2int(tmpFlkStr);
		int tmpFlkStr_rcm_score = flankString2int(tmpFlkStr_rcm);
		int tmpFlkStr_freq = flkFreqVec[tmpFlkStr_score];
		int tmpFlkStr_rcm_freq = flkFreqVec[tmpFlkStr_rcm_score];
		if(tmpFlkStr_score < tmpFlkStr_rcm_score)
		{
			flkFreq_ofs << tmpFlkStr << "\t" << tmpFlkStr_rcm << "\t" << tmpFlkStr_freq + tmpFlkStr_rcm_freq 
				<< "\t" << tmpFlkStr_freq << "\t" << tmpFlkStr_rcm_freq << endl; 
		}
		else if(tmpFlkStr_score == tmpFlkStr_rcm_score)
		{
			flkFreq_ofs << tmpFlkStr << "\tNULL\t" << tmpFlkStr_freq << "\t" << tmpFlkStr_freq << "\t0" << endl; 
		}
		else
		{}
	}

	flkFreq_ofs.close();
	SJ_ifs.close();
	return 0;
}