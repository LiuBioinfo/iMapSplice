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

void parse2num(string& tmpStr, int& num_1, int& num_2, int& num_3, int& num_4, 
	int& num_5, int& num_6, int& num_7, int& num_8, int& num_9)
{
	vector<string> tmpStrVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 9; tmp++)
	{
		int tmpTabLoc = tmpStr.find(",", startLoc + 1);
		string tmpFieldStr = tmpStr.substr(startLoc, tmpTabLoc - startLoc);
		tmpStrVec.push_back(tmpFieldStr);
		startLoc = tmpTabLoc + 1;
	}
	string tmpLastFieldStr = tmpStr.substr(startLoc);
	tmpStrVec.push_back(tmpLastFieldStr);
	num_1 = atoi(tmpStrVec[1].c_str());
	num_2 = atoi(tmpStrVec[2].c_str());
	num_3 = atoi(tmpStrVec[3].c_str());
	num_4 = atoi(tmpStrVec[4].c_str());
	num_5 = atoi(tmpStrVec[5].c_str());
	num_6 = atoi(tmpStrVec[6].c_str());
	num_7 = atoi(tmpStrVec[7].c_str());
	num_8 = atoi(tmpStrVec[8].c_str());
	num_9 = atoi(tmpStrVec[9].c_str());
}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		cout << "Executable inputLuanAlignmentImprovement outputFile totalReadNum" << endl;
		exit(1);
	}
	string totalReadNumStr = argv[3];
	int totalReadNum = atoi(totalReadNumStr.c_str());

	string inputLuanAlignmentImprovement = argv[1];
	string outputFile = argv[2];
	ifstream luanAlignmentImprovement_ifs(inputLuanAlignmentImprovement.c_str());
	//string outputFormattedAlignmentImprovement = outputFilePrefix + "alignmentImprovement_formatted.txt";
	ofstream formattedAlignmentImprovement_ofs(outputFile.c_str());
	
	vector<string> strVec;
	for(int tmp = 0; tmp < 10; tmp++)
	{
		string tmpStr;
		getline(luanAlignmentImprovement_ifs, tmpStr);
		strVec.push_back(tmpStr);
	}

	int unm_2_unm, unm_2_mulInc, unm_2_uniInc, unm_2_mulOvrWin, unm_2_uniOvrWin, unm_2_mulOvrOut, unm_2_uniOvrOut, unm_2_mulPrf, unm_2_uniPrf;
	int mulInc_2_unm, mulInc_2_mulInc, mulInc_2_uniInc, mulInc_2_mulOvrWin, mulInc_2_uniOvrWin, mulInc_2_mulOvrOut, mulInc_2_uniOvrOut, mulInc_2_mulPrf, mulInc_2_uniPrf;
	int uniInc_2_unm, uniInc_2_mulInc, uniInc_2_uniInc, uniInc_2_mulOvrWin, uniInc_2_uniOvrWin, uniInc_2_mulOvrOut, uniInc_2_uniOvrOut, uniInc_2_mulPrf, uniInc_2_uniPrf;
	int mulOvrWin_2_unm, mulOvrWin_2_mulInc, mulOvrWin_2_uniInc, mulOvrWin_2_mulOvrWin, mulOvrWin_2_uniOvrWin, mulOvrWin_2_mulOvrOut, mulOvrWin_2_uniOvrOut, mulOvrWin_2_mulPrf, mulOvrWin_2_uniPrf;
	int uniOvrWin_2_unm, uniOvrWin_2_mulInc, uniOvrWin_2_uniInc, uniOvrWin_2_mulOvrWin, uniOvrWin_2_uniOvrWin, uniOvrWin_2_mulOvrOut, uniOvrWin_2_uniOvrOut, uniOvrWin_2_mulPrf, uniOvrWin_2_uniPrf;
	int mulOvrOut_2_unm, mulOvrOut_2_mulInc, mulOvrOut_2_uniInc, mulOvrOut_2_mulOvrWin, mulOvrOut_2_uniOvrWin, mulOvrOut_2_mulOvrOut, mulOvrOut_2_uniOvrOut, mulOvrOut_2_mulPrf, mulOvrOut_2_uniPrf;
	int uniOvrOut_2_unm, uniOvrOut_2_mulInc, uniOvrOut_2_uniInc, uniOvrOut_2_mulOvrWin, uniOvrOut_2_uniOvrWin, uniOvrOut_2_mulOvrOut, uniOvrOut_2_uniOvrOut, uniOvrOut_2_mulPrf, uniOvrOut_2_uniPrf;
	int mulPrf_2_unm, mulPrf_2_mulInc, mulPrf_2_uniInc, mulPrf_2_mulOvrWin, mulPrf_2_uniOvrWin, mulPrf_2_mulOvrOut, mulPrf_2_uniOvrOut, mulPrf_2_mulPrf, mulPrf_2_uniPrf;								
	int uniPrf_2_unm, uniPrf_2_mulInc, uniPrf_2_uniInc, uniPrf_2_mulOvrWin, uniPrf_2_uniOvrWin, uniPrf_2_mulOvrOut, uniPrf_2_uniOvrOut, uniPrf_2_mulPrf, uniPrf_2_uniPrf;

	parse2num(strVec[1], unm_2_unm, unm_2_mulInc, unm_2_uniInc, unm_2_mulOvrWin, unm_2_uniOvrWin, unm_2_mulOvrOut, unm_2_uniOvrOut, unm_2_mulPrf, unm_2_uniPrf);
	parse2num(strVec[2], mulInc_2_unm, mulInc_2_mulInc, mulInc_2_uniInc, mulInc_2_mulOvrWin, mulInc_2_uniOvrWin, mulInc_2_mulOvrOut, mulInc_2_uniOvrOut, mulInc_2_mulPrf, mulInc_2_uniPrf);
	parse2num(strVec[3], uniInc_2_unm, uniInc_2_mulInc, uniInc_2_uniInc, uniInc_2_mulOvrWin, uniInc_2_uniOvrWin, uniInc_2_mulOvrOut, uniInc_2_uniOvrOut, uniInc_2_mulPrf, uniInc_2_uniPrf);
	parse2num(strVec[4], mulOvrWin_2_unm, mulOvrWin_2_mulInc, mulOvrWin_2_uniInc, mulOvrWin_2_mulOvrWin, mulOvrWin_2_uniOvrWin, mulOvrWin_2_mulOvrOut, mulOvrWin_2_uniOvrOut, mulOvrWin_2_mulPrf, mulOvrWin_2_uniPrf);
	parse2num(strVec[5], uniOvrWin_2_unm, uniOvrWin_2_mulInc, uniOvrWin_2_uniInc, uniOvrWin_2_mulOvrWin, uniOvrWin_2_uniOvrWin, uniOvrWin_2_mulOvrOut, uniOvrWin_2_uniOvrOut, uniOvrWin_2_mulPrf, uniOvrWin_2_uniPrf);
	parse2num(strVec[6], mulOvrOut_2_unm, mulOvrOut_2_mulInc, mulOvrOut_2_uniInc, mulOvrOut_2_mulOvrWin, mulOvrOut_2_uniOvrWin, mulOvrOut_2_mulOvrOut, mulOvrOut_2_uniOvrOut, mulOvrOut_2_mulPrf, mulOvrOut_2_uniPrf);
	parse2num(strVec[7], uniOvrOut_2_unm, uniOvrOut_2_mulInc, uniOvrOut_2_uniInc, uniOvrOut_2_mulOvrWin, uniOvrOut_2_uniOvrWin, uniOvrOut_2_mulOvrOut, uniOvrOut_2_uniOvrOut, uniOvrOut_2_mulPrf, uniOvrOut_2_uniPrf);
	parse2num(strVec[8], mulPrf_2_unm, mulPrf_2_mulInc, mulPrf_2_uniInc, mulPrf_2_mulOvrWin, mulPrf_2_uniOvrWin, mulPrf_2_mulOvrOut, mulPrf_2_uniOvrOut, mulPrf_2_mulPrf, mulPrf_2_uniPrf);
	parse2num(strVec[9], uniPrf_2_unm, uniPrf_2_mulInc, uniPrf_2_uniInc, uniPrf_2_mulOvrWin, uniPrf_2_uniOvrWin, uniPrf_2_mulOvrOut, uniPrf_2_uniOvrOut, uniPrf_2_mulPrf, uniPrf_2_uniPrf);

	int unm_2_mulOvr = unm_2_mulOvrWin + unm_2_mulOvrOut;
	int unm_2_uniOvr = unm_2_uniOvrWin + unm_2_uniOvrOut;
	
	int mulInc_2_mulOvr = mulInc_2_mulOvrWin + mulInc_2_mulOvrOut;
	int mulInc_2_uniOvr = mulInc_2_uniOvrWin + mulInc_2_uniOvrOut;

	int uniInc_2_mulOvr = uniInc_2_mulOvrWin + uniInc_2_mulOvrOut;
	int uniInc_2_uniOvr = uniInc_2_uniOvrWin + uniInc_2_uniOvrOut;

	int mulOvr_2_unm = mulOvrWin_2_unm + mulOvrOut_2_unm;
	int mulOvr_2_mulInc = mulOvrWin_2_mulInc + mulOvrOut_2_mulInc;
	int mulOvr_2_uniInc = mulOvrWin_2_uniInc + mulOvrOut_2_uniInc;
	int mulOvr_2_mulOvr = mulOvrWin_2_mulOvrWin + mulOvrWin_2_mulOvrOut + mulOvrOut_2_mulOvrWin + mulOvrOut_2_mulOvrOut;
	int mulOvr_2_uniOvr = mulOvrWin_2_uniOvrWin + mulOvrWin_2_uniOvrOut + mulOvrOut_2_uniOvrWin + mulOvrOut_2_uniOvrOut;
	int mulOvr_2_mulPrf = mulOvrWin_2_mulPrf + mulOvrOut_2_mulPrf;
	int mulOvr_2_uniPrf = mulOvrWin_2_uniPrf + mulOvrOut_2_uniPrf;

	int uniOvr_2_unm = uniOvrWin_2_unm + uniOvrOut_2_unm;
	int uniOvr_2_mulInc = uniOvrWin_2_mulInc + uniOvrOut_2_mulInc;
	int uniOvr_2_uniInc = uniOvrWin_2_uniInc + uniOvrOut_2_uniInc;
	int uniOvr_2_mulOvr = uniOvrWin_2_mulOvrWin + uniOvrWin_2_mulOvrOut + uniOvrOut_2_mulOvrWin + uniOvrOut_2_mulOvrOut;
	int uniOvr_2_uniOvr = uniOvrWin_2_uniOvrWin + uniOvrWin_2_uniOvrOut + uniOvrOut_2_uniOvrWin + uniOvrOut_2_uniOvrOut;
	int uniOvr_2_mulPrf = uniOvrWin_2_mulPrf + uniOvrOut_2_mulPrf;
	int uniOvr_2_uniPrf = uniOvrWin_2_uniPrf + uniOvrOut_2_uniPrf;

	int mulPrf_2_mulOvr = mulPrf_2_mulOvrWin + mulPrf_2_mulOvrOut;
	int mulPrf_2_uniOvr = mulPrf_2_uniOvrWin + mulPrf_2_uniOvrOut;

	int uniPrf_2_mulOvr = uniPrf_2_mulOvrWin + uniPrf_2_mulOvrOut;
	int uniPrf_2_uniOvr = uniPrf_2_uniOvrWin + uniPrf_2_uniOvrOut;

	int unm_2_total = unm_2_unm + unm_2_mulInc + unm_2_uniInc + unm_2_mulOvr + unm_2_uniOvr + unm_2_mulPrf + unm_2_uniPrf;
	int mulInc_2_total = mulInc_2_unm + mulInc_2_mulInc + mulInc_2_uniInc + mulInc_2_mulOvr + mulInc_2_uniOvr + mulInc_2_mulPrf + mulInc_2_uniPrf;
	int uniInc_2_total = uniInc_2_unm + uniInc_2_mulInc + uniInc_2_uniInc + uniInc_2_mulOvr + uniInc_2_uniOvr + uniInc_2_mulPrf + uniInc_2_uniPrf;
	int mulOvr_2_total = mulOvr_2_unm + mulOvr_2_mulInc + mulOvr_2_uniInc + mulOvr_2_mulOvr + mulOvr_2_uniOvr + mulOvr_2_mulPrf + mulOvr_2_uniPrf;
	int uniOvr_2_total = uniOvr_2_unm + uniOvr_2_mulInc + uniOvr_2_uniInc + uniOvr_2_mulOvr + uniOvr_2_uniOvr + uniOvr_2_mulPrf + uniOvr_2_uniPrf;
	int mulPrf_2_total = mulPrf_2_unm + mulPrf_2_mulInc + mulPrf_2_uniInc + mulPrf_2_mulOvr + mulPrf_2_uniOvr + mulPrf_2_mulPrf + mulPrf_2_uniPrf;
	int uniPrf_2_total = uniPrf_2_unm + uniPrf_2_mulInc + uniPrf_2_uniInc + uniPrf_2_mulOvr + uniPrf_2_uniOvr + uniPrf_2_mulPrf + uniPrf_2_uniPrf;

	int all_2_total = unm_2_total + mulInc_2_total + uniInc_2_total + mulOvr_2_total + uniOvr_2_total + mulPrf_2_total + uniPrf_2_total;

	int total_2_unm = unm_2_unm + mulInc_2_unm + uniInc_2_unm + mulOvr_2_unm + uniOvr_2_unm + mulPrf_2_unm + uniPrf_2_unm;
	int total_2_mulInc = unm_2_mulInc + mulInc_2_mulInc + uniInc_2_mulInc + mulOvr_2_mulInc + uniOvr_2_mulInc + mulPrf_2_mulInc + uniPrf_2_mulInc;
	int total_2_uniInc = unm_2_uniInc + mulInc_2_uniInc + uniInc_2_uniInc + mulOvr_2_uniInc + uniOvr_2_uniInc + mulPrf_2_uniInc + uniPrf_2_uniInc;
	int total_2_mulOvr = unm_2_mulOvr + mulInc_2_mulOvr + uniInc_2_mulOvr + mulOvr_2_mulOvr + uniOvr_2_mulOvr + mulPrf_2_mulOvr + uniPrf_2_mulOvr;
	int total_2_uniOvr = unm_2_uniOvr + mulInc_2_uniOvr + uniInc_2_uniOvr + mulOvr_2_uniOvr + uniOvr_2_uniOvr + mulPrf_2_uniOvr + uniPrf_2_uniOvr;
	int total_2_mulPrf = unm_2_mulPrf + mulInc_2_mulPrf + uniInc_2_mulPrf + mulOvr_2_mulPrf + uniOvr_2_mulPrf + mulPrf_2_mulPrf + uniPrf_2_mulPrf;
	int total_2_uniPrf = unm_2_uniPrf + mulInc_2_uniPrf + uniInc_2_uniPrf + mulOvr_2_uniPrf + uniOvr_2_uniPrf + mulPrf_2_uniPrf + uniPrf_2_uniPrf;	

	int total_2_all = total_2_unm + total_2_mulInc + total_2_uniInc + total_2_mulOvr + total_2_uniOvr + total_2_mulPrf + total_2_uniPrf;

	int total_2_total;
	cout << "all_2_total: " << all_2_total << endl;
	cout << "total_2_all: " << total_2_all << endl;
	if(all_2_total != total_2_all)
		cout << "all_2_total != total_2_all" << endl;
	else
		total_2_total = total_2_all;

	int unm_2_unm_toAdd = totalReadNum - total_2_all;
	int unm_2_unm_new = unm_2_unm + unm_2_unm_toAdd;
	int unm_2_total_new = unm_2_total + unm_2_unm_toAdd;
	int total_2_unm_new = total_2_unm + unm_2_unm_toAdd;
	int total_2_total_new = total_2_total + unm_2_unm_toAdd;

	formattedAlignmentImprovement_ofs << "&&Unm&Mul-Inc&Uni-Inc&Mul-Ovr&Uni-Ovr&Mul-Prf&Uni-Prf&Total\\\\" << endl;
	formattedAlignmentImprovement_ofs << "&Unm&" << unm_2_unm_new << "&" << unm_2_mulInc << "&" << unm_2_uniInc 
		<< "&" << unm_2_mulOvr << "&" << unm_2_uniOvr << "&" << unm_2_mulPrf << "&" << unm_2_uniPrf << "&" << unm_2_total_new << "\\\\" << endl;
	formattedAlignmentImprovement_ofs << "&Mul-Inc&" << mulInc_2_unm << "&" << mulInc_2_mulInc << "&" << mulInc_2_uniInc 
		<< "&" << mulInc_2_mulOvr << "&" << mulInc_2_uniOvr << "&" << mulInc_2_mulPrf << "&" << mulInc_2_uniPrf << "&" << mulInc_2_total << "\\\\" << endl;
	formattedAlignmentImprovement_ofs << "&Uni-Inc&" << uniInc_2_unm << "&" << uniInc_2_mulInc << "&" << uniInc_2_uniInc 
		<< "&" << uniInc_2_mulOvr << "&" << uniInc_2_uniOvr << "&" << uniInc_2_mulPrf << "&" << uniInc_2_uniPrf << "&" << uniInc_2_total << "\\\\" << endl;		
	formattedAlignmentImprovement_ofs << "&Mul-Ovr&" << mulOvr_2_unm << "&" << mulOvr_2_mulInc << "&" << mulOvr_2_uniInc 
		<< "&" << mulOvr_2_mulOvr << "&" << mulOvr_2_uniOvr << "&" << mulOvr_2_mulPrf << "&" << mulOvr_2_uniPrf << "&" << mulOvr_2_total << "\\\\" << endl;
	formattedAlignmentImprovement_ofs << "&Uni-Ovr&" << uniOvr_2_unm << "&" << uniOvr_2_mulInc << "&" << uniOvr_2_uniInc 
		<< "&" << uniOvr_2_mulOvr << "&" << uniOvr_2_uniOvr << "&" << uniOvr_2_mulPrf << "&" << uniOvr_2_uniPrf << "&" << uniOvr_2_total << "\\\\" << endl;	
	formattedAlignmentImprovement_ofs << "&Mul-Prf&" << mulPrf_2_unm << "&" << mulPrf_2_mulInc << "&" << mulPrf_2_uniInc 
		<< "&" << mulPrf_2_mulOvr << "&" << mulPrf_2_uniOvr << "&" << mulPrf_2_mulPrf << "&" << mulPrf_2_uniPrf << "&" << mulPrf_2_total << "\\\\" << endl;
	formattedAlignmentImprovement_ofs << "&Uni-Prf&" << uniPrf_2_unm << "&" << uniPrf_2_mulInc << "&" << uniPrf_2_uniInc 
		<< "&" << uniPrf_2_mulOvr << "&" << uniPrf_2_uniOvr << "&" << uniPrf_2_mulPrf << "&" << uniPrf_2_uniPrf << "&" << uniPrf_2_total << "\\\\" << endl;
	formattedAlignmentImprovement_ofs << "&Total&" << total_2_unm_new << "&" << total_2_mulInc << "&" << total_2_uniInc 
		<< "&" << total_2_mulOvr << "&" << total_2_uniOvr << "&" << total_2_mulPrf << "&" << total_2_uniPrf << "&" << total_2_total_new << "\\\\" << endl;
	luanAlignmentImprovement_ifs.close();
	formattedAlignmentImprovement_ofs.close();
	return 0;
}