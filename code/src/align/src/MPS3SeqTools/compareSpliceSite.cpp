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
//#include <hash_map>
#include <map>
#include <set>

#include "../general/read_block_test.h"
#include "../general/index_info.h"
#include "../general/splice_info.h"
//#include "../general/alignmentToJunc_supportNum.h"
#include "../general/fixSingleAnchorNWDP_info.h"

using namespace std;

void returnSpliceSite(string& tmpJuncStr, int& chrNameInt,
	 int& donerEndPos, int& acceptorStartPos, Index_Info* indexInfo)
{
	vector<string> juncFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 3; tmp++)
	{
		int tabLoc = tmpJuncStr.find("\t", startLoc);
		string tmpJuncField = tmpJuncStr.substr(startLoc, tabLoc-startLoc);
		juncFieldVec.push_back(tmpJuncField);
		startLoc = tabLoc + 1;
	}
	string chrNameStr = juncFieldVec[0];
	chrNameInt = indexInfo->convertStringToInt(chrNameStr);
	string donerEndPosStr = juncFieldVec[1];
	donerEndPos = atoi(donerEndPosStr.c_str());
	string acceptorStartPosStr = juncFieldVec[2];
	acceptorStartPos = atoi(acceptorStartPosStr.c_str());
}

bool returnSpliceSiteSearchResults(
	vector< set<int> >& tmpSpliceSiteSetVec,
	int offset,
	int tmpChrNameInt, int tmpSpliceSite,
	vector< set<int> >& spliceSiteSetVec_compared2_detected,
	vector< set<int> >& spliceSiteSetVec_compared2_correct,
	vector< set<int> >& spliceSiteSetVec_compared2_incorrect
	)
{
	set<int>::iterator intSetIter;
	for(int tmpPos = tmpSpliceSite-offset; tmpPos <= tmpSpliceSite+offset; tmpPos++)
	{
		intSetIter = tmpSpliceSiteSetVec[tmpChrNameInt].find(tmpPos);
		if(intSetIter != tmpSpliceSiteSetVec[tmpChrNameInt].end())
		{
			spliceSiteSetVec_compared2_detected[tmpChrNameInt].insert(tmpPos);
			spliceSiteSetVec_compared2_correct[tmpChrNameInt].insert(tmpSpliceSite);
			return true;
		}
	}
	spliceSiteSetVec_compared2_incorrect[tmpChrNameInt].insert(tmpSpliceSite);
	return false;
}

int returnSpliceSiteSearchResults(
	vector< set<int> >& donerEndPosSetVec,
	vector< set<int> >& acceptorStartPosSetVec,
	int offset, int tmpChrNameInt, 
	int tmpDonerEndPos, int tmpAcceptorStartPos,

	vector< set<int> >& donerEndPosSetVec_compared2_detected,
	vector< set<int> >& acceptorStartPosSetVec_compared2_detected,

	vector< set<int> >& donerEndPosSetVec_compared2_correct,
	vector< set<int> >& acceptorStartPosSetVec_compared2_correct,
	
	vector< set<int> >& donerEndPosSetVec_compared2_incorrect,
	vector< set<int> >& acceptorStartPosSetVec_compared2_incorrect	
	)
{
	bool tmpDonerEndPosFound_bool = returnSpliceSiteSearchResults(
		donerEndPosSetVec, offset, tmpChrNameInt, tmpDonerEndPos,
		donerEndPosSetVec_compared2_detected,
		donerEndPosSetVec_compared2_correct,
		donerEndPosSetVec_compared2_incorrect
		);
	bool tmpAcceptorStartPosFound_bool = returnSpliceSiteSearchResults(
		acceptorStartPosSetVec, offset, tmpChrNameInt, tmpAcceptorStartPos,
		acceptorStartPosSetVec_compared2_detected,
		acceptorStartPosSetVec_compared2_correct,
		acceptorStartPosSetVec_compared2_incorrect
		);
	if(tmpDonerEndPosFound_bool && tmpAcceptorStartPosFound_bool)
	{
		return 1;
	}
	else if(tmpDonerEndPosFound_bool)
	{
		return 2;
	}
	else if(tmpAcceptorStartPosFound_bool)
	{
		return 3;
	}
	else
	{
		return 0;
	}
}

int returnItemNumInSetVec(vector< set<int> >& tmpSpliceSiteSetVec)
{
	int tmpItemNum = 0;
	for(int tmp = 0; tmp < tmpSpliceSiteSetVec.size(); tmp++)
	{
		tmpItemNum = tmpItemNum + tmpSpliceSiteSetVec[tmp].size();
	}
	return tmpItemNum;
}

int main(int argc, char** argv)
{
	if(argc != 6)
	{
		cout << "Executable indexFolderPath inputGroundTruthSJ inputComparedSJ outputSJfile offset" << endl;
		exit(1);
	}
	string inputGroundTruthSJfile = argv[2];
	string inputComparedSJfile = argv[3];
	string outputFolder = argv[4];
	string mkdir = "mkdir -p " + outputFolder;
	system(mkdir.c_str());
	string outputSJfilePrefix = outputFolder + "/";
	string offsetStr = argv[5];
	int offset = atoi(offsetStr.c_str());
	string outputSJbothSpliceSitesCorrectFile = outputSJfilePrefix + "bothSpliceSitesCorrect.junc";
	string outputSJbothSpliceSitesIncorrectFile = outputSJfilePrefix + "nonSpliceSitesCorrect.junc";
	string outputSJdonerSpliceSiteCorrectFile = outputSJfilePrefix + "donerSpliceSiteCorrect.junc";
	string outputSJacceptorSpliceSiteCorrectFile = outputSJfilePrefix + "acceptorSpliceSiteCorrect.junc";
	string logFile = outputSJfilePrefix + "log.txt";
	ofstream log_ofs(logFile.c_str());

	string indexFolderPath = argv[1];
	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	cout << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);

	vector< set<int> > donerEndPosSetVec;
	vector< set<int> > acceptorStartPosSetVec;
	int chromNum = indexInfo->returnChromNum();
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpDonerEndPosSet;
		donerEndPosSetVec.push_back(tmpDonerEndPosSet);
		set<int> tmpAcceptorStartPosSet;
		acceptorStartPosSetVec.push_back(tmpAcceptorStartPosSet);
	}

	vector< set<int> > donerEndPosSetVec_compared2_detected;
	vector< set<int> > acceptorStartPosSetVec_compared2_detected;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpDonerEndPosSet;
		donerEndPosSetVec_compared2_detected.push_back(tmpDonerEndPosSet);
		set<int> tmpAcceptorStartPosSet;
		acceptorStartPosSetVec_compared2_detected.push_back(tmpAcceptorStartPosSet);
	}	

	vector< set<int> > donerEndPosSetVec_compared2_correct;
	vector< set<int> > acceptorStartPosSetVec_compared2_correct;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpDonerEndPosSet;
		donerEndPosSetVec_compared2_correct.push_back(tmpDonerEndPosSet);
		set<int> tmpAcceptorStartPosSet;
		acceptorStartPosSetVec_compared2_correct.push_back(tmpAcceptorStartPosSet);
	}	

	vector< set<int> > donerEndPosSetVec_compared2_incorrect;
	vector< set<int> > acceptorStartPosSetVec_compared2_incorrect;
	for(int tmp = 0; tmp < chromNum; tmp++)
	{
		set<int> tmpDonerEndPosSet;
		donerEndPosSetVec_compared2_incorrect.push_back(tmpDonerEndPosSet);
		set<int> tmpAcceptorStartPosSet;
		acceptorStartPosSetVec_compared2_incorrect.push_back(tmpAcceptorStartPosSet);
	}	

	int groundTruthSJ_num = 0;
	ifstream GroundTruthSJ_ifs(inputGroundTruthSJfile.c_str());
 	while(1)
	{
		if(GroundTruthSJ_ifs.eof())
			break;
		string tmpGroundTruthSJstr;
		getline(GroundTruthSJ_ifs, tmpGroundTruthSJstr);
		if(GroundTruthSJ_ifs.eof())
			break;
		if(tmpGroundTruthSJstr.substr(0, 3) != "chr")
			continue;
		int tmpGroundTruthSJchrNameInt, 
			tmpGroundTruthSJdonerEndPos, tmpGroundTruthSJacceptorStartPos;
		returnSpliceSite(tmpGroundTruthSJstr, tmpGroundTruthSJchrNameInt, 
			tmpGroundTruthSJdonerEndPos, tmpGroundTruthSJacceptorStartPos, indexInfo);
		int tmpJuncSize = tmpGroundTruthSJacceptorStartPos - tmpGroundTruthSJdonerEndPos - 1;
		if((tmpJuncSize < 50)||(tmpJuncSize > 300000))
			continue;
		groundTruthSJ_num ++;
		donerEndPosSetVec[tmpGroundTruthSJchrNameInt].insert(
			tmpGroundTruthSJdonerEndPos);
		acceptorStartPosSetVec[tmpGroundTruthSJchrNameInt].insert(
			tmpGroundTruthSJacceptorStartPos);
	}
	int groundTruth_donerSpliceSite_num = returnItemNumInSetVec(donerEndPosSetVec);
	int groundTruth_acceptorSpliceSite_num = returnItemNumInSetVec(acceptorStartPosSetVec);

	ofstream bothSpliceSitesCorrect_ofs(outputSJbothSpliceSitesCorrectFile.c_str());
	ofstream bothSpliceSitesIncorrect_ofs(outputSJbothSpliceSitesIncorrectFile.c_str());
	ofstream donerSpliceSiteCorrect_ofs(outputSJdonerSpliceSiteCorrectFile.c_str());
	ofstream acceptorSpliceSiteCorrect_ofs(outputSJacceptorSpliceSiteCorrectFile.c_str());
	
	int comparedSJ_num = 0;
	int bothSpliceSitesCorrect_num = 0;
	int bothSpliceSitesIncorrect_num = 0;
	int donerSpliceSiteCorrect_num = 0;
	int acceptorSpliceSiteCorrect_num = 0;
	ifstream comparedSJ_ifs(inputComparedSJfile.c_str());
	while(1)
	{
		if(comparedSJ_ifs.eof())
			break;
		string tmpComparedSJstr;
		getline(comparedSJ_ifs, tmpComparedSJstr);
		if(comparedSJ_ifs.eof())
			break;
		if(tmpComparedSJstr.substr(0, 3) != "chr")
			continue;
		int tmpComparedSJchrNameInt, 
			tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos;
		returnSpliceSite(tmpComparedSJstr, tmpComparedSJchrNameInt, 
			tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos, indexInfo);
		int tmpJunctionSize = tmpComparedSJacceptorStartPos - tmpComparedSJdonerEndPos - 1;
		if((tmpJunctionSize < 50)||(tmpJunctionSize > 300000))
			continue;
		int tmpSpliceSiteSearchResults = returnSpliceSiteSearchResults(
				donerEndPosSetVec, acceptorStartPosSetVec,
				offset, tmpComparedSJchrNameInt, 
				tmpComparedSJdonerEndPos, tmpComparedSJacceptorStartPos,
				donerEndPosSetVec_compared2_detected,
				acceptorStartPosSetVec_compared2_detected,
				donerEndPosSetVec_compared2_correct,
				acceptorStartPosSetVec_compared2_correct,
				donerEndPosSetVec_compared2_incorrect,
				acceptorStartPosSetVec_compared2_incorrect
				);		
		comparedSJ_num ++;
		if(tmpSpliceSiteSearchResults == 0)
		{
			bothSpliceSitesIncorrect_num ++;
			bothSpliceSitesIncorrect_ofs << tmpComparedSJstr << endl;
		}
		else if(tmpSpliceSiteSearchResults == 1)
		{
			bothSpliceSitesCorrect_num ++;
			bothSpliceSitesCorrect_ofs << tmpComparedSJstr << endl;
		}
		else if(tmpSpliceSiteSearchResults == 2)
		{
			donerSpliceSiteCorrect_num ++;
			donerSpliceSiteCorrect_ofs << tmpComparedSJstr << endl;
		}
		else if(tmpSpliceSiteSearchResults == 3)
		{
			acceptorSpliceSiteCorrect_num ++;
			acceptorSpliceSiteCorrect_ofs << tmpComparedSJstr << endl;
		}		
	}
	
	double bothSpliceSitesCorrect_num_perc = ((double)bothSpliceSitesCorrect_num / (double)comparedSJ_num) * 100;
	double bothSpliceSitesIncorrect_num_perc = ((double)bothSpliceSitesIncorrect_num / (double)comparedSJ_num) * 100;	
	double donerSpliceSiteCorrect_num_perc = ((double)donerSpliceSiteCorrect_num / (double)comparedSJ_num) * 100;
	double acceptorSpliceSiteCorrect_num_perc = ((double)acceptorSpliceSiteCorrect_num / (double)comparedSJ_num) * 100;

	int donerSpliceSite_compared2_correct_num = returnItemNumInSetVec(donerEndPosSetVec_compared2_correct);
	int donerSpliceSite_compared2_incorrect_num = returnItemNumInSetVec(donerEndPosSetVec_compared2_incorrect);
	int donerSpliceSite_compared2_num = donerSpliceSite_compared2_correct_num + donerSpliceSite_compared2_incorrect_num;
	int donerSpliceSite_compared2_detected_num = returnItemNumInSetVec(donerEndPosSetVec_compared2_detected);
	double donerSpliceSite_compared2_detected_num_perc = ((double)donerSpliceSite_compared2_detected_num / (double)groundTruth_donerSpliceSite_num) * 100;	
	double donerSpliceSite_compared2_correct_num_perc = ((double)donerSpliceSite_compared2_correct_num / (double)donerSpliceSite_compared2_num) * 100;
	double donerSpliceSite_compared2_incorrect_num_perc = ((double)donerSpliceSite_compared2_incorrect_num / (double)donerSpliceSite_compared2_num) * 100;

	int acceptorSpliceSite_compared2_correct_num = returnItemNumInSetVec(acceptorStartPosSetVec_compared2_correct);
	int acceptorSpliceSite_compared2_incorrect_num = returnItemNumInSetVec(acceptorStartPosSetVec_compared2_incorrect);
	int acceptorSpliceSite_compared2_num = acceptorSpliceSite_compared2_correct_num + acceptorSpliceSite_compared2_incorrect_num;
	int acceptorSpliceSite_compared2_detected_num = returnItemNumInSetVec(acceptorStartPosSetVec_compared2_detected);
	double acceptorSpliceSite_compared2_detected_num_perc = ((double)acceptorSpliceSite_compared2_detected_num / (double)groundTruth_acceptorSpliceSite_num) * 100;
	double acceptorSpliceSite_compared2_correct_num_perc = ((double)acceptorSpliceSite_compared2_correct_num / (double)acceptorSpliceSite_compared2_num) * 100; 
	double acceptorSpliceSite_compared2_incorrect_num_perc = ((double)acceptorSpliceSite_compared2_incorrect_num / (double)acceptorSpliceSite_compared2_num) * 100; 


	log_ofs << "groundTruth SJ #: " << groundTruthSJ_num << endl;
	log_ofs << "groundTruth doner splice site #: " << groundTruth_donerSpliceSite_num << endl;
	log_ofs << "groundTruth acceptor splice site #: " << groundTruth_acceptorSpliceSite_num << endl;
	log_ofs << endl;
	log_ofs << "compared SJ #: " << comparedSJ_num << endl;
	log_ofs << "compared doner splice site #: " << donerSpliceSite_compared2_num << endl;
	log_ofs << "compared acceptor splice site #: " << acceptorSpliceSite_compared2_num << endl;
	log_ofs << endl << "************************************************************************" << endl << endl;
	log_ofs << "bothSpliceSites Correct SJ #: " << bothSpliceSitesCorrect_num << " -- " << bothSpliceSitesCorrect_num_perc << "%" << endl;
	log_ofs << "donerSpliceSite Correct SJ #: " << donerSpliceSiteCorrect_num << " -- " << donerSpliceSiteCorrect_num_perc << "%" << endl;
	log_ofs << "acceptorSpliceSite Correct SJ #: " << acceptorSpliceSiteCorrect_num << " -- " << acceptorSpliceSiteCorrect_num_perc << "%" << endl;
	log_ofs << "bothSpliceSites Incorrect SJ #: " << bothSpliceSitesIncorrect_num << " -- " << bothSpliceSitesIncorrect_num_perc << "%" << endl;
	log_ofs << endl << "************************************************************************" << endl << endl;
	log_ofs << "donerSpliecSite detected #: " << donerSpliceSite_compared2_detected_num << " -- " << donerSpliceSite_compared2_detected_num_perc << "%" << endl;
	log_ofs << "donerSpliecSite correct #: " << donerSpliceSite_compared2_correct_num << " -- " << donerSpliceSite_compared2_correct_num_perc << "%" << endl; 
	log_ofs << "donerSpliecSite incorrect #: " << donerSpliceSite_compared2_incorrect_num << " -- " << donerSpliceSite_compared2_incorrect_num_perc << "%" << endl; 
	log_ofs << endl;
	log_ofs << "acceptorSpliecSite detected #: " << acceptorSpliceSite_compared2_detected_num << " -- " << acceptorSpliceSite_compared2_detected_num_perc << "%" << endl;
	log_ofs << "acceptorSpliecSite correct #: " << acceptorSpliceSite_compared2_correct_num << " -- " << acceptorSpliceSite_compared2_correct_num_perc << "%" << endl; 
	log_ofs << "acceptorSpliecSite incorrect #: " << acceptorSpliceSite_compared2_incorrect_num << " -- " << acceptorSpliceSite_compared2_incorrect_num_perc << "%" << endl; 	
	delete indexInfo;
	GroundTruthSJ_ifs.close();
	comparedSJ_ifs.close();
	parameter_ifs.close();
	bothSpliceSitesCorrect_ofs.close();
	bothSpliceSitesIncorrect_ofs.close();
	donerSpliceSiteCorrect_ofs.close();
	acceptorSpliceSiteCorrect_ofs.close();
	log_ofs.close();
	return 0;
}