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

#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/alignInferJunctionHash_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 7)
	{
		cout << "Executable indexInfoFolderPath outputFolder name_compared2 SJ_compared2 ";
		cout << " name_1 SJ_1 (name_2 SJ_2 ...)" << endl;
		exit(1);
	}
	int toCompareAlignerSJnum = (argc - 5)/2;
	if((toCompareAlignerSJnum * 2) != (argc - 5))
	{
		cout << "incorrect argv #" << endl;
		exit(1);
	}	

	string outputFolderStr = argv[2];
	outputFolderStr += "/";
	string mkdir= "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());	
	string logFile = outputFolderStr + "log.txt";
	ofstream log_ofs(logFile.c_str());
	string SJsupNumFile = outputFolderStr + "SJsupNum.txt";
	ofstream SJsupNum_ofs(SJsupNumFile.c_str());
	string SJsupNumFile_sharedByAll = outputFolderStr + "SJsupNum_sharedByAll.txt";
	ofstream SJsupNum_sharedByAll_ofs(SJsupNumFile_sharedByAll.c_str());
	string SJsupNumFile_sharedByAll_log10 = outputFolderStr + "SJsupNum_sharedByAll_log10.txt";
	ofstream SJsupNum_sharedByAll_log10_ofs(SJsupNumFile_sharedByAll_log10.c_str());

	cout << "start to initiate indexInfo" << endl;
	log_ofs << "start to initiate indexInfo" << endl;
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());	
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	parameter_ifs.close();

	cout << "start to initiate ground truth alignInferJuncHash" << endl;
	log_ofs << "start to initiate ground truth alignInferJuncHash" << endl;
	string groundTruthName = argv[3];
	string groundTruthSJfilePath = argv[4];
	AlignInferJunctionHash_Info* GTjuncHash 
		= new AlignInferJunctionHash_Info();
	GTjuncHash->initiateAlignInferJunctionInfo(chromNum);
	GTjuncHash->insertJuncFromJuncFile_chrNamePos_supportNum(
			groundTruthSJfilePath, indexInfo);

	cout << "start to read aligner names and SJ file paths" << endl;
	log_ofs << "start to read aligner names and SJ file paths" << endl;
	vector<string> alignerNameVec;
	vector<string> inputSJfileVec;
	vector<AlignInferJunctionHash_Info*> juncHashVec;
	for(int tmp = 5; tmp < argc; )
	{
		string tmpAlignName = argv[tmp];
		string tmpInputSJfile = argv[tmp+1];
		cout << "tmpAligner: " << endl << tmpAlignName << endl;
		cout << "tmpInputSJfile: " << endl << tmpInputSJfile << endl;		
		log_ofs << "tmpAligner: " << endl << tmpAlignName << endl;
		log_ofs << "tmpInputSJfile: " << endl << tmpInputSJfile << endl;		
		alignerNameVec.push_back(tmpAlignName);
		inputSJfileVec.push_back(tmpInputSJfile);
		log_ofs << "tmpAlignInferJuncHash initiates ... " << endl;
		AlignInferJunctionHash_Info* tmpJuncHash 
			= new AlignInferJunctionHash_Info();
		tmpJuncHash->initiateAlignInferJunctionInfo(chromNum);
		log_ofs << "insert SJs from juncFile 2 juncHash" << endl << endl;
		tmpJuncHash->insertJuncFromJuncFile_chrNamePos_supportNum(
			tmpInputSJfile, indexInfo);
		juncHashVec.push_back(tmpJuncHash);
		tmp += 2;
	}	

	cout << "start to do compared2otherAlignInferJuncHashVecAndReturnSJsupNumVec " << endl;
	log_ofs << "start to do compared2otherAlignInferJuncHashVecAndReturnSJsupNumVec " << endl;
	int min_intron_size = 50;
	int max_intron_size = 300000;
	int totalJuncNum_groundTruth = 0;
	int validJuncNum_groundTruth = 0;
	int invalidJuncNum_groundTruth = 0;
	log_ofs << "min_intron_size: " << min_intron_size << endl;
	log_ofs << "max_intron_size: " << max_intron_size << endl << endl;
	vector<int> chrNameIntVecInGroundTruthJuncHash;
	vector<int> donerSpliceSiteVecInGroundTruthJuncHash;
	vector<int> acceptorSpliceSiteVecInGroundTruthJuncHash;
	vector<int> supNumInGroundTruthJuncHash;
	vector< vector<int> > supNumVecInJuncHashVec2compare;
	vector<bool> sharedByAllOrNotBoolVec;	
	GTjuncHash->compared2otherAlignInferJuncHashVecAndReturnSJsupNumVec(
		juncHashVec, chrNameIntVecInGroundTruthJuncHash,
		donerSpliceSiteVecInGroundTruthJuncHash,
		acceptorSpliceSiteVecInGroundTruthJuncHash, supNumInGroundTruthJuncHash, 
		supNumVecInJuncHashVec2compare, sharedByAllOrNotBoolVec,
		min_intron_size, max_intron_size, totalJuncNum_groundTruth,
		validJuncNum_groundTruth, invalidJuncNum_groundTruth, indexInfo);


	log_ofs << "totalJuncNum_groundTruth: " << totalJuncNum_groundTruth << endl;
	log_ofs << "validJuncNum_groundTruth: " << validJuncNum_groundTruth << endl;
	log_ofs << "invalidJuncNum_groundTruth: " << invalidJuncNum_groundTruth << endl << endl;

	cout << "start to output SJ sup num results " << endl;
	log_ofs << "start to output SJ sup num results " << endl;
	SJsupNum_ofs << "chr\tdonerEndPos\tacceptorStartPos\t" << groundTruthName << "_raw";
	SJsupNum_sharedByAll_ofs << "chr\tdonerEndPos\tacceptorStartPos\t" << groundTruthName << "_raw";
	SJsupNum_sharedByAll_log10_ofs << "chr\tdonerEndPos\tacceptorStartPos\t" << groundTruthName << "_log10";
	for(int tmp = 0; tmp < toCompareAlignerSJnum; tmp++)
	{
		SJsupNum_ofs << "\t" << alignerNameVec[tmp] << "_raw";
		SJsupNum_sharedByAll_ofs << "\t" << alignerNameVec[tmp] << "_raw";
		SJsupNum_sharedByAll_log10_ofs << "\t" << alignerNameVec[tmp] << "_log10";
	}
	SJsupNum_ofs << endl;
	SJsupNum_sharedByAll_ofs << endl;
	SJsupNum_sharedByAll_log10_ofs << endl;

	cout << "start to output SJsupNum results for all aligners in the same files" << endl;
	log_ofs << "start to output SJsupNum results for all aligners in the same files" << endl;
	for(int tmp = 0; tmp < validJuncNum_groundTruth; tmp++)
	{
		//cout << "tmp in validJuncNum_groundTruth: " << tmp << endl;
		int tmpChrNameInt = chrNameIntVecInGroundTruthJuncHash[tmp];
		string tmpChrNameStr = indexInfo->returnChrNameStr(tmpChrNameInt);
		//cout << "tmpChrNameStr: " << tmpChrNameStr << endl;
		int tmpDonerSpliceSite = donerSpliceSiteVecInGroundTruthJuncHash[tmp];
		int tmpAcceptorSpliceSite = acceptorSpliceSiteVecInGroundTruthJuncHash[tmp];
		//cout << "tmpDonerSpliceSite: " << tmpDonerSpliceSite << endl;
		int tmpSupNum = supNumInGroundTruthJuncHash[tmp];
		//cout << "tmpSupNum: " << tmpSupNum << endl;
		bool sharedByAllOrNot_bool = sharedByAllOrNotBoolVec[tmp];
		//cout << "sharedByAllOrNot_bool: " << sharedByAllOrNot_bool << endl;
		if(sharedByAllOrNot_bool)
		{	
			double tmpSupNum_log10 = log10((double)tmpSupNum);
			//cout << "tmpSupNum_log10: " << tmpSupNum_log10 << endl;
			SJsupNum_ofs << tmpChrNameStr << "\t" << tmpDonerSpliceSite << "\t" 
				<< tmpAcceptorSpliceSite << "\t" << tmpSupNum;
			SJsupNum_sharedByAll_ofs << tmpChrNameStr << "\t" << tmpDonerSpliceSite << "\t" 
				<< tmpAcceptorSpliceSite << "\t" << tmpSupNum;
			SJsupNum_sharedByAll_log10_ofs << tmpChrNameStr << "\t" << tmpDonerSpliceSite << "\t" 
				<< tmpAcceptorSpliceSite << "\t" << tmpSupNum_log10;
			for(int tmpAligner = 0; tmpAligner < toCompareAlignerSJnum; tmpAligner++)
			{
				//cout << "tmpAligner: " << tmpAligner << endl;
				//cout << "supNumVecInJuncHashVec2compare.size(): " << supNumVecInJuncHashVec2compare.size() << endl;
				//cout << "(supNumVecInJuncHashVec2compare[tmp]).size(): " << (supNumVecInJuncHashVec2compare[tmp]).size() << endl;
				int tmpAlignerJuncSupNum = (supNumVecInJuncHashVec2compare[tmp])[tmpAligner];
				//cout << "tmpAlignerJuncSupNum: " << tmpAlignerJuncSupNum << endl;
				double tmpAlignerJuncSupNum_log10 = log10((double)tmpAlignerJuncSupNum);
				SJsupNum_ofs << "\t" << tmpAlignerJuncSupNum;
				SJsupNum_sharedByAll_ofs << "\t" << tmpAlignerJuncSupNum;
				SJsupNum_sharedByAll_log10_ofs << "\t" << tmpAlignerJuncSupNum_log10;
			}
			SJsupNum_ofs << endl;
			SJsupNum_sharedByAll_ofs << endl;
			SJsupNum_sharedByAll_log10_ofs << endl;
		}
		else
		{
			SJsupNum_ofs << tmpChrNameStr << "\t" << tmpDonerSpliceSite << "\t" 
				<< tmpAcceptorSpliceSite << "\t" << tmpSupNum;
			for(int tmpAligner = 0; tmpAligner < toCompareAlignerSJnum; tmpAligner++)
			{
				int tmpAlignerJuncSupNum = (supNumVecInJuncHashVec2compare[tmp])[tmpAligner];
				SJsupNum_ofs << "\t" << tmpAlignerJuncSupNum;
			}
			SJsupNum_ofs << endl;
		}
	}

	cout << "start to output SJsupNum results for each aligner in separate files" << endl;
	log_ofs << "start to output SJsupNum results for each aligner in separate files" << endl;
	for(int tmp = 0; tmp < toCompareAlignerSJnum; tmp++)
	{
		string tmpAlignerName = alignerNameVec[tmp];
		string tmpSJsupNumFile_tmpAligner = outputFolderStr + "SJsupNum_sharedByAll_log10." + tmpAlignerName + ".txt";
		ofstream tmpSJsupNum_tmpAligner_ofs(tmpSJsupNumFile_tmpAligner.c_str());
		for(int tmpJunc = 0; tmpJunc < validJuncNum_groundTruth; tmpJunc++)
		{
			bool sharedByAllOrNot_bool = sharedByAllOrNotBoolVec[tmpJunc];
			if(sharedByAllOrNot_bool)
			{
				int tmpChrNameInt = chrNameIntVecInGroundTruthJuncHash[tmpJunc];
				string tmpChrNameStr = indexInfo->returnChrNameStr(tmpChrNameInt);
				int tmpDonerSpliceSite = donerSpliceSiteVecInGroundTruthJuncHash[tmpJunc];
				int tmpAcceptorSpliceSite = acceptorSpliceSiteVecInGroundTruthJuncHash[tmpJunc];
				int tmpSupNum_groundTruth = supNumInGroundTruthJuncHash[tmpJunc];				
				double tmpSupNum_groundTruth_log10 = log10((double)tmpSupNum_groundTruth);
				int tmpSupNum_tmpAligner = (supNumVecInJuncHashVec2compare[tmpJunc])[tmp];
				double tmpSupNum_tmpAligner_log10 = log10((double)tmpSupNum_tmpAligner);
				tmpSJsupNum_tmpAligner_ofs << tmpChrNameStr << "\t" << tmpDonerSpliceSite << "\t"
					<< tmpAcceptorSpliceSite << "\t" << tmpSupNum_groundTruth << "\t"
					<< tmpSupNum_groundTruth_log10 << "\t" << tmpSupNum_tmpAligner_log10 << "\t"
					<< tmpAlignerName << endl;
			}
		}
		tmpSJsupNum_tmpAligner_ofs.close();
	}

	SJsupNum_sharedByAll_log10_ofs.close();
	SJsupNum_sharedByAll_ofs.close();
	SJsupNum_ofs.close();
	log_ofs.close();
	delete indexInfo;
	return 0;
}