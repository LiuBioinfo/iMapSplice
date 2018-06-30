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

#include "general/read_block_test.h"
#include "general/index_info.h"
#include "general/splice_info.h"
#include "general/alignmentToJunc_supportNum.h"

using namespace std;

int getSJdistanceFromSJstr_gt(const string& SJstr)
{
	int tabLocation1 = SJstr.find('\t', 0);
	int tabLocation2 = SJstr.find('\t', tabLocation1+1);
	int tabLocation3 = SJstr.find('\t', tabLocation2+1);	

	string donerEndPosStr = SJstr.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
	//cout << "donerEndPosStr: " << donerEndPosStr << endl;
	string acceptorStartPosStr = SJstr.substr(tabLocation2+1, tabLocation3- tabLocation2-1);
	//cout << "acceptorStartPosStr: " << acceptorStartPosStr << endl;
	int donerEndPos = atoi(donerEndPosStr.c_str());
	//cout << "donerEndPos: " << donerEndPos << endl;
	int acceptorStartPos = atoi(acceptorStartPosStr.c_str());
	//cout << "acceptorStartPos: " << acceptorStartPos << endl;
	int SJdistance = acceptorStartPos - donerEndPos - 1;	
	return SJdistance;
}

int getSJdistanceFromSJstr(const string& SJstr)
{
	int tabLocation1 = SJstr.find('\t', 0);
	int tabLocation2 = SJstr.find('\t', tabLocation1+1);
	int tabLocation3 = SJstr.find('\t', tabLocation2+1);	

	string donerEndPosStr = SJstr.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
	//cout << "donerEndPosStr: " << donerEndPosStr << endl;
	string acceptorStartPosStr = SJstr.substr(tabLocation2+1, tabLocation3- tabLocation2-1);
	//cout << "acceptorStartPosStr: " << acceptorStartPosStr << endl;
	int donerEndPos = atoi(donerEndPosStr.c_str());
	//cout << "donerEndPos: " << donerEndPos << endl;
	int acceptorStartPos = atoi(acceptorStartPosStr.c_str());
	//cout << "acceptorStartPos: " << acceptorStartPos << endl;
	int SJdistance = acceptorStartPos - donerEndPos - 1;	
	//cout << "SJdistance: " << SJdistance << endl;
	return SJdistance;
}

int getSJsupportNumFromSJstr(const string& SJstr)
{
	int tabLocation1 = SJstr.find('\t', 0);
	int tabLocation2 = SJstr.find('\t', tabLocation1+1);
	int tabLocation3 = SJstr.find('\t', tabLocation2+1);
	int tabLocation4 = SJstr.find('\t', tabLocation3+1);
	int tabLocation5 = SJstr.find('\t', tabLocation4+1);
	string supportNumStr = SJstr.substr(tabLocation4+1, tabLocation5- tabLocation4-1);
	int supportNum = atoi(supportNumStr.c_str());
	return supportNum;
}

int main(int argc, char** argv)
{
	if(argc < 9)
	{
		cout << "Executable inputIndexFile offset outputSJstatsFile inputGroundTruthSJ SJ_size_min SJ_size_max SupportNum_max aligner_1 inputSamSJ_1 ..." << endl;
		exit(1);
	}
	string indexParameterFileStr = argv[1];
	indexParameterFileStr = indexParameterFileStr + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());	
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	int chromNum = indexInfo->returnChromNum();
	parameter_ifs.close();

	string offsetStr = argv[2];
	int offset = atoi(offsetStr.c_str());	

	string outputSJstatsFile = argv[3];
	ofstream output_ofs(outputSJstatsFile.c_str());

	string inputGroundTruthSJfile = argv[4];
	ifstream inputGroundTruthSJ_ifs(inputGroundTruthSJfile.c_str());

	string SJ_size_min_str = argv[5];
	string SJ_size_max_str = argv[6];

	string supportNum_max_str = argv[7];

	int SJ_size_min = atoi(SJ_size_min_str.c_str());
	int SJ_size_max = atoi(SJ_size_max_str.c_str());
	int supportNum_max = atoi(supportNum_max_str.c_str());

	int groundTruthSJnum = 0;

	vector<string> alignerNameVec;
	vector<string> inputSJfileVec;
	cout << "*****************************************************" << endl;
	for(int tmp = 8; tmp < argc; )
	{
		string tmpAlignName = argv[tmp];
		string tmpInputSJfile = argv[tmp+1];
		cout << "tmpAligner: " << endl << tmpAlignName << endl;
		cout << "tmpInputSJfile: " << endl << tmpInputSJfile << endl;		
		alignerNameVec.push_back(tmpAlignName);
		inputSJfileVec.push_back(tmpInputSJfile);
		tmp += 2;
	}

	int* SJnum_found_supportNum_duplicate = (int*)malloc((supportNum_max+1) * (alignerNameVec.size()) * sizeof(int));
	int* SJnum_found_supportNum = (int*)malloc((supportNum_max+1) * (alignerNameVec.size()) * sizeof(int));
	int* SJnum_notFound_supportNum = (int*)malloc((supportNum_max+1) * (alignerNameVec.size()) * sizeof(int));
	
	//CHECK GROUND TRUTH sj FILE
	int tmpGroundTruthSJnum = 0;
	while(!inputGroundTruthSJ_ifs.eof())
	{
		string SJstr;
		getline(inputGroundTruthSJ_ifs, SJstr);
		if(inputGroundTruthSJ_ifs.eof())
			break;
		int tmpSJdistance = getSJdistanceFromSJstr_gt(SJstr);
		int tmpSJsupportNum = getSJsupportNumFromSJstr(SJstr);
		//getInfoFromSJstr_groundTruth(SJstr, tmpSJdistance, tmpSJsupportNum);

		if((tmpSJdistance >= SJ_size_min)&&(tmpSJdistance <= SJ_size_max))
			tmpGroundTruthSJnum ++;
	}
	inputGroundTruthSJ_ifs.close();

	string SJ_file_groundTruth_str = inputGroundTruthSJfile;
	AlignmentToJunc_supportNum_Info* juncInfo_groundTruth = new AlignmentToJunc_supportNum_Info();
	juncInfo_groundTruth->initiateAlignmentToJuncInfo(chromNum);
	int juncInfo_groundTruth_validSJnum, juncInfo_groundTruth_invalidSJnum;// juncInfo_2_validSJnum, juncInfo_2_invalidSJnum;
	string invalidJuncFile_groundTruth = SJ_file_groundTruth_str + ".invalidSJ";
	juncInfo_groundTruth->getAlignmentToJuncInfoFromSJfile(SJ_file_groundTruth_str, indexInfo,
		SJ_size_min, SJ_size_max, juncInfo_groundTruth_validSJnum, juncInfo_groundTruth_invalidSJnum,
		invalidJuncFile_groundTruth);

	cout << "*****************************************************" << endl;
	groundTruthSJnum = tmpGroundTruthSJnum;
	cout << "groundTruthSJnum: " << groundTruthSJnum << endl;
	cout << "*****************************************************" << endl;
	for(int tmpAligner = 0; tmpAligner < alignerNameVec.size(); tmpAligner++)
	{
		string tmp_SJ_file_str = inputSJfileVec[tmpAligner];

		// get alignment 2 juncInfo From SJ file
		cout << "getAlignmentToJuncInfoFromSJfile file  ..." << endl;
		string tmp_compare_file_str = tmp_SJ_file_str + ".SJcomp";
		int juncInfo_tmp_validSJnum, juncInfo_tmp_invalidSJnum;
		AlignmentToJunc_supportNum_Info* tmpJuncInfo = new AlignmentToJunc_supportNum_Info(); 
		tmpJuncInfo->initiateAlignmentToJuncInfo(chromNum);
		string invalidJuncFile_tmp = tmp_SJ_file_str + ".invalidSJ";
		tmpJuncInfo->getAlignmentToJuncInfoFromSJfile(tmp_SJ_file_str, indexInfo, 
			SJ_size_min, SJ_size_max, juncInfo_tmp_validSJnum, juncInfo_tmp_invalidSJnum,
			invalidJuncFile_tmp);
		tmpJuncInfo->comparedToOtherAlignmentToJuncInfo(juncInfo_groundTruth, offset, indexInfo,
			tmp_compare_file_str, SJ_size_min, SJ_size_max, 
			juncInfo_groundTruth_invalidSJnum, juncInfo_tmp_invalidSJnum);
		delete tmpJuncInfo;

		string tmpOutputSJfile_prefix = tmp_compare_file_str + "." + int_to_str(offset);
		string tmpInputSJfile_found_duplicate = tmpOutputSJfile_prefix + ".correctSJ";
		string tmpInputSJfile_found = tmpOutputSJfile_prefix + ".discoveredSJ";
		string tmpInputSJfile_notFound = tmpOutputSJfile_prefix + ".incorrectSJ";
		ifstream tmpInputSJfile_found_ifs_duplicate(tmpInputSJfile_found_duplicate.c_str());		
		ifstream tmpInputSJfile_found_ifs(tmpInputSJfile_found.c_str());
		ifstream tmpInputSJfile_notFound_ifs(tmpInputSJfile_notFound.c_str());
		
		//check found SJ file
		cout << "start to check discoveredSJ file ..." << endl;
		while(!tmpInputSJfile_found_ifs_duplicate.eof())
		{
			string SJstr;
			getline(tmpInputSJfile_found_ifs_duplicate, SJstr);
			//cout << "SJstr: " << SJstr << endl;
			//break;
			if(tmpInputSJfile_found_ifs_duplicate.eof())
				break;
			int tmpSJdistance = getSJdistanceFromSJstr(SJstr);	
			int tmpSJsupportNum = getSJsupportNumFromSJstr(SJstr);
			//cout << "tmpSJdistance: " << tmpSJdistance << endl;
			//cout << "tmpSJsupportNum: " << tmpSJsupportNum << endl;
			if( (tmpSJdistance >= SJ_size_min) && (tmpSJdistance <= SJ_size_max))
			{
				int tmpMax = tmpSJsupportNum;
				if(tmpMax > supportNum_max)
					tmpMax = supportNum_max;
				for(int tmpSJsupportNum_tmp = 1 + tmpAligner * (supportNum_max+1); tmpSJsupportNum_tmp <= tmpMax + tmpAligner * (supportNum_max+1); tmpSJsupportNum_tmp++)
				{				
					SJnum_found_supportNum_duplicate[tmpSJsupportNum_tmp]++;
				}
			}
		}

		cout << "start to check correct splice junctions" << endl;
		// check found groundTruth SJ file
		while(!tmpInputSJfile_found_ifs.eof())
		{
			string SJstr;
			getline(tmpInputSJfile_found_ifs, SJstr);
			if(tmpInputSJfile_found_ifs.eof())
				break;
			int tmpSJdistance = getSJdistanceFromSJstr(SJstr);	
			int tmpSJsupportNum = getSJsupportNumFromSJstr(SJstr);
			if( (tmpSJdistance >= SJ_size_min) && (tmpSJdistance <= SJ_size_max))
			{
				int tmpMax = tmpSJsupportNum;
				if(tmpMax > supportNum_max)
					tmpMax = supportNum_max;
				for(int tmpSJsupportNum_tmp = 1 + tmpAligner * (supportNum_max+1); tmpSJsupportNum_tmp <= tmpMax + tmpAligner * (supportNum_max+1); tmpSJsupportNum_tmp++)
				{	
					SJnum_found_supportNum[tmpSJsupportNum_tmp]++;
				}
			}
		}

		cout << "start to check incorrect splice junctions" << endl;

		// check unfound SJ file
		while(!tmpInputSJfile_notFound_ifs.eof())
		{
			string SJstr;
			getline(tmpInputSJfile_notFound_ifs, SJstr);
			if(tmpInputSJfile_notFound_ifs.eof())
				break;
			int tmpSJdistance = getSJdistanceFromSJstr(SJstr);
			int tmpSJsupportNum = getSJsupportNumFromSJstr(SJstr);
			if( (tmpSJdistance >= SJ_size_min) && (tmpSJdistance <= SJ_size_max))
			{
				int tmpMax = tmpSJsupportNum;
				if(tmpMax > supportNum_max)
					tmpMax = supportNum_max;
				for(int tmpSJsupportNum_tmp = 1 + tmpAligner * (supportNum_max+1); tmpSJsupportNum_tmp <= tmpMax + tmpAligner * (supportNum_max+1); tmpSJsupportNum_tmp++)
				{
					SJnum_notFound_supportNum[tmpSJsupportNum_tmp]++;
				}
			}
		}
		cout << "start to output stats ..." << endl << "*********************************************" << endl;
		string outputSJstatsFile_tmpAligner_sensitivity 
			= outputSJstatsFile + "." + alignerNameVec[tmpAligner] + ".sensitivity";
		string outputSJstatsFile_tmpAligner_consistentNum
			= outputSJstatsFile + "." + alignerNameVec[tmpAligner] + ".consistentNum";			
		string outputSJstatsFile_tmpAligner_specificity
			= outputSJstatsFile + "." + alignerNameVec[tmpAligner] + ".specificity";
		ofstream outputSJstatsFile_tmpAligner_sensitivity_ofs(outputSJstatsFile_tmpAligner_sensitivity.c_str());
		ofstream outputSJstatsFile_tmpAligner_consistentNum_ofs(outputSJstatsFile_tmpAligner_consistentNum.c_str());
		ofstream outputSJstatsFile_tmpAligner_specificity_ofs(outputSJstatsFile_tmpAligner_specificity.c_str());
		for(int tmpSupportNum = 1; tmpSupportNum <= supportNum_max; tmpSupportNum++)
		{			
			int tmpSJtotalNum_aligner = SJnum_found_supportNum_duplicate[tmpSupportNum + tmpAligner * (supportNum_max+1)] 
				+ SJnum_notFound_supportNum[tmpSupportNum + tmpAligner * (supportNum_max+1)];
			//cout << "tmpSJtotalNum_aligner: " << tmpSJtotalNum_aligner << endl << "aligner: " << alignerNameVec[tmpAligner] << endl;;
			double tmpSensitivity = ((double)SJnum_found_supportNum[tmpSupportNum + tmpAligner * (supportNum_max+1)]/((double)groundTruthSJnum))*100;
			double tmpSpecificity = ((double)SJnum_found_supportNum_duplicate[tmpSupportNum + tmpAligner * (supportNum_max+1)]/((double)tmpSJtotalNum_aligner))*100;
			output_ofs << tmpSensitivity << " " 
				<< 100.00-tmpSpecificity << " " << alignerNameVec[tmpAligner] << endl;
			outputSJstatsFile_tmpAligner_sensitivity_ofs << tmpSensitivity << endl;
			outputSJstatsFile_tmpAligner_consistentNum_ofs << SJnum_found_supportNum_duplicate[tmpSupportNum + tmpAligner * (supportNum_max+1)] << endl;
			outputSJstatsFile_tmpAligner_specificity_ofs << tmpSpecificity << endl;
		}

		outputSJstatsFile_tmpAligner_sensitivity_ofs.close();
		outputSJstatsFile_tmpAligner_consistentNum_ofs.close();
		outputSJstatsFile_tmpAligner_specificity_ofs.close();
		tmpInputSJfile_found_ifs_duplicate.close();
		tmpInputSJfile_found_ifs.close();
		tmpInputSJfile_notFound_ifs.close();
	}

	free(SJnum_found_supportNum_duplicate);
	free(SJnum_found_supportNum);
	free(SJnum_notFound_supportNum);

	return 0;
}