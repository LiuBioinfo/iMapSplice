// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ALIGNINFERJUNCTIONHASH_INFO_VEC_H
#define ALIGNINFERJUNCTIONHASH_INFO_VEC_H

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

#include "alignInfer_info.h"
#include "alignInferJunctionHash_info.h"

class AlignInferJunctionHash_Info_Vec
{
private:
	vector<AlignInferJunctionHash_Info*> alignInferJuncHashInfoVec;

	//AlignInferJunctionHash_Info alignInferJuncHashInfo_merged;
public:
	AlignInferJunctionHash_Info_Vec()
	{}

	void freeMemory()
	{
		for(int tmp = 0; tmp < alignInferJuncHashInfoVec.size(); tmp++)
		{
			delete alignInferJuncHashInfoVec[tmp];
		}
	}

	void initiateAlignInferJunctionHashInfoVec(
		int alignInferJuncHashInfoVecSize, int chromNum)
	{
		for(int tmp = 0; tmp < alignInferJuncHashInfoVecSize; tmp++)
		{
			AlignInferJunctionHash_Info* tmpAlignInferJuncHashInfo 
				= new AlignInferJunctionHash_Info();
			tmpAlignInferJuncHashInfo->initiateAlignInferJunctionInfo(chromNum); 
			alignInferJuncHashInfoVec.push_back(tmpAlignInferJuncHashInfo);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePosOnly_parallel(
		vector<string>& alignmentFileVec, Index_Info* indexInfo,
		int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePosOnly_parallel(
				alignmentFileVec[tmp], indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePos_supportNum_parallel(
		vector<string>& alignmentFileVec, Index_Info* indexInfo,
		int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePos_supportNum_parallel(
				alignmentFileVec[tmp], indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePos_supportNum_anchorSize_parallel(
		vector<string>& alignmentFileVec, Index_Info* indexInfo,
		int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePos_supportNum_anchorSize_parallel(
				alignmentFileVec[tmp], indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
		}
	}

	void insertJuncFromAlignmentFileVec_chrNamePos_supportNum_anchorSize_XM_parallel(
		vector<string>& alignmentFileVec, Index_Info* indexInfo,
		int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		for(int tmp = 0; tmp < alignmentFileVec.size(); tmp++)
		{
			cout << "tmpAlignmentFile: " << alignmentFileVec[tmp] << endl;
			insertJuncFromAlignmentFile_chrNamePos_supportNum_anchorSize_XM_parallel(
				alignmentFileVec[tmp], indexInfo, alignInferJuncHashInfoVecSize, log_ofs);
		}
	}

	void insertJuncFromAlignmentFile_chrNamePosOnly_parallel(
		string& alignmentFile, Index_Info* indexInfo, int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		int normalRecordNum = 2000000;
		vector<string> readAlignmentSAMstrVec(normalRecordNum);
		
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;
		int alignmentTotalNum;

		string tmpLineStr;
		for(tmpTurn = 0; ; tmpTurn++)
		{
			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				getline(sam_ifs, tmpLineStr);
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				readAlignmentSAMstrVec[recordNumTmp] = tmpLineStr;
			}

			alignmentTotalNum += realRecordNum;

			// log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);		
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// log_ofs << "start to process alignments, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << "start to process alignments, turn: " << tmpTurn + 1 << endl;

			omp_set_num_threads(alignInferJuncHashInfoVecSize);

			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
			{
				if((readAlignmentSAMstrVec[tmpOpenMP] == "")
					||(readAlignmentSAMstrVec[tmpOpenMP].at(0) == '@'))
					continue;
				int threadNO = omp_get_thread_num();
				//cout << "start to process: " << endl
				//	<< readAlignmentSAMstrVec[tmpOpenMP] << endl;
				alignInferJuncHashInfoVec[threadNO]->
					getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePosOnly(
						readAlignmentSAMstrVec[tmpOpenMP], indexInfo);
			}
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_chrNamePos_supportNum_parallel(
		string& alignmentFile, Index_Info* indexInfo, int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		int normalRecordNum = 2000000;
		vector<string> readAlignmentSAMstrVec(normalRecordNum);
		
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;
		int alignmentTotalNum;

		string tmpLineStr;
		for(tmpTurn = 0; ; tmpTurn++)
		{
			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				getline(sam_ifs, tmpLineStr);
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				readAlignmentSAMstrVec[recordNumTmp] = tmpLineStr;
			}

			alignmentTotalNum += realRecordNum;

			// log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);		
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// log_ofs << "start to process alignments, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << "start to process alignments, turn: " << tmpTurn + 1 << endl;

			omp_set_num_threads(alignInferJuncHashInfoVecSize);

			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
			{
				if((readAlignmentSAMstrVec[tmpOpenMP] == "")
					||(readAlignmentSAMstrVec[tmpOpenMP].at(0) == '@'))
					continue;
				int threadNO = omp_get_thread_num();
				//cout << "start to process: " << endl
				//	<< readAlignmentSAMstrVec[tmpOpenMP] << endl;
				alignInferJuncHashInfoVec[threadNO]->
					getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum(
						readAlignmentSAMstrVec[tmpOpenMP], indexInfo);
			}
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_chrNamePos_supportNum_anchorSize_parallel(
		string& alignmentFile, Index_Info* indexInfo, int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		int normalRecordNum = 2000000;
		vector<string> readAlignmentSAMstrVec(normalRecordNum);
		
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;
		int alignmentTotalNum;

		string tmpLineStr;
		for(tmpTurn = 0; ; tmpTurn++)
		{
			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				getline(sam_ifs, tmpLineStr);
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				readAlignmentSAMstrVec[recordNumTmp] = tmpLineStr;
			}

			alignmentTotalNum += realRecordNum;

			// log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);		
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// log_ofs << "start to process alignments, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << "start to process alignments, turn: " << tmpTurn + 1 << endl;

			omp_set_num_threads(alignInferJuncHashInfoVecSize);

			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
			{
				if((readAlignmentSAMstrVec[tmpOpenMP] == "")
					||(readAlignmentSAMstrVec[tmpOpenMP].at(0) == '@'))
					continue;
				int threadNO = omp_get_thread_num();
				//cout << "start to process: " << endl
				//	<< readAlignmentSAMstrVec[tmpOpenMP] << endl;
				alignInferJuncHashInfoVec[threadNO]->
					getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize(
						readAlignmentSAMstrVec[tmpOpenMP], indexInfo);
			}
		}
		sam_ifs.close();
	}

	void insertJuncFromAlignmentFile_chrNamePos_supportNum_anchorSize_XM_parallel(
		string& alignmentFile, Index_Info* indexInfo, int alignInferJuncHashInfoVecSize, ofstream& log_ofs)
	{
		ifstream sam_ifs(alignmentFile.c_str());
		int normalRecordNum = 2000000;
		vector<string> readAlignmentSAMstrVec(normalRecordNum);
		
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;
		int alignmentTotalNum;

		string tmpLineStr;
		for(tmpTurn = 0; ; tmpTurn++)
		{
			if(EndOfRecord)
				break;
			int recordNum = normalRecordNum;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "start to read SAM file, turn: " << tmpTurn + 1 << endl;
			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				getline(sam_ifs, tmpLineStr);
				if(sam_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}
				readAlignmentSAMstrVec[recordNumTmp] = tmpLineStr;
			}

			alignmentTotalNum += realRecordNum;

			// log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			// nowtime = time(NULL);
			// local = localtime(&nowtime);		
			// log_ofs << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// log_ofs << "start to process alignments, turn: " << tmpTurn + 1 << endl;
			// cout << endl << "[" << asctime(local) 
			// 	<< "finish reading SAM file, turn: " << tmpTurn + 1 << endl;
			// cout << "start to process alignments, turn: " << tmpTurn + 1 << endl;

			omp_set_num_threads(alignInferJuncHashInfoVecSize);

			#pragma omp parallel for schedule(dynamic)
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP ++)
			{
				if((readAlignmentSAMstrVec[tmpOpenMP] == "")
					||(readAlignmentSAMstrVec[tmpOpenMP].at(0) == '@'))
					continue;
				int threadNO = omp_get_thread_num();
				//cout << "start to process: " << endl << readAlignmentSAMstrVec[tmpOpenMP] << endl;
				alignInferJuncHashInfoVec[threadNO]->
					getAlignInferInfoFromSAM_InsertIntoAlignInferJunctionHash_chrNamePos_supportNum_anchorSize_XM(
						readAlignmentSAMstrVec[tmpOpenMP], indexInfo);
			}
		}
		sam_ifs.close();
	}

	void mergeAlignInferJuncHashInfoInVec2one_chrNamePosOnly(
		AlignInferJunctionHash_Info* mergedAlignInferJunctionHashInfo, Index_Info* indexInfo)
	{
		cout << "start to merge alignInferJuncHashInfoInVec ..." << endl;
		int tmpVecSize = alignInferJuncHashInfoVec.size();
		cout << "alignInferJuncHashInfoVec.size(): " << tmpVecSize << endl;
		for(int tmp = 0; tmp < tmpVecSize; tmp++)
		{
			cout << "tmpAlignInferJuncHashInfo index: " << tmp << endl;
			mergedAlignInferJunctionHashInfo->mergeWithAnotherAlignInferJuncHash_chrNamePosOnly(
				alignInferJuncHashInfoVec[tmp], indexInfo);
		}
	}

	void mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum(
		AlignInferJunctionHash_Info* mergedAlignInferJunctionHashInfo, Index_Info* indexInfo)
	{
		cout << "start to merge alignInferJuncHashInfoInVec ..." << endl;
		int tmpVecSize = alignInferJuncHashInfoVec.size();
		cout << "alignInferJuncHashInfoVec.size(): " << tmpVecSize << endl;
		for(int tmp = 0; tmp < tmpVecSize; tmp++)
		{
			cout << "tmpAlignInferJuncHashInfo index: " << tmp << endl;
			mergedAlignInferJunctionHashInfo->mergeWithAnotherAlignInferJuncHash_chrNamePos_supportNum(
				alignInferJuncHashInfoVec[tmp], indexInfo);
		}
	}

	void mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum_anchorSize(
		AlignInferJunctionHash_Info* mergedAlignInferJunctionHashInfo, Index_Info* indexInfo)
	{
		cout << "start to merge alignInferJuncHashInfoInVec ..." << endl;
		int tmpVecSize = alignInferJuncHashInfoVec.size();
		cout << "alignInferJuncHashInfoVec.size(): " << tmpVecSize << endl;
		for(int tmp = 0; tmp < tmpVecSize; tmp++)
		{
			cout << "tmpAlignInferJuncHashInfo index: " << tmp << endl;
			mergedAlignInferJunctionHashInfo->mergeWithAnotherAlignInferJuncHash_chrNamePos_supportNum_anchorSize(
				alignInferJuncHashInfoVec[tmp], indexInfo);
		}
	}

	void mergeAlignInferJuncHashInfoInVec2one_chrNamePos_supportNum_anchorSize_XM(
		AlignInferJunctionHash_Info* mergedAlignInferJunctionHashInfo, Index_Info* indexInfo)
	{
		cout << "start to merge alignInferJuncHashInfoInVec ..." << endl;
		int tmpVecSize = alignInferJuncHashInfoVec.size();
		cout << "alignInferJuncHashInfoVec.size(): " << tmpVecSize << endl;
		for(int tmp = 0; tmp < tmpVecSize; tmp++)
		{
			cout << "tmpAlignInferJuncHashInfo index: " << tmp << endl;
			mergedAlignInferJunctionHashInfo->mergeWithAnotherAlignInferJuncHash_chrNamePos_supportNum_anchorSize_XM(
				alignInferJuncHashInfoVec[tmp], indexInfo);
		}
	}
};

#endif