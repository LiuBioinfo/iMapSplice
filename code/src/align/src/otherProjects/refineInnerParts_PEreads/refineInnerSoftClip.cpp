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

#include "../../general/otherFunc.h"
#include "../../general/read_block_test.h"
#include "../../general/index_info.h"
#include "../../general/splice_info.h"
#include "../../general/alignInferJunctionHash_info.h"
#include "../../general/fixSingleAnchorNWDP_info.h"
#include "../../general/pairAlignInfo_refineInnerSoftClip.h"

using namespace std;

time_t nowtime;
struct tm *local;

int main(int argc, char** argv)
{
	if(argc != 5)
	{
		cout << "Executable <InputIndexInfo> <inputAlignInferJuncHash> <inputInnerSoftClipSAM> <outputFolder>" << endl;
		exit(1);
	}
	bool do_nwdp_bool = false;
	do_nwdp_bool = true;
	bool do_novelSJsearch_bool = false;
	do_novelSJsearch_bool = true;

	string indexFolderPath = argv[1];
	string inputAlignInferJuncHashPath = argv[2];
	string inputSAMpath = argv[3];
	string outputFolderStr = argv[4];

	outputFolderStr += "/";
	string mkdir = "mkdir -p " + outputFolderStr;
	system(mkdir.c_str());
	string log_file = outputFolderStr + "log";
	ofstream log_ofs(log_file.c_str());	
	string noInnerSoftClip_sam_file = outputFolderStr + "noInnerSoftClip.sam";
	ofstream noInnerSoftClipSam_ofs(noInnerSoftClip_sam_file.c_str());
	string innerSoftClip_unfixed_sam_file = outputFolderStr + "innerSoftClip_unfixed.sam";
	ofstream innerSoftClipUnfixedSam_ofs(innerSoftClip_unfixed_sam_file.c_str());
	string innerSoftClip_fixed_sam_file = outputFolderStr + "innerSoftClip_fixed.sam";
	ofstream innerSoftClipFixedSam_ofs(innerSoftClip_fixed_sam_file.c_str());

	indexFolderPath += "/";
	string indexStr = indexFolderPath;
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string indexParameterFileStr = indexFolderPath + "/_parameter";
	ifstream parameter_ifs(indexParameterFileStr.c_str());
	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "initiate indexInfo" << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "initiate indexInfo" << endl;
	Index_Info* indexInfo = new Index_Info(parameter_ifs);
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	cout << "indexSize: " << indexInfo->returnIndexSize() << endl;
	log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	cout << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	indexInfo->initiate();	
	indexInfo->initiateChrNameIndexArray(1000);
	cout << "finish loading chromosomes" << endl;
	log_ofs << "finish loading chromosomes" << endl;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "...";
	cout << "start to initaite alignInferJunctionHashInfo " << endl;
	log_ofs << endl << "[" << asctime(local) << "...";
	log_ofs << "start to initaite alignInferJunctionHashInfo " << endl;	
	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	int chromNum = indexInfo->returnChromNum();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(chromNum);
	cout << "start to insert alignInferJunction into alignInferJunctionHashInfo" << endl;
	log_ofs << "start to insert alignInferJunction into alignInferJunctionHashInfo" << endl;
	string inputJuncFile = inputAlignInferJuncHashPath;
	alignInferJunctionHashInfo->insertJuncFromJuncFile_chrNamePosOnly(inputJuncFile, indexInfo);
	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(chromNum);	
	alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
	int junctionNum_in_alignInferJuncHashInfo = alignInferJunctionHashInfo->returnAlignInferInfoVecSize();
	cout << "end of inserting alignInferJunction into alignInferJunctionHashInfo, junctionNum = " 
		<< junctionNum_in_alignInferJuncHashInfo << endl;
	log_ofs << "end of inserting alignInferJunction into alignInferJunctionHashInfo, junctionNum = " 
		<< junctionNum_in_alignInferJuncHashInfo << endl;

	ifstream sam_ifs(inputSAMpath.c_str());
	while(!(sam_ifs.eof()))
	{
		string samStr_for, samStr_rcm;
		getline(sam_ifs, samStr_for);
		 if(sam_ifs.eof()||(samStr_for == ""))
		 	break;
		if(samStr_for.at(0) == '@')
			continue;
		getline(sam_ifs, samStr_rcm);
		
		PairAlignInfo_RefineInnerSoftClip* tmpPairAlignInfo 
			= new PairAlignInfo_RefineInnerSoftClip();
		tmpPairAlignInfo->initiateWithSamStr(samStr_for, samStr_rcm, indexInfo);

		bool innerSoftClipExistsBool = tmpPairAlignInfo->innerSoftClip_exists_bool();
		if(!innerSoftClipExistsBool)
		{
			noInnerSoftClipSam_ofs << tmpPairAlignInfo->returnSamStr_original_for() << endl;
			noInnerSoftClipSam_ofs << tmpPairAlignInfo->returnSamStr_original_rcm() << endl;
		}
		else
		{
			//cout << endl << "************************************************" << endl; 
			//cout << "samStr_for: " << endl << samStr_for << endl;
			//cout << "samStr_rcm: " << endl << samStr_rcm << endl;
			string tmpSamStr;
			bool innerSoftClipForTailExistsBool 
				= tmpPairAlignInfo->innerSoftClip_for_tail_exists_bool();
			bool innerSoftClipForTailFixedBool;
			//cout << "innerSoftClipForTailExistsBool: " << innerSoftClipForTailExistsBool << endl;
			if(innerSoftClipForTailExistsBool)
			{
				//cout << "start to do fix_tail_for" << endl;
				tmpPairAlignInfo->fix_tail_for(indexInfo, 
					alignInferJunctionHashInfo, SJ,
					do_nwdp_bool, do_novelSJsearch_bool);
				//cout << "end of fix_tail_for" << endl;
				innerSoftClipForTailFixedBool 
					= tmpPairAlignInfo->innerSoftClip_for_tail_fixed_bool();
				//cout << "innerSoftClipForTailFixedBool " << innerSoftClipForTailFixedBool << endl;
				if(innerSoftClipForTailFixedBool)
				{
					tmpSamStr = tmpSamStr + tmpPairAlignInfo->returnSamStr_fixedTail_for(indexInfo);
				}
				else
				{
					tmpSamStr = tmpSamStr + tmpPairAlignInfo->returnSamStr_original_for();
				}
			}
			else
			{
				tmpSamStr = tmpSamStr + tmpPairAlignInfo->returnSamStr_original_for();
			}

			bool innerSoftClipRcmHeadExistsBool
				= tmpPairAlignInfo->innerSoftClip_rcm_head_exists_bool();
			bool innerSoftClipRcmHeadFixedBool;
			//cout << "innerSoftClipRcmHeadExistsBool: " << innerSoftClipRcmHeadExistsBool << endl;
			if(innerSoftClipRcmHeadExistsBool)
			{
				//cout << "start to do fix_head_rcm " << endl;
				tmpPairAlignInfo->fix_head_rcm(indexInfo, 
					alignInferJunctionHashInfo, SJ,
					do_nwdp_bool, do_novelSJsearch_bool);
				//cout << "end of fix_head_rcm " << endl;
				innerSoftClipRcmHeadFixedBool
					= tmpPairAlignInfo->innerSoftClip_rcm_head_fixed_bool();
				//cout << "innerSoftClipRcmHeadFixedBool: " << innerSoftClipRcmHeadFixedBool << endl;
				if(innerSoftClipRcmHeadFixedBool)
				{
					tmpSamStr = tmpSamStr + "\n" + tmpPairAlignInfo->returnSamStr_fixedHead_rcm(indexInfo);
				}
				else
				{
					//cout << "tmpSamStr: " << endl;
					//cout << tmpSamStr << endl;
					tmpSamStr = tmpSamStr + "\n" + tmpPairAlignInfo->returnSamStr_original_rcm();
				}
				//cout << "end of generating samStr" << endl;
				//cout << "tmpSamStr: " << tmpSamStr << endl;
			}
			else
			{
				tmpSamStr = tmpSamStr + "\n" + tmpPairAlignInfo->returnSamStr_original_rcm();
			}			

			if(((innerSoftClipForTailExistsBool)&&(!innerSoftClipForTailFixedBool))
				||((innerSoftClipRcmHeadExistsBool)&&(!innerSoftClipRcmHeadFixedBool)))
			{
				innerSoftClipUnfixedSam_ofs << tmpSamStr << endl;
			}
			else
			{
				innerSoftClipFixedSam_ofs << tmpSamStr << endl;
			}
		}
		delete tmpPairAlignInfo;
	}
	delete SJ;
	delete alignInferJunctionHashInfo;
	delete indexInfo;
	
	log_ofs.close();
	noInnerSoftClipSam_ofs.close();
	innerSoftClipFixedSam_ofs.close();
	innerSoftClipUnfixedSam_ofs.close();
	return 0;
}