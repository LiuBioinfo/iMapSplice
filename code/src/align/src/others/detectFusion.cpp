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
#include <omp.h>
#include <time.h>

//#include "switch.h"
#include "read_block_test.h"
#include "bwtmap_info.h"
#include "DoubleAnchorScore.h"
#include "sbndm.h"
#include "otherFunc.h"
#include "index_info.h"
#include "constantDefinitions.h"
#include "segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "splice_info.h"
#include "fixGapRelationParameters.h"
#include "read_info.h"
#include "seg_info.h"
#include "gap_info.h"
#include "align_info.h"
#include "spliceJunction_info.h"
#include "unmapEnd_info.h"
#include "unfixedHead.h"
#include "unfixedTail.h"
#include "sam2junc.h"
#include "detectFusion_new.h"

#define PreIndexSize 268435456

using namespace std;  

int main(int argc, char**argv)
{
    if(argc < 4)
	{
		//cout << "Executable <InputReads> <InputReads_PE> <OutputSAM> <threads_num> <Fasta_or_Fastq> <HumanOrMouse>" << endl;
		cout << "Executable <InputUnpairedAlignInfoRecords> <outputDir> <indexFilePrefix> <secondLevelIndex> <inputIncomplePairedFile>" << endl;

		exit(0); 
	}

	bool load2ndLevelIndexBool = true;
	bool detectExactFusionSiteBool = true;
	bool load2ndLevelIndexBool_compressedSize = true;

	string inputRecordsStr = argv[1];

	string indexStr = argv[3];

    string outputDirStr = argv[2];

    //string inputIncompletePairedFileStr = argv[4];

    string secondLevelIndexStr = argv[4];

    string inputIncomplePairedFile = argv[5];

   	string mkdirOutputCommand = "mkdir -p " + outputDirStr;
   	system(mkdirOutputCommand.c_str());

   	string processLogStr = outputDirStr + "/process.log";
   	ofstream log_ofs(processLogStr.c_str());

    log_ofs << "inputRecordsStr: " << inputRecordsStr << endl;
    log_ofs << "outputDir: " << outputDirStr << endl;
    log_ofs << "indexStr: " << indexStr << endl;
   	

   	indexStr += "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs);
	//log_ofs << "index: " << indexStr << endl;
	///////////////////////////////////////
 
	log_ofs << "start to load whole genome" << endl;
	char *chrom; 

	chrom = (char*)malloc((indexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (indexInfo->indexSize) * sizeof(char)); 

	indexInfo->chromString = chrom;
	log_ofs << "chromSize = " <<(indexInfo->chromString).size() << endl;
	
	log_ofs << "start to load every chromosome" << endl;
	//chromStr[0] = 
	(indexInfo->chromStr).push_back((indexInfo->chromString).substr(0, (indexInfo->chrEndPosInGenome)[0]+1));
	(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[0]+1));
	for(int tmp = 1; tmp < indexInfo->chromNum; tmp++)
	{
		//chromStr[tmp] = 
		(indexInfo->chromStr).push_back((indexInfo->chromString).substr((indexInfo->chrEndPosInGenome)[tmp-1]+2, 
			(indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));	
		(indexInfo->chromLength).push_back(((indexInfo->chrEndPosInGenome)[tmp]-(indexInfo->chrEndPosInGenome)[tmp-1]-1));
	}	
	log_ofs << "finish loading all index" << endl;


	FusionDetection_Info* fusionDetectionInfo = new FusionDetection_Info();

	fusionDetectionInfo->generateCandidateFusionMap(inputRecordsStr, indexInfo); // input file: one end unmapped file
	fusionDetectionInfo->generateCandidateFusionMap("/data/homes/lxauky/adSA_result/chrAll/result_0529/1M/output.sam.incomplete.alignInfo", 
		indexInfo);
	//fusionDectionInfo->generateCandidateFusionArrayIndexMap();

	log_ofs << "fusion info: \nsize of fusionCandidate: " << (fusionDetectionInfo->fusionCandidate).size() << endl;
		//<< "\nsize of -+: " << (fusionDetectionInfo->fusionCandidate_ReveForw).size()
		//<< "\nsize of --: " << (fusionDetectionInfo->fusionCandidate_ReveReve).size() << endl;
	
	cout << fusionDetectionInfo->getCandidateFusionMapStr() << endl;

	log_ofs << "start to generate 2 fusionEndSets " << endl;
	fusionDetectionInfo->generateCandidateFusionSetAndCandidateFusionArrayIndexMap();
	log_ofs << "finish generating 2 fusionEndSets " << endl;

	cout << fusionDetectionInfo->getFusionSetStr() << endl;

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;

	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

	vector<unsigned int*> secondLevelLcp;
	vector<unsigned int*> secondLevelUp;
	vector<unsigned int*> secondLevelDown;
	vector<unsigned int*> secondLevelNext;
	
	if(load2ndLevelIndexBool)
	{
		log_ofs << "start to load second-level index ..." << endl;
		
		int secondLevelIndexNO = 0;
		for(int tmpChrNO = 0; tmpChrNO < indexInfo->chromNum; tmpChrNO ++)
		{
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->secondLevelIndexPartsNum)[tmpChrNO]; tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;
				
				string inputIndexFileStr = secondLevelIndexStr + "/" + indexInfo->chrNameStr[tmpChrNO] + "/" 
					//+ indexInfo->chrNameStr[tmpChrNO] + "_part."
					+ tmpFileNumStr + "/";//"." + "test3_";


				string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
				ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
				string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
				ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);

				//if(!load2ndLevelIndexBool_compressedSize)
				//{
				
					string secondLevelIndexFileLcpStr = inputIndexFileStr + "lcp";	
					ifstream secondLevelLcp_file_ifs(secondLevelIndexFileLcpStr.c_str(), ios::binary);	
					string secondLevelIndexFileUpStr = inputIndexFileStr + "up";	
					ifstream secondLevelUp_file_ifs(secondLevelIndexFileUpStr.c_str(), ios::binary);
					string secondLevelIndexFileDownStr = inputIndexFileStr + "down";	
					ifstream secondLevelDown_file_ifs(secondLevelIndexFileDownStr.c_str(), ios::binary);
					string secondLevelIndexFileNextStr = inputIndexFileStr + "next";	
					ifstream secondLevelNext_file_ifs(secondLevelIndexFileNextStr.c_str(), ios::binary);
					
				//}
				//else
				//{
					string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
					ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
					string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
					ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
					string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
					ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					
				//}

				int sizeOfIndex = indexInfo->secondLevelIndexNormalSize + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
				}
				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
				{
					(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
				}

				bool No_ATGC_Bool = true;
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					char ch = tmpSecondLevelChrom[tmpMallocSpace];
					if((ch == 'A')||(ch == 'T')||(ch == 'G')||(ch == 'C'))
					{
						No_ATGC_Bool = false;
						break;
					}
				}				
				if(No_ATGC_Bool)
				{
					(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
				}	

				secondLevelChrom.push_back(tmpSecondLevelChrom);
				
				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
				secondLevelSa.push_back(tmpSecondLevelSa);

				if(!load2ndLevelIndexBool_compressedSize)
				{
					unsigned int* tmpSecondLevelLcp = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelLcp_file_ifs.read((char*)tmpSecondLevelLcp, sizeOfIndex * sizeof(unsigned int));
					secondLevelLcp.push_back(tmpSecondLevelLcp);
					
					unsigned int* tmpSecondLevelUp = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelUp_file_ifs.read((char*)tmpSecondLevelUp, sizeOfIndex * sizeof(unsigned int));
					secondLevelUp.push_back(tmpSecondLevelUp);
					
					unsigned int* tmpSecondLevelDown = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelDown_file_ifs.read((char*)tmpSecondLevelDown, sizeOfIndex * sizeof(unsigned int));
					secondLevelDown.push_back(tmpSecondLevelDown);
					
					unsigned int* tmpSecondLevelNext = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelNext_file_ifs.read((char*)tmpSecondLevelNext, sizeOfIndex * sizeof(unsigned int));
					secondLevelNext.push_back(tmpSecondLevelNext);
				}
				else
				{
					BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
					secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
					secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);
					
					unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
					secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
					secondLevelChildTab.push_back(tmpSecondLevelChildTab);

					BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
					secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
					secondLevelDetChild.push_back(tmpSecondLevelDetChild);
				}


				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();
				//if(!load2ndLevelIndexBool_compressedSize)
				//{
					secondLevelLcp_file_ifs.close();
					secondLevelUp_file_ifs.close();
					secondLevelDown_file_ifs.close();
					secondLevelNext_file_ifs.close();
				//}
				//else
				//{
					secondLevelLcpCompress_file_ifs.close();
					secondLevelChildTab_file_ifs.close();
					secondLevelDetChild_file_ifs.close();					
				//}				

				secondLevelIndexNO ++;

			}
			log_ofs << "finish loading 2nd-level index of " << indexInfo->chrNameStr[tmpChrNO] << endl; 
		}
		log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
		//loadIndex_end = clock(); 
	}

	log_ofs << "start to load incomplete paired alignments ......" << endl;

	if(detectExactFusionSiteBool)
	{
		log_ofs << "start to detect exact fusion site" << endl;
		fusionDetectionInfo->detectExactFusionFromPairedAlignFile(inputIncomplePairedFile, 
			secondLevelChrom,
			secondLevelSa,
			secondLevelLcpCompress,
			secondLevelChildTab,
			secondLevelDetChild,
			indexInfo);
	}


	return 0;
}