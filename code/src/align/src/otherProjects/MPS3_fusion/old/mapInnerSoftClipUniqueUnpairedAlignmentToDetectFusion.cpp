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
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

#include "../../general/extractUnmapAlignment2ReadFile.h"
#include "../../phase1/arrayQueue.h"
#include "../../stats_info.h"
#include "../../constantDefinitions.h"
#include "../../general/option_info.h"
#include "../../general/read_block_test.h"
#include "../../general/bwtmap_info.h"
#include "../../general/DoubleAnchorScore.h"
#include "../../general/sbndm.h"
#include "../../general/otherFunc.h"
#include "../../general/index_info.h"
#include "../../general/enhanced_suffix_array_info.h"
#include "../../general/annotation_info.h"
#include "../../phase1/repeatRegion.h"
#include "../../general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "../../general/splice_info.h"
#include "../../general/fixGapRelationParameters.h"
#include "../../general/read_info.h"
#include "../../general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "../../general/fixDoubleAnchorNWDP_info.h"
#include "../../general/fixDoubleAnchorMatch_info.h"
#include "../../general/fixDoubleAnchorInsertion_info.h"
#include "../../general/fixDoubleAnchorDeletion_info.h"
#include "../../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../../general/fixDoubleAnchorSplice_info.h"
#include "../../general/fixDoubleAnchorCirRNA_info.h"
#include "../../general/path_info.h"
#include "../../general/gap_info.h"
#include "../../general/align_info.h"
#include "../../general/peAlign_info.h"
#include "../../general/groupSeg_info.h"
#include "../../general/alignInferJunctionHash_info_vec.h"
#include "../../phase2/spliceJunctionHash_info.h"
#include "../../phase2/unmapEnd_info.h"
#include "../../phase2/unfixedHead.h"
#include "../../phase2/unfixedTail.h"
#include "../../phase2/incompleteLongHead.h"
#include "../../phase2/incompleteLongTail.h"
#include "../../phase2/sam2junc.h"
#include "../../fixHeadTail.h"
#include "../../phase2/fixOneEndUnmapped.h"
#include "../../fixPhase1.h"
#include "../../general/readSeqPreProcessing.h"
#include "../../general/headerSection_info.h"
#include "../../general/otherFunc2.h"
#include "../../general/alignmentToJunc.h"
#include "../../../general/localMapUnfixedReadEnd2DetectFusion_info.h"

using namespace std;

bool forOrRcm(int flag)
{
	if(flag & 0x10)
		return false;
	else
		return true;
}

void extractFieldVal_IH_XI(int& tmpIHint, int& tmpXIint,
	string& tmpUnpairedPEalignmentStr)
{
	int IHfieldLoc = tmpUnpairedPEalignmentStr.find("IH:i:");
	int XIfieldLoc = tmpUnpairedPEalignmentStr.find("XI:i:");
	int IHfieldLoc_nextTab = tmpUnpairedPEalignmentStr.find("\t", IHfieldLoc);
	int lineLength = tmpUnpairedPEalignmentStr.length();
	string IHintStr = tmpUnpairedPEalignmentStr.substr(IHfieldLoc + 5, IHfieldLoc_nextTab - IHfieldLoc - 5);
	string XIintStr = tmpUnpairedPEalignmentStr.substr(XIfieldLoc + 5, lineLength - XIfieldLoc - 5);
	tmpIHint = atoi(IHintStr.c_str());
	tmpXIint = atoi(XIintStr.c_str());
}

int getEndPosOfSpecificJumpCode(int startPos, vector<Jump_Code>& cigarStringJumpCodeVec, int jumpCodeIndex)
{
	int tmpEndPos = 0;
	for(int tmpIndex = 0; tmpIndex <= jumpCodeIndex; tmpIndex++)
	{
		string tmpJumpCodeType = cigarStringJumpCodeVec[tmpIndex].type;
		int tmpJumpCodeLength = cigarStringJumpCodeVec[tmpIndex].len;
		if(tmpJumpCodeType == "S")
		{}
		else if(tmpJumpCodeType == "M")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "I")
		{}
		else if(tmpJumpCodeType == "D")
			tmpEndPos += tmpJumpCodeLength;
		else if(tmpJumpCodeType == "N")
			tmpEndPos += tmpJumpCodeLength;
		else
		{
			cout << "incorrect jumpCode type" << endl;
			exit(1);
		}								
	}
	return (tmpEndPos + startPos-1);
}

void cigarString2jumpCodeVec(string& jumpCodeStr, vector<Jump_Code>& cigarStringJumpCodeVec)
{
	int tmpJumpCodeLength;
	string tmpJumpCodeType;
	int jumpCodeStartPosInCigarStr = 0;
	int jumpCodeEndPosInCigarStr;
	
	string candidateJumpCodeType = "SMNIDX";
	while(1)
	{
		jumpCodeEndPosInCigarStr = 
			jumpCodeStr.find_first_of(candidateJumpCodeType, jumpCodeStartPosInCigarStr);
		if(jumpCodeEndPosInCigarStr == jumpCodeStr.npos)
			{break;}
		else
		{
			tmpJumpCodeLength = 
				atoi((jumpCodeStr.substr(jumpCodeStartPosInCigarStr, jumpCodeEndPosInCigarStr - jumpCodeStartPosInCigarStr)).c_str());
			tmpJumpCodeType = jumpCodeStr.substr(jumpCodeEndPosInCigarStr, 1);
			cigarStringJumpCodeVec.push_back(Jump_Code(tmpJumpCodeLength, tmpJumpCodeType));
			jumpCodeStartPosInCigarStr = jumpCodeEndPosInCigarStr + 1;
		}
	}
}

void getSAMinfoFromSamStr(
	int& tmpChrNameInt, int& tmpStartMapPos, int& tmpEndMapPos,
	int& unfixedHeadLength, int& unfixedTailLength, bool& Nor_or_Rcm_bool,
	string& tmpReadSeq, string& tmpQualSeq,
	string& tmpSamStr, Index_Info* indexInfo)
{
	vector<string> tmpSamFieldVec;
	int startLoc = 0;
	for(int tmp = 0; tmp < 11; tmp++)
	{
		int tabLoc = tmpSamStr.find("\t", startLoc);
		tmpSamFieldVec.push_back(tmpSamStr.substr(startLoc, tabLoc - startLoc));
		startLoc = tabLoc + 1;
	}

	string tmpChrName = tmpSamFieldVec[2];
	tmpChrNameInt = indexInfo->convertStringToInt(tmpChrName);

	string tmpStartMapPosStr = tmpSamFieldVec[3];
	tmpStartMapPos = atoi(tmpStartMapPosStr.c_str());
	
	string cigarString = tmpSamFieldVec[5];
	vector<Jump_Code> tmpJumpCodeVec;
	cigarString2jumpCodeVec(cigarString, tmpJumpCodeVec);
	tmpEndMapPos = getEndPosOfSpecificJumpCode(tmpStartMapPos, tmpJumpCodeVec, tmpJumpCodeVec.size()-1);

	if(tmpJumpCodeVec[0].type == "S")
		unfixedHeadLength = tmpJumpCodeVec[0].len;
	else
		unfixedHeadLength = 0;	

	if(tmpJumpCodeVec[tmpJumpCodeVec.size()-1].type == "S")
		unfixedTailLength = tmpJumpCodeVec[tmpJumpCodeVec.size()-1].len;
	else
		unfixedTailLength = 0;

	int flagInt = atoi(tmpSamFieldVec[1].c_str());
	Nor_or_Rcm_bool = forOrRcm(flagInt);	

	tmpReadSeq = tmpSamFieldVec[9];
	tmpQualSeq = tmpSamFieldVec[10];
}


int main(int argc char** argv)
{
	if(argc != 7)
	{
		cout << "Executable GlobalIndex LocalIndex inputUnpairedPEalignmentPath ";
		cout << "outputFolder threads_num fasta_fastq" << endl;
		exit(1);	
	}

	string globalIndexStr = argv[1];
	string localIndexStr = argv[2];
	string unpairedPEalignmentPath = argv[3];
	ifstream unpairedPEalignment_ifs(unpairedPEalignmentPath.c_str());

	string outputFolderStr = argv[4];
	string threadsNumStr = argv[5];
	int threads_num = atoi(threadsNumStr.c_str());
	bool fasta_or_fastq_bool = true;
	string fasta_or_fastq_str = argv[5];
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta"))
		fasta_or_fastq_bool = true;
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq"))
		fasta_or_fastq_bool = false;
	else
	{
		cout << "Please set the correct format, fasta or fastq" << endl;
		exit(1);
	}
	cout << "creat folders and files ......" << endl;

	Annotation_Info* annotationInfo = new Annotation_Info();

	string mkdirOutputCommand = "mkdir -p " + outputFolderStr;
	system(mkdirOutputCommand.c_str());

	string log_ofs_str = outputFolderStr + "/process.log";
	ofstream log_ofs(log_ofs_str.c_str());
	string runtime_ofs_str = outputFolderStr + "/runtime.log";
	ofstream runtime_log_ofs(runtime_ofs_str.c_str());
	string settings_ofs_str = outputFolderStr + "/settings.log";
	ofstream settings_log_ofs(settings_ofs_str.c_str());

	string keptSAMpath = outputFolderStr + "/kept.sam";
	ofstream keptSAM_ofs(keptSAMpath.c_str());
	string fixedFusionSAMpath = outputFolderStr + "/fusion.sam";
	ofstream fusionSAM_ofs(fixedFusionSAMpath.c_str()); 

	string indexStr = globalIndexStr + "/";
	string chrom_bit_file = indexStr; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);
	string parameter_file = indexStr; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	Index_Info* indexInfo = new Index_Info(parameter_file_ifs, settings_log_ofs);

	settings_log_ofs << "index: " << indexStr << endl;
	/////////////////////////////////////// 
	log_ofs << "start to load whole genome" << endl;
	char *chrom; chrom = (char*)malloc((indexInfo->returnIndexSize()) * sizeof(char)); chrom_bit_file_ifs.read((char*)chrom, (indexInfo->returnIndexSize()) * sizeof(char)); 
	settings_log_ofs << "indexSize: " << indexInfo->returnIndexSize() << endl;
	indexInfo->readGenome(chrom);
	settings_log_ofs << "chromSize: " << indexInfo->returnChromStringLength() << endl;
	log_ofs << "start to load every chromosome" << endl;
	indexInfo->initiate();	
	log_ofs << "start to initiate chrNameIndexArray" << endl;
	indexInfo->initiateChrNameIndexArray(1000);
	log_ofs << "finish loading chromosomes" << endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////    	Load Second Level Index      ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "start to initaite alignInferJunctionHashInfo " << endl;

	AlignInferJunctionHash_Info* alignInferJunctionHashInfo = new AlignInferJunctionHash_Info();
	alignInferJunctionHashInfo->initiateAlignInferJunctionInfo(indexInfo->returnChromNum());
	cout << "start to insert alignInferJunction into alignInferJunctionHashInfo" << endl;
	string inputJuncFile = inputAlignInferJuncHashFile;
	alignInferJunctionHashInfo->insertJuncFromJuncFile_onlyChrNamePos(inputJuncFile, indexInfo);


	time_t nowtime;
	struct tm *local;

	nowtime = time(NULL);
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... start to load 2nd level index ......" << endl << endl; 
	log_ofs << endl << "[" << asctime(local) << "... load 2nd level index starts ......" << endl << endl; 

	vector<char*> secondLevelChrom;
	vector<unsigned int*> secondLevelSa;

	vector<BYTE*> secondLevelLcpCompress;
	vector<unsigned int*> secondLevelChildTab;
	vector<BYTE*> secondLevelDetChild;

	string secondLevelIndexStr = localIndexStr + "/";

	if(load2ndLevelIndexBool)
	{
		log_ofs << "start to load second-level index ..." << endl;
		
		int secondLevelIndexNO = 0;
		for(int tmpChrNO = 0; tmpChrNO < indexInfo->returnChromNum(); tmpChrNO ++)
		{
			for(int tmpSecondLevelIndexNO = 1; tmpSecondLevelIndexNO <= (indexInfo->returnSecondLevelIndexPartsNum(tmpChrNO)); tmpSecondLevelIndexNO ++)
			{
				char tmpFileNumChar[4];
				sprintf(tmpFileNumChar, "%d", tmpSecondLevelIndexNO);
				string tmpFileNumStr = tmpFileNumChar;
				
				string inputIndexFileStr = secondLevelIndexStr + "/" + 
					//indexInfo->chrNameStr[tmpChrNO] 
					indexInfo->returnChrNameStr(tmpChrNO)
					+ "/" 
					//+ indexInfo->chrNameStr[tmpChrNO] + "_part."
					+ tmpFileNumStr + "/";//"." + "test3_";


				string secondLevelIndexFileChromStr = inputIndexFileStr + "chrom"; 
				ifstream secondLevelChrom_file_ifs(secondLevelIndexFileChromStr.c_str(), ios::binary);
				string secondLevelIndexFileSaStr = inputIndexFileStr + "SA";
				ifstream secondLevelSA_file_ifs(secondLevelIndexFileSaStr.c_str(), ios::binary);

					string secondLevelIndexFileLcpCompressStr = inputIndexFileStr + "_lcpCompress";	
					ifstream secondLevelLcpCompress_file_ifs(secondLevelIndexFileLcpCompressStr.c_str(), ios::binary);	
					string secondLevelIndexFileChildTabStr = inputIndexFileStr + "childTab";	
					ifstream secondLevelChildTab_file_ifs(secondLevelIndexFileChildTabStr.c_str(), ios::binary);
					string secondLevelIndexFileDetChildStr = inputIndexFileStr + "detChild";	
					ifstream secondLevelDetChild_file_ifs(secondLevelIndexFileDetChildStr.c_str(), ios::binary);					

				int sizeOfIndex = indexInfo->returnSecondLevelIndexNormalSize() + 1;
				char* tmpSecondLevelChrom = (char*)malloc(sizeOfIndex * sizeof(char));
				for(int tmpMallocSpace = 0; tmpMallocSpace < sizeOfIndex; tmpMallocSpace++)
				{
					tmpSecondLevelChrom[tmpMallocSpace] = '0';
				}
				secondLevelChrom_file_ifs.read((char*)tmpSecondLevelChrom, sizeOfIndex * sizeof(char));
				if(tmpSecondLevelChrom[sizeOfIndex-1] != 'X')
				{
					//(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
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
					//(indexInfo->invalidSecondLevelIndexNOset).insert(secondLevelIndexNO + 1);
					indexInfo->insert2invalidSecondLevelIndexNOset(secondLevelIndexNO + 1);
				}	

				secondLevelChrom.push_back(tmpSecondLevelChrom);
				
				unsigned int* tmpSecondLevelSa = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelSA_file_ifs.read((char*)tmpSecondLevelSa, sizeOfIndex * sizeof(unsigned int));
				secondLevelSa.push_back(tmpSecondLevelSa);

				BYTE* tmpSecondLevelLcpCompress = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress_file_ifs.read((char*)tmpSecondLevelLcpCompress, sizeOfIndex * sizeof(BYTE));
				secondLevelLcpCompress.push_back(tmpSecondLevelLcpCompress);
					
				unsigned int* tmpSecondLevelChildTab = (unsigned int*)malloc(sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab_file_ifs.read((char*)tmpSecondLevelChildTab, sizeOfIndex * sizeof(unsigned int));
				secondLevelChildTab.push_back(tmpSecondLevelChildTab);

				BYTE* tmpSecondLevelDetChild = (BYTE*)malloc(sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild_file_ifs.read((char*)tmpSecondLevelDetChild, sizeOfIndex * sizeof(BYTE));
				secondLevelDetChild.push_back(tmpSecondLevelDetChild);

				secondLevelChrom_file_ifs.close();
				secondLevelSA_file_ifs.close();

				secondLevelLcpCompress_file_ifs.close();
				secondLevelChildTab_file_ifs.close();
				secondLevelDetChild_file_ifs.close();							

				secondLevelIndexNO ++;
			}
			log_ofs << "finish loading 2nd-level index of " << indexInfo->returnChrNameStr(tmpChrNO) << endl; 
		}
		log_ofs << "finish loading ALL 2nd-level index !" << endl;
		log_ofs << indexInfo->getInvalidSecondLevelIndexNOstr() << endl;
	}

	while(!unpairedPEalignment_ifs.eof())
	{
		string tmpUnpairedPEalignmentStr;
		getline(unpairedPEalignment_ifs, tmpUnpairedPEalignmentStr);
		if((unpairedPEalignment_ifs.eof())||(tmpUnpairedPEalignmentStr == ""))
			break;
		int tmpIHint, tmpXIint;
		vector<string> tmpOriSamStrVec;
		vector<string> tmpEndSamStrVec_1;
		vector<string> tmpEndSamStrVec_2;
		extractFieldVal_IH_XI(tmpIHint, tmpXIint, tmpUnpairedPEalignmentStr);
		tmpOriSamStrVec.push_back(tmpUnpairedPEalignmentStr);
		if(tmpIHint == 0)
		{}
		else
		{
			tmpEndSamStrVec_1.push_back(tmpUnpairedPEalignmentStr);
			for(int tmp = 1; tmp < tmpIHint; tmp++)
			{
				string tmpStr;
				getline(unpairedPEalignment_ifs, tmpStr);
				tmpEndSamStrVec_1.push_back(tmpStr);
				tmpOriSamStrVec.push_back(tmpStr);
			}
		}
		if(tmpXIint == 0)
		{
			string tmpStr;
			getline(unpairedPEalignment_ifs, tmpStr);
			tmpOriSamStrVec.push_back(tmpStr);
		}
		else
		{
			for(int tmp = 0; tmp < tmpXIint; tmp++)
			{
				string tmpStr;
				getline(unpairedPEalignment_ifs, tmpStr);
				tmpEndSamStrVec_2.push_back(tmpStr);
				tmpOriSamStrVec.push_back(tmpStr);
			}
		}

		if((tmpEndSamStrVec_1.size() == 1)&&(tmpEndSamStrVec_2.size() == 1))
		{
			int tmpChrNameInt_1, tmpChrNameInt_2, tmpChrStartMapPos_1, 
				tmpChrEndMapPos_1, tmpChrStartMapPos_2, tmpChrEndMapPos_2, 
				unfixedHeadLength_1, unfixedTailLength_1,
				unfixedHeadLength_2, unfixedTailLength_2;
			bool forOrRevBool_1, forOrRevBool_2;
			string tmpReadSeq_1, tmpReadSeq_2, tmpQualSeq_1, tmpQualSeq_2;

			getSAMinfoFromSamStr(tmpChrNameInt_1, tmpChrStartMapPos_1, 
				tmpChrEndMapPos_1, unfixedHeadLength_1, unfixedTailLength_1, 
				forOrRevBool_1, tmpReadSeq_1, tmpQualSeq_1, tmpEndSamStrVec_1[0], indexInfo);

			getSAMinfoFromSamStr(tmpChrNameInt_2, tmpChrStartMapPos_2, 
				tmpChrEndMapPos_2, unfixedHeadLength_2, unfixedTailLength_2, 
				forOrRevBool_2, tmpReadSeq_2, tmpQualSeq_2, tmpEndSamStrVec_2[0], indexInfo);

			if(forOrRevBool_1 && forOrRevBool_2)
			{
				if((unfixedHeadLength_1 > 0)||(unfixedHeadLength_2 > 0))
				{
					bool firstCaseFixedBool = false;
					bool secondCaseFixedBool = false;
					string fixedFirstCaseSAMstr, fixedSecondCaseSAMstr;
					if(unfixedHeadLength_1 > 0) // try Nor_2, Rcm_1/2 + Nor_1/2
					{
						string unfixedReadSubSeq // rcm of unfixedHead_1 
							= covertStringToReverseComplement(
								tmpReadSeq_1.substr(0,unfixedHeadLength_1));
						int leftMostMapPos // startPos of Nor_2
							= tmpChrStartMapPos_2;
						int leftMostMapPos_chrNameInt
							= tmpChrNameInt_2;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						firstCaseFixedBool = tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithLeftMostMapPos(
							fixedFirstCaseSAMstr, unfixedReadSubSeq, leftMostMapPos, 
							leftMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo;
					}
					if(unfixedHeadLength_2 > 0)
					{
						string unfixedReadSubSeq // rcm of unfixedHead_1 
							= covertStringToReverseComplement(
								tmpReadSeq_2.substr(0,unfixedHeadLength_2));
						int leftMostMapPos // startPos of Nor_2
							= tmpChrStartMapPos_1;
						int leftMostMapPos_chrNameInt
							= tmpChrNameInt_1;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						secondCaseFixedBool = tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithLeftMostMapPos(
							fixedSecondCaseSAMstr, unfixedReadSubSeq, leftMostMapPos, 
							leftMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo;
					}
					if((!firstCaseFixedBool)&&(!secondCaseFixedBool))
					{
						for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
							keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;						
					}
				}
				else
				{
					for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
						keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
				}
			}
			else if(forOrRevBool_1 && (!forOrRevBool_2))
			{
				if((unfixedTailLength_1 > 0)||(unfixedHeadLength_2 > 0))
				{
					bool firstCaseFixedBool = false;
					bool secondCaseFixedBool = false;
					string fixedFirstCaseSAMstr, fixedSecondCaseSAMstr;
					if(unfixedTailLength_1 > 0)
					{
						int readLength_1 = tmpReadSeq_1.length();
						string unfixedReadSubSeq = tmpReadSeq_1.substr(
							readLength_1 - unfixedTailLength_1, unfixedTailLength_1);
						int rightMostMapPos = tmpChrEndMapPos_2;
						int rightMostMapPos_chrNameInt = tmpChrNameInt_2;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						firstCaseFixedBool = tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithRightMostMapPos(
							fixedFirstCaseSAMstr, unfixedReadSubSeq, rightMostMapPos,
							rightMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo;
					}
					if(unfixedHeadLength_2 > 0)
					{
						string unfixedReadSubSeq = tmpReadSeq_2.substr(0,unfixedHeadLength_2);
						int leftMostMapPos = tmpChrStartMapPos_1;
						int leftMostMapPos_chrNameInt = tmpChrNameInt_1;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						secondCaseFixedBool = tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithLeftMostMapPos(
							fixedSecondCaseSAMstr, unfixedReadSubSeq, leftMostMapPos,
							leftMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo;
					}
					if((!firstCaseFixedBool)&&(!secondCaseFixedBool))
					{
						for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
							keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
					}					
				}
				else
				{
					for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
						keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
				}
			}
			else if((!forOrRevBool_1) && forOrRevBool_2) // Nor_2 rev_1
			{
				if((unfixedHeadLength_1 > 0)||(unfixedTailLength_2 > 0))
				{
					bool firstCaseFixedBool = false;
					bool secondCaseFixedBool = false;
					string fixedFirstCaseSAMstr, fixedSecondCaseSAMstr;
					if(unfixedHeadLength_1 > 0)
					{
						string unfixedReadSubSeq = tmpReadSeq_1.substr(0, unfixedHeadLength_1);
						int leftMostMapPos = tmpChrStartMapPos_2;
						int leftMostMapPos_chrNameInt = tmpChrNameInt_2;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						firstCaseFixedBool = tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithLeftMostMapPos(
							fixedFirstCaseSAMstr, unfixedReadSubSeq, leftMostMapPos,
							leftMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo
					}
					if(unfixedTailLength_2 > 0)
					{
						int readLength_2 = tmpReadSeq_2.length();
						string unfixedReadSubSeq = tmpReadSeq_2.substr(
							readLength_2 - unfixedTailLength_2, unfixedTailLength_2);
						int rightMostMapPos = tmpChrEndMapPos_1;
						int rightMostMapPos_chrNameInt = tmpChrNameInt_1;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						secondCaseFixedBool = tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithRightMostMapPos(
							fixedSecondCaseSAMstr, unfixedReadSubSeq, rightMostMapPos,
							rightMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo;
					}
					if((!firstCaseFixedBool)&&(!secondCaseFixedBool))
					{
						for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
							keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
					}						
				}
				else
				{
					for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
						keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
				}
			}
			else // rev_1 && rev_2
			{
				if((unfixedTailLength_1 > 0)||(unfixedTailLength_2 > 0))
				{
					bool firstCaseFixedBool = false;
					bool secondCaseFixedBool = false;
					string fixedFirstCaseSAMstr, fixedSecondCaseSAMstr;
					if(unfixedTailLength_1 > 0)
					{
						int readLength_1 = tmpReadSeq_1.length();
						string unfixedReadSubSeq = covertStringToReverseComplement(
							tmpReadSeq_1.substr(readLength_1 - unfixedTailLength_1, unfixedTailLength_1));
						int rightMostMapPos = tmpChrEndMapPos_2;
						int rightMostMapPos_chrNameInt = tmpChrNameInt_2;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						firstCaseFixedBool 
							= tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithRightMostMapPos(
								fixedFirstCaseSAMstr, unfixedReadSubSeq, rightMostMapPos,
								rightMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo
					}
					if(unfixedTailLength_2 > 0)
					{
						int readLength_2 = tmpReadSeq_2.length();
						string unfixedReadSubSeq = covertStringToReverseComplement(
							tmpReadSeq_2.substr(readLength_2 - unfixedTailLength_2, unfixedTailLength_2));
						int rightMostMapPos = tmpChrEndMapPos_1;
						int rightMostMapPos_chrNameInt = tmpChrNameInt_2;
						LocalMapUnfixedReadEnd2DetectFusion_Info* tmpLocalMapUnfixedReadEnd2DetectFusionInfo
							= new LocalMapUnfixedReadEnd2DetectFusion_Info(); 
						secondCaseFixedBool 
							= tmpLocalMapUnfixedReadEnd2DetectFusionInfo->mapSeq2genomeWithRightMostMapPos(
								fixedSecondCaseSAMstr, unfixedReadSubSeq, rightMostMapPos,
								rightMostMapPos_chrNameInt);
						delete tmpLocalMapUnfixedReadEnd2DetectFusionInfo;
					}
					if((!firstCaseFixedBool)&&(!secondCaseFixedBool))
					{
						for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
							keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
					}
				}
				else
				{
					for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
						keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
				}
			}
		}
		else
		{
			for(int tmp = 0; tmp < tmpOriSamStrVec.size(); tmp++)
			{
				keptSAM_ofs << tmpOriSamStrVec[tmp] << endl;
			}
		}
	}
	unpairedPEalignment_ifs.close();

	return 0;
}