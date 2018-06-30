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


#include "../general/extractUnmapAlignment2ReadFile.h"
#include "../phase1/arrayQueue.h"
#include "../stats_info.h"
#include "../constantDefinitions.h"
#include "../general/option_info.h"
#include "../general/read_block_test.h"
#include "../general/bwtmap_info.h"
#include "../general/DoubleAnchorScore.h"
#include "../general/sbndm.h"
#include "../general/otherFunc.h"
#include "../general/index_info.h"
#include "../general/enhanced_suffix_array_info.h"
#include "../general/annotation_info.h"
#include "../phase1/repeatRegion.h"
#include "../general/segmentMapping.h"
//#include "segmentMapping_secondLevel.h"
#include "../general/splice_info.h"
#include "../general/fixGapRelationParameters.h"
#include "../general/read_info.h"
#include "../general/seg_info.h"
//#include "general/fixDoubleAnchor_annotation_info.h"
#include "../general/fixDoubleAnchorNWDP_info.h"
#include "../general/fixDoubleAnchorMatch_info.h"
#include "../general/fixDoubleAnchorInsertion_info.h"
#include "../general/fixDoubleAnchorDeletion_info.h"
#include "../general/fixDoubleAnchorSplice_complicate_info.h"
#include "../general/fixDoubleAnchorSplice_info.h"
#include "../general/fixDoubleAnchorCirRNA_info.h"
#include "../general/path_info.h"
#include "../general/gap_info.h"
#include "../general/align_info.h"
#include "../general/peAlign_info.h"
#include "../general/groupSeg_info.h"
#include "../general/alignInferJunctionHash_info_vec.h"
#include "../phase2/spliceJunctionHash_info.h"
#include "../phase2/unmapEnd_info.h"
#include "../phase2/unfixedHead.h"
#include "../phase2/unfixedTail.h"
#include "../phase2/incompleteLongHead.h"
#include "../phase2/incompleteLongTail.h"
#include "../phase2/sam2junc.h"
#include "../fixHeadTail.h"
#include "../phase2/fixOneEndUnmapped.h"
#include "../fixPhase1.h"
#include "../general/readSeqPreProcessing.h"
#include "../general/headerSection_info.h"
#include "../general/otherFunc2.h"
#include "../general/alignmentToJunc.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cout << "Exetuable GlobalIndex LocalIndex inputAlignmentInfo_incomplete";
		cout << " inputAlignInferJuncHashFile outputFolder thread_num fasta_or_fastq SE_or_PE_bool" << endl;
		exit(1);
	}	
	string globalIndexStr = argv[1];
	string localIndexStr = argv[2];
	string inputAlignInfoStr_incomplete = argv[3];
	//string inputSJstr = argv[4];
	string inputAlignInferJuncHashFile = argv[4];
	string outputFolderStr = argv[5];
	string threadNumStr = argv[6];
	string fasta_or_fastq_str = argv[7];

	cout << "to initiate ......" << endl;

	int normalRecordNum_fixHeadTail = 1;//500000;
	bool load2ndLevelIndexBool = true;
	
	bool SE_or_PE_bool = false;
	string SE_or_PE_str = argv[8];
	if(SE_or_PE_str == "SE")
		SE_or_PE_bool = true;
	else if(SE_or_PE_str == "PE")
		SE_or_PE_bool = false;
	else
	{
		cout << "incorrect SE / PE setttings" << endl;
		exit(1);
	}

	//cout << "to initiate 2 ......" << endl;

	bool annotation_provided_bool = false;
	string annotation_file_path = "";
	bool checkQualSeqForReadSegSeq = false;
	bool Do_annotation_only_bool = false;
	bool DoSam2JuncBool = true;
	bool DoRemappingOnUnfixedHeadTailAlignmentBool = true;
	//cout << "to initiate 3 ......" << endl;
	bool Do_extendHeadTail_fixHeadTail = true;
	//Do_extendHeadTail_fixHeadTail = false;
	bool Do_fixHeadTail_remapping = true;
	//Do_fixHeadTail_remapping = false;
	bool Do_fixHeadTail_greedyMapping = true;
	//Do_fixHeadTail_greedyMapping = false;
	bool Do_fixHeadTail_remappingAndTargetMapping = true;
	//Do_fixHeadTail_remappingAndTargetMapping = false;
	bool Do_fixHeadTail_remappingAgain = true;
	Do_fixHeadTail_remappingAgain = false;

	//cout << "threadNumStr: " << threadNumStr << endl;

	int threads_num = atoi(threadNumStr.c_str());
	bool fasta_or_fastq_bool;
	if((fasta_or_fastq_str == "fasta")||(fasta_or_fastq_str == "Fasta")||(fasta_or_fastq_str == "FASTA"))
	{
		fasta_or_fastq_bool = true;
	}
	else if((fasta_or_fastq_str == "fastq")||(fasta_or_fastq_str == "Fastq")||(fasta_or_fastq_str == "FASTQ"))
	{
		fasta_or_fastq_bool = false;
	}
	else
	{
		cout << "input file format invalid" << endl;
		exit(1);
	}

	cout << "creat folders and files ......" << endl;

	Annotation_Info* annotationInfo = new Annotation_Info();

	string mkdirOutputCommand = "mkdir -p " + outputFolderStr;
	system(mkdirOutputCommand.c_str());

	string outputDirStr = outputFolderStr + "/";
	string outputDirStr_phase2 = outputDirStr + "phase2_output";
	string mkdirOutputCommand_phase2 = "mkdir -p " + outputDirStr_phase2;
	system(mkdirOutputCommand_phase2.c_str());

	string log_ofs_str = outputDirStr + "/process.log";
	ofstream log_ofs(log_ofs_str.c_str());
	string runtime_ofs_str = outputDirStr + "/runtime.log";
	ofstream runtime_log_ofs(runtime_ofs_str.c_str());
	string settings_ofs_str = outputDirStr + "/settings.log";
	ofstream settings_log_ofs(settings_ofs_str.c_str());

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
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////   Do REMAPPING On unfixed head/tail Reads    ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	nowtime = time(NULL);
	//struct tm *local;
	local = localtime(&nowtime);
	cout << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 
	runtime_log_ofs << endl << "[" << asctime(local) << "... fixing unfixed-head/tail reads starts ......" << endl << endl ; 

	string OutputSamFile_fixHeadTail_complete_pair = outputDirStr + "/phase2_output/fixHeadTail_complete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_pair_ofs(OutputSamFile_fixHeadTail_complete_pair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_pair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_pair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_pair_ofs(OutputSamFile_fixHeadTail_incomplete_pair.c_str());

	string OutputSamFile_fixHeadTail_complete_unpair = outputDirStr + "/phase2_output/fixHeadTail_complete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_complete_unpair_ofs(OutputSamFile_fixHeadTail_complete_unpair.c_str());

	string OutputSamFile_fixHeadTail_incomplete_unpair = outputDirStr + "/phase2_output/fixHeadTail_incomplete_unpair.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_unpair_ofs(OutputSamFile_fixHeadTail_incomplete_unpair.c_str());	

	string OutputSamFile_fixHeadTail_pair_lowScore = outputDirStr + "/phase2_output/fixHeadTail_pair_lowScore.sam";
	ofstream OutputSamFile_fixHeadTail_pair_lowScore_ofs(OutputSamFile_fixHeadTail_pair_lowScore.c_str());	

	string OutputSamFile_fixHeadTail_complete_SE = outputDirStr + "/phase2_output/fixHeadTail_complete_SE.sam";
	ofstream OutputSamFile_fixHeadTail_complete_SE_ofs(OutputSamFile_fixHeadTail_complete_SE.c_str());

	string OutputSamFile_fixHeadTail_incomplete_SE = outputDirStr + "/phase2_output/fixHeadTail_incomplete_SE.sam";
	ofstream OutputSamFile_fixHeadTail_incomplete_SE_ofs(OutputSamFile_fixHeadTail_incomplete_SE.c_str());

	string OutputSamFile_fixHeadTail_lowScore_SE = outputDirStr + "/phase2_output/fixHeadTail_lowScore_SE.sam";
	ofstream OutputSamFile_fixHeadTail_lowScore_SE_ofs(OutputSamFile_fixHeadTail_lowScore_SE.c_str());

	// start to read splice junction
	int junctionNum = 0;
	int junctionNum_in_alignInferJuncHashInfo = 0;
	int junctionNum_in_annotation = 0;

	SJhash_Info* SJ = new SJhash_Info();
	SJ->initiateAreaAndStringHash(indexInfo->returnChromNum());	
	//string tmpIntermediateJunctionFile = inputSJstr;
	//string InputSpliceJunction = tmpIntermediateJunctionFile;

	////////  take annotation as reference for remapping //////////////
	//InputSpliceJunction = "/data/homes/lxauky/GroundTruth2sam/resutls/sim1_test1/simulated_reads_test1_twoEnds.sam.noRandom.junc";
	// settings_log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	// log_ofs << "InputSpliceJunction: " << InputSpliceJunction << endl; 
	// cout << "InputSpliceJunction: " << InputSpliceJunction << endl;
	///////////////////////////////////////////////////////////////////

	//if(annotation_provided_bool && Do_annotation_only_bool)
	//	InputSpliceJunction = annotation_file_path;	
	
	if(DoRemappingOnUnfixedHeadTailAlignmentBool)
	{
    	log_ofs << "start to build spliceJunction Hash" << endl;
    	cout << "start to build spliceJunction Hash" << endl;
    	bool spliceJunctionHashExists = true;

		string entryString;
		int tabLocation1, tabLocation2, tabLocation3, tabLocation4, tabLocation5;
		char entry[500];
		int chrInt;
		int spliceStartPos;
		int spliceEndPos;
		string chrIntString;
		string spliceStartPosString;
		string spliceEndPosString;

		/////////////////////////////////////////////////////////////////////////////
		//////////////////////  string hash /////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////

		if(!Do_annotation_only_bool)
		{	

			// loading SJs generated from aligned reads.
			/*FILE *fp_spliceJunction = fopen(InputSpliceJunction.c_str(), "r");
			fgets(entry, sizeof(entry), fp_spliceJunction);
			while(!feof(fp_spliceJunction))
			{
				fgets(entry, sizeof(entry), fp_spliceJunction);
				if(feof(fp_spliceJunction))
					break;
				junctionNum ++;
				entryString = entry;
				tabLocation1 = entryString.find('\t', 0);
				tabLocation2 = entryString.find('\t', tabLocation1+1);
				tabLocation3 = entryString.find('\t', tabLocation2+1);
				chrIntString = entryString.substr(0, tabLocation1);
				spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
				spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
				chrInt = indexInfo->convertStringToInt(chrIntString);
				spliceStartPos = atoi(spliceStartPosString.c_str());
				spliceEndPos = atoi(spliceEndPosString.c_str());	
				SJ->insert2AreaAndStringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
			}
			fclose(fp_spliceJunction);*/
			alignInferJunctionHashInfo->convert2SJhashInfo(SJ, indexInfo);
			junctionNum_in_alignInferJuncHashInfo = alignInferJunctionHashInfo->returnAlignInferInfoVecSize();
		}
		// log_ofs << "after inserting SJs generated from alignments, junctionNum = " << junctionNum << endl;
		// log_ofs << "start to insert SJs generated from annotation (if provided) " << endl;
		// settings_log_ofs << "after inserting SJs generated from alignments, junctionNum = " << junctionNum << endl;
		// settings_log_ofs << "start to insert SJs generated from annotation (if provided) " << endl;
		// cout << "after inserting SJs generated from alignments, junctionNum = " << junctionNum << endl;
		// cout << "start to insert SJs generated from annotation (if provided) " << endl;
		if(annotation_provided_bool)
		{	
			// loading SJs in annotation file (if provided)
			FILE *fp_annotatedSJ_file = fopen(annotation_file_path.c_str(), "r");
			fgets(entry, sizeof(entry), fp_annotatedSJ_file);
			while(!feof(fp_annotatedSJ_file))
			{
				fgets(entry, sizeof(entry), fp_annotatedSJ_file);
				if(feof(fp_annotatedSJ_file))
					break;
				junctionNum_in_annotation ++;
				entryString = entry;
				tabLocation1 = entryString.find('\t', 0);
				tabLocation2 = entryString.find('\t', tabLocation1+1);
				tabLocation3 = entryString.find('\t', tabLocation2+1);
				chrIntString = entryString.substr(0, tabLocation1);
				spliceStartPosString = entryString.substr(tabLocation1+1, tabLocation2-tabLocation1-1);
				spliceEndPosString = entryString.substr(tabLocation2+1, tabLocation3-tabLocation2-1);
				chrInt = indexInfo->convertStringToInt(chrIntString);
				spliceStartPos = atoi(spliceStartPosString.c_str());
				spliceEndPos = atoi(spliceEndPosString.c_str());	
				SJ->insert2AreaAndStringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
			}
			fclose(fp_annotatedSJ_file);
		}
		junctionNum = junctionNum_in_alignInferJuncHashInfo + junctionNum_in_annotation;
		if(junctionNum == 0)
		{
			spliceJunctionHashExists = false;
		}

		log_ofs << "After inserting SJs generated from alignments and annotation" << endl;
		log_ofs << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
		log_ofs << "\tjunctionNum in annotation = " << junctionNum_in_annotation << endl;
		log_ofs << "finish building spliceJunction Hash" << endl;		
		log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		settings_log_ofs << "After inserting SJs generated from alignments and annotation" << endl;
		settings_log_ofs << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
		settings_log_ofs << "\tjunctionNum in annotation = " << junctionNum_in_annotation << endl;
		settings_log_ofs << "finish building spliceJunction Hash" << endl;		
		settings_log_ofs << "start doing remapping on unfixed head/tail alignments" << endl;
		cout << "After inserting SJs generated from alignments and annotation" << endl;
		cout << "\tjunctionNum in alignments = " << junctionNum_in_alignInferJuncHashInfo << endl;
		cout << "\tjunctionNum in annotation = " << junctionNum_in_annotation << endl;
		cout << "finish building spliceJunction Hash" << endl;		
		cout << "start doing remapping on unfixed head/tail alignments" << endl;

		string tmpAlignIncompletePair = inputAlignInfoStr_incomplete;
		string headTailSoftClippingFile = tmpAlignIncompletePair;// + ".all";

		ifstream inputUnfixedHeadTailRecord_ifs(headTailSoftClippingFile.c_str());

		int normalRecordNum = normalRecordNum_fixHeadTail; //1000000;
		bool EndOfRecord = false;
		int tmpTurn = 0;
		int realRecordNum;// = normalRecordNum;

		for(tmpTurn = 0; /*tmpTurn < TurnNum*/; tmpTurn++)
		{
			vector<string> line1StrVec(normalRecordNum);
			vector<string> line2StrVec(normalRecordNum);
			vector<string> line3StrVec(normalRecordNum);
			vector<string> line4StrVec(normalRecordNum);
			vector<string> line5StrVec(normalRecordNum);
			vector<string> line6StrVec(normalRecordNum);
			vector<string> line7StrVec(normalRecordNum);
			vector<string> line8StrVec(normalRecordNum);
			vector<string> line9StrVec(normalRecordNum);
			vector<string> line10StrVec(normalRecordNum);
			vector<string> peAlignSamVec_complete_pair(normalRecordNum);
			vector<string> peAlignSamVec_incomplete_pair(normalRecordNum);
			vector<string> peAlignSamVec_complete_unpair(normalRecordNum);
			vector<string> peAlignSamVec_incomplete_unpair(normalRecordNum);
			vector<string> peAlignSamVec_pair_lowScore(normalRecordNum);

			vector<string> seAlignSamVec_complete(normalRecordNum);
			vector<string> seAlignSamVec_incomplete(normalRecordNum);
			vector<string> seAlignSamVec_lowScore(normalRecordNum);
			//vector<string> SeAlignSamStrVec_inComplete(normalRecordNum);

			if(EndOfRecord)
				break;

			int recordNum = normalRecordNum;

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;
			cout << endl << "[" << asctime(local)
				<< "start to read Head/Tail file record, turn: " << tmpTurn+1 << endl;

			realRecordNum = normalRecordNum;

			for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
			{
				string line1, line2, line3, line4, line5, line6, line7, 
					line8, line9, line10, line11;

				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}				
				getline(inputUnfixedHeadTailRecord_ifs, line11);
				if(inputUnfixedHeadTailRecord_ifs.eof())
				{
					realRecordNum = recordNumTmp;
					EndOfRecord = true;
					break;
				}		
				getline(inputUnfixedHeadTailRecord_ifs, line1);
				getline(inputUnfixedHeadTailRecord_ifs, line2);
				getline(inputUnfixedHeadTailRecord_ifs, line3);
				getline(inputUnfixedHeadTailRecord_ifs, line4);
				getline(inputUnfixedHeadTailRecord_ifs, line5);
				getline(inputUnfixedHeadTailRecord_ifs, line6);
				getline(inputUnfixedHeadTailRecord_ifs, line7);
				getline(inputUnfixedHeadTailRecord_ifs, line8);
				getline(inputUnfixedHeadTailRecord_ifs, line9);
				getline(inputUnfixedHeadTailRecord_ifs, line10);
				//getline(inputUnfixedHeadTailRecord_ifs, line11);

				line1StrVec[recordNumTmp] = line1;
				line2StrVec[recordNumTmp] = line2;
				line3StrVec[recordNumTmp] = line3;
				line4StrVec[recordNumTmp] = line4;
				line5StrVec[recordNumTmp] = line5;
				line6StrVec[recordNumTmp] = line6;
				line7StrVec[recordNumTmp] = line7;
				line8StrVec[recordNumTmp] = line8;
				line9StrVec[recordNumTmp] = line9;
				line10StrVec[recordNumTmp] = line10;
			}

			runtime_log_ofs << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			//cout << "realRecordNum: " << realRecordNum << " turn: " << tmpTurn+1 << endl;
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;					
			cout << endl << "[" << asctime(local)
				<< "finish reading Head/Tail records file, turn: " << tmpTurn+1 << endl;
			cout << endl << "[" << asctime(local)
				<< "start to fix Head/Tail, turn: " << tmpTurn+1 << endl;					

			omp_set_num_threads(threads_num);
			//omp_set_num_threads(1);
			#pragma omp parallel for			
			for(int tmpOpenMP = 0; tmpOpenMP < realRecordNum; tmpOpenMP++)
			{
				////////////////  parse long head reads record after 1-mapping process  ///////////////////////////////////////
				int threadNO = omp_get_thread_num();

				PE_Read_Info peReadInfo;// = new PE_Read_Info();
				PE_Read_Alignment_Info* peAlignInfo = new PE_Read_Alignment_Info();
				peAlignInfo->generatePeReadInfoAndPeAlignInfo_toFixIncompleteAlignment_getline(
					line1StrVec[tmpOpenMP], line2StrVec[tmpOpenMP], line3StrVec[tmpOpenMP],
					line4StrVec[tmpOpenMP], line5StrVec[tmpOpenMP], line6StrVec[tmpOpenMP],
					line7StrVec[tmpOpenMP], line8StrVec[tmpOpenMP], //line9StrVec[tmpOpenMP],
					line9StrVec[tmpOpenMP], line10StrVec[tmpOpenMP], peReadInfo, 
					indexInfo, fasta_or_fastq_bool, SE_or_PE_bool);		

				//cout << "tmpReadName: " << endl << line1StrVec[tmpOpenMP] << endl;
				//cout << endl << endl << "readName_1: " << peReadInfo.returnReadName_1() << endl;
				//cout << "readName_2: " << peReadInfo.returnReadName_2() << endl;
				
				#ifdef MAP_INFO
				cout << endl << endl << "readName_1: " << peReadInfo.returnReadName_1() << endl;
				cout << "readName_2: " << peReadInfo.returnReadName_2() << endl;
				cout << "start fixHeadTail: " << endl;
				cout << "PeAlignInfo:" << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				#endif

				FixHeadTailInfo* fixHeadTailInfo = new FixHeadTailInfo();
				if(Do_fixHeadTail_remapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);	
				}
				#ifdef MAP_INFO
				cout << "after remapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				#endif
				
				if(Do_fixHeadTail_greedyMapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_greedyMappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail,
						annotation_provided_bool, Do_annotation_only_bool, annotationInfo,
						MAX_SPLICE_DISTANCE_PHASE2, checkQualSeqForReadSegSeq, SE_or_PE_bool
						);	
				}				
				#ifdef MAP_INFO
				cout << "after greedyMapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;
				#endif
				
				if(Do_fixHeadTail_remappingAndTargetMapping)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingAndTargetMapping(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail, checkQualSeqForReadSegSeq, SE_or_PE_bool);	
				}
				#ifdef MAP_INFO
				cout << "after remapping And Target Mapping:" << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;				
				#endif
				if(Do_fixHeadTail_remappingAgain)
				{
					fixHeadTailInfo->fixHeadTail_areaAndStringHash_new_remappingOnly(
						peReadInfo, peAlignInfo, SJ, 
						secondLevelChrom,
						secondLevelSa,
						secondLevelLcpCompress,
						secondLevelChildTab,
						secondLevelDetChild,
						spliceJunctionHashExists,
						indexInfo, Do_extendHeadTail_fixHeadTail, SE_or_PE_bool);
				}
				
				fixHeadTailInfo->fixHeadTail_extend2end_finalStepForAligner(peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);
				fixHeadTailInfo->fixHeadTail_extend2end_fixIndel(peReadInfo, peAlignInfo, indexInfo, SE_or_PE_bool);	
				#ifdef MAP_INFO
				cout << "after extending: " << endl;
				cout << "PeAlignInfo: " << endl << peAlignInfo->returnPeAlignInfoStr() << endl;								
				#endif
				// remove duplicate mismatch
				peAlignInfo->removeDuplicateMismatch(SE_or_PE_bool);

				if(SE_or_PE_bool)
					peAlignInfo->chooseBestAlignment_final_SE();
				else
					peAlignInfo->chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
				//cout << "seAlignVec_final: " << endl << peAlignInfo->returnSeAlignVec_final() << endl;
				if(SE_or_PE_bool)
				{
					seAlignSamVec_complete[tmpOpenMP] = "";
					seAlignSamVec_incomplete[tmpOpenMP] = "";
					seAlignSamVec_lowScore[tmpOpenMP] = "";					

					bool allFinalAlignComplete_bool 
						= peAlignInfo->allFinalAlignmentComplete_SE();
					bool unique_bool = peAlignInfo->checkUniqueOrMulti_SE();
				    //statsInfo->increExistingAlignNum_phase1_SE(threadNO, allFinalAlignComplete_bool, unique_bool)
				    if(allFinalAlignComplete_bool)
				    {	    	
				    	bool completeAlignmentScore_tooLow 
				    		= peAlignInfo->alignScoreTooLow_bool_SE(peReadInfo);
				    	if(completeAlignmentScore_tooLow)
				    		seAlignSamVec_lowScore[tmpOpenMP] 
				    			= peAlignInfo->getSAMformatForUnmapped_SE(
				    				peReadInfo, fasta_or_fastq_bool);
					   	else
				    		seAlignSamVec_complete[tmpOpenMP] 
				    			= peAlignInfo->getSAMformatForFinalAlignment_SE(
				    				peReadInfo, fasta_or_fastq_bool);
				    }
				    else
				    	seAlignSamVec_incomplete[tmpOpenMP] 
				    		= peAlignInfo->getSAMformatForFinalAlignment_SE(
				    			peReadInfo, fasta_or_fastq_bool);
				}				
				else
				{	
					bool pairExistsBool = peAlignInfo->finalPairExistsBool();									
					peAlignSamVec_complete_pair[tmpOpenMP] = "";
					peAlignSamVec_incomplete_pair[tmpOpenMP] = "";
					peAlignSamVec_complete_unpair[tmpOpenMP] = "";
					peAlignSamVec_incomplete_unpair[tmpOpenMP] = "";
					peAlignSamVec_pair_lowScore[tmpOpenMP] = "";

					if(pairExistsBool) // some pair exists, all completed, print out paired SAM info
					{
						bool allFinalPairAlignmentCompleteBool = peAlignInfo->allAlignmentInFinalPairCompleted();	
						bool unique_bool = peAlignInfo->checkUniqueOrMulti();
						bool align_lowScore_bool = peAlignInfo->alignPairScoreTooLow_bool(peReadInfo);
						if(align_lowScore_bool)
							peAlignSamVec_pair_lowScore[tmpOpenMP] 
								= peAlignInfo->getSAMformatForBothEndsUnmapped(peReadInfo, fasta_or_fastq_bool);
						else if(allFinalPairAlignmentCompleteBool)
							peAlignSamVec_complete_pair[tmpOpenMP] 
								= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(
								peReadInfo, fasta_or_fastq_bool);
						else
							peAlignSamVec_incomplete_pair[tmpOpenMP] 
								= peAlignInfo->getSAMformatForFinalPair_secondaryOrNot(
									peReadInfo, fasta_or_fastq_bool);
					}
					else //if((!pairExistsBool) && (allAlignmentCompleteBool)) // no pair exists, all complete, print out original SAM info
					{
						bool allUnpairAlignmentCompleteBool = peAlignInfo->allUnpairedAlignmentCompleted();
						//statsInfo->increUnpairedNum_fixHeadTail(threadNO, allUnpairAlignmentCompleteBool);
						if(allUnpairAlignmentCompleteBool)
							peAlignSamVec_complete_unpair[tmpOpenMP] = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
								peReadInfo, fasta_or_fastq_bool);
						else
							peAlignSamVec_incomplete_unpair[tmpOpenMP] = peAlignInfo->getSAMformatForUnpairedAlignments_secondaryOrNot(
								peReadInfo, fasta_or_fastq_bool);	
					}
				}

				delete fixHeadTailInfo;
				peAlignInfo->memoryFree();
				delete peAlignInfo;
			}
			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "start to output ... turn: " << tmpTurn+1 << endl;
			cout << endl << "[" << asctime(local)
				<< "finish fixing Head/Tail, turn: " << tmpTurn+1 << endl;// << endl;
			cout << endl << "[" << asctime(local)
				<< "start to output ... turn: " << tmpTurn+1 << endl;
			if(SE_or_PE_bool)
			{
				for(int tmp = 0; tmp < realRecordNum; tmp++)
				{
					if(seAlignSamVec_complete[tmp] != "")
					{
						OutputSamFile_fixHeadTail_complete_SE_ofs << seAlignSamVec_complete[tmp] << endl;
					}					
					if(seAlignSamVec_incomplete[tmp] != "")
					{
						OutputSamFile_fixHeadTail_incomplete_SE_ofs << seAlignSamVec_incomplete[tmp] << endl;
					}
					if(seAlignSamVec_lowScore[tmp] != "")
					{
						OutputSamFile_fixHeadTail_lowScore_SE_ofs << seAlignSamVec_lowScore[tmp] << endl;
					}
				}				
			}
			else
			{	
				for(int tmp = 0; tmp < realRecordNum; tmp++)
				{
					if(peAlignSamVec_complete_pair[tmp] != "")
					{
						OutputSamFile_fixHeadTail_complete_pair_ofs << peAlignSamVec_complete_pair[tmp] << endl;
					}
					if(peAlignSamVec_incomplete_pair[tmp] != "")
					{	
						OutputSamFile_fixHeadTail_incomplete_pair_ofs << peAlignSamVec_incomplete_pair[tmp] << endl;
					}
					if(peAlignSamVec_complete_unpair[tmp] != "")
					{
						OutputSamFile_fixHeadTail_complete_unpair_ofs << peAlignSamVec_complete_unpair[tmp] << endl;
					}
					if(peAlignSamVec_incomplete_unpair[tmp] != "")	
					{
						OutputSamFile_fixHeadTail_incomplete_unpair_ofs << peAlignSamVec_incomplete_unpair[tmp] << endl;
					}
					if(peAlignSamVec_pair_lowScore[tmp] != "")
					{
						OutputSamFile_fixHeadTail_pair_lowScore_ofs << peAlignSamVec_pair_lowScore[tmp] << endl;
					}
				}		
			}		

			nowtime = time(NULL);
			local = localtime(&nowtime);
			runtime_log_ofs << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;
			cout << endl << "[" << asctime(local)
				<< "finish output, turn: " << tmpTurn+1 << endl << endl;
		}
		inputUnfixedHeadTailRecord_ifs.close();
	}
	delete(SJ);

	OutputSamFile_fixHeadTail_complete_pair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_pair_ofs.close();
	OutputSamFile_fixHeadTail_complete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_unpair_ofs.close();
	OutputSamFile_fixHeadTail_pair_lowScore_ofs.close();

	OutputSamFile_fixHeadTail_complete_SE_ofs.close();
	OutputSamFile_fixHeadTail_incomplete_SE_ofs.close();
	OutputSamFile_fixHeadTail_lowScore_SE_ofs.close();
	
	nowtime = time(NULL);
	local = localtime(&nowtime);

	cout << endl << "[" << asctime(local) 
		<< "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	log_ofs << endl << "[" << asctime(local) 
		<< "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  
	runtime_log_ofs << endl << "[" << asctime(local) 
		<< "... fixing unfixed-head/tail reads ends ......" << endl << endl ;  


	return 0;
}