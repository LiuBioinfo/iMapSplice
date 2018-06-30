// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef PHASE1_PARALLELPROCESSES_H
#define PHASE1_PARALLELPROCESSES_H

using namespace std;

void io_stage_phase1(ifstream& input_ifs_1,
	ifstream& input_ifs_2,
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,
	int inputReadNumInBatchArray_phase1,
	int inputTimePerc_phase1,
	int outputTimePerc_phase1,
	ofstream& log_ofs,
	InputReadPreProcess* readPreProcessInfo,
	ofstream& tmpAlignCompleteRead_ofs,
	ofstream& tmpAlignIncompletePair_ofs,
	ofstream& tmpAlignOneEndUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_lowScore_ofs,
	ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
	ofstream& tmpAlignIncompletePair_SAM_ofs,
	ofstream& input_log_ofs,
	ofstream& output_log_ofs,
	bool fasta_or_fastq_bool,
	bool SE_or_PE_bool
	)
{
	int tmpBatchIndex_input = 0;
	int tmpBatchIndex_output = 0;

	int tmpThread = omp_get_thread_num();
	log_ofs << "input thread: " << tmpThread << endl;
	
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	input_log_ofs << endl << "[" << asctime(local) << "... input of Phase1 starts ......" << endl << endl;  

	string readNameStr_1, readNameStr_2, readSeqStr_1, readSeqStr_2;
		
	getline(input_ifs_1, readNameStr_1);
	getline(input_ifs_1, readSeqStr_1);
	getline(input_ifs_2, readNameStr_2);
	getline(input_ifs_2, readSeqStr_2);

	if(fasta_or_fastq_bool)
	{	
		readArrayQueue->initiateWith1stRead(readNameStr_1, readNameStr_2,
			readSeqStr_1, readSeqStr_2, input_log_ofs);
		tmpBatchIndex_input ++;
	}
	else
	{
		string readCommentStr_1, readCommentStr_2, readQualSeq_1, readQualSeq_2;
		getline(input_ifs_1, readCommentStr_1);
		getline(input_ifs_1, readQualSeq_1);
		getline(input_ifs_2, readCommentStr_2);
		getline(input_ifs_2, readQualSeq_2);
		readArrayQueue->initiateWith1stRead_fq(readNameStr_1, readNameStr_2,
			readSeqStr_1, readSeqStr_2, readQualSeq_1, readQualSeq_2, input_log_ofs);		
		tmpBatchIndex_input ++;
	}

	bool input_stage_end_bool = false;
	bool output_stage_end_bool = false;
	while(1)
	{
		if(!input_stage_end_bool)
		{
			// read fa/fq file and insert the readArray into readArrayList
			for(int tmp = 0; tmp < inputTimePerc_phase1; tmp++)
			{
				if((input_ifs_1.eof())||(input_ifs_2.eof()))
				{
					time_t nowtime;
					nowtime = time(NULL);
					struct tm *local;
					local = localtime(&nowtime);
					input_log_ofs << endl << "[" << asctime(local) 
						<< "... end of input ......" << endl << endl;  
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;
				}
				string tmpReadName_1, tmpReadName_2;
				getline(input_ifs_1, tmpReadName_1);
				getline(input_ifs_2, tmpReadName_2); 
				if((input_ifs_1.eof())||(input_ifs_2.eof()))
				{
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;
				}
				string tmpReadSeq_1, tmpReadSeq_2;
				getline(input_ifs_1, tmpReadSeq_1);
				getline(input_ifs_2, tmpReadSeq_2);

				if(fasta_or_fastq_bool)
				{
					readArrayQueue->getSeqFromInputFile(
						tmpReadName_1.substr(1), 
						tmpReadName_2.substr(1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2),
						tmpReadSeq_1,
						tmpReadSeq_2,
						inputReadNumInBatchArray_phase1, input_log_ofs, tmpBatchIndex_input);
				}
				else
				{
					string tmpCommentStr_1, tmpCommentStr_2, tmpReadQualSeq_1, tmpReadQualSeq_2;
					getline(input_ifs_1, tmpCommentStr_1);
					getline(input_ifs_1, tmpReadQualSeq_1);
					getline(input_ifs_2, tmpCommentStr_2);
					getline(input_ifs_2, tmpReadQualSeq_2);
					readArrayQueue->getSeqFromInputFile_fq(
						tmpReadName_1.substr(1), 
						tmpReadName_2.substr(1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2),
						tmpReadSeq_1,
						tmpReadSeq_2,						
						tmpReadQualSeq_1,
						tmpReadQualSeq_2,
						inputReadNumInBatchArray_phase1, input_log_ofs, tmpBatchIndex_input);
				}
			}
		}
		if(!output_stage_end_bool)
		{	
			for(int tmp = 0; tmp < outputTimePerc_phase1; tmp++)
			{		
				if(resultArrayQueue->atLeast3Node())
				{
					time_t nowtime;
					nowtime = time(NULL);
					struct tm *local;
					local = localtime(&nowtime);
					tmpBatchIndex_output ++;
					output_log_ofs << endl << "[" << asctime(local) 
						<< "... output regular array to file starts ......" << endl;  					
					output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
					resultArrayQueue->outputFrontResultArray(
						tmpAlignCompleteRead_ofs,
						tmpAlignIncompletePair_ofs,
						tmpAlignOneEndUnmapped_ofs,
						tmpAlignBothEndsUnmapped_ofs,
						tmpAlignBothEndsUnmapped_lowScore_ofs,
						tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
						tmpAlignIncompletePair_SAM_ofs);
					resultArrayQueue->popFromResultQueue();
				}
				else if(resultArrayQueue->only2Node())
				{
					if(endOfProcessing_bool)
					{	
						time_t nowtime;
						nowtime = time(NULL);
						struct tm *local;
						local = localtime(&nowtime);
						tmpBatchIndex_output ++;
						output_log_ofs << endl << "[" << asctime(local) 
							<< "... output last array to file starts ......" << endl;
						output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
						resultArrayQueue->outputFrontResultArray(
							tmpAlignCompleteRead_ofs,
							tmpAlignIncompletePair_ofs,
							tmpAlignOneEndUnmapped_ofs,
							tmpAlignBothEndsUnmapped_ofs,
							tmpAlignBothEndsUnmapped_lowScore_ofs,
							tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
							tmpAlignIncompletePair_SAM_ofs);
						resultArrayQueue->popFromResultQueue();
					}
					else
						continue;
				}
				else if(resultArrayQueue->resultQueueEmpty())
				{
					if(endOfProcessing_bool)
					{
						time_t nowtime;
						nowtime = time(NULL);
						struct tm *local;
						local = localtime(&nowtime);
						output_log_ofs << endl << "[" << asctime(local) 
							<< "... end of output ......" << endl << endl;  
						output_stage_end_bool = true;
						break;	
					}
				}
				else
				{
					//cout << "other cases in output stage" << endl;
					continue;
				}
			}
		}
		if(input_stage_end_bool && output_stage_end_bool)
			break;
	}
	log_ofs << "end of input/output " << endl;
}

/*
void io_stage_phase1_new(ifstream& input_ifs_1,
	ifstream& input_ifs_2,
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,
	int inputReadNumInBatchArray_phase1,
	int inputTimePerc_phase1,
	int outputTimePerc_phase1,
	ofstream& log_ofs,
	InputReadPreProcess* readPreProcessInfo,
	ofstream& tmpAlignCompleteRead_ofs,
	ofstream& tmpAlignIncompletePair_ofs,
	ofstream& tmpAlignOneEndUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_ofs,
	ofstream& tmpAlignBothEndsUnmapped_lowScore_ofs,
	ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
	ofstream& tmpAlignIncompletePair_SAM_ofs,
	ofstream& input_log_ofs,
	ofstream& output_log_ofs,
	bool fasta_or_fastq_bool,
	bool SE_or_PE_bool
	)
{
	int tmpBatchIndex_input = 0;
	int tmpBatchIndex_output = 0;

	int tmpThread = omp_get_thread_num();
	log_ofs << "input thread: " << tmpThread << endl;
	
	time_t nowtime;
	nowtime = time(NULL);
	struct tm *local;
	local = localtime(&nowtime);
	input_log_ofs << endl << "[" << asctime(local) << "... input of Phase1 starts ......" << endl << endl;  

	string readNameStr_1, readNameStr_2, readSeqStr_1, readSeqStr_2;
		
	getline(input_ifs_1, readNameStr_1);
	getline(input_ifs_1, readSeqStr_1);
	getline(input_ifs_2, readNameStr_2);
	getline(input_ifs_2, readSeqStr_2);

	if(fasta_or_fastq_bool)
	{	
		readArrayQueue->initiateWith1stRead(readNameStr_1, readNameStr_2,
			readSeqStr_1, readSeqStr_2, input_log_ofs);
		tmpBatchIndex_input ++;
	}
	else
	{
		string readCommentStr_1, readCommentStr_2, readQualSeq_1, readQualSeq_2;
		getline(input_ifs_1, readCommentStr_1);
		getline(input_ifs_1, readQualSeq_1);
		getline(input_ifs_2, readCommentStr_2);
		getline(input_ifs_2, readQualSeq_2);
		readArrayQueue->initiateWith1stRead_fq(readNameStr_1, readNameStr_2,
			readSeqStr_1, readSeqStr_2, readQualSeq_1, readQualSeq_2, input_log_ofs);		
		tmpBatchIndex_input ++;
	}

	bool input_stage_end_bool = false;
	bool output_stage_end_bool = false;
	while(1)
	{
		if(!input_stage_end_bool)
		{
			// read fa/fq file and insert the readArray into readArrayList
			for(int tmp = 0; tmp < inputTimePerc_phase1; tmp++)
			{
				if((input_ifs_1.eof())||(input_ifs_2.eof()))
				{
					time_t nowtime;
					nowtime = time(NULL);
					struct tm *local;
					local = localtime(&nowtime);
					input_log_ofs << endl << "[" << asctime(local) 
						<< "... end of input ......" << endl << endl;  
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;
				}
				string tmpReadName_1, tmpReadName_2;
				getline(input_ifs_1, tmpReadName_1);
				getline(input_ifs_2, tmpReadName_2); 
				if((input_ifs_1.eof())||(input_ifs_2.eof()))
				{
					input_stage_end_bool = true;
					endOfFile_bool = true;
					break;
				}
				string tmpReadSeq_1, tmpReadSeq_2;
				getline(input_ifs_1, tmpReadSeq_1);
				getline(input_ifs_2, tmpReadSeq_2);

				if(fasta_or_fastq_bool)
				{
					readArrayQueue->getSeqFromInputFile(
						tmpReadName_1.substr(1), 
						tmpReadName_2.substr(1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2),
						tmpReadSeq_1,
						tmpReadSeq_2,
						inputReadNumInBatchArray_phase1, input_log_ofs, tmpBatchIndex_input);
				}
				else
				{
					string tmpCommentStr_1, tmpCommentStr_2, tmpReadQualSeq_1, tmpReadQualSeq_2;
					getline(input_ifs_1, tmpCommentStr_1);
					getline(input_ifs_1, tmpReadQualSeq_1);
					getline(input_ifs_2, tmpCommentStr_2);
					getline(input_ifs_2, tmpReadQualSeq_2);
					readArrayQueue->getSeqFromInputFile_fq(
						tmpReadName_1.substr(1), 
						tmpReadName_2.substr(1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_1),
						//readPreProcessInfo->upperCaseReadSeq(tmpReadSeq_2),
						tmpReadSeq_1,
						tmpReadSeq_2,						
						tmpReadQualSeq_1,
						tmpReadQualSeq_2,
						inputReadNumInBatchArray_phase1, input_log_ofs, tmpBatchIndex_input);
				}
			}
		}
		if(!output_stage_end_bool)
		{	
			for(int tmp = 0; tmp < outputTimePerc_phase1; tmp++)
			{		
				if(resultArrayQueue->atLeast3Node())
				{
					time_t nowtime;
					nowtime = time(NULL);
					struct tm *local;
					local = localtime(&nowtime);
					tmpBatchIndex_output ++;
					output_log_ofs << endl << "[" << asctime(local) 
						<< "... output regular array to file starts ......" << endl;  					
					output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
					resultArrayQueue->outputFrontResultArray(
						tmpAlignCompleteRead_ofs,
						tmpAlignIncompletePair_ofs,
						tmpAlignOneEndUnmapped_ofs,
						tmpAlignBothEndsUnmapped_ofs,
						tmpAlignBothEndsUnmapped_lowScore_ofs,
						tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
						tmpAlignIncompletePair_SAM_ofs);
					resultArrayQueue->popFromResultQueue();
				}
				else if(resultArrayQueue->only2Node())
				{
					if(endOfProcessing_bool)
					{	
						time_t nowtime;
						nowtime = time(NULL);
						struct tm *local;
						local = localtime(&nowtime);
						tmpBatchIndex_output ++;
						output_log_ofs << endl << "[" << asctime(local) 
							<< "... output last array to file starts ......" << endl;
						output_log_ofs << "tmpBatchIndex: " << tmpBatchIndex_output << endl << endl;
						resultArrayQueue->outputFrontResultArray(
							tmpAlignCompleteRead_ofs,
							tmpAlignIncompletePair_ofs,
							tmpAlignOneEndUnmapped_ofs,
							tmpAlignBothEndsUnmapped_ofs,
							tmpAlignBothEndsUnmapped_lowScore_ofs,
							tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
							tmpAlignIncompletePair_SAM_ofs);
						resultArrayQueue->popFromResultQueue();
					}
					else
						continue;
				}
				else if(resultArrayQueue->resultQueueEmpty())
				{
					if(endOfProcessing_bool)
					{
						time_t nowtime;
						nowtime = time(NULL);
						struct tm *local;
						local = localtime(&nowtime);
						output_log_ofs << endl << "[" << asctime(local) 
							<< "... end of output ......" << endl << endl;  
						output_stage_end_bool = true;
						break;	
					}
				}
				else
				{
					//cout << "other cases in output stage" << endl;
					continue;
				}
			}
		}
		if(input_stage_end_bool && output_stage_end_bool)
			break;
	}
	log_ofs << "end of input/output " << endl;
}*/

void process_stage_phase1(
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,	
	int threadNumForProcess, 

	unsigned int* sa, 
	BYTE* lcpCompress,
	unsigned int* childTab, 
	char* chrom,
	BYTE* verifyChild, 
	Index_Info* indexInfo,
	int* preIndexMapLengthArray, 
	unsigned int* preIndexIntervalStartArray,
	unsigned int* preIndexIntervalEndArray,
	vector<RepeatRegion_Info*>& repeatRegionInfoVec,
	bool Do_cirRNA, 
	bool Do_extendHeadTail_phase1,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,			
	bool outputDirectlyBool_Phase1Only,
	bool Do_Phase1_Only,	
	Stats_Info* statsInfo, 
	bool fasta_or_fastq_bool,
	ofstream& mapping_log_ofs,//, vector<ofstream*>& mapping_log_ofs_vec
	bool checkQualSeqForReadSegSeq,
	bool SE_or_PE_bool
	)
{
	int tmpBatchIndex = 0;
	while(1)
	{
		int tmpBatchArraySize = 0;
		if(readArrayQueue->atLeast3Node())
		{
			mapping_log_ofs << endl << "at least 3 node" << endl;
			tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
		}
		else if(readArrayQueue->only2Node())
		{
			//mapping_log_ofs << "only 2 node" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... only 2 node -- end of file" << endl;
				tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
			}
			else
			{
				//mapping_log_ofs << "... only 2 node -- not end of file" << endl;
				continue;
			}
		}
		else if(readArrayQueue->readQueueEmpty())
		{
			//mapping_log_ofs << "readQueueEmpty" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... readQueueEmpty -- end of file" << endl;
				break;
			}
			else
			{
				//mapping_log_ofs << "... readQueueEmpty -- not end of file" << endl;
				continue;
			}
		}
		else
		{
			//mapping_log_ofs << "other cases" << endl;
			continue;
		}
		time_t nowtime;
		nowtime = time(NULL);
		struct tm *local;
		local = localtime(&nowtime);
		mapping_log_ofs << "*******************************************************************************" 
			<< endl << "[" << asctime(local) 
			<< "... mapping a batch of reads starts ......" << endl;  
		tmpBatchIndex ++;
		mapping_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl;
		mapping_log_ofs << "threadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "tmpBatchArraySize: " << tmpBatchArraySize << endl;			

		Result_Array* tmpResultArray = new Result_Array(tmpBatchArraySize);
		omp_set_num_threads(threadNumForProcess);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < tmpBatchArraySize; tmpOpenMP++)
		{	
			int threadNO = omp_get_thread_num();
			
			string tmpReadName_1 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_1(tmpOpenMP);
			string tmpReadSeq_1 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_1(tmpOpenMP);
			string tmpReadName_2 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_2(tmpOpenMP);
			string tmpReadSeq_2 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_2(tmpOpenMP);	
			
			string tmpReadQualSeq_1 = "";
			string tmpReadQualSeq_2 = "";

			if(!fasta_or_fastq_bool)
			{	
				tmpReadQualSeq_1 = readArrayQueue->returnFrontNodeReadQualSeq_1(tmpOpenMP);
				tmpReadQualSeq_2 = readArrayQueue->returnFrontNodeReadQualSeq_2(tmpOpenMP);
			}

 			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(tmpReadName_1, tmpReadName_2, tmpReadSeq_1, tmpReadSeq_2,
				tmpReadQualSeq_1, tmpReadQualSeq_2, fasta_or_fastq_bool, SE_or_PE_bool);	
    		FixPhase1Info fixPhase1Info;
			fixPhase1Info.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);						
			
			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, 
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, 
				MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);

			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);
			
			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.initiatePeAlignInfo(
				fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, 
				indexInfo, SE_or_PE_bool);		

			if(Do_Phase1_Only)
			{
				peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			}
			else
			{	
				peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(),
					readInfo.returnReadSeqLength_2());
			}
  			
			peAlignInfo.output_phase1(
				outputDirectlyBool_Phase1Only,
				Do_Phase1_Only,				
				tmpResultArray,
				repeatRegionInfoVec,
				readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo,
				fasta_or_fastq_bool);

			fixPhase1Info.memoryFree();
			peAlignInfo.memoryFree();

		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		mapping_log_ofs << endl << "*******************************************************************************"
			<< endl << "[" << asctime(local) 
			<< "... mapping a batch of reads ends ......" << endl << endl << endl; 		
		resultArrayQueue->pushBack2ResultArrayQueue(tmpResultArray);
		readArrayQueue->popFromReadQueue();
	}
	endOfProcessing_bool = true;
}

void process_stage_phase1_mpsI(
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,	
	int threadNumForProcess,
	unsigned int* sa, 
	BYTE* lcpCompress,
	unsigned int* childTab, 
	char* chrom,
	BYTE* verifyChild, 
	Index_Info* indexInfo,
	int* preIndexMapLengthArray, 
	unsigned int* preIndexIntervalStartArray,
	unsigned int* preIndexIntervalEndArray,
	vector<RepeatRegion_Info*>& repeatRegionInfoVec,
	bool Do_cirRNA, 
	bool Do_extendHeadTail_phase1,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,			
	bool outputDirectlyBool_Phase1Only,
	bool Do_Phase1_Only,	
	Stats_Info* statsInfo, 
	bool fasta_or_fastq_bool,
	ofstream& mapping_log_ofs,//, vector<ofstream*>& mapping_log_ofs_vec
	bool checkQualSeqForReadSegSeq,
	bool SE_or_PE_bool,
	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq)
{
	int tmpBatchIndex = 0;
	while(1)
	{
		int tmpBatchArraySize = 0;
		if(readArrayQueue->atLeast3Node())
		{
			mapping_log_ofs << endl << "at least 3 node" << endl;
			tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
		}
		else if(readArrayQueue->only2Node())
		{
			//mapping_log_ofs << "only 2 node" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... only 2 node -- end of file" << endl;
				tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
			}
			else
			{
				//mapping_log_ofs << "... only 2 node -- not end of file" << endl;
				continue;
			}
		}
		else if(readArrayQueue->readQueueEmpty())
		{
			//mapping_log_ofs << "readQueueEmpty" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... readQueueEmpty -- end of file" << endl;
				break;
			}
			else
			{
				//mapping_log_ofs << "... readQueueEmpty -- not end of file" << endl;
				continue;
			}
		}
		else
		{
			//mapping_log_ofs << "other cases" << endl;
			continue;
		}
		time_t nowtime;
		nowtime = time(NULL);
		struct tm *local;
		local = localtime(&nowtime);
		mapping_log_ofs << "*******************************************************************************" 
			<< endl << "[" << asctime(local) 
			<< "... mapping a batch of reads starts ......" << endl;  
		tmpBatchIndex ++;
		mapping_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl;
		mapping_log_ofs << "threadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "tmpBatchArraySize: " << tmpBatchArraySize << endl;			

		Result_Array* tmpResultArray = new Result_Array(tmpBatchArraySize);
		omp_set_num_threads(threadNumForProcess);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < tmpBatchArraySize; tmpOpenMP++)
		{	
			int threadNO = omp_get_thread_num();
			
			string tmpReadName_1 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_1(tmpOpenMP);
			string tmpReadSeq_1 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_1(tmpOpenMP);
			string tmpReadName_2 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_2(tmpOpenMP);
			string tmpReadSeq_2 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_2(tmpOpenMP);	
			string tmpReadQualSeq_1 = "";
			string tmpReadQualSeq_2 = "";
			if(!fasta_or_fastq_bool)
			{	
				tmpReadQualSeq_1 = readArrayQueue->returnFrontNodeReadQualSeq_1(tmpOpenMP);
				tmpReadQualSeq_2 = readArrayQueue->returnFrontNodeReadQualSeq_2(tmpOpenMP);
			}

 			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(tmpReadName_1, tmpReadName_2, tmpReadSeq_1, tmpReadSeq_2,
				tmpReadQualSeq_1, tmpReadQualSeq_2, fasta_or_fastq_bool, SE_or_PE_bool);	
    		FixPhase1Info fixPhase1Info;
			fixPhase1Info.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);						

			#ifdef PERSONALIZED_CHR_SEQ
			fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			#endif

			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);

			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);
			
			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.initiatePeAlignInfo(fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo, SE_or_PE_bool);		

			if(Do_Phase1_Only)
			{
				peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			}
			else
			{	
				peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(),
					readInfo.returnReadSeqLength_2());
			}
			int overlapLength_top2candiAlignment = peAlignInfo.getMaxOverlapLengthInAlignmentPair(
				readInfo.returnReadSeqLength_1(), readInfo.returnReadSeqLength_2());  			
			peAlignInfo.output_phase1(outputDirectlyBool_Phase1Only, Do_Phase1_Only, tmpResultArray, repeatRegionInfoVec, 
				readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo, fasta_or_fastq_bool, overlapLength_top2candiAlignment);

			fixPhase1Info.memoryFree();
			peAlignInfo.memoryFree();
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		mapping_log_ofs << endl << "*******************************************************************************"
			<< endl << "[" << asctime(local) << "... mapping a batch of reads ends ......" << endl << endl << endl; 		
		resultArrayQueue->pushBack2ResultArrayQueue(tmpResultArray);
		readArrayQueue->popFromReadQueue();
	}
	endOfProcessing_bool = true;
}


void process_stage_phase1_main(
	Read_Array_Queue* readArrayQueue,
	Result_Array_Queue* resultArrayQueue,
	bool& endOfFile_bool,
	bool& endOfProcessing_bool,	
	int threadNumForProcess,
	unsigned int* sa, 
	BYTE* lcpCompress,
	unsigned int* childTab, 
	char* chrom,
	BYTE* verifyChild, 
	Index_Info* indexInfo,
	int* preIndexMapLengthArray, 
	unsigned int* preIndexIntervalStartArray,
	unsigned int* preIndexIntervalEndArray,
	vector<RepeatRegion_Info*>& repeatRegionInfoVec,
	bool Do_cirRNA, 
	bool Do_extendHeadTail_phase1,
	bool annotation_provided_bool, 
	bool Do_annotation_only_bool, 
	Annotation_Info* annotationInfo,			
	bool outputDirectlyBool_Phase1Only,
	bool Do_Phase1_Only,	
	Stats_Info* statsInfo, 
	bool fasta_or_fastq_bool,
	ofstream& mapping_log_ofs,//, vector<ofstream*>& mapping_log_ofs_vec
	bool checkQualSeqForReadSegSeq,
	bool SE_or_PE_bool,
	unsigned int* sa_SNP, 
	BYTE* lcpCompress_SNP,
	unsigned int* childTab_SNP, 
	char* chrom_SNP,
	BYTE* verifyChild_SNP, 
	Index_Info* indexInfo_SNP,
	int SNPlocInSyntheticSNPseq)
{
	int tmpBatchIndex = 0;
	while(1)
	{
		int tmpBatchArraySize = 0;
		if(readArrayQueue->atLeast3Node())
		{
			mapping_log_ofs << endl << "at least 3 node" << endl;
			tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
		}
		else if(readArrayQueue->only2Node())
		{
			//mapping_log_ofs << "only 2 node" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... only 2 node -- end of file" << endl;
				tmpBatchArraySize = readArrayQueue->returnFrontNodeSize();
			}
			else
			{
				//mapping_log_ofs << "... only 2 node -- not end of file" << endl;
				continue;
			}
		}
		else if(readArrayQueue->readQueueEmpty())
		{
			//mapping_log_ofs << "readQueueEmpty" << endl;
			if(endOfFile_bool)
			{
				mapping_log_ofs << endl << "... readQueueEmpty -- end of file" << endl;
				break;
			}
			else
			{
				//mapping_log_ofs << "... readQueueEmpty -- not end of file" << endl;
				continue;
			}
		}
		else
		{
			//mapping_log_ofs << "other cases" << endl;
			continue;
		}
		time_t nowtime;
		nowtime = time(NULL);
		struct tm *local;
		local = localtime(&nowtime);
		mapping_log_ofs << "*******************************************************************************" 
			<< endl << "[" << asctime(local) 
			<< "... mapping a batch of reads starts ......" << endl;  
		tmpBatchIndex ++;
		mapping_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl;
		mapping_log_ofs << "threadNum: " << threadNumForProcess << endl;	
		mapping_log_ofs << "tmpBatchArraySize: " << tmpBatchArraySize << endl;			

		Result_Array* tmpResultArray = new Result_Array(tmpBatchArraySize);
		omp_set_num_threads(threadNumForProcess);
		#pragma omp parallel for schedule(dynamic)
		for(int tmpOpenMP = 0; tmpOpenMP < tmpBatchArraySize; tmpOpenMP++)
		{	
			int threadNO = omp_get_thread_num();
			
			string tmpReadName_1 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_1(tmpOpenMP);
			string tmpReadSeq_1 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_1(tmpOpenMP);
			string tmpReadName_2 = //"seq.1/1";
				readArrayQueue->returnFrontNodeReadName_2(tmpOpenMP);
			string tmpReadSeq_2 = //"AACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTATAACGTGCATCGACTGCTAT";
				readArrayQueue->returnFrontNodeReadSeq_2(tmpOpenMP);	
			string tmpReadQualSeq_1 = "";
			string tmpReadQualSeq_2 = "";
			if(!fasta_or_fastq_bool)
			{	
				tmpReadQualSeq_1 = readArrayQueue->returnFrontNodeReadQualSeq_1(tmpOpenMP);
				tmpReadQualSeq_2 = readArrayQueue->returnFrontNodeReadQualSeq_2(tmpOpenMP);
			}

 			PE_Read_Info readInfo; //= new PE_Read_Info();
			readInfo.initiateReadInfo(tmpReadName_1, tmpReadName_2, tmpReadSeq_1, tmpReadSeq_2,
				tmpReadQualSeq_1, tmpReadQualSeq_2, fasta_or_fastq_bool, SE_or_PE_bool);	
    		FixPhase1Info fixPhase1Info;
			fixPhase1Info.fixPhase1_segInfo(sa, lcpCompress, childTab, chrom, verifyChild, indexInfo, 
				preIndexMapLengthArray, preIndexIntervalStartArray, preIndexIntervalEndArray,
				readInfo, repeatRegionInfoVec[threadNO], checkQualSeqForReadSegSeq, SE_or_PE_bool);						

			#ifdef PERSONALIZED_CHR_SEQ
			# ifdef VARY_SNP_MER
			fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo_varySNPmer(sa_SNP, lcpCompress_SNP, 
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo);
			# else
			fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			# endif

			# ifdef VARY_SNP_MER
			fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo_varySNPmer(sa_SNP, lcpCompress_SNP, 
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo);
			# else
			fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
				childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			# endif
			#endif

			// #ifdef PERSONALIZED_CHR_SEQ
			// fixPhase1Info.fixPhase1_map2syntheticSNPtransSeq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
			// 	childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			// fixPhase1Info.fixPhase1_targetMap2syntheticSNPseq_updateOriSegInfo(sa_SNP, lcpCompress_SNP, 
			// 	childTab_SNP, chrom_SNP, verifyChild_SNP, indexInfo_SNP, readInfo, indexInfo, SNPlocInSyntheticSNPseq);
			// #endif

			fixPhase1Info.fixPhase1_pathInfo(Do_cirRNA, indexInfo, readInfo, annotation_provided_bool, 
				Do_annotation_only_bool, annotationInfo, MAX_SPLICE_DISTANCE_PHASE1, SE_or_PE_bool);

			fixPhase1Info.fixPhase1_gapInfo(readInfo, indexInfo, Do_cirRNA, Do_extendHeadTail_phase1,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo, SE_or_PE_bool);
			
			PE_Read_Alignment_Info peAlignInfo;
			peAlignInfo.initiatePeAlignInfo(fixPhase1Info.pathInfo_Nor1, fixPhase1Info.pathInfo_Rcm1, 
				fixPhase1Info.pathInfo_Nor2, fixPhase1Info.pathInfo_Rcm2, indexInfo, SE_or_PE_bool);		

			if(Do_Phase1_Only)
			{
				if(SE_or_PE_bool)
					peAlignInfo.chooseBestAlignment_final_SE();
				else
					peAlignInfo.chooseBestAlignment_selectAllIfMultiBest_filterOutNoncanonical_SJpenalty();
			}
			else
			{
				if(SE_or_PE_bool)
					peAlignInfo.alignmentFilter_fixPhase1_SE(readInfo.returnReadLength_SE());
				else
					peAlignInfo.alignmentFilter_fixPhase1_SJpenalty(readInfo.returnReadSeqLength_1(),
						readInfo.returnReadSeqLength_2());
			}
			int overlapLength_top2candiAlignment = peAlignInfo.getMaxOverlapLengthInAlignmentPair(
				readInfo.returnReadSeqLength_1(), readInfo.returnReadSeqLength_2());  			
			peAlignInfo.output_phase1(outputDirectlyBool_Phase1Only, Do_Phase1_Only, tmpResultArray, repeatRegionInfoVec, 
				readInfo, indexInfo, tmpOpenMP, threadNO, statsInfo, fasta_or_fastq_bool, overlapLength_top2candiAlignment);

			fixPhase1Info.memoryFree();
			peAlignInfo.memoryFree();
		}
		nowtime = time(NULL);
		local = localtime(&nowtime);
		mapping_log_ofs << endl << "*******************************************************************************"
			<< endl << "[" << asctime(local) << "... mapping a batch of reads ends ......" << endl << endl << endl; 		
		resultArrayQueue->pushBack2ResultArrayQueue(tmpResultArray);
		readArrayQueue->popFromReadQueue();
	}
	endOfProcessing_bool = true;
}


#endif