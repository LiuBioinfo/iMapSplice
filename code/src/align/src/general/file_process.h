// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FILE_PROCESS_H
#define FILE_PROCESS_H

#include <string>
#include <string.h>
#include "readSeqPreProcessing.h"
//#include "splice_info.h"

using namespace std;

class FileProcessInfo
{
public:
	FileProcessInfo()
	{}

	void output_fixHeadTail(int realRecordNum,
		vector<string>& peAlignSamVec_complete_pair,
		vector<string>& peAlignSamVec_incomplete_pair,
		vector<string>& peAlignSamVec_complete_unpair,
		vector<string>& peAlignSamVec_incomplete_unpair,
		vector<string>& peAlignSamVec_complete_pair_alignInfo,
		vector<string>& peAlignSamVec_incomplete_pair_alignInfo,
		ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
		ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
		ofstream& OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs,
		bool outputAlignInfoAndSamForAllPairedAlignmentBool)
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

				if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				{
					if(peAlignSamVec_complete_pair_alignInfo[tmp] != "")
						OutputSamFile_fixHeadTail_complete_pair_alignInfo_ofs << peAlignSamVec_complete_pair_alignInfo[tmp] << endl;
					if(peAlignSamVec_incomplete_pair_alignInfo[tmp] != "")
						OutputSamFile_fixHeadTail_incomplete_pair_alignInfo_ofs << peAlignSamVec_incomplete_pair_alignInfo[tmp] << endl;					
				}
			}		
	}

	void output_phase1(int realRecordNum,
		vector<string>& PeAlignSamStrVec_complete,
		vector<string>& PeAlignInfoStrVec_inCompletePair,
		vector<string>& PeAlignInfoStrVec_oneEndUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped,
		vector<string>& PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion,
		vector<string>& PeAlignSamStrVec_inCompletePair,
		vector<string>& PeAlignInfoStrVec_completePaired,
		ofstream& tmpAlignCompleteRead_ofs,
		ofstream& tmpAlignIncompletePair_ofs,
		ofstream& tmpAlignOneEndUnmapped_ofs,
		ofstream& tmpAlignBothEndsUnmapped_ofs,
		ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
		ofstream& tmpAlignIncompletePair_SAM_ofs,
		ofstream& tmpAlignCompleteRead_alignInfo_ofs,
		bool outputAlignInfoAndSamForAllPairedAlignmentBool)
	{
			for(int tmp = 0; tmp < realRecordNum; tmp++)
			{	
				if(PeAlignSamStrVec_complete[tmp] != "")
				{
					tmpAlignCompleteRead_ofs << PeAlignSamStrVec_complete[tmp] << endl;
				}			

				if(PeAlignInfoStrVec_inCompletePair[tmp] != "")
				{
					tmpAlignIncompletePair_ofs << PeAlignInfoStrVec_inCompletePair[tmp] << endl;
				}
				
				if(PeAlignInfoStrVec_oneEndUnmapped[tmp] != "")
				{
					tmpAlignOneEndUnmapped_ofs << PeAlignInfoStrVec_oneEndUnmapped[tmp] << endl;
				}
				
				if(PeAlignSamStrVec_bothEndsUnmapped[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_ofs << PeAlignSamStrVec_bothEndsUnmapped[tmp] << endl;
				}
				if(PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs << PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] << endl;
				}

				if(PeAlignSamStrVec_inCompletePair[tmp] != "")
				{
					tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;
				}

				if(outputAlignInfoAndSamForAllPairedAlignmentBool)
				{
					if(PeAlignInfoStrVec_completePaired[tmp] != "")
					{
						tmpAlignCompleteRead_alignInfo_ofs << PeAlignInfoStrVec_completePaired[tmp] << endl;
					}
				}

			}			
	}

	void generateMultiThreadOutputFile(const string& tmpAlignCompleteRead,
		const string& tmpAlignCompleteRead_alignInfo,
		const string& tmpAlignOneEndUnmapped,
		const string& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile,
		const string& tmpAlignBothEndsUnmapped,
		const string& tmpAlignIncompletePair,
		const string& tmpAlignIncompletePair_SAM,
		//const string& tmpAlignCompleteRead_multi,
		//const string& tmpAlignCompleteRead_unique,
		vector< ofstream* >& tmpAlignCompleteRead_ofs_vec,
		vector< ofstream* >& tmpAlignOneEndUnmapped_ofs_vec,
		vector< ofstream* >& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs_vec,
		vector< ofstream* >& tmpAlignBothEndsUnmapped_ofs_vec,
		vector< ofstream* >& tmpAlignIncompletePair_ofs_vec,
		vector< ofstream* >& tmpAlignIncompletePair_SAM_ofs_vec,
		vector< ofstream* >& tmpAlignCompleteRead_alignInfo_ofs_vec,
		//#ifdef CHECK_MULTI
		//vector< ofstream* >& tmpAlignCompleteRead_ofs_vec_multi,
		//vector< ofstream* >& tmpAlignCompleteRead_ofs_vec_unique, 		
		int threads_num)
	{
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			string tmp_tmpAlignCompleteRead = tmpAlignCompleteRead + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignCompleteRead_ofs = new ofstream(tmp_tmpAlignCompleteRead.c_str());
			tmpAlignCompleteRead_ofs_vec.push_back(tmp_tmpAlignCompleteRead_ofs);

			string tmp_tmpAlignOneEndUnmapped = tmpAlignOneEndUnmapped + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignOneEndUnmapped_ofs = new ofstream(tmp_tmpAlignOneEndUnmapped.c_str());
			tmpAlignOneEndUnmapped_ofs_vec.push_back(tmp_tmpAlignOneEndUnmapped_ofs);

			string tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile = tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs = new ofstream(tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile.c_str());
			tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs_vec.push_back(tmp_tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs);

			string tmp_tmpAlignBothEndsUnmapped = tmpAlignBothEndsUnmapped + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignBothEndsUnmapped_ofs = new ofstream(tmp_tmpAlignBothEndsUnmapped.c_str());
			tmpAlignBothEndsUnmapped_ofs_vec.push_back(tmp_tmpAlignBothEndsUnmapped_ofs);

			string tmp_tmpAlignIncompletePair = tmpAlignIncompletePair + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignIncompletePair_ofs = new ofstream(tmp_tmpAlignIncompletePair.c_str());
			tmpAlignIncompletePair_ofs_vec.push_back(tmp_tmpAlignIncompletePair_ofs);

			string tmp_tmpAlignIncompletePair_SAM = tmpAlignIncompletePair_SAM + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignIncompletePair_SAM_ofs = new ofstream(tmp_tmpAlignIncompletePair_SAM.c_str());
			tmpAlignIncompletePair_SAM_ofs_vec.push_back(tmp_tmpAlignIncompletePair_SAM_ofs);

			string tmp_tmpAlignCompleteRead_alignInfo = tmpAlignCompleteRead_alignInfo + "." + int_to_str(tmp);
			ofstream* tmp_tmpAlignCompleteRead_alignInfo_ofs = new ofstream(tmp_tmpAlignCompleteRead_alignInfo.c_str());
			tmpAlignCompleteRead_alignInfo_ofs_vec.push_back(tmp_tmpAlignCompleteRead_alignInfo_ofs);
		}
	}

	bool inputRead_phase1(
		ifstream& inputRead_ifs, ifstream& inputRead_PE_ifs,
		int& realRecordNum, 
		//int recordNumTmp,
		int recordNum,
		vector<string>& readName1Vec,
		vector<string>& readSeq1Vec,
		vector<string>& readQualSeq1Vec,
		vector<string>& readName2Vec,
		vector<string>& readSeq2Vec,
		vector<string>& readQualSeq2Vec,
		int& readLengthMax_tmp,
		bool InputAsFastq,
		InputReadPreProcess* readPreProcessInfo)
	{
		for(int recordNumTmp = 0; recordNumTmp < recordNum; recordNumTmp++)
		{		
	    	if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
	    	{
				realRecordNum = recordNumTmp;
				//EndOfRecord = true;
				//break;    			
	    		return true;
	    	}

	    	string line1;
	    	getline(inputRead_ifs, line1); // readName_1

	    	if((inputRead_ifs.eof())||(inputRead_PE_ifs.eof()))
	    	{
				realRecordNum = recordNumTmp;
				//EndOfRecord = true;
				//break;    			
	    		return true;
	    	}

	       	readName1Vec[recordNumTmp] = line1.substr(1);
	    	string line2;
	    	getline(inputRead_ifs, line2); // readSeq_1
	    	string line2_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2);		
	    	readSeq1Vec[recordNumTmp] = line2_afterProcess;    	
		
			int readLength_1 = line2.length();
			if(readLength_1 > readLengthMax_tmp)
				readLengthMax_tmp = readLength_1;		
	   		string line3, line4;
	    	if(InputAsFastq)
	    	{
	    		getline(inputRead_ifs, line3);
	   			getline(inputRead_ifs, line4);
	    		readQualSeq1Vec[recordNumTmp] = line4;
	   		}
			string line1_PE, line2_PE;   
	    	getline(inputRead_PE_ifs, line1_PE); // readName_2
	    	readName2Vec[recordNumTmp] = line1_PE.substr(1);
	   		getline(inputRead_PE_ifs, line2_PE);

			string line2_PE_afterProcess = readPreProcessInfo->upperCaseReadSeq(line2_PE);
			readSeq2Vec[recordNumTmp] = line2_PE_afterProcess;   		
	    		
	    	int readLength_2 = line2_PE.length();
			if(readLength_2 > readLengthMax_tmp)
	    		readLengthMax_tmp = readLength_2;
	    	string line3_PE, line4_PE;	
	    	if(InputAsFastq)
	    	{
	    		getline(inputRead_PE_ifs, line3_PE);
	    		getline(inputRead_PE_ifs, line4_PE);
	    		readQualSeq2Vec[recordNumTmp] = line4_PE;
	    	}
	    }
	    return false;
	}

	void mergePhase1OutputFiles(int threads_num, const string& tmpAlignCompleteRead,
		const string& tmpAlignOneEndUnmapped, const string& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile,
		const string& tmpAlignBothEndsUnmapped, const string& tmpAlignIncompletePair_SAM,
		const string& tmpAlignCompleteRead_alignInfo, //const string& tmpAlignCompleteRead_multi,
		//const string& tmpAlignCompleteRead_unique, 
		ofstream& log_ofs)
	{
		time_t nowtime_start;
		nowtime_start = time(NULL);
		struct tm *local_start;
		local_start = localtime(&nowtime_start);
		log_ofs << endl << "[" << asctime(local_start) << "start to merge phase1 output files ..." << endl;
		string cat_cmd_completeRead = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_completeRead = cat_cmd_completeRead + " " + tmpAlignCompleteRead + "." + int_to_str(tmp); 
		} 
		cat_cmd_completeRead = cat_cmd_completeRead + " > " + tmpAlignCompleteRead;
		system(cat_cmd_completeRead.c_str());

		string cat_cmd_oneEndUnmapped = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_oneEndUnmapped = cat_cmd_oneEndUnmapped + " " + tmpAlignOneEndUnmapped + "." + int_to_str(tmp);
		}
		cat_cmd_oneEndUnmapped = cat_cmd_oneEndUnmapped + " > " + tmpAlignOneEndUnmapped;
		system(cat_cmd_oneEndUnmapped.c_str());

		string cat_cmd_mappedToRepeatRegion = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_mappedToRepeatRegion = cat_cmd_mappedToRepeatRegion + " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile + "." + int_to_str(tmp);
		}
		cat_cmd_mappedToRepeatRegion = cat_cmd_mappedToRepeatRegion + " > " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile;
		system(cat_cmd_mappedToRepeatRegion.c_str());

		string cat_cmd_bothEndsUnmapped = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_bothEndsUnmapped = cat_cmd_bothEndsUnmapped + " " + tmpAlignBothEndsUnmapped + "." + int_to_str(tmp);
		}
		cat_cmd_bothEndsUnmapped = cat_cmd_bothEndsUnmapped + " > " + tmpAlignBothEndsUnmapped; 
		system(cat_cmd_bothEndsUnmapped.c_str());
		/*
		string cat_cmd_tmpAlignIncompletePair = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " " + tmpAlignIncompletePair + "." + int_to_str(tmp);
		}
		cat_cmd_tmpAlignIncompletePair = cat_cmd_tmpAlignIncompletePair + " > " + tmpAlignIncompletePair;
		system(cat_cmd_tmpAlignIncompletePair.c_str());
		*/
		string cat_cmd_tmpAlignIncompletePair_SAM = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_tmpAlignIncompletePair_SAM = cat_cmd_tmpAlignIncompletePair_SAM + " " + tmpAlignIncompletePair_SAM + "." + int_to_str(tmp);
		}
		cat_cmd_tmpAlignIncompletePair_SAM = cat_cmd_tmpAlignIncompletePair_SAM + " > " + tmpAlignIncompletePair_SAM;
		system(cat_cmd_tmpAlignIncompletePair_SAM.c_str());

		string cat_cmd_tmpAlignCompleteRead_alignInfo = "cat";
		for(int tmp = 1; tmp <= threads_num; tmp++)
		{
			cat_cmd_tmpAlignCompleteRead_alignInfo = cat_cmd_tmpAlignCompleteRead_alignInfo + " " + tmpAlignCompleteRead_alignInfo + "." + int_to_str(tmp);
		}
		cat_cmd_tmpAlignCompleteRead_alignInfo = cat_cmd_tmpAlignCompleteRead_alignInfo + " > " + tmpAlignCompleteRead_alignInfo;
		system(cat_cmd_tmpAlignCompleteRead_alignInfo.c_str());

		time_t nowtime;
		nowtime = time(NULL);
		struct tm *local;
		local = localtime(&nowtime);
		log_ofs << endl << "[" << asctime(local) << "finish merging files ...." << endl;		
	}

	void mergeAllFinalFiles(bool Do_Phase1_Only, const string& tmpHeadSectionInfo, 
		const string& tmpAlignCompleteRead, const string& tmpAlignIncompletePair_SAM,
		const string& tmpAlignOneEndUnmapped, const string& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile,
		const string& tmpAlignBothEndsUnmapped, const string& OutputSamFile_oneEndMapped,
		const string& OutputSamFile_fixHeadTail_complete_pair, const string& OutputSamFile_fixHeadTail_incomplete_pair,
		const string& OutputSamFile_fixHeadTail_complete_unpair, const string& OutputSamFile_fixHeadTail_incomplete_unpair,
		const string& OutputSamFile_oneEndMapped_unpairComplete,  
		const string& outputDirStr
		)
	{
		string finalOutputSam = outputDirStr + "/output.sam";
		if(Do_Phase1_Only)
		{
			string cat_cmd = "cat "
				+ tmpHeadSectionInfo
				+ " " + tmpAlignCompleteRead  
				//+ " " + tmpAlignOneEndUnmapped
				+ " " + tmpAlignIncompletePair_SAM 
				+ " " + tmpAlignOneEndUnmapped
				+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
				+ " " + tmpAlignBothEndsUnmapped
				+ " > " + finalOutputSam;
			system(cat_cmd.c_str()); 
		}
		else
		{
			string cat_cmd = "cat " 
				+ tmpHeadSectionInfo
				+ " " + tmpAlignCompleteRead 
				+ " " + OutputSamFile_oneEndMapped 
				+ " " + OutputSamFile_fixHeadTail_complete_pair
				+ " " + OutputSamFile_fixHeadTail_incomplete_pair
				+ " " + OutputSamFile_fixHeadTail_complete_unpair
				+ " " + OutputSamFile_fixHeadTail_incomplete_unpair
				+ " " + OutputSamFile_oneEndMapped_unpairComplete
				+ " " + tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile
				+ " " + tmpAlignBothEndsUnmapped 
				+ " > " + finalOutputSam;
			system(cat_cmd.c_str()); 		
		}
	}

	void removeAllTemporaryFiles(const string& outputDirStr)
	{
		string phase1_dir_path = outputDirStr + "/phase1_output";
		string phase2_dir_path = outputDirStr + "/phase2_output";
		string remove_phase1_output_cmd = "rm -r " + phase1_dir_path;
		string remove_phase2_output_cmd = "rm -r " + phase2_dir_path; 

		system(remove_phase1_output_cmd.c_str());
		system(remove_phase2_output_cmd.c_str());
	}
};

#endif