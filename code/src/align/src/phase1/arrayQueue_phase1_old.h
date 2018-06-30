// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ARRAYQUEUE_PHASE1_H
#define ARRAYQUEUE_PHASE1_H
	
#include <stdio.h>
#include <stdlib.h>
#include <string>

//#define TotalReadNumInReadArray 500

using namespace std;

class Read_Array
{
private:
	vector<string> inputReadName_1;
	vector<string> inputReadSeq_1;
	vector<string> inputReadQualSeq_1;

	vector<string> inputReadName_2;
	vector<string> inputReadSeq_2;
	vector<string> inputReadQualSeq_2;
public:
	Read_Array* nextReadNode;
	Read_Array()
	{}

	void pushBack2ReadArray(const string& tmpReadName_1,
		const string& tmpReadName_2,
		const string& tmpReadSeq_1,
		const string& tmpReadSeq_2)
	{
		inputReadName_1.push_back(tmpReadName_1);
		inputReadSeq_1.push_back(tmpReadSeq_1);

		inputReadName_2.push_back(tmpReadName_2);
		inputReadSeq_2.push_back(tmpReadSeq_2);		
	}

	void pushBack2ReadArray_fq(const string& tmpReadName_1,
		const string& tmpReadName_2,
		const string& tmpReadSeq_1,
		const string& tmpReadSeq_2,
		const string& tmpReadQualSeq_1,
		const string& tmpReadQualSeq_2)
	{
		inputReadName_1.push_back(tmpReadName_1);
		inputReadSeq_1.push_back(tmpReadSeq_1);
		inputReadQualSeq_1.push_back(tmpReadQualSeq_1);

		inputReadName_2.push_back(tmpReadName_2);
		inputReadSeq_2.push_back(tmpReadSeq_2);		
		inputReadQualSeq_2.push_back(tmpReadQualSeq_2);
	}

	string returnReadSeq_1(int index)
	{
		return inputReadSeq_1[index];
	}

	string returnReadSeq_2(int index)
	{
		return inputReadSeq_2[index];
	}

	string returnReadQualSeq_1(int index)
	{
		return inputReadQualSeq_1[index];
	}

	string returnReadQualSeq_2(int index)
	{
		return inputReadQualSeq_2[index];
	}

	string returnReadName_1(int index)
	{
		return inputReadName_1[index];
	}

	string returnReadName_2(int index)
	{
		return inputReadName_2[index];
	}

	int returnSize()
	{
		return inputReadName_1.size();
	}
};

class Read_Array_Queue
{
private:
public:
	Read_Array* frontNode;
	Read_Array* tailNode;
	int queueSize;
	/*Read_Array_Queue()
	{
		queueSize = 0;
	}*/
	
	Read_Array_Queue()
	{
		queueSize = 0;
		Read_Array* p = new Read_Array();
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}

		p->nextReadNode = NULL;
		frontNode = p;
		tailNode = p;
	}

	void initiate()
	{
		queueSize = 0;
		Read_Array* p = new Read_Array();
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}

		p->nextReadNode = NULL;
		frontNode = p;
		tailNode = p;
	}

	void initiateWith1stRead(const string& tmpReadName_1,
		const string& tmpReadName_2,
		const string& tmpReadSeq_1,
		const string& tmpReadSeq_2, ofstream& input_log_ofs)
	{
		queueSize ++;
		Read_Array* p = new Read_Array();
		p->nextReadNode = NULL;
		tailNode->nextReadNode = p;
		tailNode = p;			
		tailNode->pushBack2ReadArray(tmpReadName_1,
			tmpReadName_2, tmpReadSeq_1, tmpReadSeq_2);		
		
		time_t nowtime;
		nowtime = time(NULL);
		struct tm *local;
		local = localtime(&nowtime);		
		input_log_ofs << endl << "[" << asctime(local) 
			<< "... input 2 new array starts ......" << endl;  
		input_log_ofs << "tmpBatchIndex: 1" << endl << endl;
	}

	void initiateWith1stRead_fq(const string& tmpReadName_1,
		const string& tmpReadName_2,
		const string& tmpReadSeq_1,
		const string& tmpReadSeq_2,
		const string& tmpReadQualSeq_1,
		const string& tmpReadQualSeq_2, ofstream& input_log_ofs)
	{
		queueSize ++;
		Read_Array* p = new Read_Array();
		p->nextReadNode = NULL;
		tailNode->nextReadNode = p;
		tailNode = p;			
		tailNode->pushBack2ReadArray_fq(tmpReadName_1,
			tmpReadName_2, tmpReadSeq_1, tmpReadSeq_2,
			tmpReadQualSeq_1, tmpReadQualSeq_2);		
		
		time_t nowtime;
		nowtime = time(NULL);
		struct tm *local;
		local = localtime(&nowtime);
		input_log_ofs << endl << "[" << asctime(local) 
			<< "... input 2 new array starts ......" << endl;  
		input_log_ofs << "tmpBatchIndex: 1" << endl << endl;
	}

	bool atLeast3Node()
	{
		if(frontNode == tailNode)
			return false;
		else if((frontNode->nextReadNode) == tailNode)
			return false;
		else
			return true;
	}

	bool only2Node()
	{
		if((frontNode->nextReadNode) == tailNode)
			return true;
		else 
			return false;
	}

	bool readQueueEmpty()
	{
		return (frontNode == tailNode);
	}

	void popFromReadQueue()
	{
		if(frontNode == tailNode)
		{
			//cout << "The queue is empty" << endl;
			exit(1);
		}
		else
		{
			Read_Array* p = frontNode->nextReadNode;
			frontNode->nextReadNode = p->nextReadNode;
			if(tailNode == p) // when only one node exits
				tailNode = frontNode;
			delete p; p = NULL;
		}
	}

	void getSeqFromInputFile(
		const string& tmpReadName_1,
		const string& tmpReadName_2,
		const string& tmpReadSeq_1,
		const string& tmpReadSeq_2,
		int totalReadNumInReadArray, ofstream& input_log_ofs, int& tmpBatchIndex)
	{
		if(tailNode->returnSize() < totalReadNumInReadArray)
		{
			tailNode->pushBack2ReadArray(
				tmpReadName_1, tmpReadName_2,
				tmpReadSeq_1, tmpReadSeq_2);
		}
		else
		{		
			time_t nowtime;
			nowtime = time(NULL);
			struct tm *local;
			local = localtime(&nowtime);
			tmpBatchIndex ++;
			input_log_ofs << endl << "[" << asctime(local) 
				<< "... input 2 new array starts ......" << endl;  
			input_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl << endl;
			//cout << "add a new readNode" << endl;
			Read_Array* p = new Read_Array();
			p->nextReadNode = NULL;
			tailNode->nextReadNode = p;
			tailNode = p;			
			tailNode->pushBack2ReadArray(
				tmpReadName_1, tmpReadName_2,
				tmpReadSeq_1, tmpReadSeq_2);
		}
	}

	void getSeqFromInputFile_fq(
		const string& tmpReadName_1,
		const string& tmpReadName_2,
		const string& tmpReadSeq_1,
		const string& tmpReadSeq_2,
		const string& tmpReadQualSeq_1,
		const string& tmpReadQualSeq_2,
		int totalReadNumInReadArray, ofstream& input_log_ofs, int& tmpBatchIndex)
	{
		if(tailNode->returnSize() < totalReadNumInReadArray)
		{
			tailNode->pushBack2ReadArray_fq(
				tmpReadName_1, tmpReadName_2,
				tmpReadSeq_1, tmpReadSeq_2,
				tmpReadQualSeq_1, tmpReadQualSeq_2);
		}
		else
		{		
			time_t nowtime;
			nowtime = time(NULL);
			struct tm *local;
			local = localtime(&nowtime);
			tmpBatchIndex ++;
			input_log_ofs << endl << "[" << asctime(local) 
				<< "... input 2 new array starts ......" << endl;  
			input_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl << endl;
			//cout << "add a new readNode" << endl;
			Read_Array* p = new Read_Array();
			p->nextReadNode = NULL;
			tailNode->nextReadNode = p;
			tailNode = p;			
			tailNode->pushBack2ReadArray_fq(
				tmpReadName_1, tmpReadName_2,
				tmpReadSeq_1, tmpReadSeq_2,
				tmpReadQualSeq_1, tmpReadQualSeq_2);
		}
	}

	string returnFrontNodeReadQualSeq_1(int index)
	{
		return ((frontNode->nextReadNode)->returnReadQualSeq_1(index));
	}
	string returnFrontNodeReadQualSeq_2(int index)
	{
		return ((frontNode->nextReadNode)->returnReadQualSeq_2(index));
	}

	string returnFrontNodeReadSeq_1(int index)
	{
		return ((frontNode->nextReadNode)->returnReadSeq_1(index));
	}
	string returnFrontNodeReadSeq_2(int index)
	{
		return ((frontNode->nextReadNode)->returnReadSeq_2(index));
	}

	string returnFrontNodeReadName_1(int index)
	{
		return ((frontNode->nextReadNode)->returnReadName_1(index));
	}
	string returnFrontNodeReadName_2(int index)
	{
		return ((frontNode->nextReadNode)->returnReadName_2(index));
	}
	int returnFrontNodeSize()
	{
		return ((frontNode->nextReadNode)->returnSize());
	}
};

class Result_Array
{
private:
	vector<string> PeAlignSamStrVec_complete;
	vector<string> PeAlignInfoStrVec_inCompletePair;
	vector<string> PeAlignInfoStrVec_oneEndUnmapped;
	vector<string> PeAlignSamStrVec_bothEndsUnmapped;
	vector<string> PeAlignSamStrVec_bothEndsUnmapped_lowScore;
	vector<string> PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion;
	vector<string> PeAlignSamStrVec_inCompletePair;
	vector<string> PeAlignInfoStrVec_completePaired;	
public:
	Result_Array* nextResultNode;
	Result_Array(int size)
	{
		for(int tmp = 0; tmp < size; tmp++)
		{
			PeAlignSamStrVec_complete.push_back("");
			PeAlignInfoStrVec_inCompletePair.push_back("");
			PeAlignInfoStrVec_oneEndUnmapped.push_back("");
			PeAlignSamStrVec_bothEndsUnmapped.push_back("");
			PeAlignSamStrVec_bothEndsUnmapped_lowScore.push_back("");
			PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion.push_back("");
			PeAlignSamStrVec_inCompletePair.push_back("");
			PeAlignInfoStrVec_completePaired.push_back("");	
		}
	}

	void insert_PeAlignSamStrVec_complete(const string& tmpStr, int index)
	{
		PeAlignSamStrVec_complete[index] = tmpStr;
	}
	void insert_PeAlignInfoStrVec_inCompletePair(const string& tmpStr, int index)
	{
		PeAlignInfoStrVec_inCompletePair[index] = tmpStr;
	}
	void insert_PeAlignInfoStrVec_oneEndUnmapped(const string& tmpStr, int index)
	{
		PeAlignInfoStrVec_oneEndUnmapped[index] = tmpStr;
	}
	void insert_PeAlignSamStrVec_bothEndsUnmapped(const string& tmpStr, int index)
	{
		PeAlignSamStrVec_bothEndsUnmapped[index] = tmpStr;
	}
	void insert_PeAlignSamStrVec_bothEndsUnmapped_lowScore(const string& tmpStr, int index)
	{
		PeAlignSamStrVec_bothEndsUnmapped_lowScore[index] = tmpStr;
	}
	void insert_PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion(const string& tmpStr, int index)
	{
		PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[index] = tmpStr;
	}
	void insert_PeAlignSamStrVec_inCompletePair(const string& tmpStr, int index)
	{
		PeAlignSamStrVec_inCompletePair[index] = tmpStr;
	}
	void insert_PeAlignInfoStrVec_completePaired(const string& tmpStr, int index)
	{
		PeAlignInfoStrVec_completePaired[index] = tmpStr;
	}

	int returnSize()
	{
		return PeAlignSamStrVec_complete.size();
	}

	void outputResultArray(
		ofstream& tmpAlignCompleteRead_ofs,
		ofstream& tmpAlignIncompletePair_ofs,
		ofstream& tmpAlignOneEndUnmapped_ofs,
		ofstream& tmpAlignBothEndsUnmapped_ofs,
		ofstream& tmpAlignBothEndsUnmapped_lowScore_ofs,
		ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
		ofstream& tmpAlignIncompletePair_SAM_ofs)
	{
			for(int tmp = 0; tmp < PeAlignSamStrVec_complete.size(); tmp++)
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
				if(PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_lowScore_ofs << PeAlignSamStrVec_bothEndsUnmapped_lowScore[tmp] << endl;
				}
				if(PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] != "")
				{
					tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs << PeAlignSamStrVec_bothEndsUnmapped_mappedToRepeatRegion[tmp] << endl;
				}
				if(PeAlignSamStrVec_inCompletePair[tmp] != "")
				{
					tmpAlignIncompletePair_SAM_ofs << PeAlignSamStrVec_inCompletePair[tmp] << endl;
				}
			}		
	}
};

class Result_Array_Queue
{
private:

public:
	Result_Array* frontNode;
	Result_Array* tailNode;
	Result_Array_Queue()
	{
		Result_Array* p = new Result_Array(0);
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}

		p->nextResultNode = NULL;
		frontNode = p;
		tailNode = p;
	}	

	void pushBack2ResultArrayQueue(Result_Array* tmpResultArray)
	{
		tmpResultArray->nextResultNode = NULL;
		tailNode->nextResultNode = tmpResultArray;		
		tailNode = tmpResultArray;
	}

	void popFromResultQueue()
	{
		if(frontNode == tailNode)
		{
			//cout << "The queue is empty" << endl;
			exit(1);
		}
		else
		{
			Result_Array* p = frontNode->nextResultNode;
			frontNode->nextResultNode = p->nextResultNode;
			if(tailNode == p) // when only one node exits
				tailNode = frontNode;
			delete p; p = NULL;
		}
	}

	bool atLeast3Node()
	{
		if(frontNode == tailNode)
			return false;
		else if((frontNode->nextResultNode) == tailNode)
			return false;
		else
			return true;
	}

	bool only2Node()
	{
		if((frontNode->nextResultNode) == tailNode)
			return true;
		else 
			return false;
	}

	bool resultQueueEmpty()
	{
		return (frontNode == tailNode);
	}

	void outputFrontResultArray(
		ofstream& tmpAlignCompleteRead_ofs,
		ofstream& tmpAlignIncompletePair_ofs,
		ofstream& tmpAlignOneEndUnmapped_ofs,
		ofstream& tmpAlignBothEndsUnmapped_ofs,
		ofstream& tmpAlignBothEndsUnmapped_lowScore_ofs,
		ofstream& tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
		ofstream& tmpAlignIncompletePair_SAM_ofs)
	{
		(frontNode->nextResultNode)->outputResultArray(
			tmpAlignCompleteRead_ofs,
			tmpAlignIncompletePair_ofs,
			tmpAlignOneEndUnmapped_ofs,
			tmpAlignBothEndsUnmapped_ofs,
			tmpAlignBothEndsUnmapped_lowScore_ofs,
			tmpAlignBothEndsUnmapped_mappedToRepeatRegionFile_ofs,
			tmpAlignIncompletePair_SAM_ofs);
	}
};

#endif