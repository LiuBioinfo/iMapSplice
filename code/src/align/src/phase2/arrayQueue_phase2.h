// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef ARRAYQUEUE_PHASE2_H
#define ARRAYQUEUE_PHASE2_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

class AlignInfoInput_Array
{
private:
	vector<string> inputReadName_alignNum_1;
	vector<string> inputReadSeq_1;
	vector<string> inputReadQualSeq_1;
	vector<string> inputReadName_alignNum_2;
	vector<string> inputReadSeq_2;
	vector<string> inputReadQualSeq_2;
	vector<string> inputAlignInfo_Nor1;
	vector<string> inputAlignInfo_Rcm1;
	vector<string> inputAlignInfo_Nor2;
	vector<string> inputAlignInfo_Rcm2;	
public:
	AlignInfoInput_Array* nextAlignInfoNode;
	AlignInfoInput_Array()
	{}

	void pushBack2AlignInfoInputArray(
		const string& tmpReadName_alignNum_1,
		const string& tmpReadSeq_1,
		const string& tmpReadQualSeq_1,
		const string& tmpReadName_alignNum_2,
		const string& tmpReadSeq_2,
		const string& tmpReadQualSeq_2,
		const string& tmpAlignInfo_Nor1,
		const string& tmpAlignInfo_Rcm1,
		const string& tmpAlignInfo_Nor2,
		const string& tmpAlignInfo_Rcm2
		)
	{
		inputReadName_alignNum_1.push_back(tmpReadName_alignNum_1);
		inputReadSeq_1.push_back(tmpReadSeq_1);
		inputReadQualSeq_1.push_back(tmpReadQualSeq_1);
		inputReadName_alignNum_2.push_back(tmpReadName_alignNum_2);
		inputReadSeq_2.push_back(tmpReadSeq_2);
		inputReadQualSeq_2.push_back(tmpReadQualSeq_2);
		inputAlignInfo_Nor1.push_back(tmpAlignInfo_Nor1);
		inputAlignInfo_Rcm1.push_back(tmpAlignInfo_Rcm1);
		inputAlignInfo_Nor2.push_back(tmpAlignInfo_Nor2);
		inputAlignInfo_Rcm2.push_back(tmpAlignInfo_Rcm2);		
	}

	string returnAlignInfo_Nor1(int index)
	{
		return inputAlignInfo_Nor1[index];
	}
	string returnAlignInfo_Rcm1(int index)
	{
		return inputAlignInfo_Rcm1[index];
	}	
	string returnAlignInfo_Nor2(int index)
	{
		return inputAlignInfo_Nor2[index];
	}
	string returnAlignInfo_Rcm2(int index)
	{
		return inputAlignInfo_Rcm2[index];
	}	
	string returnReadName_alignNum_1(int index)
	{
		return inputReadName_alignNum_1[index];
	}
	string returnReadNameAlignNum_1(int index)
	{
		return inputReadName_alignNum_1[index];
	}
	string returnReadSeq_1(int index)
	{
		return inputReadSeq_1[index];
	}
	string returnReadQualSeq_1(int index)
	{
		return inputReadQualSeq_1[index];
	}

	string returnReadName_alignNum_2(int index)
	{
		return inputReadName_alignNum_2[index];
	}
	string returnReadNameAlignNum_2(int index)
	{
		return inputReadName_alignNum_2[index];
	}
	string returnReadSeq_2(int index)
	{
		return inputReadSeq_2[index];
	}
	string returnReadQualSeq_2(int index)
	{
		return inputReadQualSeq_2[index];
	}
	int returnSize()
	{
		return inputReadName_alignNum_1.size();
	}
};

class AlignInfoInput_Array_Queue
{
private:

public:
	AlignInfoInput_Array* frontNode;
	AlignInfoInput_Array* tailNode;
	int queueSize;

	void free()
	{
		frontNode = NULL;
		tailNode = NULL;
	}	

	AlignInfoInput_Array_Queue()
	{
		queueSize = 0;
		AlignInfoInput_Array* p = new AlignInfoInput_Array();
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}
		p->nextAlignInfoNode = NULL;
		frontNode = p;
		tailNode = p;		
	}

	void initiate()
	{
		queueSize = 0;
		AlignInfoInput_Array* p = new AlignInfoInput_Array();
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}

		p->nextAlignInfoNode = NULL;
		frontNode = p;
		tailNode = p;
	}	

	void initiateWith1stAlignInfo(
		const string& tmpReadName_alignNum_1,
		const string& tmpReadSeq_1,
		const string& tmpReadQualSeq_1,
		const string& tmpReadName_alignNum_2,
		const string& tmpReadSeq_2,
		const string& tmpReadQualSeq_2,
		const string& tmpAlignInfo_Nor1,
		const string& tmpAlignInfo_Rcm1,
		const string& tmpAlignInfo_Nor2,
		const string& tmpAlignInfo_Rcm2,
		ofstream& input_log_ofs)
	{
		queueSize ++;
		AlignInfoInput_Array* p = new AlignInfoInput_Array();
		p->nextAlignInfoNode = NULL;
		tailNode->nextAlignInfoNode = p;
		tailNode = p;			
		tailNode->pushBack2AlignInfoInputArray(
			tmpReadName_alignNum_1,
			tmpReadSeq_1,
			tmpReadQualSeq_1,
			tmpReadName_alignNum_2,
			tmpReadSeq_2,
			tmpReadQualSeq_2,
			tmpAlignInfo_Nor1,
			tmpAlignInfo_Rcm1,
			tmpAlignInfo_Nor2,
			tmpAlignInfo_Rcm2);

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
		else if((frontNode->nextAlignInfoNode) == tailNode)
			return false;
		else
			return true;
	}

	bool only2Node()
	{
		if((frontNode->nextAlignInfoNode) == tailNode)
			return true;
		else 
			return false;
	}

	bool inputQueueEmpty()
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
			AlignInfoInput_Array* p = frontNode->nextAlignInfoNode;
			frontNode->nextAlignInfoNode = p->nextAlignInfoNode;
			if(tailNode == p) // when only one node exits
				tailNode = frontNode;
			delete p; p = NULL;
		}
	}	

	void getAlignInfoFromInputFile(
		const string& tmpReadName_alignNum_1,
		const string& tmpReadSeq_1,
		const string& tmpReadQualSeq_1,
		const string& tmpReadName_alignNum_2,
		const string& tmpReadSeq_2,
		const string& tmpReadQualSeq_2,
		const string& tmpAlignInfo_Nor1,
		const string& tmpAlignInfo_Rcm1,
		const string& tmpAlignInfo_Nor2,
		const string& tmpAlignInfo_Rcm2,
		int totalAlignInfoNumInAlignInfoInputArray,
		ofstream& input_log_ofs, int& tmpBatchIndex
		)
	{
		if(tailNode->returnSize() < totalAlignInfoNumInAlignInfoInputArray)
		{
			tailNode->pushBack2AlignInfoInputArray(
				tmpReadName_alignNum_1,
				tmpReadSeq_1,
				tmpReadQualSeq_1,
				tmpReadName_alignNum_2,
				tmpReadSeq_2,
				tmpReadQualSeq_2,
				tmpAlignInfo_Nor1,
				tmpAlignInfo_Rcm1,
				tmpAlignInfo_Nor2,
				tmpAlignInfo_Rcm2);
		}	
		else
		{
			time_t nowtime;
			nowtime = time(NULL);
			struct tm *local;
			local = localtime(&nowtime);
			tmpBatchIndex ++;
			input_log_ofs << endl << "[" << asctime(local) 
				<< "... input 2 new alignInfoInputArray starts ......" << endl;  			
			input_log_ofs << "tmpBatchIndex: " << tmpBatchIndex << endl << endl;
			AlignInfoInput_Array* p = new AlignInfoInput_Array();
			p->nextAlignInfoNode = NULL;
			tailNode->nextAlignInfoNode = p;
			tailNode = p;
			tailNode->pushBack2AlignInfoInputArray(
				tmpReadName_alignNum_1,
				tmpReadSeq_1,
				tmpReadQualSeq_1,
				tmpReadName_alignNum_2,
				tmpReadSeq_2,
				tmpReadQualSeq_2,
				tmpAlignInfo_Nor1,
				tmpAlignInfo_Rcm1,
				tmpAlignInfo_Nor2,
				tmpAlignInfo_Rcm2);
		}
	}

	string returnFrontNodeAlignInfo_Nor1(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnAlignInfo_Nor1(index));
	}
	string returnFrontNodeAlignInfo_Rcm1(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnAlignInfo_Rcm1(index));
	}	
	string returnFrontNodeAlignInfo_Nor2(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnAlignInfo_Nor2(index));
	}
	string returnFrontNodeAlignInfo_Rcm2(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnAlignInfo_Rcm2(index));
	}	

	string returnFrontNodeReadQualSeq_1(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadQualSeq_1(index));
	}
	string returnFrontNodeReadQualSeq_2(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadQualSeq_2(index));
	}

	string returnFrontNodeReadSeq_1(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadSeq_1(index));
	}
	string returnFrontNodeReadSeq_2(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadSeq_2(index));
	}

	string returnFrontNodeReadName_alignNum_1(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadName_alignNum_1(index));
	}
	string returnFrontNodeReadName_alignNum_2(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadName_alignNum_2(index));
	}
	string returnFrontNodeReadNameAlignNum_1(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadName_alignNum_1(index));
	}
	string returnFrontNodeReadNameAlignNum_2(int index)
	{
		return ((frontNode->nextAlignInfoNode)->returnReadName_alignNum_2(index));
	}	
	int returnFrontNodeSize()
	{
		return ((frontNode->nextAlignInfoNode)->returnSize());
	}	
};

class Result_FixOneEndUnmapped_Array
{
private:
	vector<string> peAlignInfoVec_fixUnpair;
	vector<string> peAlignSamVec_fixUnpair;
	vector<string> peAlignSamVec_unpair_fixUnpair;
	//vector<string> peAlignInfoVec_pair_complete_fixUnpair;
	vector<string> peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair;
public:
	Result_FixOneEndUnmapped_Array* nextResultFixOneEndUnmappedArrayNode;
	
	// void free()
	// {
	// 	nextResultFixOneEndUnmappedArrayNode = NULL;
	// }

	Result_FixOneEndUnmapped_Array(int size)
	{
		for(int tmp = 0; tmp < size; tmp++)
		{
			peAlignInfoVec_fixUnpair.push_back("");
			peAlignSamVec_fixUnpair.push_back("");
			peAlignSamVec_unpair_fixUnpair.push_back("");
			//peAlignInfoVec_pair_complete_fixUnpair.push_back("");
			peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair.push_back("");
		}
	}

	void insert_PeAlignInfoVec_fixUnpair(const string& tmpStr, int index)
	{
		peAlignInfoVec_fixUnpair[index] = tmpStr;
	}

	void insert_PeAlignSamVec_fixUnpair(const string& tmpStr, int index)
	{
		peAlignSamVec_fixUnpair[index] = tmpStr;
	}

	void insert_PeAlignSamVec_unpair_fixUnpair(const string& tmpStr, int index)
	{
		peAlignSamVec_unpair_fixUnpair[index] = tmpStr;
	}

	// void insert_PeAlignInfoVec_pair_complete_fixUnpair(const string& tmpStr, int index)
	// {
	// 	peAlignInfoVec_pair_complete_fixUnpair[index] = tmpStr;
	// }

	void insert_PeAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair(const string& tmpStr, int index)
	{
		peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[index] = tmpStr;
	}

	int returnSize()
	{
		return peAlignInfoVec_fixUnpair.size();
	}

	void outputResultArray(
		ofstream& OutputSamFile_oneEndMapped_ofs,
		ofstream& tmpAlignIncompletePair_ofs,
		ofstream& OutputSamFile_oneEndMapped_unpair_ofs,
		ofstream& OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs
		)
	{
		for(int tmp = 0; tmp < peAlignInfoVec_fixUnpair.size(); tmp++)
		{
			if(peAlignSamVec_fixUnpair[tmp] != "")
			{
				OutputSamFile_oneEndMapped_ofs << peAlignSamVec_fixUnpair[tmp] << endl;
			}
			if(peAlignInfoVec_fixUnpair[tmp] != "")
			{
				tmpAlignIncompletePair_ofs << peAlignInfoVec_fixUnpair[tmp] << endl;
			}
			if(peAlignSamVec_unpair_fixUnpair[tmp] != "")
			{
				OutputSamFile_oneEndMapped_unpair_ofs << peAlignSamVec_unpair_fixUnpair[tmp] << endl;
			}
			if(peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] != "")
			{
				OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs << peAlignSamVec_bothEndsUnmapped_lowScore_fixUnpair[tmp] << endl;
			}
		}			
	}
};

class Result_FixOneEndUnmapped_Array_Queue
{
private:

public:
	Result_FixOneEndUnmapped_Array* frontNode;
	Result_FixOneEndUnmapped_Array* tailNode;

	void free()
	{
		frontNode = NULL;
		tailNode = NULL;
	}

	Result_FixOneEndUnmapped_Array_Queue()
	{
		Result_FixOneEndUnmapped_Array* p = new Result_FixOneEndUnmapped_Array(0);
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}

		p->nextResultFixOneEndUnmappedArrayNode = NULL;
		frontNode = p;
		tailNode = p;
	}

	void pushBack2ResultArrayQueue(Result_FixOneEndUnmapped_Array* tmpResultArray)
	{
		tmpResultArray->nextResultFixOneEndUnmappedArrayNode = NULL;
		tailNode->nextResultFixOneEndUnmappedArrayNode = tmpResultArray;		
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
			Result_FixOneEndUnmapped_Array* p = frontNode->nextResultFixOneEndUnmappedArrayNode;
			frontNode->nextResultFixOneEndUnmappedArrayNode = p->nextResultFixOneEndUnmappedArrayNode;
			if(tailNode == p) // when only one node exits
				tailNode = frontNode;
			delete p; p = NULL;
		}
	}

	bool atLeast3Node()
	{
		if(frontNode == tailNode)
			return false;
		else if((frontNode->nextResultFixOneEndUnmappedArrayNode) == tailNode)
			return false;
		else
			return true;
	}

	bool only2Node()
	{
		if((frontNode->nextResultFixOneEndUnmappedArrayNode) == tailNode)
			return true;
		else 
			return false;
	}

	bool resultQueueEmpty()
	{
		return (frontNode == tailNode);
	}

	void outputFrontResultArray(
		ofstream& OutputSamFile_oneEndMapped_ofs,
		ofstream& tmpAlignIncompletePair_ofs,
		ofstream& OutputSamFile_oneEndMapped_unpair_ofs,
		ofstream& OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs)
	{
		(frontNode->nextResultFixOneEndUnmappedArrayNode)->outputResultArray(
			OutputSamFile_oneEndMapped_ofs,
			tmpAlignIncompletePair_ofs,
			OutputSamFile_oneEndMapped_unpair_ofs,
			OutputSamFile_oneEndMapped_bothEndsUnmapped_lowScore_ofs);
	}	
};

class Result_FixHeadTail_Array
{
private:
	vector<string> peAlignSamVec_complete_pair;
	vector<string> peAlignSamVec_incomplete_pair;
	vector<string> peAlignSamVec_complete_unpair;
	vector<string> peAlignSamVec_incomplete_unpair;
	vector<string> peAlignSamVec_pair_lowScore;
public:
	Result_FixHeadTail_Array* nextResultFixHeadTailArrayNode;
	Result_FixHeadTail_Array(int size)
	{
		for(int tmp = 0; tmp < size; tmp++)
		{
			peAlignSamVec_complete_pair.push_back("");
			peAlignSamVec_incomplete_pair.push_back("");
			peAlignSamVec_complete_unpair.push_back("");
			peAlignSamVec_incomplete_unpair.push_back("");
			peAlignSamVec_pair_lowScore.push_back("");
		}
	}
	void insert_peAlignSamVec_complete_pair(const string& tmpStr, int index)
	{
		peAlignSamVec_complete_pair[index] = tmpStr;
	}
	void insert_peAlignSamVec_incomplete_pair(const string& tmpStr, int index)
	{
		peAlignSamVec_incomplete_pair[index] = tmpStr;
	}
	void insert_peAlignSamVec_complete_unpair(const string& tmpStr, int index)
	{
		peAlignSamVec_complete_unpair[index] = tmpStr;
	}
	void insert_peAlignSamVec_incomplete_unpair(const string& tmpStr, int index)
	{
		peAlignSamVec_incomplete_unpair[index] = tmpStr;
	}
	void insert_peAlignSamVec_pair_lowScore(const string& tmpStr, int index)
	{
		peAlignSamVec_pair_lowScore[index] = tmpStr;
	}
	int returnSize()
	{
		return peAlignSamVec_complete_pair.size();
	}
	void outputResultArray(
		ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
		ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
		ofstream& OutputSamFile_fixHeadTail_pair_lowScore_ofs
		)
	{
		int size = this->returnSize();
		for(int tmp = 0; tmp < size; tmp++)
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
};

class Result_FixHeadTail_Array_Queue
{
private:

public:
	Result_FixHeadTail_Array* frontNode;
	Result_FixHeadTail_Array* tailNode;

	void free()
	{
		frontNode = NULL;
		tailNode = NULL;
	}

	Result_FixHeadTail_Array_Queue()
	{
		Result_FixHeadTail_Array* p = new Result_FixHeadTail_Array(0);
		if(NULL == p)
		{
			cout << "failed to malloc a node" << endl;
		}

		p->nextResultFixHeadTailArrayNode = NULL;
		frontNode = p;
		tailNode = p;
	}

	void pushBack2ResultArrayQueue(Result_FixHeadTail_Array* tmpResultArray)
	{
		tmpResultArray->nextResultFixHeadTailArrayNode = NULL;
		tailNode->nextResultFixHeadTailArrayNode = tmpResultArray;		
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
			Result_FixHeadTail_Array* p = frontNode->nextResultFixHeadTailArrayNode;
			frontNode->nextResultFixHeadTailArrayNode = p->nextResultFixHeadTailArrayNode;
			if(tailNode == p) // when only one node exits
				tailNode = frontNode;
			delete p; p = NULL;
		}
	}

	bool atLeast3Node()
	{
		if(frontNode == tailNode)
			return false;
		else if((frontNode->nextResultFixHeadTailArrayNode) == tailNode)
			return false;
		else
			return true;
	}

	bool only2Node()
	{
		if((frontNode->nextResultFixHeadTailArrayNode) == tailNode)
			return true;
		else 
			return false;
	}

	bool resultQueueEmpty()
	{
		return (frontNode == tailNode);
	}

	void outputFrontResultArray(
		ofstream& OutputSamFile_fixHeadTail_complete_pair_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_pair_ofs,
		ofstream& OutputSamFile_fixHeadTail_complete_unpair_ofs,
		ofstream& OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
		ofstream& OutputSamFile_fixHeadTail_pair_lowScore_ofs)
	{
		(frontNode->nextResultFixHeadTailArrayNode)->outputResultArray(
			OutputSamFile_fixHeadTail_complete_pair_ofs,
			OutputSamFile_fixHeadTail_incomplete_pair_ofs,
			OutputSamFile_fixHeadTail_complete_unpair_ofs,
			OutputSamFile_fixHeadTail_incomplete_unpair_ofs,
			OutputSamFile_fixHeadTail_pair_lowScore_ofs);
	}
};


#endif