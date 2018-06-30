// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORDELETION_INFO_H
#define FIXDOUBLEANCHORDELETION_INFO_H

//#include "nw_DP.h"
//#include "fixDoubleAnchorNWDP_info.h"

class FixDoubleAnchor_Deletion_Info
{
private:
	// ***************  final results ...*************************//
	int best_deletion_length;
	int best_deletion_prefix_match_length;
	int best_deletion_mismatch;

	vector<int> best_deletion_mismatchPosVec;
	vector<char> best_deletion_mismatchCharVec;


	//*****************   inter results .... **********************//
	int max_extension_forward_doner;
	int max_extension_backward_acceptor;

	vector<int> donerMapCumulativeMismatchNumVec; 
	// donerMapCumulativeMismatchNumVec[tmp] = Y; 
	// Y is mismatch # if readSeq.substr(toFixSeqLocInRead_start-1, tmp+1) 
	// or toFixSeq.substr(0, tmp+1)  mapped to chromSeq.substr(toFixSeqMapPos_doner - 1, tmp+1); 
	vector<int> donerMapMismatchPosVec;
	vector<char> donerMapMismatchCharVec;

	vector<int> acceptorMapCumulativeMismatchNumVec;
	// acceptorMapCumulativeMismatchNumVec[tmp] = Y;
	// Y is mismatch # if readSeq.substr(toFixSeqLocInRead_end - 1 - tmp, tmp+1) 
	// or toFixSeq.substr(toFixSeqLength - 1- tmp, tmp+1) mapped to chromSeq.substr(toFixSeqMapPos_end - tmp - 1, tmp+1);
	vector<int> acceptorMapMismatchPosVec;
	vector<char> acceptorMapMismatchCharVec;

	vector< pair<int, int > > deletionSiteVec; // < deletionSite(prefix_match_length), mismatch >

	int mismatchPos_interval_min;
public:
	int returnBestDeletionMismatchPosVecSize()
	{
		return best_deletion_mismatchPosVec.size();
	}

	int returnBestDeletionMismatchCharVecSize()
	{
		return best_deletion_mismatchCharVec.size();
	}

	int returnBestDeletionMismatchPos(int index)
	{
		return best_deletion_mismatchPosVec[index];
	}

	char returnBestDeletionMismatchChar(int index)
	{
		return best_deletion_mismatchCharVec[index];
	}

	int return_best_deletion_length()
	{
		return best_deletion_length;
	}
	int return_best_deletion_prefix_match_length()
	{
		return best_deletion_prefix_match_length;
	}
	int return_best_deletion_mismatch()
	{
		return best_deletion_mismatch;
	}

	int getMinMismatchInterval(vector<int>& best_deletion_mismatchPosVec_sorted)
	{
		for(int tmp = 0; tmp < best_deletion_mismatchPosVec.size(); tmp++)
		{
			best_deletion_mismatchPosVec_sorted.push_back(best_deletion_mismatchPosVec[tmp]);
		}				
		sort(best_deletion_mismatchPosVec_sorted.begin(), best_deletion_mismatchPosVec_sorted.end());		
		int tmpMinMismatchPosInterval = 100;
		for(int tmp = 0; tmp < best_deletion_mismatchPosVec_sorted.size()-1; tmp++)
		{
			int tmpIndex_1 = tmp;
			int tmpIndex_2 = tmp + 1;
			int tmpMismatchPos_1 = best_deletion_mismatchPosVec_sorted[tmpIndex_1];
			int tmpMismatchPos_2 = best_deletion_mismatchPosVec_sorted[tmpIndex_2];
			int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
			if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
				tmpMinMismatchPosInterval = tmpMismatchPosInterval;
		}
		return tmpMinMismatchPosInterval;
	}

		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> best_deletion_mismatchPosVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(best_deletion_mismatchPosVec_sorted);
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < best_deletion_mismatchPosVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = best_deletion_mismatchPosVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = best_deletion_mismatchPosVec_sorted[tmpIndex_2];
				int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
				if(tmpMismatchPosInterval < minimumMismatchGap)
				{
					//cout << "some thing wrong in getMinMismatchInterval_2" << endl;
					exit(1);
				}
				else if(tmpMismatchPosInterval == minimumMismatchGap)
				{	
					if(!checkedMinimumGap_bool)
						checkedMinimumGap_bool = true;
					else
						return tmpMismatchPosInterval;
				}
				else
				{
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval_2)
						tmpMinMismatchPosInterval_2 = tmpMismatchPosInterval;
				}
			}
			return tmpMinMismatchPosInterval_2;
		}

	bool filterInterResults()
	{
		int mismatchPosVecSize = best_deletion_mismatchPosVec.size();
		if(mismatchPosVecSize >= 3)
		{
			int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
			if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
				return false;
		}
		return true;
	}


	FixDoubleAnchor_Deletion_Info()
	{
		best_deletion_length = -1;
		best_deletion_prefix_match_length = -1;
		best_deletion_mismatch = -1;
		mismatchPos_interval_min = 3;
	}

	void copyMismatchPos2TargetVec(vector<int>& targetVec)
	{	
		int vecSize = this->returnBestDeletionMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = this->returnBestDeletionMismatchPos(tmp);
			targetVec.push_back(tmpMismatchPos);	
		}	
	}

	void copyMismatchChar2TargetVec(vector<char>& targetVec)
	{
		int vecSize = this->returnBestDeletionMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = this->returnBestDeletionMismatchChar(tmp);
			targetVec.push_back(tmpMismatchPos);	
		}	
	}

	void generateBestDeletionMismatchVec(int startLocInRead)
	{
		//if(STORE_MISMATCH_CHA)
			this->generateBestDeletionMismatchVec_Pos_Char(startLocInRead);
		//else
		//	this->generateBestDeletionMismatchVec_Pos(startLocInRead);
	}

	// void generateBestDeletionMismatchVec_Pos(int startLocInRead)
	// {
	// 	//int best_deletion_suffix_match_length = toFixSeqLength - best_deletion_prefix_match_length
	// 	for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
	// 	{
	// 		int tmpMismatchPos = donerMapMismatchPosVec[tmp];
	// 		if(tmpMismatchPos <= best_deletion_prefix_match_length + startLocInRead - 1)
	// 		{
	// 			best_deletion_mismatchPosVec.push_back(tmpMismatchPos);
	// 		}
	// 	}

	// 	for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
	// 	{
	// 		int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
	// 		if(tmpMismatchPos > best_deletion_prefix_match_length + startLocInRead - 1)
	// 		{
	// 			best_deletion_mismatchPosVec.push_back(tmpMismatchPos);
	// 		}
	// 	}			
	// }

	void generateBestDeletionMismatchVec_Pos_Char(int startLocInRead)
	{
		// cout << "start to generateBestDeletionMismatchVec ..." << endl;
		// cout << "best_deletion_prefix_match_length: " << best_deletion_prefix_match_length << endl;
		// cout << "startLocInRead: " << startLocInRead << endl;
		// cout << "tmpMismatchPos_doner: " << endl;
		for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = donerMapMismatchPosVec[tmp];
			if(tmpMismatchPos <= best_deletion_prefix_match_length + startLocInRead - 1)
			{
				//cout << " " << tmpMismatchPos  << ",";
				best_deletion_mismatchPosVec.push_back(tmpMismatchPos);
				best_deletion_mismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
			}
		}
		//cout << endl;
		//cout << "tmpMismatchPos_acceptor: ";
		for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
			if(tmpMismatchPos > best_deletion_prefix_match_length + startLocInRead - 1)
			{
				//cout << " " << tmpMismatchPos  << ",";
				best_deletion_mismatchPosVec.push_back(tmpMismatchPos);
				best_deletion_mismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
			}
		}		
	}

	void generateBestDeletionMismatchVec_Pos_Char_copyMismatchInfoToTargetVec(
		int startLocInRead, vector<int>& targetMismatchPosVec, vector<char>& targetMismatchCharVec)
	{
		for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = donerMapMismatchPosVec[tmp];
			if(tmpMismatchPos <= best_deletion_prefix_match_length + startLocInRead - 1)
			{
				targetMismatchPosVec.push_back(tmpMismatchPos);
				targetMismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
			}
		}

		for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
			if(tmpMismatchPos > best_deletion_prefix_match_length + startLocInRead - 1)
			{
				targetMismatchPosVec.push_back(tmpMismatchPos);
				targetMismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
			}
		}		
	}

	// void nw_DP(
	// 	int toFix_read_start, int toFix_read_end,
	// 	int toFix_chrom_start, int toFix_chrom_end,
	// 	const string& readSeq_inProcess,
	// 	Index_Info* indexInfo, int chrom_name_int,
	// 	int max_allowed_mismatchNum)
	// {
	// 	string readSeqToProcess 
	// 		= readSeq_inProcess.substr(toFix_read_start-1, toFix_read_end - toFix_read_start + 1);
	// 	string chromSeqToProcess 
	// 		= indexInfo->returnChromStrSubstr(chrom_name_int, toFix_chrom_start, toFix_chrom_end - toFix_chrom_start + 1);
	// }	

	bool detectBestDeletion_lessMismatch(
		int toFix_read_start, int toFix_read_end,
		int toFix_chrom_start, int toFix_chrom_end,
		const string& readSeq_inProcess,
		Index_Info* indexInfo, int chrom_name_int,
		int max_allowed_mismatchNum)
	{
		// cout << "**********************\n start to fix deletion ...." << endl;
		// cout << "toFix_read_start: " << toFix_read_start << endl;
		// cout << "toFix_read_end: " << toFix_read_end << endl;
		// cout << "toFix_chrom_start: " << toFix_chrom_start << endl;
		// cout << "toFix_chrom_end: " << toFix_chrom_end << endl;
		int toFix_read_length = toFix_read_end - toFix_read_start + 1;
		int toFix_chrom_length = toFix_chrom_end - toFix_chrom_start + 1;
		best_deletion_length = toFix_chrom_length - toFix_read_length;
		//cout << "best_deletion_length: " << best_deletion_length << endl;

		if(best_deletion_length <= 0)
			return false;
		this->scanGenomeAndReadSeq(
			toFix_read_start, toFix_read_end,
			toFix_chrom_start, toFix_chrom_end, readSeq_inProcess,
			indexInfo, chrom_name_int, max_allowed_mismatchNum);		
		this->generateDeletionSiteVec(
			toFix_read_start, toFix_read_end,
			toFix_chrom_start, toFix_chrom_end, 
			indexInfo, chrom_name_int, max_allowed_mismatchNum);	

		bool fixDeletionBool = this->selectBestDeletionSite(max_allowed_mismatchNum);
		if(fixDeletionBool)
		{
			this->generateBestDeletionMismatchVec(toFix_read_start);
			bool filter_bool = this->filterInterResults();
			if(filter_bool)
			{
				return true;
			}
			else 
				return false;
		}
		else
			return false;
		//cout << "end of fixDeletion ....." << endl;
		//cout << "******************************" << endl;
		//return fixDeletionBool; 
	}

	void generateDeletionSiteVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum)
	{
		//cout << "start to generateDeletionSiteVec: " << endl;
		int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
		for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner; 
			tmp_forwared_doner ++ )
		{
			int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
			if(tmp_backward_acceptor > max_extension_backward_acceptor)
				continue;
			int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
			int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
			int tmp_doner_mismatch 
				= donerMapCumulativeMismatchNumVec[index_donerMapCumulativeMismatchNumVec];
			int tmp_acceptor_mismatch 
				= acceptorMapCumulativeMismatchNumVec[index_acceptorMapCumulativeMismatchNumVec];
			int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;
			if(tmp_mismatch_sum <= max_allowed_mismatchNum)
			{
				deletionSiteVec.push_back(pair<int,int>(tmp_forwared_doner, tmp_mismatch_sum));
				//cout << "; tmp_forwared_doner: " << tmp_forwared_doner << " tmp_mismatch_sum: " << tmp_mismatch_sum;
			}
		}
		//cout << endl;
	}

	bool selectBestDeletionSite(int max_allowed_mismatchNum)
	{
		bool bestSite_bool = false;
		int tmp_mismatch_min = 20;
		for(int tmp = 0; tmp < deletionSiteVec.size(); tmp++)
		{
			if(deletionSiteVec[tmp].second <= tmp_mismatch_min)
			{
				tmp_mismatch_min = deletionSiteVec[tmp].second;
				best_deletion_prefix_match_length = deletionSiteVec[tmp].first;				
			}
		}
		if(tmp_mismatch_min <= max_allowed_mismatchNum)
		{
			best_deletion_mismatch = tmp_mismatch_min;
			return true;
		}
		else
			return false;
	}

	void scanGenomeAndReadSeq(
		int toFix_read_start, int toFix_read_end,
		int toFix_chrom_start, int toFix_chrom_end,
		const string& readSeq_inProcess,
		Index_Info* indexInfo, int chrNameInt,
		int max_allowed_mismatchNum)
	{
		int toFix_read_length = toFix_read_end - toFix_read_start + 1;
		int toFix_chrom_length = toFix_chrom_end - toFix_chrom_start + 1;
		
		max_extension_forward_doner = toFix_read_length;
		max_extension_backward_acceptor = toFix_read_length;		

		int tmpDonerMapCumulativeMismatchNum = 0;
		donerMapCumulativeMismatchNumVec.push_back(0);
		for(int tmp = 0; tmp < toFix_read_length; tmp++)
		{
			char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFix_chrom_start + tmp);
			if(readSeq_inProcess.at(toFix_read_start - 1 + tmp) 
				!= charInRef )
			{
				tmpDonerMapCumulativeMismatchNum ++;
				//if(STORE_MISMATCH_POS)
				//{
					donerMapMismatchPosVec.push_back(tmp+ toFix_read_start);
					//if(STORE_MISMATCH_CHA)
					//{
						donerMapMismatchCharVec.push_back(charInRef);
					//}
				//}
			}
			if(tmpDonerMapCumulativeMismatchNum > max_allowed_mismatchNum)
			{
				max_extension_forward_doner = tmp;
				break;
			}
			donerMapCumulativeMismatchNumVec.push_back(tmpDonerMapCumulativeMismatchNum);
			//cout << "tmp: " << tmp << " tmpDonerMapCumulativeMismatchNum: " << tmpDonerMapCumulativeMismatchNum << endl;
		}

		int tmpAcceptorMapCumulativeMismatchNum = 0;
		acceptorMapCumulativeMismatchNumVec.push_back(0);
		for(int tmp = 0; tmp < toFix_read_length; tmp++)
		{
			char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFix_chrom_end - tmp);
			if(readSeq_inProcess.at(toFix_read_end - tmp - 1) 
				!= charInRef )
			{
				tmpAcceptorMapCumulativeMismatchNum ++;
				//if(STORE_MISMATCH_POS)
				//{
					acceptorMapMismatchPosVec.push_back( toFix_read_end - tmp);
					//if(STORE_MISMATCH_CHA)
					//{
						acceptorMapMismatchCharVec.push_back(charInRef);
					//}
				//}
			}
			if(tmpAcceptorMapCumulativeMismatchNum > max_allowed_mismatchNum)
			{
				max_extension_backward_acceptor = tmp;
				break;
			}
			acceptorMapCumulativeMismatchNumVec.push_back(tmpAcceptorMapCumulativeMismatchNum);
			//cout << "tmp: " << tmp << "tmpAcceptorMapCumulativeMismatchNum: " << tmpAcceptorMapCumulativeMismatchNum << endl;
		}				
	}
};

#endif