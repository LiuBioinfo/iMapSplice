// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORINSERTION_INFO_H
#define FIXDOUBLEANCHORINSERTION_INFO_H

//#include "nw_DP.h"

class FixDoubleAnchor_Insertion_Info
{
private:
	//*******************  Final results ******************//	
	int best_insertion_length;
	int best_insertion_prefix_match_length;
	int best_insertion_mismatch;
	int best_insertion_suffix_match_length;

	vector<int> best_insertion_mismatchPosVec;
	vector<char> best_insertion_mismatchCharVec;


	//*******************  Inter results ******************//
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

	vector< pair<int, int > > insertionSiteVec; // < insertionSite(prefix_match_length), mismatch >	

	int mismatchPos_interval_min;
public:
	FixDoubleAnchor_Insertion_Info()
	{
		best_insertion_length = -1;
		best_insertion_prefix_match_length = -1;
		best_insertion_mismatch = -1 ;
		mismatchPos_interval_min = 3;
	}
	int return_best_insertion_length()
	{
		return best_insertion_length;
	}
	int return_best_insertion_prefix_match_length()
	{
		return best_insertion_prefix_match_length;
	}
	int return_best_insertion_suffix_match_length()
	{
		return best_insertion_suffix_match_length;
	}
	int return_best_insertion_mismatch()
	{
		return best_insertion_mismatch;
	}

	int returnBestInsertionMismatchPosVecSize()
	{
		return best_insertion_mismatchPosVec.size();
	}
	int returnBestInsertionMismatchCharVecSize()
	{
		return best_insertion_mismatchCharVec.size();
	}

	int returnBestInsertionMismatchPos(int index)
	{
		return best_insertion_mismatchPosVec[index];
	}
	char returnBestInsertionMismatchChar(int index)
	{
		return best_insertion_mismatchCharVec[index];
	}

	int getMinMismatchInterval(vector<int>& best_insertion_mismatchPosVec_sorted)
	{
		for(int tmp = 0; tmp < best_insertion_mismatchPosVec.size(); tmp++)
		{
			best_insertion_mismatchPosVec_sorted.push_back(best_insertion_mismatchPosVec[tmp]);
		}				
		sort(best_insertion_mismatchPosVec_sorted.begin(), best_insertion_mismatchPosVec_sorted.end());		

		int tmpMinMismatchPosInterval = 100;
		for(int tmp = 0; tmp < best_insertion_mismatchPosVec_sorted.size()-1; tmp++)
		{
			int tmpIndex_1 = tmp;
			int tmpIndex_2 = tmp + 1;
			int tmpMismatchPos_1 = best_insertion_mismatchPosVec_sorted[tmpIndex_1];
			int tmpMismatchPos_2 = best_insertion_mismatchPosVec_sorted[tmpIndex_2];
			int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
			if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
				tmpMinMismatchPosInterval = tmpMismatchPosInterval;
		}
		return tmpMinMismatchPosInterval;
	}
		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> best_insertion_mismatchPosVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(best_insertion_mismatchPosVec_sorted);
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < best_insertion_mismatchPosVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = best_insertion_mismatchPosVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = best_insertion_mismatchPosVec_sorted[tmpIndex_2];
				int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
				if(tmpMismatchPosInterval < minimumMismatchGap)
				{
					cout << "some thing wrong in getMinMismatchInterval_2" << endl;
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
		int mismatchPosVecSize = best_insertion_mismatchPosVec.size();
		if(mismatchPosVecSize >= 3)
		{
			int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
			if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
				return false;
		}
		return true;
	}

	void copyMismatchPos2TargetVec(vector<int>& targetVec)
	{
		int vecSize = this->returnBestInsertionMismatchPosVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			int tmpMismatchPos = this->returnBestInsertionMismatchPos(tmp);
			targetVec.push_back(tmpMismatchPos);	
		}	
	}
	void copyMismatchChar2TargetVec(vector<char>& targetVec)
	{
		int vecSize = this->returnBestInsertionMismatchCharVecSize();
		for(int tmp = 0; tmp < vecSize; tmp++)
		{
			char tmpMismatchChar = this->returnBestInsertionMismatchChar(tmp);
			targetVec.push_back(tmpMismatchChar);
		}
	}


	void generateBestInsertionMismatchVec(int startLocInRead)
	{
		//if(STORE_MISMATCH_CHA)
			this->generateBestInsertionMismatchVec_Pos_Char(startLocInRead);
		//else
		//	this->generateBestInsertionMismatchVec_Pos(startLocInRead);
	}

	void generateBestInsertionMismatchVec_Pos(int startLocInRead)
	{
		//int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
		//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
		for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = donerMapMismatchPosVec[tmp];
			if(tmpMismatchPos <= best_insertion_prefix_match_length + startLocInRead - 1)
			{
				best_insertion_mismatchPosVec.push_back(tmpMismatchPos);
			}
		}

		for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
			if(tmpMismatchPos > (best_insertion_prefix_match_length + best_insertion_length + startLocInRead - 1))
			{
				best_insertion_mismatchPosVec.push_back(tmpMismatchPos);
			}
		}		
	}

	void generateBestInsertionMismatchVec_Pos_Char(int startLocInRead)
	{
		//int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
		//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
		for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = donerMapMismatchPosVec[tmp];
			if(tmpMismatchPos <= best_insertion_prefix_match_length + startLocInRead - 1)
			{
				best_insertion_mismatchPosVec.push_back(tmpMismatchPos);
				best_insertion_mismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
			}
		}

		for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
			if(tmpMismatchPos > (best_insertion_prefix_match_length + best_insertion_length + startLocInRead - 1))
			{
				best_insertion_mismatchPosVec.push_back(tmpMismatchPos);
				best_insertion_mismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
			}
		}		
	}

	void generateBestInsertionMismatchVec_Pos_Char_copyMismatchInfoToTargetVec(
		int startLocInRead, vector<int>& tmpGapMismatchPosVec, vector<char>& tmpGapMismatchCharVec)
	{
		//int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
		//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
		for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = donerMapMismatchPosVec[tmp];
			if(tmpMismatchPos <= best_insertion_prefix_match_length + startLocInRead - 1)
			{
				tmpGapMismatchPosVec.push_back(tmpMismatchPos);
				tmpGapMismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
			}
		}

		for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
			if(tmpMismatchPos > (best_insertion_prefix_match_length + best_insertion_length + startLocInRead - 1))
			{
				tmpGapMismatchPosVec.push_back(tmpMismatchPos);
				tmpGapMismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
			}
		}		
	}

	string returnInsertionSiteVecStr()
	{
		string tmpStr = "insertionSiteVec: \n";
		for(int tmp = 0; tmp < insertionSiteVec.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(insertionSiteVec[tmp].first) + " " + int_to_str(insertionSiteVec[tmp].second) + ", ";
		}
		tmpStr += "\n";
		return tmpStr;
	}

	string returnDonerAndAcceptorMapCumulativeMismatchNumVecStr()
	{
		string tmpStr = "donerMapCumulativeMismatchNumVec:\n";
		for(int tmp = 0; tmp < donerMapCumulativeMismatchNumVec.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(donerMapCumulativeMismatchNumVec[tmp]) + ", ";
		}
		tmpStr = tmpStr + "\nacceptorMapCumulativeMismatchNumVec:\n";
		for(int tmp = 0; tmp < acceptorMapCumulativeMismatchNumVec.size(); tmp++)
		{
			tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(acceptorMapCumulativeMismatchNumVec[tmp]) + ", ";
		}
		tmpStr += "\n";
		return tmpStr;
	}



	bool detectBestInsertion_lessMismatch(
		int toFix_read_start, int toFix_read_end,
		int toFix_chrom_start, int toFix_chrom_end,
		const string& readSeq_inProcess,
		Index_Info* indexInfo, int chrom_name_int,
		int max_allowed_mismatchNum)
	{
		int toFix_read_length = toFix_read_end - toFix_read_start + 1;
		int toFix_chrom_length = toFix_chrom_end - toFix_chrom_start + 1;
		best_insertion_length = toFix_read_length - toFix_chrom_length;
		//cout << "start to detectBestInsertion_lessMismatch ..." << endl;
		//cout << "toFix_read_start: " << toFix_read_start << endl;
		//cout << "toFix_read_end: " << toFix_read_end << endl;
		//cout << "toFix_read_length: " << toFix_read_length << endl << endl;
		//cout << "toFix_chrom_start: " << toFix_chrom_start << endl;
		//cout << "toFix_chrom_end: " << toFix_chrom_end << endl;
		//cout << "toFix_chrom_length: " << toFix_chrom_length << endl << endl;

		if(best_insertion_length <= 0)
			return false;
		
		//cout << "start to scan genome ..." << endl;

		this->scanGenomeAndReadSeq(
			toFix_read_start, toFix_read_end,
			toFix_chrom_start, toFix_chrom_end, readSeq_inProcess,
			indexInfo, chrom_name_int, max_allowed_mismatchNum
			);		

		//cout << "scan genome results: " << endl << endl;
		//cout << "max_extension_forward_doner: " << max_extension_forward_doner << endl << endl;
		//cout << "max_extension_backward_acceptor: " << max_extension_backward_acceptor << endl << endl;
		//cout << this->returnDonerAndAcceptorMapCumulativeMismatchNumVecStr() << endl;
		//cout << "**************************" << endl << "start to generateInsertionSiteVec ...." << endl;	
		this->generateInsertionSiteVec(
			toFix_read_start, toFix_read_end,
			toFix_chrom_start, toFix_chrom_end, 
			indexInfo, chrom_name_int, max_allowed_mismatchNum);	

		//cout << "generated insertion vec: " << endl;
		//cout << this->returnInsertionSiteVecStr() << endl;
		//cout << "finish generating all insertion site ..." << endl;
		//cout << "********************************" << endl << "start to select best insertion site ...." << endl;
		bool fixInsertionBool = this->selectBestInsertionSite(max_allowed_mismatchNum, toFix_read_length);
		//cout << "finish selecting best insertion site ? : " << fixInsertionBool << endl;	
		if(fixInsertionBool)
		{
			this->generateBestInsertionMismatchVec(toFix_read_start);
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
		//return fixInsertionBool;
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
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		for(int tmp = 0; tmp < toFix_read_length; tmp++)
		{
			int tmpBasePosInChr = toFix_chrom_start + tmp;
			if(tmpBasePosInChr > tmpChromLength)
				break;
			char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFix_chrom_start + tmp);
			if(readSeq_inProcess.at(toFix_read_start - 1 + tmp) 
				!= charInRef )
			{
				tmpDonerMapCumulativeMismatchNum ++;
				donerMapMismatchPosVec.push_back(tmp+toFix_read_start);
				donerMapMismatchCharVec.push_back(charInRef);			
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
			if(toFix_chrom_end - tmp < 1)
				break;
			char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFix_chrom_end - tmp);
			if(readSeq_inProcess.at(toFix_read_end - tmp - 1) 
				!= charInRef )
			{
				tmpAcceptorMapCumulativeMismatchNum ++;
				acceptorMapMismatchPosVec.push_back(toFix_read_end - tmp);
				acceptorMapMismatchCharVec.push_back(charInRef);		
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

	void generateInsertionSiteVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum)
	{
		int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
		int toFixChromSeqLength = toFixSeqMapPos_acceptor - toFixSeqMapPos_doner + 1;
		best_insertion_length = toFixSeqLength - toFixChromSeqLength;

		for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner; 
			tmp_forwared_doner ++ )
		{
			int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner - best_insertion_length;
			if((tmp_backward_acceptor > max_extension_backward_acceptor)||(tmp_backward_acceptor < 0))
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
				insertionSiteVec.push_back(pair<int,int>(tmp_forwared_doner, tmp_mismatch_sum));
			}
		}
	}

	bool selectBestInsertionSite(int max_allowed_mismatchNum, int toFix_read_length)
	{
		bool bestSite_bool = false;
		int tmp_mismatch_min = 20;
		for(int tmp = 0; tmp < insertionSiteVec.size(); tmp++)
		{
			if(insertionSiteVec[tmp].second <= tmp_mismatch_min)
			{
				tmp_mismatch_min = insertionSiteVec[tmp].second;
				best_insertion_prefix_match_length = insertionSiteVec[tmp].first;				
			}
		}
		if(tmp_mismatch_min <= max_allowed_mismatchNum)
		{
			best_insertion_mismatch = tmp_mismatch_min;
			best_insertion_suffix_match_length = toFix_read_length - best_insertion_prefix_match_length - best_insertion_length;
			return true;
		}
		else
			return false;
	}
};

#endif