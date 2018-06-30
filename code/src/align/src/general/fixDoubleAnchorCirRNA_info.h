// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORCIRRNA_INFO_H
#define FIXDOUBLEANCHORCIRRNA_INFO_H

class FixDoubleAnchor_CirRNA_Info
{
private:
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

	int max_extension_forward_doner;
	int max_extension_backward_acceptor;

	vector< pair< int, pair<int, int> > > spliceSiteVec; // < < splice_site = forwared_doner_# , mismatch # >, flank_string_case >

	int bestSpliceSiteIndex;
	vector<int> bestSpliceMismatchPosVec;
	vector<char> bestSpliceMismatchCharVec;

	int mismatchPos_interval_min;	
public:
	FixDoubleAnchor_CirRNA_Info()
	{
		bestSpliceSiteIndex = -1;
		mismatchPos_interval_min = 3;
		//indel_length_max = INDEL_BESIDE_SJ_LEN_MAX;
	}	

		int returnBestSpliceMismatchPosVecSize()
		{
			return bestSpliceMismatchPosVec.size();
		}
		int returnBestSpliceMismatchCharVecSize()
		{
			return bestSpliceMismatchCharVec.size();
		}

		int returnMismatchPosInRead(int index_bestSpliceMismatchPosVec)
		{
			return bestSpliceMismatchPosVec[index_bestSpliceMismatchPosVec];
		}
		char returnMismatchChar(int index_bestSpliceMismatchCharVec)
		{
			return bestSpliceMismatchCharVec[index_bestSpliceMismatchCharVec];
		}

		int returnBestSplice_flankString()
		{
			int flankString_case = spliceSiteVec[bestSpliceSiteIndex].second.second;
			return flankString_case;
		}
		bool returnBestSplice_canonicalOrNot()
		{
			return this->SJ_canonicalOrNot(this->returnBestSplice_flankString());
		}
		int returnBestCirRNAspliceSite_mismatchNum()
		{
			int bestSJ_mismatchNum = spliceSiteVec[bestSpliceSiteIndex].second.first;
			return bestSJ_mismatchNum;
		}

		int returnBestCirRNAspliceSite_prefixMatchLength()
		{
			return spliceSiteVec[bestSpliceSiteIndex].first;
		}

		int returnFlankStringCase(const string& flank_string)   
		{
			if(flank_string == "ATAC")
			{
				return 1;
			}
			else if(flank_string == "GTAT")
			{
				return 2;
			}
			else if(flank_string == "CTGC")
			{
				return 3;
			}
			else if(flank_string == "GCAG")
			{
				return 4;
			}
			else if(flank_string == "GTAG")
			{
				return 5;
			}
			else if(flank_string == "CTAC")
			{
				return 6;
			}
			else
			{
				return 0;
			}
		}

		int getMinMismatchInterval(vector<int>& bestSpliceMismatchPosVec_sorted)
		{
				// cout << "before sorting the vector<int> ..." << endl;
				// for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
				// {
				// 	cout << bestSpliceMismatchPosVec[tmp] << ",";
				// }
				// cout << endl;
				// sort(bestSpliceMismatchPosVec.begin(), bestSpliceMismatchPosVec.end());
				// cout << "after sorting ..." << endl;
				// for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
				// {
				// 	cout << bestSpliceMismatchPosVec[tmp] << ",";
				// }
				// cout << endl;
				for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
				{
					bestSpliceMismatchPosVec_sorted.push_back(bestSpliceMismatchPosVec[tmp]);
				}				
				sort(bestSpliceMismatchPosVec_sorted.begin(), bestSpliceMismatchPosVec_sorted.end());

				int tmpMinMismatchPosInterval = 100;
				for(int tmp = 0; tmp < bestSpliceMismatchPosVec_sorted.size()-1; tmp++)
				{
					int tmpIndex_1 = tmp;
					int tmpIndex_2 = tmp + 1;
					int tmpMismatchPos_1 = bestSpliceMismatchPosVec_sorted[tmpIndex_1];
					int tmpMismatchPos_2 = bestSpliceMismatchPosVec_sorted[tmpIndex_2];
					int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
						tmpMinMismatchPosInterval = tmpMismatchPosInterval;
				}
				return tmpMinMismatchPosInterval;
		}

		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> bestSpliceMismatchPosVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(bestSpliceMismatchPosVec_sorted);
			//cout << "minimumMismatchGap: " << minimumMismatchGap << endl;
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < bestSpliceMismatchPosVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = bestSpliceMismatchPosVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = bestSpliceMismatchPosVec_sorted[tmpIndex_2];
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

		bool filterInterResults(int maxMismatchNum_semi_noncanonical,
			int bestSpliceSiteIndex)
		{
			//cout << "start to filter interResults ..." << endl;
			int mismatchPosVecSize = bestSpliceMismatchPosVec.size();
			//cout << "mismatchPosVecSize: " << mismatchPosVecSize << endl;
			if(mismatchPosVecSize >= 3)
			{
				int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
				//cout << "tmpMinMismatchPosInterval_2: " << tmpMinMismatchPosInterval << endl;
				if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
					return false;
			}

			int tmpSJcase = (spliceSiteVec[bestSpliceSiteIndex].second).second;
			if(tmpSJcase < 5)
			{
				if(mismatchPosVecSize > maxMismatchNum_semi_noncanonical)
					return false;
			}
			return true;
		}

		double returnSJpenalty(int flank_string_case)
		{
			if((flank_string_case == 5)||(flank_string_case == 6)) 
				return CANONICAL_SJ_PENALTY;
			else if((flank_string_case >= 1)&&(flank_string_case <= 4))
				return SEMICANONICAL_SJ_PENALTY;
			else if(flank_string_case == 0)
				return NONCANONICAL_SJ_PENALTY;
			else
			{
				cout << "incorrect flank_string_case in FixDoubleAnchorSplice_Info.h" << endl;
				exit(1);
			}
		}

		
		bool SJ_canonicalOrNot(int flank_string_case)
		{
			if((flank_string_case == 5)||(flank_string_case == 6)) 
				return true;
			else
				return false;
		}

		bool SJ_semicanonicalOrNot(int flank_string_case)
		{
			if((flank_string_case >= 1)&&(flank_string_case <= 4)) 
				return true;
			else
				return false;
		}

		bool SJ_noncanonicalOrNot(int flank_string_case)
		{
			if(flank_string_case == 0) 
				return true;
			else
				return false;
		}


	void scanGenomeAndReadSeq(
		int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
		int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
		const string& readSeq_inProcess, 
		Index_Info* indexInfo, int chrNameInt, 
		int max_allowed_mismatchNum //int* mismatchNum_detected, int* prefixMatch_length
		)	
	{
		//cout << endl << "start to scan scanGenomeAndReadSeq " << endl << endl;
		int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
		max_extension_forward_doner = toFixSeqLength;
		max_extension_backward_acceptor = toFixSeqLength;

		int tmpDonerMapCumulativeMismatchNum = 0;
		donerMapCumulativeMismatchNumVec.push_back(0);
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		for(int tmp = 0; tmp < toFixSeqLength; tmp++)
		{
			int tmpBasePosInChr = toFixSeqMapPos_doner + tmp;
			if(tmpBasePosInChr > tmpChromLength)
				break;
			char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_doner + tmp);
			if(readSeq_inProcess.at(toFixSeqLocInRead_start - 1 + tmp) 
				!= charInRef )
			{
				tmpDonerMapCumulativeMismatchNum ++;
				donerMapMismatchPosVec.push_back(tmp+toFixSeqLocInRead_start);
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
		for(int tmp = 0; tmp < toFixSeqLength; tmp++)
		{
			if(toFixSeqMapPos_acceptor - tmp < 1)
				break;
			char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_acceptor - tmp);
			if(readSeq_inProcess.at(toFixSeqLocInRead_end - tmp - 1) 
				!= charInRef )
			{
				tmpAcceptorMapCumulativeMismatchNum ++;
				acceptorMapMismatchPosVec.push_back(toFixSeqLocInRead_end - tmp);
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
		return;
		//int tmpSpliceSite = 0;
	}

	void generateSpliceSiteVec_noAnnotationProvided(
		int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
		int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
		Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum)
	{
		//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
		int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
		int tmpChromLength = indexInfo->returnChromLength(chrNameInt);
		for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner; tmp_forwared_doner ++ )
		{
			int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
			if(tmp_backward_acceptor > max_extension_backward_acceptor)
				continue;
			int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
			int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
			int tmp_doner_mismatch = donerMapCumulativeMismatchNumVec[index_donerMapCumulativeMismatchNumVec];
			int tmp_acceptor_mismatch = acceptorMapCumulativeMismatchNumVec[index_acceptorMapCumulativeMismatchNumVec];
			int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;

			if(tmp_mismatch_sum <= max_allowed_mismatchNum)
			{
				int tmpChromPos_doner = toFixSeqMapPos_doner - 1 + tmp_forwared_doner;
				int tmpChromPos_acceptor = toFixSeqMapPos_acceptor + 1 - tmp_backward_acceptor;
				if((tmpChromPos_doner + 2 > tmpChromLength)||(tmpChromPos_acceptor - 2 < 1))
					continue;
				string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
				int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
				spliceSiteVec.push_back(pair< int, pair<int, int> >(
				tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
			}
		}
		return;
	}

	int selectBestSpliceSite() 
	// penalty = mismatchNum + SJ_penalty(canon=0, semi=1.5, noncanon=2);
	{
		double currentBestSJ_penalty = 100000.0;
		int currentBestSJ_index = -1;
		for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
		{
			double tmpMismatchNum = (double)spliceSiteVec[tmp].second.first;
			double tmpSJpenalty = this->returnSJpenalty(spliceSiteVec[tmp].second.second);
			double tmpPenalty = tmpMismatchNum + tmpSJpenalty;
			if(tmpPenalty < currentBestSJ_penalty)
			{
				currentBestSJ_penalty = tmpPenalty;
				currentBestSJ_index = tmp;
			}
		}			
		return currentBestSJ_index;
	}

	void generateBestSpliceMismatchVec_Pos_Char(int toFixSeqLength, int startLocInRead)
	{
		int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
		//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
		for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = donerMapMismatchPosVec[tmp];
			if(tmpMismatchPos <= bestSplice_doner_length + startLocInRead - 1)
			{
				bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
				bestSpliceMismatchCharVec.push_back(donerMapMismatchCharVec[tmp]);
			}
		}
		for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
		{
			int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
			if(tmpMismatchPos > bestSplice_doner_length + startLocInRead - 1)
			{
				bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
				bestSpliceMismatchCharVec.push_back(acceptorMapMismatchCharVec[tmp]);
			}
		}
	}

	bool detectBestCirRNAspliceSite(//_prefer_canonical_lessMismatch
		int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
		int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
		const string& readSeq_inProcess, 
		Index_Info* indexInfo, int chrNameInt, 
		int max_allowed_mismatchNum)
	{
		//cout << "start to detect best splice site " << endl;
		this->scanGenomeAndReadSeq(
			toFixSeqLocInRead_start, toFixSeqLocInRead_end,
			toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, 
			readSeq_inProcess, indexInfo, chrNameInt, 
			max_allowed_mismatchNum //int* mismatchNum_detected, int* prefixMatch_length
			);			

		//cout << "scanGenomeResults: " << endl 
		//	<< this->returnDonerAndAcceptorMapCumulativeMismatchNumVecStr() << endl;
		this->generateSpliceSiteVec_noAnnotationProvided(
			toFixSeqLocInRead_start, toFixSeqLocInRead_end,
			toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
			indexInfo, chrNameInt, max_allowed_mismatchNum);
		//cout << "candi splice site: " << endl 
		//	<< this->returnSpliceSiteVecStr() << endl; 
		//bestSpliceSiteIndex = this->selectBestSpliceSite_prefer_canonical_lessMismatch();
		bestSpliceSiteIndex = this->selectBestSpliceSite();
		//cout << "bestSpliceSiteIndex: " << bestSpliceSiteIndex << endl;
		//int LengthOfSeqPerMismatchAllowed_semi_noncanonical = 15;
		//int maxMismatchNum_semi_noncanonical = readSeq_inProcess.length() / LengthOfSeqPerMismatchAllowed_semi_noncanonical;
		int maxMismatchNum_semi_noncanonical = max_allowed_mismatchNum;
		if(bestSpliceSiteIndex >= 0)
		{
			this->generateBestSpliceMismatchVec_Pos_Char(toFixSeqLocInRead_end-toFixSeqLocInRead_start+1, toFixSeqLocInRead_start);
			bool filter_bool = filterInterResults(maxMismatchNum_semi_noncanonical, bestSpliceSiteIndex);
			if(filter_bool)
			{
				return true;
			}
			else
				return false;
		}
		else
			return false;
	}
};
#endif