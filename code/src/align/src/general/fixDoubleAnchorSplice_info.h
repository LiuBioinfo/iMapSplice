// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORSPLICE_INFO_H
#define FIXDOUBLEANCHORSPLICE_INFO_H

//#include "fixDoubleAnchorSplice_complicate_info.h"
#include "nw_DP.h"

class FixDoubleAnchor_Splice_Info
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
		void copyMismatchPos2TargetVec(vector<int>& targetMismatchPosVec)
		{
			for(int tmp = 0; tmp < bestSpliceMismatchPosVec.size(); tmp++)
			{
				targetMismatchPosVec.push_back(bestSpliceMismatchPosVec[tmp]);
			}
		}
		void copyMismatchChar2TargetVec(vector<char>& targetMismatchCharVec)
		{
			for(int tmp = 0; tmp < bestSpliceMismatchCharVec.size(); tmp++)
			{
				targetMismatchCharVec.push_back(bestSpliceMismatchCharVec[tmp]);
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

		FixDoubleAnchor_Splice_Info()
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

		void generateBestSpliceMismatchVec(int toFixSeqLength, int startLocInRead)
		{
			if(STORE_MISMATCH_CHA)
				this->generateBestSpliceMismatchVec_Pos_Char(toFixSeqLength, startLocInRead);
			else
				this->generateBestSpliceMismatchVec_Pos(toFixSeqLength, startLocInRead);
		}

		void generateBestSpliceMismatchVec_Pos(int toFixSeqLength, int startLocInRead)
		{
			int bestSplice_doner_length = spliceSiteVec[bestSpliceSiteIndex].first;
			//int bestSplice_acceptor_length = toFixSeqLength - bestSplice_doner_length;
			for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = donerMapMismatchPosVec[tmp];
				if(tmpMismatchPos <= bestSplice_doner_length + startLocInRead - 1)
				{
					bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
				}
			}

			for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
			{
				int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
				if(tmpMismatchPos > bestSplice_doner_length + startLocInRead - 1)
				{
					bestSpliceMismatchPosVec.push_back(tmpMismatchPos);
				}
			}
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

		string returnSpliceSiteVecStr()
		{
			string tmpStr = "spliceSiteVec:\n";
			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				tmpStr = tmpStr + "doner_length: " + int_to_str(spliceSiteVec[tmp].first) + " mismatch#: " 
					+ int_to_str((spliceSiteVec[tmp].second).first) + " flank_string_case: " 
					+ int_to_str((spliceSiteVec[tmp].second).second) + "\n";
			}
			return tmpStr;
		}
		int returnComparedPenalty(int mismatchPenalty,
					int semiCanonical_penalty, int nonCanonical_penalty)
		{
			int tmpMismatchNum_bestSplice = this->returnBestSplice_mismatchNum();
			int tmpFlankStringCase_bestSplice = this->returnBestSplice_flankString();
			int tmpMismatch_penalty = tmpMismatchNum_bestSplice *  mismatchPenalty;
			int tmpSpliceTypePenalty;
			if(tmpFlankStringCase_bestSplice >= 5)
			{
				tmpSpliceTypePenalty = 0;
			}
			else if(tmpFlankStringCase_bestSplice == 0)
			{
				tmpSpliceTypePenalty = nonCanonical_penalty;
			}
			else
				tmpSpliceTypePenalty = semiCanonical_penalty;
			return tmpMismatch_penalty + tmpSpliceTypePenalty;
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
		int returnBestSplice_mismatchNum()
		{
			int bestSJ_mismatchNum = spliceSiteVec[bestSpliceSiteIndex].second.first;
			return bestSJ_mismatchNum;
		}

		int returnBestSplice_prefixMatchLength()
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

		bool fixSpliceResultConfident()
		{
			if((this->returnBestSplice_mismatchNum() == 0)
				||((this->returnBestSplice_canonicalOrNot())&&(this->returnBestSplice_mismatchNum() == 1)))
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		/*bool detectBestSpliceSite_canonicalOnly_lessMismatch(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
				int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
				const string& readSeq_inProcess, 
				Index_Info* indexInfo, int chrNameInt, 
				int max_allowed_mismatchNum,
				bool annotation_provided_bool, bool Do_annotation_only_bool,
				Annotation_Info* annotationInfo //int* mismatchNum_detected, int* prefixMatch_length
				)	
		{
			//cout << "start detectBestSpliceSite_prefer_canonical_lessMismatch function ..." << endl;
			//cout << "toFixSeqLocInRead_start: " << toFixSeqLocInRead_start << endl 
			//	<< "toFixSeqLocInRead_end: " << toFixSeqLocInRead_end << endl
			//	<< "toFixSeqMapPos_doner: " << toFixSeqMapPos_doner << endl 
			//	<< "toFixSeqMapPos_acceptor: " << toFixSeqMapPos_acceptor << endl
			//	<< "chrNameInt: " << chrNameInt << endl
			//	<< "max_allowed_mismatchNum: " << max_allowed_mismatchNum << endl; 
			this->scanGenomeAndReadSeq(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, 
				readSeq_inProcess, indexInfo, chrNameInt, 
				max_allowed_mismatchNum //int* mismatchNum_detected, int* prefixMatch_length
				);	
			this->generateSpliceSiteVec(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				indexInfo, chrNameInt, max_allowed_mismatchNum,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
			if(annotation_provided_bool && Do_annotation_only_bool)
			{
				bestSpliceSiteIndex = this->selectBestSpliceSite_prefer_canonical_lessMismatch();
				//this->selectBestSpliceSite_canonicalOnly_lessMismatch();
			}
			else
			{
				bestSpliceSiteIndex = this->selectBestSpliceSite_canonicalOnly_lessMismatch();
			}

			if(bestSpliceSiteIndex >= 0)
				return true;
			else
				return false;
		}*/

		bool detectBestSpliceSite_withoutAnnotation(//_prefer_canonical_lessMismatch
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
				this->generateBestSpliceMismatchVec(toFixSeqLocInRead_end-toFixSeqLocInRead_start+1, toFixSeqLocInRead_start);
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
		
		bool detectBestSpliceSite(//_prefer_canonical_lessMismatch
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt, 
			int max_allowed_mismatchNum,
			bool annotation_provided_bool, 
			bool Do_annotation_only_bool, Annotation_Info* annotationInfo)	
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
			
			this->generateSpliceSiteVec(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				indexInfo, chrNameInt, max_allowed_mismatchNum,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);

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
				this->generateBestSpliceMismatchVec(toFixSeqLocInRead_end-toFixSeqLocInRead_start+1, toFixSeqLocInRead_start);
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
		/*
		int selectBestSpliceSite_prefer_canonical_lessMismatch_new()  // canonical - semicanonical - noncanonical lessmismatch
		{
			int currentBestSJ_index_canonical = -1;
			int currentBestSJ_index_semicanonical = -1;
			int currentBestSJ_index_noncanonical = -1;

			int currentBestSJ_mismatchNum_canonical = 1000000;
			int currentBestSJ_mismatchNum_semicanonical = 1000000;
			int currentBestSJ_mismatchNum_noncanonical = 1000000;

			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				bool tmpSJ_canonical = this->SJ_canonicalOrNot(spliceSiteVec[tmp].second.second);
				if(tmpSJ_canonical)
				{
					int currentMismatchNum = spliceSiteVec[tmp].second.first;
					if(currentMismatchNum < currentBestSJ_mismatchNum_canonical)
					{
						currentBestSJ_index_canonical = tmp;
						currentBestSJ_mismatchNum_canonical = currentMismatchNum;						
					}
				}
				else
				{}
			}

			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				bool tmpSJ_semicanonical = this->SJ_semicanonicalOrNot(spliceSiteVec[tmp].second.second);
				if(tmpSJ_semicanonical)
				{
					int currentMismatchNum = spliceSiteVec[tmp].second.first;
					if(currentMismatchNum < currentBestSJ_mismatchNum_semicanonical)
					{
						currentBestSJ_index_semicanonical = tmp;
						currentBestSJ_mismatchNum_semicanonical = currentMismatchNum;						
					}
				}
				else
				{}
			}

			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				bool tmpSJ_noncanonical = this->SJ_noncanonicalOrNot(spliceSiteVec[tmp].second.second);
				if(tmpSJ_noncanonical)
				{
					int currentMismatchNum = spliceSiteVec[tmp].second.first;
					if(currentMismatchNum < currentBestSJ_mismatchNum_noncanonical)
					{
						currentBestSJ_index_noncanonical = tmp;
						currentBestSJ_mismatchNum_noncanonical = currentMismatchNum;						
					}
				}
				else
				{}
			}

			if(currentBestSJ_index_canonical >= 0)
			{
				return currentBestSJ_index_canonical;
			}
			else if(currentBestSJ_index_semicanonical >= 0)
			{
				return currentBestSJ_index_semicanonical;
			}
			else
			{
				return currentBestSJ_index_noncanonical;
			}
		}*/

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
		/*
		int selectBestSpliceSite_prefer_canonical_lessMismatch()
		{
			//cout << "start to select best" << endl;
			//string currentBestSJ_flankString;
			bool currentBestSJ_canonical = false;
			int currentBestSJ_mismatchNum = 1000000;
			int currentBestSJ_index = -1;
			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				//int tmpFlankStringCase = spliceSiteVec[tmp].second.second;
				bool tmpSJ_canonical = this->SJ_canonicalOrNot(spliceSiteVec[tmp].second.second);
				int currentMismatchNum = spliceSiteVec[tmp].second.first;
				if(tmpSJ_canonical)
				{
					if(currentBestSJ_canonical)
					{
						if(currentMismatchNum < currentBestSJ_mismatchNum)
						{
							currentBestSJ_index = tmp;
							currentBestSJ_mismatchNum = currentMismatchNum;
						}
						else
						{
							//continue;
						}
					}
					else
					{
						currentBestSJ_index = tmp;
						currentBestSJ_mismatchNum = currentMismatchNum;
						currentBestSJ_canonical = true;
					}					
				}
				else // tmp_SJ -- not canonical
				{
					if(currentBestSJ_canonical)
					{
						//continue;
					}
					else
					{
						if(currentMismatchNum < currentBestSJ_mismatchNum)
						{
							currentBestSJ_index = tmp;
							currentBestSJ_mismatchNum = currentMismatchNum;
						}
						else if(currentMismatchNum == currentBestSJ_mismatchNum)
						{
							if(this->SJ_semicanonicalOrNot(spliceSiteVec[tmp].second.second))
							{
								currentBestSJ_index = tmp;
								currentBestSJ_mismatchNum = currentMismatchNum;
							}
							else
							{}
						}
						else
						{}						
					}
				}
			}
			//cout << "currentBestSJ_index: " << currentBestSJ_index << endl;
			return currentBestSJ_index;
		}

		int selectBestSpliceSite_canonicalOnly_lessMismatch()
		{
			//cout << "start to select best" << endl;
			//string currentBestSJ_flankString;
			bool currentBestSJ_canonical = false;
			int currentBestSJ_mismatchNum = 1000000;
			int currentBestSJ_index = -1;
			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp++)
			{
				//int tmpFlankStringCase = spliceSiteVec[tmp].second.second;
				bool tmpSJ_canonical = this->SJ_canonicalOrNot(spliceSiteVec[tmp].second.second);
				int currentMismatchNum = spliceSiteVec[tmp].second.first;
				if(tmpSJ_canonical)
				{
					if(currentMismatchNum < currentBestSJ_mismatchNum)
					{
						currentBestSJ_index = tmp;
						currentBestSJ_mismatchNum = currentMismatchNum;
					}				
				}
				else // tmp_SJ -- not canonical
				{}
			}
			//cout << "currentBestSJ_index: " << currentBestSJ_index << endl;
			return currentBestSJ_index;
		}*/

		void generateSpliceSiteVec_noAnnotationProvided(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum)
		{
			//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;

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
					string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
					int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
					spliceSiteVec.push_back(pair< int, pair<int, int> >(
						tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
				}
			}
			return;
		}


		void generateSpliceSiteVec_doAnnotationOnly(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum, 
			Annotation_Info* annotationInfo)
		{
			//cout << endl << "start to generateSpliceSiteVec: " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;

			for(int tmp_forwared_doner = 0; tmp_forwared_doner <= max_extension_forward_doner; tmp_forwared_doner ++ )
			{
				int tmp_backward_acceptor = toFixSeqLength - tmp_forwared_doner;
				if(tmp_backward_acceptor > max_extension_backward_acceptor)
					continue;

				// search in annotated SJs
				int tmpChromPos_doner = toFixSeqMapPos_doner - 1 + tmp_forwared_doner;
				int tmpChromPos_acceptor = toFixSeqMapPos_acceptor + 1 - tmp_backward_acceptor;				
				bool foundInAnnotation_bool = annotationInfo->SJfoundInAnnotation(chrNameInt, tmpChromPos_doner,
					tmpChromPos_acceptor);
				if(!foundInAnnotation_bool)
					continue;

				int index_donerMapCumulativeMismatchNumVec = tmp_forwared_doner;
				int index_acceptorMapCumulativeMismatchNumVec = tmp_backward_acceptor;
				int tmp_doner_mismatch = donerMapCumulativeMismatchNumVec[index_donerMapCumulativeMismatchNumVec];
				int tmp_acceptor_mismatch = acceptorMapCumulativeMismatchNumVec[index_acceptorMapCumulativeMismatchNumVec];
				int tmp_mismatch_sum = tmp_doner_mismatch + tmp_acceptor_mismatch;

				if(tmp_mismatch_sum <= max_allowed_mismatchNum)
				{
					string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
					int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
					spliceSiteVec.push_back(pair< int, pair<int, int> >(
						tmp_forwared_doner, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
				}
			}
			return;
		}


		void generateSpliceSiteVec(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			Index_Info* indexInfo, int chrNameInt, int max_allowed_mismatchNum, 
			bool annotation_provided_bool, bool Do_annotation_only_bool, 
			Annotation_Info* annotationInfo)
		{
			if(annotation_provided_bool)
			{
				if(Do_annotation_only_bool)
				{
					this->generateSpliceSiteVec_doAnnotationOnly(
						toFixSeqLocInRead_start, toFixSeqLocInRead_end,
						toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
						indexInfo, chrNameInt, max_allowed_mismatchNum, annotationInfo);
				}
				else
				{}
			}	
			else
			{
				this->generateSpliceSiteVec_noAnnotationProvided(
					toFixSeqLocInRead_start, toFixSeqLocInRead_end,
					toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
					indexInfo, chrNameInt, max_allowed_mismatchNum);
			}
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
			for(int tmp = 0; tmp < toFixSeqLength; tmp++)
			{
				char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_doner + tmp);
				if(readSeq_inProcess.at(toFixSeqLocInRead_start - 1 + tmp) 
					!= charInRef )
				{
					tmpDonerMapCumulativeMismatchNum ++;
					if(STORE_MISMATCH_POS)
					{
						donerMapMismatchPosVec.push_back(tmp+toFixSeqLocInRead_start);
						if(STORE_MISMATCH_CHA)
						{
							donerMapMismatchCharVec.push_back(charInRef);
						}
					}	
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
				char charInRef = indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_acceptor - tmp);
				if(readSeq_inProcess.at(toFixSeqLocInRead_end - tmp - 1) 
					!= charInRef )
				{
					tmpAcceptorMapCumulativeMismatchNum ++;
					if(STORE_MISMATCH_POS)
					{
						acceptorMapMismatchPosVec.push_back(toFixSeqLocInRead_end - tmp);
						if(STORE_MISMATCH_CHA)
						{
							acceptorMapMismatchCharVec.push_back(charInRef);
						}
					}
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
};

#endif