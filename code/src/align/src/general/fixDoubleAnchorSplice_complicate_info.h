// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHORSPLICE_COMPLICATE_INFO_H
#define FIXDOUBLEANCHORSPLICE_COMPLICATE_INFO_H

#include "nw_DP.h"
#include "fixDoubleAnchorInsertion_info.h"
#include "fixDoubleAnchorDeletion_info.h"
#include "fixDoubleAnchorSplice_info.h"

class FixDoubleAnchor_Splice_Complicate_Info
{
private:
		int bestComplicateSJ_index;

		int prefix_match_length_best;
		
		int first_jumpCode_length_best;
		string first_jumpCode_type_best;
		 
		int mid_match_length_best;
		
		int second_jumpCode_length_best;
		string second_jumpCode_type_best;
		
		int suffix_match_length_best;
			
		int mismatch_best;
			
		vector<int> GT_splice_doner_vec; // < prefix_match_length in chrom >
		vector<int> AG_splice_acceptor_vec; // < suffix_match_length in chrom >

		vector<int> CT_splice_doner_vec; // < prefix_match_length in chrom >
		vector<int> AC_splice_acceptor_vec; // < suffix_match_length in chrom >
		
		vector< pair < pair<int, int>, pair<int, int> > > validate_paired_SJ_vec; // < prefix_match_length, suffix_match_length >, <indel_len, mismatch#>

		vector< pair< pair< int, pair<int,string> >, pair< int, pair<int,string> > > > candi_complicate_SJ_vec;
		// <prefix_match_length, <first_jumpCode_length, first_jumpCode_type> >, <mid_match_length, <second_jumpCode_length, second_jumpCode_type> >
		vector< vector<int> > candi_complicate_SJ_vec_mismatchPos;
		vector< vector<char> > candi_complicate_SJ_vec_mismatchChar;
		vector< pair<int, int> > candi_complicate_SJ_penalty_vec; // <indel_length, mismatchNum>

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


		//int max_extension_forward_doner;
		//int max_extension_backward_acceptor;

		int deletion_penalty_weight;
		int insertion_penalty_weight;
		int mismatch_penalty_weight;

		vector<int> bestComplicatedSpliceMismatchPosVec;
		vector<char> bestComplicatedSpliceMismatchCharVec;

		int penalty_best;
		//int penalty_score_max;

		int mismatchPos_interval_min;
public:
		void copyMismatchPos2TargetVec(vector<int>& targetMismatchPosVec)
		{
			for(int tmp = 0; tmp < bestComplicatedSpliceMismatchPosVec.size(); tmp++)
			{
				targetMismatchPosVec.push_back(bestComplicatedSpliceMismatchPosVec[tmp]);
			}
		}
		void copyMismatchChar2TargetVec(vector<char>& targetMismatchCharVec)
		{
			for(int tmp = 0; tmp < bestComplicatedSpliceMismatchCharVec.size(); tmp++)
			{
				targetMismatchCharVec.push_back(bestComplicatedSpliceMismatchCharVec[tmp]);
			}
		}

		int returnBestComplicateSpliceMismatchPosVecSize()
		{
			return bestComplicatedSpliceMismatchPosVec.size();
		}
		int returnBestComplicateSpliceMismatchCharVecSize()
		{
			return bestComplicatedSpliceMismatchCharVec.size();
		}

		int returnMismatchPosInRead(int index_bestSpliceMismatchPosVec)
		{
			return bestComplicatedSpliceMismatchPosVec[index_bestSpliceMismatchPosVec];
		}
		char returnMismatchChar(int index_bestSpliceMismatchCharVec)
		{
			return bestComplicatedSpliceMismatchCharVec[index_bestSpliceMismatchCharVec];
		}

		int return_prefix_match_length_best()
		{
			return prefix_match_length_best;
		}
		int return_first_jumpCode_length_best()
		{
			return first_jumpCode_length_best;
		}
		string return_first_jumpCode_type_best()
		{
			return first_jumpCode_type_best;
		}
		int return_mid_match_length_best()
		{
			return mid_match_length_best;
		}
		int return_second_jumpCode_length_best()
		{
			return second_jumpCode_length_best;
		}
		string return_second_jumpCode_type_best()
		{
			return second_jumpCode_type_best;
		}
		int return_suffix_match_length_best()
		{
			return suffix_match_length_best;
		}
		int return_mismatch_bestComplicatedSplice()
		{
			return mismatch_best;
		}
		bool selectBestComplicateSJ_lessPenalty(int toFix_read_start, int toFix_read_end, int max_allowed_mismatchNum)
		{
			int toFix_read_length = toFix_read_end - toFix_read_start + 1;
			bestComplicateSJ_index = -1;
			int tmp_best_complicate_SJ_penalty = 100;
			for(int tmp = 0; tmp < candi_complicate_SJ_penalty_vec.size(); tmp ++)
			{
				int tmp_penalty_score 
					= this->returnPenalty(candi_complicate_SJ_penalty_vec[tmp].first, candi_complicate_SJ_penalty_vec[tmp].second);
				if(tmp_penalty_score <= tmp_best_complicate_SJ_penalty)
				{
					bestComplicateSJ_index = tmp;
					tmp_best_complicate_SJ_penalty = tmp_penalty_score;
				}
			}
			if(bestComplicateSJ_index >= 0)
			{
				prefix_match_length_best = (candi_complicate_SJ_vec[bestComplicateSJ_index].first).first;
		
				first_jumpCode_length_best = ((candi_complicate_SJ_vec[bestComplicateSJ_index].first).second).first;
				first_jumpCode_type_best = ((candi_complicate_SJ_vec[bestComplicateSJ_index].first).second).second;
				 
				mid_match_length_best = (candi_complicate_SJ_vec[bestComplicateSJ_index].second).first;
		
				second_jumpCode_length_best = ((candi_complicate_SJ_vec[bestComplicateSJ_index].second).second).first;
				second_jumpCode_type_best = ((candi_complicate_SJ_vec[bestComplicateSJ_index].second).second).second;
		
				int tmp_indel_length = candi_complicate_SJ_penalty_vec[bestComplicateSJ_index].first;
				if(tmp_indel_length < 0) // insertion
					suffix_match_length_best = toFix_read_length - prefix_match_length_best
						- mid_match_length_best + tmp_indel_length;
				else
					suffix_match_length_best = toFix_read_length - prefix_match_length_best
						- mid_match_length_best;

				penalty_best = tmp_best_complicate_SJ_penalty;
				mismatch_best = candi_complicate_SJ_penalty_vec[bestComplicateSJ_index].second;
				
				if(penalty_best > max_allowed_mismatchNum * mismatch_penalty_weight)
				{
					return false;
				}			
				else
					return true;
			}
			else
			{
				return false;
			}
		}
		void generateBestComplicateSJMismatchVec(int toFixSeqStartLocInRead, int toFixSeqEndLocInRead)
		{
			if(STORE_MISMATCH_CHA)
				this->generateBestComplicateSJMismatchVec_Pos_Char(toFixSeqStartLocInRead, toFixSeqEndLocInRead);
			else
				this->generateBestComplicateSJMismatchVec_Pos(toFixSeqStartLocInRead, toFixSeqEndLocInRead);
		}
		void generateBestComplicateSJMismatchVec_Pos(int toFixSeqStartLocInRead, int toFixSeqEndLocInRead)
		{

			if(second_jumpCode_type_best == "N") // indel in the left part
			{
				// push back indel mismatch pos vec
				for(int tmp = 0; tmp < candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index].size(); tmp++)
				{
					bestComplicatedSpliceMismatchPosVec.push_back((candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index])[tmp]);
				}

				// push back SJ right part mismatch pos vec;
				int validMismatchPos_acceptor_start = toFixSeqEndLocInRead - suffix_match_length_best + 1;
				for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
				{
					int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
					if(tmpMismatchPos >= validMismatchPos_acceptor_start)
					{
						bestComplicatedSpliceMismatchPosVec.push_back(tmpMismatchPos);
					}
				}
			}
			else if(first_jumpCode_type_best == "N")// indel in the right part
			{
				// push back SJ left part mismatch pos vec;
				int validMismatchPos_doner_end = toFixSeqStartLocInRead + prefix_match_length_best - 1;
				for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
				{
					int tmpMismatchPos = donerMapMismatchPosVec[tmp];
					if(tmpMismatchPos <= validMismatchPos_doner_end)
					{
						bestComplicatedSpliceMismatchPosVec.push_back(tmpMismatchPos);
					}
				}

				// push back indel mismatch pos vec;
				for(int tmp = 0; tmp < candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index].size(); tmp++)
				{
					bestComplicatedSpliceMismatchPosVec.push_back((candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index])[tmp]);
				}
			}
			else
			{
				cout << "error in generateBestComplicateSJMismatchVec_Pos in FixDoubleAnchor_Splice_Complicate_Info.h" << endl;
				exit(1);
			}
		}
		void generateBestComplicateSJMismatchVec_Pos_Char(int toFixSeqStartLocInRead, int toFixSeqEndLocInRead)
		{
			if(second_jumpCode_type_best == "N") // indel in the left part
			{
				// push back indel mismatch pos vec
				for(int tmp = 0; tmp < candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index].size(); tmp++)
				{
					bestComplicatedSpliceMismatchPosVec.push_back((candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index])[tmp]);
					bestComplicatedSpliceMismatchCharVec.push_back((candi_complicate_SJ_vec_mismatchChar[bestComplicateSJ_index])[tmp]);
				}

				// push back SJ right part mismatch pos vec;
				int validMismatchPos_acceptor_start = toFixSeqEndLocInRead - suffix_match_length_best + 1;
				for(int tmp = 0; tmp < acceptorMapMismatchPosVec.size(); tmp++)
				{
					int tmpMismatchPos = acceptorMapMismatchPosVec[tmp];
					char tmpMismatchChar = acceptorMapMismatchCharVec[tmp];
					if(tmpMismatchPos >= validMismatchPos_acceptor_start)
					{
						bestComplicatedSpliceMismatchPosVec.push_back(tmpMismatchPos);
						bestComplicatedSpliceMismatchCharVec.push_back(tmpMismatchChar);
					}
				}
			}
			else if(first_jumpCode_type_best == "N")// indel in the right part
			{
				// push back SJ left part mismatch pos vec;
				int validMismatchPos_doner_end = toFixSeqStartLocInRead + prefix_match_length_best - 1;
				for(int tmp = 0; tmp < donerMapMismatchPosVec.size(); tmp++)
				{
					int tmpMismatchPos = donerMapMismatchPosVec[tmp];
					char tmpMismatchChar = donerMapMismatchCharVec[tmp];
					if(tmpMismatchPos <= validMismatchPos_doner_end)
					{
						bestComplicatedSpliceMismatchPosVec.push_back(tmpMismatchPos);
						bestComplicatedSpliceMismatchCharVec.push_back(tmpMismatchChar);
					}
				}

				// push back indel mismatch pos vec;
				for(int tmp = 0; tmp < candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index].size(); tmp++)
				{
					bestComplicatedSpliceMismatchPosVec.push_back((candi_complicate_SJ_vec_mismatchPos[bestComplicateSJ_index])[tmp]);
					bestComplicatedSpliceMismatchCharVec.push_back((candi_complicate_SJ_vec_mismatchChar[bestComplicateSJ_index])[tmp]);
				}
			}
			else
			{
				cout << "error in generateBestComplicateSJMismatchVec_Pos_Char in FixDoubleAnchor_Splice_Complicate_Info.h" << endl;
				exit(1);
			}
		}
		int returnPenalty(int indel_length, int mismatch)
		{
			if(indel_length < 0)
			{
				return ( (0-indel_length) * insertion_penalty_weight + mismatch * mismatch_penalty_weight);
			}
			else
				return ( indel_length * deletion_penalty_weight + mismatch * mismatch_penalty_weight);
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
		string returnSpliceSignalVecStr()
		{
			string tmpStr = "spliceSignalVec:\nGT: ";
			for(int tmp = 0; tmp < GT_splice_doner_vec.size(); tmp ++)
			{
				tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(GT_splice_doner_vec[tmp]) + ", ";
			}
			tmpStr += "\nAG: ";
			for(int tmp = 0; tmp < AG_splice_acceptor_vec.size(); tmp ++)
			{
				tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(AG_splice_acceptor_vec[tmp]) + ", ";
			}
			tmpStr += "\nCT: ";
			for(int tmp = 0; tmp < CT_splice_doner_vec.size(); tmp ++)
			{
				tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(CT_splice_doner_vec[tmp]) + ", ";
			}
			tmpStr += "\nAC: ";
			for(int tmp = 0; tmp < AC_splice_acceptor_vec.size(); tmp ++)
			{
				tmpStr = tmpStr + int_to_str(tmp) + " " + int_to_str(AC_splice_acceptor_vec[tmp]) + ", ";
			}		
			return tmpStr;				
		}
		string returnValidatePairedSJvecStr()
		{
			string tmpStr = "valid Paired SJ vec: \n";
			for(int tmp = 0; tmp < validate_paired_SJ_vec.size(); tmp++)
			{
				int a = (validate_paired_SJ_vec[tmp].first).first;
				int b = (validate_paired_SJ_vec[tmp].first).second;
				int c = (validate_paired_SJ_vec[tmp].second).first;
				int d = (validate_paired_SJ_vec[tmp].second).second;				
				
				tmpStr = tmpStr + int_to_str(a) + " " + int_to_str(b) + " "
					+ int_to_str(c) + " " + int_to_str(d) + " \n";
			}
			return tmpStr;
		}
		string returnCandidateComplicatedSJvecStr()
		{
			string tmpStr = "candidateComplicateSJvec:\n";
			for(int tmp = 0; tmp < candi_complicate_SJ_vec.size(); tmp ++)
			{
				int a = (candi_complicate_SJ_vec[tmp].first).first;
				int b = ((candi_complicate_SJ_vec[tmp].first).second).first;
				string c = ((candi_complicate_SJ_vec[tmp].first).second).second;
				int d = (candi_complicate_SJ_vec[tmp].second).first;
				int e = ((candi_complicate_SJ_vec[tmp].second).second).first;
				string f = ((candi_complicate_SJ_vec[tmp].second).second).second;
				tmpStr = tmpStr + int_to_str(a) + " " + int_to_str(b) + " "
					+ c + " " + int_to_str(d) + " " + int_to_str(e) + " " + f + "\n";
			}
			return tmpStr;
		}

		void add2CandiComplicateSJvecMismatchPosVec_ins(FixDoubleAnchor_Insertion_Info* fixInsInfo)
		{
			vector<int> newMismatchPosVec;
			int vecSize = fixInsInfo->returnBestInsertionMismatchPosVecSize();
			for(int tmp = 0; tmp < vecSize; tmp++)
			{
				int tmpMismatchPos = fixInsInfo->returnBestInsertionMismatchPos(tmp);
				newMismatchPosVec.push_back(tmpMismatchPos);	
			}	
			candi_complicate_SJ_vec_mismatchPos.push_back(newMismatchPosVec);
		}
		void add2CandiComplicatedSJvecMismatchCharVec_ins(FixDoubleAnchor_Insertion_Info* fixInsInfo)
		{
			vector<char> newMismatchCharVec;
			int vecSize = fixInsInfo->returnBestInsertionMismatchCharVecSize();
			for(int tmp = 0; tmp < vecSize; tmp++)
			{
				char tmpMismatchChar = fixInsInfo->returnBestInsertionMismatchChar(tmp);
				newMismatchCharVec.push_back(tmpMismatchChar);
			}			
			candi_complicate_SJ_vec_mismatchChar.push_back(newMismatchCharVec);
		}

		void add2CandiComplicateSJvecMismatchPosVec_del(FixDoubleAnchor_Deletion_Info* fixDelInfo)
		{
			vector<int> newMismatchPosVec;
			int vecSize = fixDelInfo->returnBestDeletionMismatchPosVecSize();
			for(int tmp = 0; tmp < vecSize; tmp++)
			{
				int tmpMismatchPos = fixDelInfo->returnBestDeletionMismatchPos(tmp);
				newMismatchPosVec.push_back(tmpMismatchPos);	
			}	
			candi_complicate_SJ_vec_mismatchPos.push_back(newMismatchPosVec);
		}
		void add2CandiComplicatedSJvecMismatchCharVec_del(FixDoubleAnchor_Deletion_Info* fixDelInfo)
		{
			vector<char> newMismatchCharVec;
			int vecSize = fixDelInfo->returnBestDeletionMismatchCharVecSize();
			for(int tmp = 0; tmp < vecSize; tmp++)
			{
				char tmpMismatchChar = fixDelInfo->returnBestDeletionMismatchChar(tmp);
				newMismatchCharVec.push_back(tmpMismatchChar);
			}			
			candi_complicate_SJ_vec_mismatchChar.push_back(newMismatchCharVec);
		}

		void get_candi_complicate_SJ_vec(int toFix_read_start, int toFix_read_end,
			int toFix_chrom_start, int toFix_chrom_end,
			const string& readSeq_inProcess, Index_Info* indexInfo,
			int chrom_name_int, int max_allowed_mismatchNum
			)
		{
			for(int tmp = 0; tmp < validate_paired_SJ_vec.size(); tmp++)
			{
				int tmp_toFix_read_seq_length = toFix_read_end - toFix_read_start + 1;

				int tmp_prefix_match_length = (validate_paired_SJ_vec[tmp].first).first;
				int tmp_suffix_match_length = (validate_paired_SJ_vec[tmp].first).second;

				int tmp_indel_length = (validate_paired_SJ_vec[tmp].second).first;
				int tmp_mismatch_num = (validate_paired_SJ_vec[tmp].second).second;
				int tmp_doner_side_mismatch_num = donerMapCumulativeMismatchNumVec[tmp_prefix_match_length];
				int tmp_acceptor_side_mismatch_num = acceptorMapCumulativeMismatchNumVec[tmp_suffix_match_length];

				int tmp_SJ_start_inChr = toFix_chrom_start + tmp_prefix_match_length;
				int tmp_SJ_end_inChr = toFix_chrom_end - tmp_suffix_match_length;
				int tmp_SJ_distance = tmp_SJ_end_inChr - tmp_SJ_start_inChr + 1; 
				if(tmp_indel_length < 0)  // insertion found
				{
					if(tmp_mismatch_num == tmp_acceptor_side_mismatch_num) // insertion in 1st part
					{
						FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info(); 
						int tmp_toFix_insertion_read_start = toFix_read_start;
						int tmp_new_prefix_length = tmp_toFix_read_seq_length - tmp_suffix_match_length;
						int tmp_toFix_insertion_read_end = toFix_read_start + tmp_new_prefix_length - 1;
						int tmp_toFix_insertion_chrom_start = toFix_chrom_start;
						int tmp_toFix_insertion_chrom_end = tmp_toFix_insertion_chrom_start + tmp_prefix_match_length - 1;
						int tmp_max_allowed_mismatchNum = max_allowed_mismatchNum - tmp_mismatch_num;

						bool fixIns_bool = insInfo->detectBestInsertion_lessMismatch(
							tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
							tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, readSeq_inProcess, indexInfo,
							chrom_name_int, tmp_max_allowed_mismatchNum);
						if( fixIns_bool && (insInfo->return_best_insertion_suffix_match_length()>0))
						{
							int tmp_best_insertion_prefix_match_length = insInfo->return_best_insertion_prefix_match_length();
							int tmp_best_insertion_length = insInfo->return_best_insertion_length();
							int tmp_best_insertion_suffix_match_length = insInfo->return_best_insertion_suffix_match_length();
							candi_complicate_SJ_vec.push_back(pair< pair<int, pair<int,string> >, pair<int, pair<int,string> > >(
								pair<int, pair<int,string> > (tmp_best_insertion_prefix_match_length, pair<int,string>(tmp_best_insertion_length, "I")),
								pair<int, pair<int,string> > (tmp_best_insertion_suffix_match_length, pair<int,string>(tmp_SJ_distance, "N"))
								));
							if(STORE_MISMATCH_POS)
							{
								//insInfo->generateBestInsertionMismatchVec(tmp_toFix_insertion_read_start);
								this->add2CandiComplicateSJvecMismatchPosVec_ins(insInfo);
								if(STORE_MISMATCH_CHA)
								{
									this->add2CandiComplicatedSJvecMismatchCharVec_ins(insInfo);
								}
							}
							candi_complicate_SJ_penalty_vec.push_back(pair<int,int>((0-tmp_best_insertion_length), insInfo->return_best_insertion_mismatch() + tmp_acceptor_side_mismatch_num));

						}
						delete insInfo;
					}
					if(tmp_mismatch_num == tmp_doner_side_mismatch_num) // insertion in 2nd part
					{
						FixDoubleAnchor_Insertion_Info* insInfo = new FixDoubleAnchor_Insertion_Info(); 
						int tmp_toFix_insertion_read_end = toFix_read_end;
						int tmp_new_suffix_length = tmp_toFix_read_seq_length - tmp_prefix_match_length;
						int tmp_toFix_insertion_read_start = tmp_toFix_insertion_read_end - tmp_new_suffix_length + 1;
						int tmp_toFix_insertion_chrom_end = toFix_chrom_end;
						int tmp_toFix_insertion_chrom_start = toFix_chrom_end - tmp_suffix_match_length + 1;
						int tmp_max_allowed_mismatchNum = max_allowed_mismatchNum - tmp_mismatch_num;

						bool fixIns_bool = insInfo->detectBestInsertion_lessMismatch(
							tmp_toFix_insertion_read_start, tmp_toFix_insertion_read_end,
							tmp_toFix_insertion_chrom_start, tmp_toFix_insertion_chrom_end, readSeq_inProcess, indexInfo,
							chrom_name_int, tmp_max_allowed_mismatchNum);
						if(fixIns_bool && (insInfo->return_best_insertion_prefix_match_length() > 0))
						{
							int tmp_best_insertion_prefix_match_length = insInfo->return_best_insertion_prefix_match_length();
							int tmp_best_insertion_length = insInfo->return_best_insertion_length();
							int tmp_best_insertion_suffix_match_length = insInfo->return_best_insertion_suffix_match_length();
							candi_complicate_SJ_vec.push_back(pair< pair<int, pair<int,string> >, pair<int, pair<int,string> > >(
								pair<int, pair<int,string> > (tmp_prefix_match_length, pair<int,string>(tmp_SJ_distance, "N")),
								pair<int, pair<int,string> > (tmp_best_insertion_prefix_match_length, pair<int,string>(tmp_best_insertion_length, "I"))
								));
							if(STORE_MISMATCH_POS)
							{
								//insInfo->generateBestInsertionMismatchVec(tmp_toFix_insertion_read_start);
								this->add2CandiComplicateSJvecMismatchPosVec_ins(insInfo);
								if(STORE_MISMATCH_CHA)
								{
									this->add2CandiComplicatedSJvecMismatchCharVec_ins(insInfo);
								}
							}
							candi_complicate_SJ_penalty_vec.push_back(pair<int,int>((0-tmp_best_insertion_length), insInfo->return_best_insertion_mismatch() + tmp_doner_side_mismatch_num));

						}
						delete insInfo;
					}
				}
				else if(tmp_indel_length > 0) // deletion found
				{
					if(tmp_mismatch_num == tmp_acceptor_side_mismatch_num) // deletion in 1st part
					{
						FixDoubleAnchor_Deletion_Info* delInfo = new FixDoubleAnchor_Deletion_Info();
						int tmp_toFix_deletion_read_start = toFix_read_start;
						int tmp_new_prefix_length = tmp_toFix_read_seq_length - tmp_suffix_match_length;
						int tmp_toFix_deletion_read_end = toFix_read_start + tmp_new_prefix_length - 1;
						int tmp_toFix_deletion_chrom_start = toFix_chrom_start;
						int tmp_toFix_deletion_chrom_end = tmp_toFix_deletion_chrom_start + tmp_prefix_match_length - 1;
						int tmp_max_allowed_mismatchNum = max_allowed_mismatchNum - tmp_mismatch_num;
						//string tmp_readSeq_inProcess = readSeq_inProcess.substr()
						bool fixDel_bool = delInfo->detectBestDeletion_lessMismatch(tmp_toFix_deletion_read_start, tmp_toFix_deletion_read_end,
							tmp_toFix_deletion_chrom_start, tmp_toFix_deletion_chrom_end, readSeq_inProcess, indexInfo,
							chrom_name_int, tmp_max_allowed_mismatchNum);
						if(fixDel_bool && (tmp_new_prefix_length - delInfo->return_best_deletion_prefix_match_length() > 0))
						{
							int tmp_best_deletion_prefix_match_length = delInfo->return_best_deletion_prefix_match_length();
							int tmp_best_deletion_length = delInfo->return_best_deletion_length();
							int tmp_best_deletion_suffix_match_length = tmp_new_prefix_length - tmp_best_deletion_prefix_match_length;
							candi_complicate_SJ_vec.push_back(pair< pair<int, pair<int,string> >, pair<int, pair<int,string> > >(
								pair<int, pair<int,string> > (tmp_best_deletion_prefix_match_length, pair<int,string>(tmp_best_deletion_length, "D")),
								pair<int, pair<int,string> > (tmp_best_deletion_suffix_match_length, pair<int,string>(tmp_SJ_distance, "N"))
								));
							if(STORE_MISMATCH_POS)
							{
								//delInfo->generateBestDeletionMismatchVec(tmp_toFix_deletion_read_start);
								this->add2CandiComplicateSJvecMismatchPosVec_del(delInfo);
								if(STORE_MISMATCH_CHA)
								{
									this->add2CandiComplicatedSJvecMismatchCharVec_del(delInfo);
								}
							}
							candi_complicate_SJ_penalty_vec.push_back(pair<int,int>(tmp_best_deletion_length, delInfo->return_best_deletion_mismatch() + tmp_acceptor_side_mismatch_num));
						}
						delete delInfo;
					}
					if(tmp_mismatch_num == tmp_doner_side_mismatch_num) // deletion in 2nd part
					{
						FixDoubleAnchor_Deletion_Info* delInfo = new FixDoubleAnchor_Deletion_Info();
						int tmp_toFix_deletion_read_end = toFix_read_end;
						int tmp_new_suffix_length = tmp_toFix_read_seq_length - tmp_prefix_match_length;
						int tmp_toFix_deletion_read_start = tmp_toFix_deletion_read_end - tmp_new_suffix_length + 1;
						int tmp_toFix_deletion_chrom_end = toFix_chrom_end;
						int tmp_toFix_deletion_chrom_start = toFix_chrom_end - tmp_suffix_match_length + 1;
						int tmp_max_allowed_mismatchNum = max_allowed_mismatchNum - tmp_mismatch_num;

						bool fixDel_bool = delInfo->detectBestDeletion_lessMismatch(tmp_toFix_deletion_read_start, tmp_toFix_deletion_read_end,
							tmp_toFix_deletion_chrom_start, tmp_toFix_deletion_chrom_end, readSeq_inProcess, indexInfo,
							chrom_name_int, tmp_max_allowed_mismatchNum);
						if(fixDel_bool && (delInfo->return_best_deletion_prefix_match_length() > 0))
						{
							int tmp_best_deletion_prefix_match_length = delInfo->return_best_deletion_prefix_match_length();
							int tmp_best_deletion_length = delInfo->return_best_deletion_length();
							int tmp_best_deletion_suffix_match_length = tmp_new_suffix_length - tmp_best_deletion_prefix_match_length;

							candi_complicate_SJ_vec.push_back(pair< pair<int, pair<int,string> >, pair<int, pair<int,string> > >(
								pair<int, pair<int,string> > (tmp_prefix_match_length, pair<int,string>(tmp_SJ_distance, "N")),
								pair<int, pair<int,string> > (tmp_best_deletion_prefix_match_length, pair<int,string>(tmp_best_deletion_length, "D"))
								));
							if(STORE_MISMATCH_POS)
							{
								//delInfo->generateBestDeletionMismatchVec(tmp_toFix_deletion_read_start);
								this->add2CandiComplicateSJvecMismatchPosVec_del(delInfo);
								if(STORE_MISMATCH_CHA)
								{
									this->add2CandiComplicatedSJvecMismatchCharVec_del(delInfo);
								}
							}
							candi_complicate_SJ_penalty_vec.push_back(pair<int,int>(tmp_best_deletion_length, delInfo->return_best_deletion_mismatch() + tmp_doner_side_mismatch_num));
						}
						delete delInfo;
					}
				}
				else
				{}
			}
		}

		void getValidateSpliceSignalPair_noAnnotationProvided(int max_allowed_mismatchNum, 
			int toFix_read_start, int toFix_read_end)
		{
			int toFixSeqLength = toFix_read_end - toFix_read_start + 1;
			for(int tmp1 = 0; tmp1 < GT_splice_doner_vec.size(); tmp1 ++)
			{
				for(int tmp2 = 0; tmp2 < AG_splice_acceptor_vec.size(); tmp2 ++)
				{
					int prefix_match_length = GT_splice_doner_vec[tmp1];
					int suffix_match_length = AG_splice_acceptor_vec[tmp2];
					int doner_side_mismatch = donerMapCumulativeMismatchNumVec[prefix_match_length];
					int acceptor_side_mismatch = acceptorMapCumulativeMismatchNumVec[suffix_match_length];
					int tmp_mismatch = doner_side_mismatch;					
					if(doner_side_mismatch >= acceptor_side_mismatch)
						tmp_mismatch = acceptor_side_mismatch;						
					if(tmp_mismatch <= max_allowed_mismatchNum)
					{
						int tmp_indel_length = this->return_indel_length(
							prefix_match_length, suffix_match_length, toFixSeqLength);
						if(tmp_indel_length < 0)
						{
							int tmp_ins_length = 0 - tmp_indel_length;
							if(tmp_ins_length <= MAX_INSERTION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 
						}
						else
						{
							if(tmp_indel_length <= MAX_DELETION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 							
						}
					}
				}
			}
			for(int tmp1 = 0; tmp1 < CT_splice_doner_vec.size(); tmp1 ++)
			{
				for(int tmp2 = 0; tmp2 < AC_splice_acceptor_vec.size(); tmp2 ++)
				{
					int prefix_match_length = CT_splice_doner_vec[tmp1];
					int suffix_match_length = AC_splice_acceptor_vec[tmp2];
					int doner_side_mismatch = donerMapCumulativeMismatchNumVec[prefix_match_length];
					int acceptor_side_mismatch = acceptorMapCumulativeMismatchNumVec[suffix_match_length];
					int tmp_mismatch = doner_side_mismatch;					
					if(doner_side_mismatch >= acceptor_side_mismatch)
						tmp_mismatch = acceptor_side_mismatch;						
					if(tmp_mismatch <= max_allowed_mismatchNum)
					{
						int tmp_indel_length = this->return_indel_length(
							prefix_match_length, suffix_match_length, toFixSeqLength);
						if(tmp_indel_length < 0)
						{
							int tmp_ins_length = 0 - tmp_indel_length;
							if(tmp_ins_length <= MAX_INSERTION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 
						}
						else
						{
							if(tmp_indel_length <= MAX_DELETION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 							
						}
					}
				}
			}
		}

		void getValidateSpliceSignalPair_doAnnotationOnly(int chrNameInt, int max_allowed_mismatchNum, 
			int toFix_read_start, int toFix_read_end,
			int toFix_chrom_start, int toFix_chrom_end, Annotation_Info* annotationInfo)
		{
			int toFixSeqLength = toFix_read_end - toFix_read_start + 1;
			for(int tmp1 = 0; tmp1 < GT_splice_doner_vec.size(); tmp1 ++)
			{
				for(int tmp2 = 0; tmp2 < AG_splice_acceptor_vec.size(); tmp2 ++)
				{
					int prefix_match_length = GT_splice_doner_vec[tmp1];
					int suffix_match_length = AG_splice_acceptor_vec[tmp2];

					// search in annotated SJs
					int tmpChromPos_doner = toFix_chrom_start - 1 + prefix_match_length;
					int tmpChromPos_acceptor = toFix_chrom_end + 1 - suffix_match_length;				
					bool foundInAnnotation_bool = annotationInfo->SJfoundInAnnotation(chrNameInt,
						tmpChromPos_doner, tmpChromPos_acceptor);
					if(!foundInAnnotation_bool)
						continue;

					int doner_side_mismatch = donerMapCumulativeMismatchNumVec[prefix_match_length];
					int acceptor_side_mismatch = acceptorMapCumulativeMismatchNumVec[suffix_match_length];
					int tmp_mismatch = doner_side_mismatch;					
					if(doner_side_mismatch >= acceptor_side_mismatch)
						tmp_mismatch = acceptor_side_mismatch;						
					if(tmp_mismatch <= max_allowed_mismatchNum)
					{
						int tmp_indel_length = this->return_indel_length(
							prefix_match_length, suffix_match_length, toFixSeqLength);
						if(tmp_indel_length < 0)
						{
							int tmp_ins_length = 0 - tmp_indel_length;
							if(tmp_ins_length <= MAX_INSERTION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 
						}
						else
						{
							if(tmp_indel_length <= MAX_DELETION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 							
						}
					}
				}
			}
			for(int tmp1 = 0; tmp1 < CT_splice_doner_vec.size(); tmp1 ++)
			{
				for(int tmp2 = 0; tmp2 < AC_splice_acceptor_vec.size(); tmp2 ++)
				{
					int prefix_match_length = CT_splice_doner_vec[tmp1];
					int suffix_match_length = AC_splice_acceptor_vec[tmp2];

					// search in annotated SJs
					int tmpChromPos_doner = toFix_chrom_start - 1 + prefix_match_length;
					int tmpChromPos_acceptor = toFix_chrom_end + 1 - suffix_match_length;				
					bool foundInAnnotation_bool = annotationInfo->SJfoundInAnnotation(chrNameInt,
						tmpChromPos_doner, tmpChromPos_acceptor);
					if(!foundInAnnotation_bool)
						continue;

					int doner_side_mismatch = donerMapCumulativeMismatchNumVec[prefix_match_length];
					int acceptor_side_mismatch = acceptorMapCumulativeMismatchNumVec[suffix_match_length];
					int tmp_mismatch = doner_side_mismatch;					
					if(doner_side_mismatch >= acceptor_side_mismatch)
						tmp_mismatch = acceptor_side_mismatch;						
					if(tmp_mismatch <= max_allowed_mismatchNum)
					{
						int tmp_indel_length = this->return_indel_length(
							prefix_match_length, suffix_match_length, toFixSeqLength);
						if(tmp_indel_length < 0)
						{
							int tmp_ins_length = 0 - tmp_indel_length;
							if(tmp_ins_length <= MAX_INSERTION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 
						}
						else
						{
							if(tmp_indel_length <= MAX_DELETION_LENGTH)
								validate_paired_SJ_vec.push_back(pair< pair<int,int>, pair<int,int> > (
									pair<int,int>(prefix_match_length, suffix_match_length), pair<int,int> (tmp_indel_length, tmp_mismatch))); 							
						}
					}
				}
			}
		}

		void getValidateSpliceSignalPair(int chrNameInt, int max_allowed_mismatchNum,
			int toFix_read_start, int toFix_read_end,
			int toFix_chrom_start, int toFix_chrom_end, 
			bool annotation_provided_bool, bool Do_annotation_only_bool,
			Annotation_Info* annotationInfo)
		{
			if(annotation_provided_bool)
			{
				if(Do_annotation_only_bool)
				{
					this->getValidateSpliceSignalPair_doAnnotationOnly(chrNameInt, 
						max_allowed_mismatchNum, 
						toFix_read_start, toFix_read_end, toFix_chrom_start, toFix_chrom_end,
						annotationInfo);
				}
				else
				{}
			}
			else
			{
				this->getValidateSpliceSignalPair_noAnnotationProvided(max_allowed_mismatchNum, 
					toFix_read_start, toFix_read_end);
			}
		}

		int return_indel_length(int prefix_match_length, int suffix_match_length, int toFix_seq_length)
		{
			return (prefix_match_length + suffix_match_length - toFix_seq_length);
		}

		FixDoubleAnchor_Splice_Complicate_Info()
		{
			deletion_penalty_weight = 3;
			insertion_penalty_weight = 4;
			mismatch_penalty_weight = 2;			
		
			bestComplicateSJ_index = -1;

			mismatchPos_interval_min = 3;
		}

		int getMinMismatchInterval(vector<int>& bestComplicatedSpliceMismatchPosVec_sorted)
		{
				for(int tmp = 0; tmp < bestComplicatedSpliceMismatchPosVec.size(); tmp++)
				{
					bestComplicatedSpliceMismatchPosVec_sorted.push_back(
						bestComplicatedSpliceMismatchPosVec[tmp]);
				}				
				sort(bestComplicatedSpliceMismatchPosVec_sorted.begin(), 
					bestComplicatedSpliceMismatchPosVec_sorted.end());

				int tmpMinMismatchPosInterval = 100;
				for(int tmp = 0; tmp < bestComplicatedSpliceMismatchPosVec_sorted.size()-1; tmp++)
				{
					int tmpIndex_1 = tmp;
					int tmpIndex_2 = tmp + 1;
					int tmpMismatchPos_1 = bestComplicatedSpliceMismatchPosVec_sorted[tmpIndex_1];
					int tmpMismatchPos_2 = bestComplicatedSpliceMismatchPosVec_sorted[tmpIndex_2];
					int tmpMismatchPosInterval = tmpMismatchPos_2 - tmpMismatchPos_1 - 1;
					if(tmpMismatchPosInterval < tmpMinMismatchPosInterval)
						tmpMinMismatchPosInterval = tmpMismatchPosInterval;
				}
				return tmpMinMismatchPosInterval;
		}
		int getMinMismatchInterval_2()// min mismatch interval beside the minimum one
		{
			vector<int> bestComplicatedSpliceMismatchPosVec_sorted;
			int minimumMismatchGap = this->getMinMismatchInterval(bestComplicatedSpliceMismatchPosVec_sorted);
			int tmpMinMismatchPosInterval_2 = 100;
			bool checkedMinimumGap_bool = false;
			for(int tmp = 0; tmp < bestComplicatedSpliceMismatchPosVec_sorted.size()-1; tmp++)
			{
				int tmpIndex_1 = tmp;
				int tmpIndex_2 = tmp + 1;
				int tmpMismatchPos_1 = bestComplicatedSpliceMismatchPosVec_sorted[tmpIndex_1];
				int tmpMismatchPos_2 = bestComplicatedSpliceMismatchPosVec_sorted[tmpIndex_2];
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
			int mismatchPosVecSize = bestComplicatedSpliceMismatchPosVec.size();
			if(mismatchPosVecSize >= 3)
			{
				int tmpMinMismatchPosInterval = this->getMinMismatchInterval_2();
				if(tmpMinMismatchPosInterval < mismatchPos_interval_min)
					return false;
			}

			// int tmpSJcase = (spliceSiteVec[bestSpliceSiteIndex].second).second;
			// if(tmpSJcase < 5)
			// {
			// 	if(mismatchPosVecSize > maxMismatchNum_semi_noncanonical)
			// 		return false;
			// }
			return true;
		}

		bool detectComplicateSplice(int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, Index_Info* indexInfo, int chrNameInt, 
			int max_allowed_mismatchNum, bool annotation_provided_bool, bool Do_annotation_only_bool,
			Annotation_Info* annotationInfo)	
		{
			//cout << "start to detectComplicateSplice " << endl;
			this->scanGenomeAndReadSeq(toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				readSeq_inProcess, indexInfo, chrNameInt);
			//cout << "scan genome results: " << endl << this->returnDonerAndAcceptorMapCumulativeMismatchNumVecStr() << endl << endl;
			this->searchSpliceSignal(toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, indexInfo, chrNameInt);
			//cout << "searchSpliceSignal: " << endl << this->returnSpliceSignalVecStr() << endl << endl;
			this->getValidateSpliceSignalPair(chrNameInt, max_allowed_mismatchNum, toFixSeqLocInRead_start, 
				toFixSeqLocInRead_end, toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				annotation_provided_bool, Do_annotation_only_bool, annotationInfo);
			//cout << "validateSplciePair: " << endl << this->returnValidatePairedSJvecStr() << endl << endl;
			this->get_candi_complicate_SJ_vec(toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				readSeq_inProcess, indexInfo, chrNameInt, max_allowed_mismatchNum);
			//cout << "complicatedSJvec: " << endl << this->returnCandidateComplicatedSJvecStr() << endl << endl;
			bool selectBest_bool = this->selectBestComplicateSJ_lessPenalty(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end, max_allowed_mismatchNum);
			
			if(selectBest_bool)
			{
				this->generateBestComplicateSJMismatchVec(toFixSeqLocInRead_start, toFixSeqLocInRead_end);
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
			//return selectBest_bool;
		}

		bool detectComplicateSplice_withoutAnnotation(int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, Index_Info* indexInfo, int chrNameInt, 
			int max_allowed_mismatchNum)	
		{
			//cout << "start to detectComplicateSplice " << endl;
			this->scanGenomeAndReadSeq(toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				readSeq_inProcess, indexInfo, chrNameInt);
			//cout << "scan genome results: " << endl << this->returnDonerAndAcceptorMapCumulativeMismatchNumVecStr() << endl << endl;
			this->searchSpliceSignal(toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, indexInfo, chrNameInt);
			//cout << "searchSpliceSignal: " << endl << this->returnSpliceSignalVecStr() << endl << endl;
			this->getValidateSpliceSignalPair_noAnnotationProvided(max_allowed_mismatchNum, 
				toFixSeqLocInRead_start, toFixSeqLocInRead_end);
			//cout << "validateSplciePair: " << endl << this->returnValidatePairedSJvecStr() << endl << endl;
			this->get_candi_complicate_SJ_vec(toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				readSeq_inProcess, indexInfo, chrNameInt, max_allowed_mismatchNum);
			//cout << "complicatedSJvec: " << endl << this->returnCandidateComplicatedSJvecStr() << endl << endl;
			bool selectBest_bool = this->selectBestComplicateSJ_lessPenalty(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end, max_allowed_mismatchNum);
			
			if(selectBest_bool)
			{
				this->generateBestComplicateSJMismatchVec(toFixSeqLocInRead_start, toFixSeqLocInRead_end);
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
			//return selectBest_bool;
		}
		
		bool compared2SpliceInfo_new(bool spliceInfoFixed_bool, FixDoubleAnchor_Splice_Info* fixSpliceInfo)
		{
			/*bool betterThanFixedSJ_bool;
			if(!spliceInfoFixed_bool)
			{	
				betterThanFixedSJ_bool = true;
			}
			else
			{
				int currentComplicatedSJ_penalty_best = penalty_best;
				int spliceInfo_penalty = fixSpliceInfo->returnComparedPenalty(mismatch_penalty_weight,
					deletion_penalty_weight*2, insertion_penalty_weight*2); // mismatchPenalty, semiSJ_penalty, non-can_penalty
				if(currentComplicatedSJ_penalty_best < spliceInfo_penalty)
					betterThanFixedSJ_bool = true;
				else 
					betterThanFixedSJ_bool = false;
			}
			return betterThanFixedSJ_bool;*/
			if(!spliceInfoFixed_bool)
			{
				return true;
			}
			if(this->return_mismatch_bestComplicatedSplice() < fixSpliceInfo->returnBestSplice_mismatchNum())
				return true;
			else 
				return false;
		}

		bool compared2SpliceInfo(bool spliceInfoFixed_bool, FixDoubleAnchor_Splice_Info* fixSpliceInfo)
		{
			bool betterThanFixedSJ_bool;
			if(!spliceInfoFixed_bool)
			{	
				betterThanFixedSJ_bool = true;
			}
			else
			{
				int currentComplicatedSJ_penalty_best = penalty_best;
				int spliceInfo_penalty = fixSpliceInfo->returnComparedPenalty(mismatch_penalty_weight,
					deletion_penalty_weight, insertion_penalty_weight); // mismatchPenalty, semiSJ_penalty, non-can_penalty
				if(currentComplicatedSJ_penalty_best < spliceInfo_penalty)
					betterThanFixedSJ_bool = true;
				else 
					betterThanFixedSJ_bool = false;
			}
			return betterThanFixedSJ_bool;
		}



		void scanGenomeAndReadSeq(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt 
			)
		// generate donerMapCumulativeMismatchNumVec & acceptorMapCumulativeMismatchNumVec
		{
			//cout << "start to scan scanGenomeAndReadSeq " << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;

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
				donerMapCumulativeMismatchNumVec.push_back(tmpDonerMapCumulativeMismatchNum);
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
				acceptorMapCumulativeMismatchNumVec.push_back(tmpAcceptorMapCumulativeMismatchNum);
			}
			return;
		}

		void searchSpliceSignal(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			//const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt 
			) // generate vector of GT_splice_doner_vec ......
		{
			int readSeq_length = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			string doner_chrom_str = indexInfo->returnChromStrSubstr(chrNameInt, toFixSeqMapPos_doner, readSeq_length);
			string acceptor_chrom_str = indexInfo->returnChromStrSubstr(chrNameInt,  toFixSeqMapPos_acceptor - readSeq_length + 1, readSeq_length); 
			
			int tmp_search_start_loc = 0;
			while(1)
			{
				int tmp_loc = doner_chrom_str.find("GT", tmp_search_start_loc);
				if(tmp_loc == string::npos)
					break;
				GT_splice_doner_vec.push_back(tmp_loc);
				tmp_search_start_loc = tmp_loc + 1;
			}
			tmp_search_start_loc = 0;
			while(1)
			{
				int tmp_loc = doner_chrom_str.find("CT", tmp_search_start_loc);
				if(tmp_loc == string::npos)
					break;
				CT_splice_doner_vec.push_back(tmp_loc);
				tmp_search_start_loc = tmp_loc + 1;
			}
			tmp_search_start_loc = 0;
			while(1)
			{
				int tmp_loc = acceptor_chrom_str.find("AG", tmp_search_start_loc);
				if(tmp_loc == string::npos)
					break;
				AG_splice_acceptor_vec.push_back(readSeq_length - tmp_loc - 2);
				tmp_search_start_loc = tmp_loc + 1;
			}
			tmp_search_start_loc = 0;
			while(1)
			{
				int tmp_loc = acceptor_chrom_str.find("AC", tmp_search_start_loc);
				if(tmp_loc == string::npos)
					break;
				AC_splice_acceptor_vec.push_back(readSeq_length - tmp_loc - 2);
				tmp_search_start_loc = tmp_loc + 1;
			}			
		}

};

#endif