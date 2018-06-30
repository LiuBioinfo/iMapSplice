// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef FIXDOUBLEANCHOR_ANNOTATION_INFO_H
#define FIXDOUBLEANCHOR_ANNOTATION_INFO_H

#include "annotation_info.h"


class FixDoubleAnchor_Annotation_Info
{
private:
		vector<int> donerMapCumulativeMismatchNumVec; 
		// donerMapCumulativeMismatchNumVec[tmp] = Y; 
		// Y is mismatch # if readSeq.substr(toFixSeqLocInRead_start-1, tmp+1) 
		// or toFixSeq.substr(0, tmp+1)  mapped to chromSeq.substr(toFixSeqMapPos_doner - 1, tmp+1); 
				
		vector<int> acceptorMapCumulativeMismatchNumVec;
		// acceptorMapCumulativeMismatchNumVec[tmp] = Y;
		// Y is mismatch # if readSeq.substr(toFixSeqLocInRead_end - 1 - tmp, tmp+1) 
		// or toFixSeq.substr(toFixSeqLength - 1- tmp, tmp+1) mapped to chromSeq.substr(toFixSeqMapPos_end - tmp - 1, tmp+1);

		int max_extension_forward_doner;
		int max_extension_backward_acceptor;

		vector< pair< int, pair<int, int> > > spliceSiteVec; // < < splice_site = forwared_doner_# , mismatch # >, flank_string_case >

		int bestSpliceSiteIndex;
public:
		FixDoubleAnchor_Annotation_Info()
		{
			bestSpliceSiteIndex = -1;
		}

		int selectBestSpliceSite()
		{
			int currentBestSJ_mismatchNum = 1000000;
			int currentBestSJ_index = -1;
			for(int tmp = 0; tmp < spliceSiteVec.size(); tmp ++)
			{
				int currentMismatchNum = spliceSiteVec[tmp].second.first;
				if(currentMismatchNum < currentBestSJ_mismatchNum)
				{
					currentBestSJ_index = tmp;
					currentBestSJ_mismatchNum = currentMismatchNum;
				}
			}
			return currentBestSJ_index;
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

		void scanGenomeAndReadSeq(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt, 
			int max_allowed_mismatchNum)	
		{
			//cout << endl << "start to scan scanGenomeAndReadSeq " << endl << endl;
			int toFixSeqLength = toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1;
			max_extension_forward_doner = toFixSeqLength;
			max_extension_backward_acceptor = toFixSeqLength;

			int tmpDonerMapCumulativeMismatchNum = 0;
			donerMapCumulativeMismatchNumVec.push_back(0);
			for(int tmp = 0; tmp < toFixSeqLength; tmp++)
			{
				if(readSeq_inProcess.at(toFixSeqLocInRead_start - 1 + tmp) 
					!= indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_doner + tmp) )
					tmpDonerMapCumulativeMismatchNum ++;
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
				if(readSeq_inProcess.at(toFixSeqLocInRead_end - tmp - 1) 
					!= indexInfo->returnOneBaseCharInGenome(chrNameInt, toFixSeqMapPos_acceptor - tmp) )
					tmpAcceptorMapCumulativeMismatchNum ++;
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

		bool detectSpliceSite_byAnnotation(
			int toFixSeqLocInRead_start, int toFixSeqLocInRead_end,
			int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor, 
			const string& readSeq_inProcess, 
			Index_Info* indexInfo, int chrNameInt, 
			int max_allowed_mismatchNum, Annotation_Info* annotationInfo)
		{
			this->scanGenomeAndReadSeq(
				toFixSeqLocInRead_start, toFixSeqLocInRead_end,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, 
				readSeq_inProcess, indexInfo, chrNameInt, 
				max_allowed_mismatchNum //int* mismatchNum_detected, int* prefixMatch_length
				);	

			vector<int> candi_doner_length_vec;
			annotationInfo->generateSpliceSite_fixDoubleAnchor_SpliceInfo(
				chrNameInt, toFixSeqMapPos_doner, toFixSeqMapPos_acceptor, 
				(toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1),
				indexInfo, candi_doner_length_vec);

			this->getFinalSpliceSite(candi_doner_length_vec, 
				toFixSeqLocInRead_end - toFixSeqLocInRead_start + 1, 
				max_allowed_mismatchNum,
				toFixSeqMapPos_doner, toFixSeqMapPos_acceptor,
				chrNameInt, indexInfo);

			bestSpliceSiteIndex = this->selectBestSpliceSite();
			if(bestSpliceSiteIndex >= 0)
				return true;
			else
				return false;
		}		

		void getFinalSpliceSite(vector<int>& candi_doner_length_vec, int toFixSeqLength,
			int max_allowed_mismatchNum, int toFixSeqMapPos_doner, int toFixSeqMapPos_acceptor,
			int chrNameInt, Index_Info* indexInfo)
		{
			for(int tmp = 0; tmp < candi_doner_length_vec.size(); tmp ++)
			{
				int tmp_doner_length = candi_doner_length_vec[tmp];
				int tmp_acceptor_length = toFixSeqLength - tmp_doner_length;
				int tmp_mismatch_sum = donerMapCumulativeMismatchNumVec[tmp_doner_length]
					+ acceptorMapCumulativeMismatchNumVec[tmp_doner_length];
				if(tmp_mismatch_sum <= max_allowed_mismatchNum)
				{
					int tmpChromPos_doner = toFixSeqMapPos_doner + tmp_doner_length - 1;
					int tmpChromPos_acceptor = toFixSeqMapPos_acceptor - tmp_acceptor_length + 1;
					string tmp_flank_string = indexInfo->returnFlankString(chrNameInt, tmpChromPos_doner, tmpChromPos_acceptor);
					int tmp_flank_string_case = this->returnFlankStringCase(tmp_flank_string); 
					spliceSiteVec.push_back(pair< int, pair<int, int> >(
						tmp_doner_length, pair<int, int> (tmp_mismatch_sum, tmp_flank_string_case)));
				}
			}
			return;
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
};

#endif