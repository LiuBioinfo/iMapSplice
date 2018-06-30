// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#ifndef SPLICEJUNCTIONHASH_INFO_H
#define SPLICEJUNCTIONHASH_INFO_H

#include <string>
#include <string.h>
#include <map>
#include <set>

using namespace std;

typedef map<string, set<int> > SpliceEndStrHash; // SJ hash (in 2nd hash, key = anchor string)
typedef map<int, SpliceEndStrHash> SplicePosHash;
typedef SplicePosHash::iterator SplicePosHashIter; 
typedef SpliceEndStrHash::iterator SpliceEndStrHashIter;

typedef map<int, set<int> > SJintHash;  // SJ hash (only 1 hash, value is all possible splice junction positions in the other side)
typedef SJintHash::iterator SJintHashIter;

typedef map<int, set<int> > SJareaHash; //( areaNO = pos/100 ) intermediate hash to directly get all possible splice junctions near a certain position
typedef SJareaHash::iterator SJareaHashIter;


class SJhash_Info
{
public: 
	vector<SplicePosHash> spliceJunctionNormal; // size = chromNum in index_info file
	vector<SplicePosHash> spliceJunctionReverse; // size = chromNum in index_info file

	int anchorStringLength;

	//vector<SJintHash> SJintHashNormal;
	//vector<SJintHash> SJintHashReverse;

	int areaSize;

	vector<SJareaHash> SJstartPosAreaHash;
	vector<SJareaHash> SJendPosAreaHash;

	SJhash_Info()
	{
		areaSize = 1000;
		anchorStringLength = 3;
	}

	void searchForCandiSJwithinAreas(int toFoundSJchrNameInt, 
		int toFoundSJdonerSite_leftMost, int toFoundSJdonerSite_rightMost,
		int toFoundSJacceptorSite_leftMost, int toFoundSJacceptorSite_rightMost,
		vector< pair<int,int> >& candiSJposPairVec)
	{
		// search for acceptor positions
		int areaNOmin_acceptorSite = (int)(toFoundSJacceptorSite_leftMost/areaSize);
		int areaNOmax_acceptorSite = (int)(toFoundSJacceptorSite_rightMost/areaSize);
		vector<int> tmpSJacceptorSiteVec;
		for(int tmpArea = areaNOmin_acceptorSite; tmpArea <= areaNOmax_acceptorSite; tmpArea ++)
		{
			SJareaHashIter tmpSJareaHashIter 
				= ((SJendPosAreaHash)[toFoundSJchrNameInt]).find(tmpArea);
			if(tmpSJareaHashIter != ((SJendPosAreaHash)[toFoundSJchrNameInt]).end())
			{
				for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
					intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
				{
					int tmpSJacceptorPos = (*intSetIter);
					if( (tmpSJacceptorPos >= toFoundSJacceptorSite_leftMost) 
						&& (tmpSJacceptorPos <= toFoundSJacceptorSite_rightMost) )
						tmpSJacceptorSiteVec.push_back(tmpSJacceptorPos);
					else
					{}
				}
			}
			else
			{}			
		}		

		for(int tmp = 0; tmp < tmpSJacceptorSiteVec.size(); tmp++)
		{
			int tmpSJacceptorSite = tmpSJacceptorSiteVec[tmp];
			vector<int> tmpSJdonerSiteVec;
			this->returnAlterDonerSpliceSiteVecWithAcceptorStartPos(
				toFoundSJchrNameInt, tmpSJacceptorSite, tmpSJdonerSiteVec);
			for(int tmp2 = 0; tmp2 < tmpSJdonerSiteVec.size(); tmp2++)
			{
				int tmpSJdonerSite = tmpSJdonerSiteVec[tmp2];
				if((tmpSJdonerSite > tmpSJacceptorSite)
					&&(tmpSJdonerSite >= toFoundSJdonerSite_leftMost)
					&&(tmpSJdonerSite <= toFoundSJacceptorSite_rightMost))
				candiSJposPairVec.push_back(pair<int,int>(tmpSJdonerSite, tmpSJacceptorSite));
			}
		}
	}

	bool searchSJdonerEndPos(int toFoundSJchrNameInt,
		int toFoundSJdonerEndPos, vector<int>& foundSJcorrespondingAcceptorStartPosVec)
	{
		SplicePosHashIter tmpPosHashIter
			= (spliceJunctionNormal)[toFoundSJchrNameInt].find(toFoundSJdonerEndPos);
		if(tmpPosHashIter == (spliceJunctionNormal)[toFoundSJchrNameInt].end())
			return false;
		//cout << "donerEndPos found !" << endl;
		SpliceEndStrHashIter tmpEndStrHashIter = (tmpPosHashIter->second).find("*");
		if(tmpEndStrHashIter != (tmpPosHashIter->second).end())
		{
			for(set<int>::iterator tmpIntSetIter = (tmpEndStrHashIter->second).begin();
				tmpIntSetIter != (tmpEndStrHashIter->second).end(); tmpIntSetIter++)
			{
				foundSJcorrespondingAcceptorStartPosVec.push_back(*tmpIntSetIter);
			}
		}
		else
		{
			cout << "error! * should be found in endStrHash" << endl;
			exit(1);
		}
		return true;
	}

	// bool searchSJacceptorStartPos(int toFoundSJacceptorStartPos, int& foundSJcorrespondingDonerEndPos)
	// {

	// }

	int returnLeftMostSJdonnerEndPos(int tmpChrNameInt, int tmpStartChrPos, int tmpEndChrPos)
	{
		int areaNOmin = (int)(tmpStartChrPos/areaSize);
		int areaNOmax = (int)(tmpEndChrPos/areaSize);

		bool someSJdonnerEndPosFoundAlready = false;
		vector<int> candiSJdonnerEndPosVec;

		for(int tmpArea = areaNOmin; tmpArea <= areaNOmax; tmpArea ++)
		{
			if(!someSJdonnerEndPosFoundAlready)
			{
				SJareaHashIter tmpSJareaHashIter
					= (SJstartPosAreaHash[tmpChrNameInt]).find(tmpArea);
				if(tmpSJareaHashIter != (SJstartPosAreaHash[tmpChrNameInt]).end())
				{
					for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
						intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
					{
						int tmpSJdonerPos = (*intSetIter);
						if( (tmpSJdonerPos >= tmpStartChrPos) && (tmpSJdonerPos <= tmpEndChrPos) )
						{
							candiSJdonnerEndPosVec.push_back(tmpSJdonerPos);
							someSJdonnerEndPosFoundAlready = true;
						}
					}
				}
				else
				{}
			}
		}
		if(someSJdonnerEndPosFoundAlready)
		{
			int tmpLeftMostSJdonnerEndPos = candiSJdonnerEndPosVec[0];
			if(candiSJdonnerEndPosVec.size() > 1)
			{	
				for(int tmp = 1; tmp < candiSJdonnerEndPosVec.size(); tmp++)
				{
					int tmpCandiSJdonnerEndPos = candiSJdonnerEndPosVec[tmp];
					if(tmpCandiSJdonnerEndPos < tmpLeftMostSJdonnerEndPos)
					{
						tmpLeftMostSJdonnerEndPos = tmpCandiSJdonnerEndPos;
					}
				}	
			}
			return tmpLeftMostSJdonnerEndPos;
		}
		else
			return -1;
	}

	int returnRightMostSJacceptorStartPos(int tmpChrNameInt, int tmpStartChrPos, int tmpEndChrPos)
	{
		int areaNOmin = (int)(tmpStartChrPos/areaSize);
		int areaNOmax = (int)(tmpEndChrPos/areaSize);

		bool someSJacceptorStartPosFoundAlready = false;
		vector<int> candiSJacceptorStartPosVec;

		for(int tmpArea = areaNOmax; tmpArea <= areaNOmin; tmpArea --)
		{
			if(!someSJacceptorStartPosFoundAlready)
			{
				SJareaHashIter tmpSJareaHashIter
					= (SJendPosAreaHash[tmpChrNameInt]).find(tmpArea);
				if(tmpSJareaHashIter != (SJendPosAreaHash[tmpChrNameInt]).end())
				{
					for(set<int>::iterator intSetIter = (tmpSJareaHashIter->second).begin();
						intSetIter != (tmpSJareaHashIter->second).end(); intSetIter ++)
					{
						int tmpSJacceptorStartPos = (*intSetIter);
						if( (tmpSJacceptorStartPos >= tmpStartChrPos) && (tmpSJacceptorStartPos <= tmpEndChrPos) )
						{
							candiSJacceptorStartPosVec.push_back(tmpSJacceptorStartPos);
							someSJacceptorStartPosFoundAlready = true;
						}
					}
				}
				else
				{}
			}
		}
		if(someSJacceptorStartPosFoundAlready)
		{
			int tmpRightMostSJacceptorStartPos = candiSJacceptorStartPosVec[0];
			if(candiSJacceptorStartPosVec.size() > 1)
			{	
				for(int tmp = 1; tmp < candiSJacceptorStartPosVec.size(); tmp++)
				{
					int tmpCandiSJacceptorStartPos = candiSJacceptorStartPosVec[tmp];
					if(tmpCandiSJacceptorStartPos > tmpRightMostSJacceptorStartPos)
					{
						tmpRightMostSJacceptorStartPos = tmpCandiSJacceptorStartPos;
					}
				}	
			}
			return tmpRightMostSJacceptorStartPos;
		}
		else
			return -1;
	}


	//////////////////////////////////////////////////////////////
	/////////////////////  SJintHash   ///////////////////////////
	//////////////////////////////////////////////////////////////
	/*
	void insert2SJintHash(
		int chrInt, int spliceStartPos, int spliceEndPos)
	{

		SJintHashIter foundIntHashIter;
		// insert to SJintHashNormal
		foundIntHashIter = SJintHashNormal[chrInt].find(spliceStartPos);
		if(foundIntHashIter == SJintHashNormal[chrInt].end())
		{
			set<int> newPosSet; 
			newPosSet.insert(spliceEndPos);
			SJintHashNormal[chrInt].insert(pair<int, set<int> > (spliceStartPos, newPosSet));
		}
		else
		{
			if((foundIntHashIter->second).find(spliceEndPos) == (foundIntHashIter->second).end())
			{
				(foundIntHashIter->second).insert(spliceEndPos);
			}
			else
			{}
		}

		//insert to SJintHashReverse
		foundIntHashIter = SJintHashReverse[chrInt].find(spliceEndPos);
		if(foundIntHashIter == SJintHashReverse[chrInt].end())
		{
			set<int> newPosSet; 
			newPosSet.insert(spliceStartPos);
			SJintHashReverse[chrInt].insert(pair<int, set<int> > (spliceEndPos, newPosSet));
		}
		else
		{
			if((foundIntHashIter->second).find(spliceStartPos) == (foundIntHashIter->second).end())
			{
				(foundIntHashIter->second).insert(spliceStartPos);
			}
			else
			{}
		}		
	}
	
	void initiateSJintHash(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SJintHash newSJintHash;
			SJintHashNormal.push_back(newSJintHash);	
		}
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SJintHash newSJintHash;
			SJintHashReverse.push_back(newSJintHash);	
		}		
	}
	*/
	//////////////////////////////////////////////////////////////////////////
	//  1. AreadHash: area - SJendSite hash    ///////////////////////////////
	//  2. StringHash: SJendSite - SJotherEndString - SJotherEndSite hash   //
	//////////////////////////////////////////////////////////////////////////

	void returnAlterDonerSpliceSiteVecWithAcceptorStartPos(
		int chrNameInt, int tmpAcceptorSpliceSite, 
		vector<int>& tmpAlterDonerSpliceSitePosVec)
	{
		set<int> tmpAlterDonerSpliceSitePosSet;

		SplicePosHash::iterator tmpSplicePosHashIter = spliceJunctionReverse[chrNameInt].find(tmpAcceptorSpliceSite);
		if(tmpSplicePosHashIter != spliceJunctionReverse[chrNameInt].end())
		{
			SpliceEndStrHash::iterator tmpSpliceEndStrHashIter;
			for(tmpSpliceEndStrHashIter = (tmpSplicePosHashIter->second).begin(); 
				tmpSpliceEndStrHashIter != (tmpSplicePosHashIter->second).end(); tmpSpliceEndStrHashIter++)
			{
				for(set<int>::iterator tmpIntSetIter = (tmpSpliceEndStrHashIter->second).begin();
					tmpIntSetIter != (tmpSpliceEndStrHashIter->second).end(); tmpIntSetIter++)
				{
					int tmpIntSetValue = (*tmpIntSetIter);
					tmpAlterDonerSpliceSitePosSet.insert(tmpIntSetValue);
				}
			}
		}

		for(set<int>::iterator tmpIntSetIter = tmpAlterDonerSpliceSitePosSet.begin();
			tmpIntSetIter != tmpAlterDonerSpliceSitePosSet.end(); tmpIntSetIter++)
		{
			tmpAlterDonerSpliceSitePosVec.push_back(*tmpIntSetIter);
		}		
	}

	void returnAlterAcceptorSpliceSiteVecWithDonerEndPos(
		int chrNameInt, int tmpDonerSpliceSite, 
		vector<int>& tmpAlterAcceptorSpliceSitePosVec)
	{
		set<int> tmpAlterAcceptorSpliceSitePosSet;

		SplicePosHash::iterator tmpSplicePosHashIter = spliceJunctionNormal[chrNameInt].find(tmpDonerSpliceSite);
		if(tmpSplicePosHashIter != spliceJunctionNormal[chrNameInt].end())
		{
			SpliceEndStrHash::iterator tmpSpliceEndStrHashIter;
			for(tmpSpliceEndStrHashIter = (tmpSplicePosHashIter->second).begin(); 
				tmpSpliceEndStrHashIter != (tmpSplicePosHashIter->second).end(); tmpSpliceEndStrHashIter++)
			{
				for(set<int>::iterator tmpIntSetIter = (tmpSpliceEndStrHashIter->second).begin();
					tmpIntSetIter != (tmpSpliceEndStrHashIter->second).end(); tmpIntSetIter++)
				{
					int tmpIntSetValue = (*tmpIntSetIter);
					tmpAlterAcceptorSpliceSitePosSet.insert(tmpIntSetValue);
				}
			}
		}

		for(set<int>::iterator tmpIntSetIter = tmpAlterAcceptorSpliceSitePosSet.begin();
			tmpIntSetIter != tmpAlterAcceptorSpliceSitePosSet.end(); tmpIntSetIter++)
		{
			tmpAlterAcceptorSpliceSitePosVec.push_back(*tmpIntSetIter);
		}
	}

	void initiateSpliceJunctionAreaHash(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			//SplicePosHash newSplicePosHash;
			//spliceJunctionNormal.push_back(newSplicePosHash);
			SJareaHash newSJareaHash;
			SJstartPosAreaHash.push_back(newSJareaHash);

		}
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			//SplicePosHash newSplicePosHash;
			//spliceJunctionReverse.push_back(newSplicePosHash);
			SJareaHash newSJareaHash;
			SJendPosAreaHash.push_back(newSJareaHash);			
		}
	}

	void initiateSpliceJunctionStringHash(int chromNum)
	{
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SplicePosHash newSplicePosHash;
			spliceJunctionNormal.push_back(newSplicePosHash);
		}
		for(int tmp = 0; tmp < chromNum; tmp++)
		{
			SplicePosHash newSplicePosHash;
			spliceJunctionReverse.push_back(newSplicePosHash);
		}
	}

	void initiateAreaAndStringHash(int chromNum)
	{
		this->initiateSpliceJunctionStringHash(chromNum);
		this->initiateSpliceJunctionAreaHash(chromNum);
	}

	void insert2AreaHash(int chrInt,
		int spliceStartPos, int spliceEndPos)
	{
		int spliceStartPosAreaNO = (int)(spliceStartPos/areaSize);
		int spliceEndPosAreaNO = (int)(spliceEndPos/areaSize);

		SJareaHashIter foundAreaHashIter;
		
		// insert to SJstartPosAreaHash
		foundAreaHashIter = SJstartPosAreaHash[chrInt].find(spliceStartPosAreaNO);
		if(foundAreaHashIter == SJstartPosAreaHash[chrInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(spliceStartPos);
			SJstartPosAreaHash[chrInt].insert(pair<int, set<int> > (spliceStartPosAreaNO, newPosSet));
		}
		else
		{
			if( (foundAreaHashIter->second).find(spliceStartPos) == (foundAreaHashIter->second).end() )
			{
				(foundAreaHashIter->second).insert(spliceStartPos);
			}
			else
			{}
		}

		// insert to SJendPosAreaHash
		foundAreaHashIter = SJendPosAreaHash[chrInt].find(spliceEndPosAreaNO);
		if(foundAreaHashIter == SJendPosAreaHash[chrInt].end())
		{
			set<int> newPosSet;
			newPosSet.insert(spliceEndPos);
			SJendPosAreaHash[chrInt].insert(pair<int, set<int> > (spliceEndPosAreaNO, newPosSet));
		}
		else
		{
			if( (foundAreaHashIter->second).find(spliceEndPos)  == (foundAreaHashIter->second).end() )
			{
				(foundAreaHashIter->second).insert(spliceEndPos);
			}
			else
			{}
		}
	}

	void insert2StringHash(int chrInt,
		int spliceStartPos, int spliceEndPos, Index_Info* indexInfo)
	{
		//string anchorString_doner = (indexInfo->chromStr)[chrInt].substr(spliceStartPos - anchorStringLength, anchorStringLength);
		string anchorString_doner = indexInfo->returnChromStrSubstr(chrInt, spliceStartPos - anchorStringLength + 1, anchorStringLength);
		//string anchorString_acceptor = (indexInfo->chromStr)[chrInt].substr(spliceEndPos - 1, anchorStringLength);
		string anchorString_acceptor = indexInfo->returnChromStrSubstr(chrInt, spliceEndPos, anchorStringLength);

		//insert to spliceJunctionNormal
		SplicePosHashIter foundPosHashIter;

		foundPosHashIter = spliceJunctionNormal[chrInt].find(spliceStartPos);
		if(foundPosHashIter != spliceJunctionNormal[chrInt].end())
		{
			SpliceEndStrHashIter foundEndStrHashIter 
				= (foundPosHashIter->second).find(anchorString_acceptor);
			if(foundEndStrHashIter != (foundPosHashIter->second).end())
			{
				set<int>::iterator intSetIter = (foundEndStrHashIter->second).find(spliceEndPos);
				if(intSetIter != (foundEndStrHashIter->second).end())
				{}  
				else
				{
					(foundEndStrHashIter->second).insert(spliceEndPos);
				}
			}
			else
			{
				set<int> newAcceptorPosSet;
				newAcceptorPosSet.insert(spliceEndPos);
				(foundPosHashIter->second).insert( 
					pair<string, set<int> > (anchorString_acceptor, newAcceptorPosSet) );
				//(foundPosHashIter->second).insert( pairn<string, set<int> > ("*", newDonerPosSet) );
			}

			SpliceEndStrHashIter foundEndStrHashIter_star = (foundPosHashIter->second).find("*");
			if(foundEndStrHashIter_star != (foundPosHashIter->second).end() )
			{
				set<int>::iterator intSetIter 
					= (foundEndStrHashIter_star->second).find(spliceEndPos);
				if(intSetIter != (foundEndStrHashIter_star->second).end())
				{}
				else
				{
					(foundEndStrHashIter_star->second).insert(spliceEndPos);
				}
			}
			else
			{
				cout << "key value = * cannot be found !!! in spliceJunction_info.h" << endl;
			}

		}
		else
		{
			set<int> newSpliceEndPosSet;
			newSpliceEndPosSet.insert(spliceEndPos);

			SpliceEndStrHash newSpliceEndStrHash;
			newSpliceEndStrHash.insert( 
				pair<string, set<int> > (anchorString_acceptor, newSpliceEndPosSet) );
			newSpliceEndStrHash.insert( 
				pair<string, set<int> > ("*", newSpliceEndPosSet) );

			spliceJunctionNormal[chrInt].insert(
				pair<int, SpliceEndStrHash> (spliceStartPos, newSpliceEndStrHash) );
		}

		foundPosHashIter = spliceJunctionReverse[chrInt].find(spliceEndPos);
		if(foundPosHashIter != spliceJunctionReverse[chrInt].end())
		{
			SpliceEndStrHashIter foundEndStrHashIter 
				= (foundPosHashIter->second).find(anchorString_doner);
			if(foundEndStrHashIter != (foundPosHashIter->second).end())
			{
				set<int>::iterator intSetIter = (foundEndStrHashIter->second).find(spliceStartPos);
				if(intSetIter != (foundEndStrHashIter->second).end())
				{}
				else
				{
					(foundEndStrHashIter->second).insert(spliceStartPos);
				}
			}
			else
			{
				set<int> newDonerPosSet;
				newDonerPosSet.insert(spliceStartPos);
				(foundPosHashIter->second).insert( 
					pair<string, set<int> > (anchorString_doner, newDonerPosSet) );
			}

			SpliceEndStrHashIter foundEndStrHashIter_star = (foundPosHashIter->second).find("*");
			if(foundEndStrHashIter_star != (foundPosHashIter->second).end() )
			{
				set<int>::iterator intSetIter 
					= (foundEndStrHashIter_star->second).find(spliceStartPos);
				if(intSetIter != (foundEndStrHashIter_star->second).end() )
				{}
				else
				{
					(foundEndStrHashIter_star->second).insert(spliceStartPos);
				} 
			}
			else
			{
				cout << "key value = * cannot be found !!! in spliceJunction_info.h" << endl;
			}

		}
		else
		{
			set<int> newSpliceStartPosSet;
			newSpliceStartPosSet.insert(spliceStartPos);

			SpliceEndStrHash newSpliceEndStrHash;
			newSpliceEndStrHash.insert( 
				pair<string, set<int> > (anchorString_doner, newSpliceStartPosSet) );
			newSpliceEndStrHash.insert(
				pair<string, set<int> > ("*", newSpliceStartPosSet) );
		
			spliceJunctionReverse[chrInt].insert(
				pair<int, SpliceEndStrHash> (spliceEndPos, newSpliceEndStrHash) );
		}
	}

	void insert2AreaAndStringHash(int chrInt,
		int spliceStartPos, int spliceEndPos, Index_Info* indexInfo)
	{
		this->insert2AreaHash(chrInt, spliceStartPos, spliceEndPos);
		this->insert2StringHash(chrInt, spliceStartPos, spliceEndPos, indexInfo);
	}
	
};

#endif