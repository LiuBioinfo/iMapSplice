// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
//#include <hash_map>
#include <map>
#include <set>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
//#include "csapp.h"
#include <sys/stat.h>
#include <errno.h>

using namespace std;

class Jump_Code
{
public:
    int len; // length fo current feature
    string type; //feature type, M/N/I/D/S

    Jump_Code(int _len, string _type)
    {
        len=_len;
        type=_type;
    }
};

void  dpm_init( int * F, char * traceback, int L1, int L2, int d )
{
        F[0] =  0 ;
        traceback[0] = 'n' ;
 
        int i=0, j=0;
 
        for( j = 1; j <= L1; j++ )
        {
                F[j] =  -j * d ;
                traceback[j] =  '-' ;
        }
        for( i = 1; i <= L2; i++ )
        {
                F[i*(L1+1)] =  -i * d ;
                traceback[i*(L1+1)] =  '|' ;
        }
}

int  max_val( int f1, int f2, int f3, char& ptr )
{
        int  max = 0 ;
 
        if( f1 >= f2 && f1 >= f3 )  
        {
                max = f1 ;
                ptr = '|' ;
        }
        else if( f2 > f3 )              
        {
                max = f2 ;
                ptr = '\\' ;
        }
        else
        {
                max = f3 ;
                ptr = '-' ;
        }
        return  max ;   
}

int max_int(int f1, int f2, int f3)
{
        int  max = 0 ;
 
        if( f1 >= f2 && f1 >= f3 )  
        {
                max = f1 ;
               // ptr = '|' ;
        }
        else if( f2 > f3 )              
        {
                max = f2 ;
                //ptr = '\\' ;
        }
        else
        {
                max = f3 ;
                //ptr = '-' ;
        }
        return  max ;     
}

char max_char(int f1, int f2, int f3)
{
        //int  max = 0 ;
        char ptr;
        if( f1 >= f2 && f1 >= f3 )  
        {
                //max = f1 ;
                ptr = '|' ;
        }
        else if( f2 > f3 )              
        {
                //max = f2 ;
                ptr = '\\' ;
        }
        else
        {
                //max = f3 ;
                ptr = '-' ;
        }
        return  ptr ;     
}

void  print_matrix( int * F, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();
 
        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << "   ";
        }
        cout << "\n  ";
 
        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << F[ (i)*(L1+1) + j ] << " ";
                }
                cout << endl;
        }
}
 
 
void  print_traceback( char * traceback, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();
 
        cout << "    ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << " ";
        }
        cout << "\n  ";
 
        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout << traceback[ (i)*(L1+1) + j ] << " ";
                }
                cout << endl;
        }
}

int nw_align(                  // Needleman-Wunsch algorithm
              int*     F,
              char*    traceback,
              string     seq_1,
              string     seq_2,
              //string&    seq_1_al,
              //string&    seq_2_al,
              int d_ins,         // Gap penalty
              int d_del,
              int d_match,
              int d_misma,    
              vector<string>& jumpCodeVec
            )
{
        //cout << "start to do nw_align ...." << endl;
        int k = 0, x = 0, y = 0;
        int fU, fD, fL ;
        //char ptr; 
        char nuc ;
        int i = 0, j = 0;
 
        const int  a = d_match;   //2 Match
        const int  b = d_misma;   //-1 Mismatch
 
        const int  s[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
                                   { b, a, b, b },
                                   { b, b, a, b },
                                   { b, b, b, a } } ;
 
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();
        //print_matrix( F, seq_1, seq_2 );
        //print_traceback(traceback, seq_1, seq_2 );  

        //cout << endl << "start to generate F and traceback matrix ...." << endl;
      
        for( i = 1; i <= L2; i++ )
        {
                for( j = 1; j <= L1; j++ )
                {
                        nuc = seq_1[ j-1 ] ;
 
                        switch( nuc )
                        {
                                case 'A':  x = 0 ;  break ;
                                case 'C':  x = 1 ;  break ;
                                case 'G':  x = 2 ;  break ;
                                case 'T':  x = 3 ;
                        }
 
                        nuc = seq_2[ i-1 ] ;
 
                        switch( nuc )
                        {
                                case 'A':  y = 0 ;  break ;
                                case 'C':  y = 1 ;  break ;
                                case 'G':  y = 2 ;  break ;
                                case 'T':  y = 3 ;
                        }
 
                        fU = F[ (i-1)*(L1+1) + j ] - d_del;
                        fD = F[ (i-1)*(L1+1) + j-1 ] + s[ x ][ y ] ;
                        fL = F[ (i)*(L1+1) + j-1 ] - d_ins;
                        //F[ (i)*(L2+1) + j ] = max_val( fU, fD, fL, ptr ) ;

                        F[ (i)*(L1+1) + j ] = max_int(fU, fD, fL);
                        traceback[ (i)*(L1+1) + j ] =  max_char(fU, fD, fL);
                        
                }
        }
        i-- ; j-- ;
        //print_matrix( F, seq_1, seq_2 );
        //print_traceback(traceback, seq_1, seq_2 );        
        //cout << "start to trace back" << endl;
        while( i > 0 || j > 0 )
        {
                //cout << "i,j: " << i << "," << j << endl;
                //cout << "traceback[i][j]: " << traceback[ (i)*(L1+1) + j ] << endl;
                switch( traceback[ (i)*(L1+1) + j ] )
                {
                        case '|' :      //seq_1_al += '-' ; 
                                        //seq_2_al += seq_2[ i-1 ] ;
                                        jumpCodeVec.push_back("D");
                                        i-- ;
                                        break ;
 
                        case '\\':      //seq_1_al += seq_1[ j-1 ] ; 
                                        //seq_2_al += seq_2[ i-1 ] ;
                                        if(seq_1[j-1] == seq_2[i-1])
                                            jumpCodeVec.push_back("M");
                                        else
                                            jumpCodeVec.push_back("X"); 
                                        i-- ;  j-- ;
                                        break ;
 
                        case '-' :      //seq_1_al += seq_1[ j-1 ] ;
                                        jumpCodeVec.push_back("I");
                                        //seq_2_al += '-' ;
                                        j--;
                                        break;
                        default:
                                        cout << "i,j: " << i << "," << j << endl;
                                        break;
                }
                k++ ;
        }

        return 0;
}

void getFinalNWDPresults(vector<string>& iniJumpCodeVec, vector<Jump_Code>& finalJumpCodeVec,
    vector<int>& mismatchPosVec, vector<char>& mismatchCharVec,
    const string& seq_read)
{
    int iniJumpCodeVecSize = iniJumpCodeVec.size();
    int lastSeqPos = 0;
    vector<string> interJumpCodeVec;
    for(int tmp = iniJumpCodeVecSize-1; tmp >= 0; tmp--)
    {
        string tmpInterJumpCodeType = iniJumpCodeVec[tmp];
        if(tmpInterJumpCodeType == "X")
        {
            lastSeqPos ++;
            mismatchPosVec.push_back(lastSeqPos);
            mismatchCharVec.push_back(seq_read.at(lastSeqPos-1));
            interJumpCodeVec.push_back("M");
        }
        else if((tmpInterJumpCodeType == "M")||(tmpInterJumpCodeType == "I"))
        {
            lastSeqPos ++;
            interJumpCodeVec.push_back(tmpInterJumpCodeType);
        }
        else
        {
            interJumpCodeVec.push_back(tmpInterJumpCodeType);
        }
    }
    int interJumpCodeVecSize = interJumpCodeVec.size();
    int lastJumpCodeLen = 1;
    string lastJumpCodeType = interJumpCodeVec[0];
    for(int tmp = 1; tmp < interJumpCodeVecSize; tmp++)
    {
        string tmpJumpCodeType = interJumpCodeVec[tmp];
        if(tmpJumpCodeType == lastJumpCodeType)
        {
            lastJumpCodeLen ++;
        }
        else
        {
            Jump_Code tmpJumpCode(lastJumpCodeLen, lastJumpCodeType);
            finalJumpCodeVec.push_back(tmpJumpCode);
            lastJumpCodeLen = 1;
            lastJumpCodeType = tmpJumpCodeType;
        }
    }
    Jump_Code tmpJumpCode(lastJumpCodeLen, lastJumpCodeType);
    finalJumpCodeVec.push_back(tmpJumpCode);
}



int nw(                                                          
        string       seq_1,          /*  Needleman-Wunsch   */
        string       seq_2          /*  algorithm for      */
        //string&      seq_1_al,       /*  global alignment   */
        //string&      seq_2_al,       /*  of nt sequence.    */
        //bool         prm
      )
{
        int d = 2;
        int d_ins = d;                 /* gap penalty */
        int d_del = d;
        int d_match = 2;
        int d_misma = -1;

        int  L1 = seq_1.length();
        int  L2 = seq_2.length();
        
        int* F = new int [(L1+1)*(L2+1)];
        char* traceback = new char [(L1+1)*(L2+1)];
        dpm_init( F, traceback, L1, L2, d );

        //cout << "start to create alignment" << endl;
        // Create alignment
        vector<string> jumpCodeVec;
        nw_align( F, traceback, seq_1, seq_2, d_ins, d_del, d_match, d_misma, jumpCodeVec);
 
        int jumpCodeVecSize = jumpCodeVec.size();
        cout << "JumpCode: ";// << endl;
        for(int tmp = jumpCodeVecSize-1; tmp >= 0; tmp--)
            cout << jumpCodeVec[tmp];
        cout << endl;

        vector<Jump_Code> finalJumpCodeVec;
        vector<int> finalMismatchPosVec;
        vector<char> finalMismatchCharVec;

        getFinalNWDPresults(jumpCodeVec, finalJumpCodeVec, 
            finalMismatchPosVec, finalMismatchCharVec, seq_1);

        cout << "Final JumpCode: ";
        for(int tmp = 0; tmp < finalJumpCodeVec.size(); tmp++)
            cout << finalJumpCodeVec[tmp].len << finalJumpCodeVec[tmp].type;
        cout << endl;
        cout << "Mismatch Pos: ";
        for(int tmp = 0; tmp < finalMismatchPosVec.size(); tmp++)
            cout << finalMismatchPosVec[tmp] << ",";
        cout << endl;
        cout << "Mismatch Char: ";
        for(int tmp = 0; tmp < finalMismatchCharVec.size(); tmp++)
            cout << finalMismatchCharVec[tmp] << ",";
        cout << endl;
        delete F;
        delete traceback;
        return  0 ;
}


int main(int argc, char** argv)
{
    if(argc != 3)
    {
        cout << "Executable seq_1 seq_2" << endl;
        exit(1);
    }
    string seq_1 = argv[1];
    string seq_2 = argv[2];

    cout << "seq_1: " << seq_1 << endl;
    cout << "seq_2: " << seq_2 << endl;

    string seq_1_align, seq_2_align;

    nw(seq_1, seq_2);
    
    string seq_1_align_final, seq_2_align_final;
    seq_1_align_final.assign(seq_1_align.rbegin(), seq_1_align.rend());
    seq_2_align_final.assign(seq_2_align.rbegin(), seq_2_align.rend());
    //cout << "seq_1_align: " << endl << seq_1_align_final << endl;
    //cout << "seq_2_align: " << endl << seq_2_align_final << endl;
    return 0;
}