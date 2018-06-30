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

using namespace std;   

#define Range 6 
#define NULL_NUM 4294967290 

void build_up_down(unsigned int *lcptab, unsigned int *up, unsigned int *down, unsigned int n)
{
    unsigned int lastIndex = NULL_NUM;	
	stack<unsigned int> up_down;
	up_down.push(0);
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		while (lcptab[i] < lcptab[up_down.top()])
		{	lastIndex = up_down.top();
			up_down.pop();
			if((lcptab[i] <= lcptab[up_down.top()]) && (lcptab[up_down.top()] != lcptab[lastIndex]))
				down[up_down.top()] = lastIndex;
		}
		// now lcptab[i] >= lcptab[up_down.top()] holds
		if(lastIndex != NULL_NUM)
		{
			up[i] = lastIndex;
			lastIndex = NULL_NUM;
		}
		up_down.push(i);	
	}
	return;
}

void build_next(unsigned int *lcptab, unsigned int *next, unsigned int n)
{
	stack<unsigned int> nextIndex;
	unsigned int lastIndex;
	unsigned int j;
	nextIndex.push(0);
	for(j = 0; j < n; j++)
	{
		while(lcptab[j] < lcptab[nextIndex.top()])
			nextIndex.pop();
		if(lcptab[j] == lcptab[nextIndex.top()])
		{
			lastIndex = nextIndex.top();
			nextIndex.pop();
			next[lastIndex] = j;
		}
		nextIndex.push(j);
	}	
	return;
}

void build_lcp(unsigned int *r, unsigned int *sa, unsigned int *lcp, unsigned int *rank, unsigned int n)
{

	unsigned int i, j; 
	int k=0;
	for (i = 0; i < n; i++) rank[sa[i]] = i;
		cout << "stop10_new" << endl;
	for (i = 0; i < n; lcp[rank[i++]] = k) 
	for (k?k--:0, (rank[i] == 0)?(j=0):(j=sa[rank[i]-1]); r[i+k] == r[j+k]; k++);
	lcp[0] = 0;
		cout << "stop11" << endl; 
	return;
}

unsigned int cmp(unsigned int *r,unsigned int a,unsigned int b,unsigned int l)
{return r[a]==r[b]&&r[a+l]==r[b+l];}  

void da(unsigned int *r,unsigned int *sa,unsigned int n,unsigned int m)
{   
	unsigned int maxn = n;

    unsigned int *wa = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *wb = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *wv = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int *ws = (unsigned int*)malloc(maxn * sizeof(unsigned int));
    unsigned int i,j,p,*x=wa,*y=wb,*t;
   
    for(i=0;i<m;i++) ws[i]=0;

    for(i=0;i<n;i++) ws[x[i]=r[i]]++; 

	for(i=1;i<m;i++) ws[i]+=ws[i-1];

    for(i=n-1;i>=0;i--) 
    {
    	sa[--ws[x[i]]]=i; 
		if(i == 0)
			break;
	}

    for(j=1,p=1;p<n;j*=2,m=p)
    {
        for(p=0,i=n-j;i<n;i++) y[p++]=i;  
        for(i=0;i<n;i++) if(sa[i]>=j) y[p++]=sa[i]-j; 
        for(i=0;i<n;i++) wv[i]=x[y[i]];  
        for(i=0;i<m;i++) ws[i]=0;
        for(i=0;i<n;i++) ws[wv[i]]++;
        for(i=1;i<m;i++) ws[i]+=ws[i-1];
        for(i=n-1;i>=0;i--) 
        {
        	sa[--ws[wv[i]]]=y[i];  
			if(i == 0)
				break;
		}
		for(t=x,x=y,y=t,p=1,x[sa[0]]=0,i=1;i<n;i++)
        x[sa[i]]=cmp(y,sa[i-1],sa[i],j)?p-1:p++; 
    }
    return;
}

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		cout << "Executable <InputChrom> <Output_name_of_lcp-up-down-next> <StringLength>" << endl;
		// StringLength = chromEndPosInGenome[last_one] + 1, chromEndPosInGenome can be found in parameter file.
		exit(0);
	}
	unsigned int MAX = atoi(argv[3]) + 1;
	cout << "MAX = " << MAX <<endl;
	char* chrom_file = argv[1];
	char* output_name = argv[2];
	string SA_file = argv[2]; SA_file.append("_SA"); ofstream SA_file_ofs(SA_file.c_str(),ios::binary); 
	string lcp_file = argv[2]; lcp_file.append("_lcp"); ofstream lcp_file_ofs(lcp_file.c_str(),ios::binary);
	string up_file = argv[2]; up_file.append("_up"); ofstream up_file_ofs(up_file.c_str(),ios::binary);
	string down_file = argv[2]; down_file.append("_down"); ofstream down_file_ofs(down_file.c_str(),ios::binary);
	string next_file = argv[2]; next_file.append("_next"); ofstream next_file_ofs(next_file.c_str(),ios::binary);
	string chrom_bit_file = argv[2]; chrom_bit_file.append("_chrom"); ofstream chrom_bit_file_ofs(chrom_bit_file.c_str(),ios::binary);
	
  	FILE *fp = fopen(chrom_file, "r");
    unsigned int *r = (unsigned int*)malloc(MAX * sizeof(unsigned int));
	char ch;
	char head[100];
	char base[Range] = {'X','A','C','G','T','N'};
	char *chrom = (char*)malloc(MAX * sizeof(char));	
	
	//fgets(head,sizeof(head),fp);
	//cout <<"head = " << head << endl;
	unsigned int chrom_base_num = 0;
	while((ch = fgetc(fp)) != EOF)
	{   
		//printf("ch = %c\n",ch); 
		if((ch == 'A')||(ch == 'a')) {chrom[chrom_base_num] = 'A'; r[chrom_base_num] = 1; chrom_base_num++;}
		else if((ch == 'C')||(ch == 'c')) {chrom[chrom_base_num] = 'C'; r[chrom_base_num] = 2; chrom_base_num++;}
		else if((ch == 'G')||(ch == 'g')) {chrom[chrom_base_num] = 'G'; r[chrom_base_num] = 3; chrom_base_num++;}
		else if((ch == 'T')||(ch == 't')) {chrom[chrom_base_num] = 'T'; r[chrom_base_num] = 4; chrom_base_num++;}
		else if((ch == 'N')||(ch == 'n')) {chrom[chrom_base_num] = 'N'; r[chrom_base_num] = 5; chrom_base_num++;}
		else if((ch == 'X')) {chrom[chrom_base_num] = 'X'; r[chrom_base_num] = 6; chrom_base_num++;}
		else if((ch == '\t')||(ch == '\n')) {continue;}
		else {printf("\n illegal input is '%c'",ch); break;}
	}
	cout << "the number of bases in Chromo is "<< chrom_base_num << endl;
	cout << "chrom is ready" << endl;
	r[MAX-1] = Range;
	chrom[MAX-1] = 'X';

	//cout << "MAX: " << MAX << endl << "chrom[MAX-1]: " << chrom[MAX-1] << endl << "chrom[MAX-2]: " 
	//	<< chrom[MAX-2] << endl << "chrom[MAX-3]: " << chrom[MAX-3] << endl
	//	<< "chrom[MAX-4]: " << chrom[MAX-4] << endl;

    unsigned int *sa = (unsigned int*)malloc(MAX * sizeof(unsigned int));
   
	cout << "start to build SA array" << endl;
	da(r,sa,MAX,Range+1);

	cout << "SA is ready" << endl;	
    //////////////////////////////////////////////////
    //unsigned int rank[MAX]={0}, lcp[MAX]={0}, up[MAX] = {0}, down[MAX] = {0}, next[MAX] = {0};  
    unsigned int *rank = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *lcp = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *up = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *down = (unsigned int*)malloc(MAX * sizeof(unsigned int));
    unsigned int *next = (unsigned int*)malloc(MAX * sizeof(unsigned int));

    cout << "start to build LCP array"<< endl;
	build_lcp(r, sa, lcp, rank, MAX); //build lcp array

	cout << "lcp is ready" << endl;	
	cout << "start to build up & down array" << endl;
	build_up_down(lcp, up, down, MAX); // build up and down array
	cout << "up & down are ready" << endl;
	cout << "start to build next array" << endl;		
	build_next(lcp, next, MAX); //build next array
	cout << "next is ready" << endl;	

	//FILE *fp_bit_sa	
	cout << "start to output Index to file "<< endl;
	SA_file_ofs.write((const char*) sa, MAX * sizeof(unsigned int));
	lcp_file_ofs.write((const char*) lcp, MAX * sizeof(unsigned int));
	up_file_ofs.write((const char*) up, MAX * sizeof(unsigned int));
	down_file_ofs.write((const char*) down, MAX * sizeof(unsigned int));
	next_file_ofs.write((const char*) next, MAX * sizeof(unsigned int));
	chrom_bit_file_ofs.write((const char*) chrom, MAX * sizeof(char));
    return 0;
}