// This file is a part of iMapSplice. Please refer to LICENSE.TXT for the LICENSE
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <iostream>
#include <fstream>

using namespace std;

int main (int argc, char *argv[])
{
	int fdin, fdout;
	char *src, *dst;
	struct stat statbuf;

 if (argc != 3)
	cout << "usage: a.out <fromfile> <tofile>" << endl;

 /* open the input file */
 if ((fdin = open (argv[1], O_RDONLY)) < 0)
	cout << "can't open " << argv[1] << "for reading" << endl;

 /* open/create the output file */
 if ((fdout = open (argv[2], O_RDWR | O_CREAT | O_TRUNC, FILE_MODE)) < 0)
	cout << "can't create" << argv[2] << " for writing" << endl;

 /* find size of input file */
 if (fstat (fdin,&statbuf) < 0)
	cout << "fstat error" << endl;

 /* go to the location corresponding to the last byte */
 if (lseek (fdout, statbuf.st_size - 1, SEEK_SET) == -1)
	cout << "lseek error" << endl;
 
 /* write a dummy byte at the last location */
 if (write (fdout, "", 1) != 1)
	cout << "write error" << endl;;

 /* mmap the input file */
 if ((src = mmap (0, statbuf.st_size, PROT_READ, MAP_SHARED, fdin, 0))
   == (caddr_t) -1)
	cout << "mmap error for input" << endl;

 /* mmap the output file */
 if ((dst = mmap (0, statbuf.st_size, PROT_READ | PROT_WRITE,
   MAP_SHARED, fdout, 0)) == (caddr_t) -1)
	cout << "mmap error for output" << endl;

 /* this copies the input file to the output file */
 memcpy (dst, src, statbuf.st_size);

} /* main */

/*
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <iostream>
#include <fstream>
//#include "csapp.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> // mmap() is defined in this header 
#include <fcntl.h>
using namespace std;

typedef unsigned char BYTE;

void mmapcopy(int fd, int size, void* dest)
{
	cout << "start to mmapCopy" << endl;

	char* start_addr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
	
	//cout << errno << endl;
	cout << "finish mmap" << endl;
		
	if (start_addr == (void *)-1)
	{
		cout << "mmap failed" << endl;
		fprintf(stderr, "mmap: %s\n", strerror(errno));
	}
	
	cout << "before memory copy " << endl;
	//memcpy(dest, start_addr, size);
	cout << "finish memory copy " << endl;
	//munmap(start_addr, size);
	return;
}


int main(int argc, char **argv)
{
	
	if (argc != 2)
	{
		printf("error.\n");
		exit(0);
	}

	string indexStr = argv[1];
 	indexStr += "/";

	string SA_file = indexStr + "_SA";
	string lcpCompress_file = indexStr + "_lcpCompress";
	string childTab_file = indexStr + "_childTab";
	string verifyChild_file = indexStr + "_detChild";

	int indexSize = 2654911540;

    unsigned int *sa; sa = (unsigned int*)malloc((indexSize) * sizeof(unsigned int)); 
	unsigned int *childTab; childTab = (unsigned int*)malloc((indexSize) * sizeof(unsigned int)); 
	
	BYTE *lcpCompress; lcpCompress = (BYTE*)malloc((indexSize) * sizeof(BYTE)); 
	BYTE *verifyChild; verifyChild = (BYTE*)malloc((indexSize) * sizeof(BYTE)); 	


	cout << endl << "start to load SA " << endl;
	struct stat statFile_sa;
	int fd_sa = open(SA_file.c_str(), O_RDWR, 0);
	if (fd_sa < 0)
		fprintf(stderr, "open: %s\n", strerror(errno));  
	fstat(fd_sa, &statFile_sa);
	cout << "stat_SA.st_size: " << statFile_sa.st_size << endl;	
	mmapcopy(fd_sa, statFile_sa.st_size, sa);
	//close(fd_sa);

	cout << "SA[0]: " << endl; 
	cout << sa[0] << endl;
	cout << "SA[1]: " << endl;
	cout << sa[1] << endl;
	cout << "SA[100]: " << sa[10] << endl;

	cout << endl << "start to load childTab " << endl;
	struct stat statFile_childTab;
	int fd_childTab = open(childTab_file.c_str(), O_RDWR, 0);
	if (fd_childTab < 0)
		fprintf(stderr, "open: %s\n", strerror(errno));  
	fstat(fd_childTab, &statFile_childTab);
	cout << "stat_childTab.st_size: " << statFile_childTab.st_size << endl;	
	mmapcopy(fd_childTab, statFile_childTab.st_size, childTab);
	close(fd_childTab);

	cout << endl << "start to load lcpCompress " << endl;
	struct stat statFile_lcpCompress;
	int fd_lcpCompress = open(lcpCompress_file.c_str(), O_RDWR, 0);
	if (fd_lcpCompress < 0)
		fprintf(stderr, "open: %s\n", strerror(errno));   
	fstat(fd_lcpCompress, &statFile_lcpCompress);
	cout << "stat_lcpCompress.st_size: " << statFile_lcpCompress.st_size << endl;	
	mmapcopy(fd_lcpCompress, statFile_lcpCompress.st_size, lcpCompress);

	cout << endl << "start to load verifyChild " << endl;
	struct stat statFile_verifyChild;
	int fd_verifyChild = open(verifyChild_file.c_str(), O_RDWR, 0);
	if (fd_verifyChild < 0)
		fprintf(stderr, "open: %s\n", strerror(errno));  
	fstat(fd_verifyChild, &statFile_verifyChild);
	cout << "stat_verifyChild.st_size: " << statFile_verifyChild.st_size << endl;	
	mmapcopy(fd_verifyChild, statFile_verifyChild.st_size, verifyChild);


	return 0;
}*/