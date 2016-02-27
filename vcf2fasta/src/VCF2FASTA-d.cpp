#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <iostream>

#define LINE_WIDTH 70
#define BUFFER_SIZE 512000


#define NO -1
#define SNP 0
#define INS 1
#define DEL 2
#define SV  3
#define OTHER  4


char getNextSymbol(FILE *file) {
	char ch;
	while(true) {
		ch = getc(file);
		if (ch != '\n')
			break;
	}
	return ch;
	
}

void writeChar(char ch, char * buf, int &bufPos, int &lineCnt, char * filename)
{
    FILE *file;
	if(bufPos == BUFFER_SIZE)
	{
        
        file = fopen(filename, "a");
		fwrite(buf, 1, BUFFER_SIZE, file);
		fclose(file);
        bufPos = 0;
	}
	
	buf[bufPos]=ch;
	bufPos++;
	if(bufPos == BUFFER_SIZE)
	{
        file = fopen(filename, "a");
		fwrite(buf, 1, BUFFER_SIZE, file);
        fclose(file);
		bufPos = 0;
	}
	lineCnt++;
	if (lineCnt == LINE_WIDTH) {
		lineCnt=0;
		buf[bufPos] = '\n';
		bufPos++;
	}
}

void writeChar(int ch, char * buf, int &bufPos, int &lineCnt, char * filename)
{
    FILE * file;
	if(bufPos == BUFFER_SIZE)
	{
        file = fopen(filename, "a");
		fwrite(buf, 1, BUFFER_SIZE, file);
        fclose(file);
		bufPos = 0;
	}
	
	buf[bufPos]=(255&ch);
	bufPos++;
	if(bufPos == BUFFER_SIZE)
	{
        file = fopen(filename, "a");
		fwrite(buf, 1, BUFFER_SIZE, file);
        fclose(file);
		bufPos = 0;
	}
	lineCnt++;
	if (lineCnt == LINE_WIDTH) {
		lineCnt=0;
		buf[bufPos] = '\n';
		bufPos++;
	}
}


struct variant
{
    bool isInsertion;
    bool doubleInsertion;
    char * alt;
};


//input arg:  <vcf> <ref> [optional: startPosOfRef]
int main (int argc, char * const argv[]) {
    
    variant * vt;
	
	int startPos = 1;
    if(argc >= 4)
        startPos = atoi(argv[3]);
    
    
    FILE * vcf;	//input file with vcf
	FILE * ref; //input file with reference seuence
	
    char  name[1000], refC[1000000],  altList[1000000], alt[1000000], filter[20], info[1000], format[100], data[200], qual[20], endPos[20];
    size_t tmpPt;
    
    printf("input vcf: %s\n", argv[1]);
	vcf = fopen(argv[1], "r");
	if (vcf == NULL) {
		printf("Cannot open %s\n", argv[1]);
		exit(8);
	}
	
    printf("input reference: %s\n", argv[2]);
	ref = fopen(argv[2], "r");
	if (ref == NULL) {
		printf("Cannot open %s\n", argv[2]);
		exit(8);
	}
    
    do {
		fscanf(vcf, "%s", data);
		
	} while (strncmp(data, "#CHROM", 6) != 0);
	fscanf(vcf, "%s %s %s %s %s %s %s %s", data, data, data, data, data, data, data, data);
    //count no samples
    tmpPt = ftell(vcf);
    int noGen = 0;
    while (getc(vcf) != '\n')
    {
        fscanf(vcf, "%s", data);
        noGen ++;
    }
    fseek( vcf, tmpPt, SEEK_SET);
    
    
	
	int  i, linePos[2][noGen], j =0, samePos=0;
	unsigned long deleteSeq[2][noGen];
    unsigned long insertSeq[2][noGen];
    unsigned long prevIns[2];
	char chrom[15];
	char * buffer;
	buffer = new char[BUFFER_SIZE];
	
	char ** buf[2];
	buf[0] = new char*[noGen];
	buf[1] = new char*[noGen];
	for (i=0; i < noGen; i++)
	{
		buf[0][i] = new char[BUFFER_SIZE];
		buf[1][i] = new char[BUFFER_SIZE];
	}
	
	int bufferPos = 0 , bufPos[2][noGen], gt_pos;
	int pos,  currentPos = 0;
	char  *ptr;
    
	
	bool  doubleInsertion=false;
	
	int delLen, currChar, noAlt, currVar, gt[2], a;
	char* p;
    p = new char[200];
    char* r;
    r = new char[200];
    
	int  chRef;
	unsigned short int prevVariant = NO;
    unsigned long int refSize, logTreshold;
    
    
    
	char * resName[2][noGen];
	//resultFileName = new char[strlen(argv[1])+ 20];
    
	for (i=0; i < noGen; i++)
	{
		
		linePos[0][i]=0;
		linePos[1][i]=0;
		deleteSeq[0][i]=0;
		deleteSeq[1][i]=0;
        insertSeq[0][i]=0;
        insertSeq[1][i]=0;
		bufPos[0][i]=0;
		bufPos[1][i]=0;
		
		
		fscanf(vcf, "%s", data);
		resName[0][i] = (char *) malloc ((strlen(argv[1])+ 20)*sizeof(char));
        resName[1][i] = (char *) malloc ((strlen(argv[1])+ 20)*sizeof(char));
		
		resName[0][i][0]='\0';
		strcat(resName[0][i], argv[1]);
		strcat(resName[0][i], ".");
		strcat(resName[0][i], data);
		strcat(resName[0][i], ".1.fa");
		strcat(resName[0][i], "\0");
				
		
		resName[1][i][0]='\0';
		strcat(resName[1][i], argv[1]);
		strcat(resName[1][i], ".");
		strcat(resName[1][i], data);
		strcat(resName[1][i], ".2.fa");
		strcat(resName[1][i], "\0");
/*		resultSeq[1][i] = fopen(resultFileName, "w");
		if (resultSeq[1][i] == NULL) {
			printf("Cannot open %s\n", resultFileName);
			exit(8);
		}*/
	}
    
    
    char header[1000];
    i=0;
    chRef=getc(ref);
	while (chRef != '\n') {
		header[i++]=chRef;
        chRef=getc(ref);
	}
    header[i++]='\n';
    header[i]='\0';
    
    FILE * file;
    for (i=0; i < noGen; i++)
	{
        file = fopen(resName[0][i], "w");
		fprintf(file, "%s", header);
        fclose(file);
        file = fopen(resName[1][i], "w");
		fprintf(file, "%s", header);
        fclose(file);
    }
    
    
    tmpPt = ftell(ref);
    fseek( ref, 0, SEEK_END);
	refSize = ftell(ref);
    fseek( ref, tmpPt, SEEK_SET);
    logTreshold = startPos+ refSize/5;
    
    
	fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
    p = strtok(format, ":");
    gt_pos = 0;
    while (strcmp(p,"GT") != 0) {
        p = strtok(NULL, ":");
        gt_pos++;
    }
    
	
    printf("Processing...\n");
    
    
	currentPos = startPos;
	
	while (currentPos < pos){
		
		chRef = getNextSymbol(ref);
        
        if(bufferPos == BUFFER_SIZE)
		{
			for (i=0; i < noGen; i++)
			{
                file = fopen(resName[0][i], "a");
				fwrite(buffer, 1, BUFFER_SIZE, file);
                fclose(file);
                file = fopen(resName[1][i], "a");
				fwrite(buffer, 1, BUFFER_SIZE, file);
                fclose(file);
				
			}
			bufferPos = 0;
		}
		
		
        buffer[bufferPos]=(255&chRef);
        bufferPos++;
        if(bufferPos == BUFFER_SIZE)
        {
            for (i=0; i < noGen; i++)
            {
                file = fopen(resName[0][i], "a");
				fwrite(buffer, 1, BUFFER_SIZE, file);
                fclose(file);
                file = fopen(resName[1][i], "a");
				fwrite(buffer, 1, BUFFER_SIZE, file);
                fclose(file);
            }
            bufferPos = 0;
        }
        j++;
        if (j == LINE_WIDTH) {
            j=0;
            buffer[bufferPos]= '\n';
            bufferPos++;
        }
		
		currentPos++;
	}
	for (i=0; i < noGen; i++)
	{
        file = fopen(resName[0][i], "a");
		fwrite(buffer, 1, bufferPos, file);
        fclose(file);
        file = fopen(resName[1][i], "a");
		fwrite(buffer, 1, bufferPos, file);
        fclose(file);
	}
	for (i=0; i < noGen; i++)
	{
		linePos[0][i] = j;
		linePos[1][i] = j;
	}
	
    
    bool variantOK;
    
    while (!feof(vcf))
	{
        variantOK = true;
        if(currentPos > logTreshold)
        {
            printf("Processed %.f%% of reference sequence\n", (float)(logTreshold-startPos)/(float)(refSize)*100);
            logTreshold = logTreshold + refSize/5;
        }
        
        
        currChar = 0;
        noAlt = 1;
        while(altList[currChar] != '\0')
        {
            if(altList[currChar] == ',')
                noAlt++;
            currChar++;
        };
        
        vt = new variant[noAlt];
        
        currChar = 0;
        currVar = 0;
        while (altList[currChar] != '\0')
        {
            a = 0;
            while (altList[currChar] != ',' && altList[currChar] != '\0')
            {
                alt[a] = altList[currChar];
                a++;
                currChar++;
            }
            alt[a] = '\0';
            vt[currVar].alt = new char[strlen(alt)+2];
            memcpy(vt[currVar].alt , alt, strlen(alt));
            vt[currVar].alt[strlen(alt)]='\0';
            
            if (altList[currChar] == ',')
                currChar++;
            currVar++;
        }
        
        
		while (currentPos < pos) {
			
            
			chRef = getNextSymbol(ref);
			
			for (i=0; i < noGen; i++)
			{
				if(deleteSeq[0][i] > 0)
				{
					deleteSeq[0][i]--;
				}
				else
				{
					writeChar(chRef, buf[0][i], bufPos[0][i], linePos[0][i], resName[0][i]);
				}
				
				if(deleteSeq[1][i] > 0)
				{
					deleteSeq[1][i]--;
				}
				else
				{
					writeChar(chRef, buf[1][i], bufPos[1][i], linePos[1][i], resName[1][i]);
				}
			}
			currentPos++;
		}
		
		if(!samePos)
			chRef = getNextSymbol(ref);
		
        if (strcmp("PASS", filter) != 0 && strcmp(".", filter) != 0)
        {
            printf("*** WARNING ***\n");
            printf("\tPosition: %d. Variant's FILTER value is neither 'PASS' nor '.': %s\n", pos, filter);
            printf("\tVariant will be processed anyway. If it should not, it should be removed from the VCF file\n");
            
        }
        if (strstr(info, "VT=SNP") != NULL)   //VT is a SNP
        {
            if(samePos)
            {
                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                printf("\tPrevious variant at this position (%d) was: SNP and current variant at this position is SNP!\n", pos);
                
                while (getc(vcf) != '\n')
                    ;
                fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                p = strtok(format, ":");
                gt_pos = 0;
                while (strcmp(p,"GT") != 0) {
                    p = strtok(NULL, ":");
                    gt_pos++;
                }
                samePos = 0;
                if(pos == currentPos)
                    samePos = 1;
                else
                    currentPos++;
                
                continue;
            }
            prevVariant = SNP;
            
            
            if (strlen(refC) != 1)
            {
                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                printf("\tPosition: %d. It is a SNP and REF length is not 1\n", pos);
                
                while (getc(vcf) != '\n')
                    ;
                fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                p = strtok(format, ":");
                gt_pos = 0;
                while (strcmp(p,"GT") != 0) {
                    p = strtok(NULL, ":");
                    gt_pos++;
                }
                samePos = 0;
                if(pos == currentPos)
                    samePos = 1;
                else
                    currentPos++;
                continue;
                
            }
            
            if((255&chRef) != refC[0])
            {
                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                printf("\tPosition: %d. Char in reference sequence and in vcf file don't match.\n", pos);
                
                while (getc(vcf) != '\n')
                    ;
                fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                p = strtok(format, ":");
                gt_pos = 0;
                while (strcmp(p,"GT") != 0) {
                    p = strtok(NULL, ":");
                    gt_pos++;
                }
                samePos = 0;
                if(pos == currentPos)
                    samePos = 1;
                else
                    currentPos++;
                continue;
                
            }
            
            for (a = 0; a < noAlt; a++)
            {
                if (strlen(vt[a].alt) != 1)
                {
                    printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                    printf("\tPosition: %d. It is a SNP and ALT length is not 1.\n", pos);
                    
                    while (getc(vcf) != '\n')
                        ;
                    fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                    p = strtok(format, ":");
                    gt_pos = 0;
                    while (strcmp(p,"GT") != 0) {
                        p = strtok(NULL, ":");
                        gt_pos++;
                    }
                    samePos = 0;
                    if(pos == currentPos)
                        samePos = 1;
                    else
                        currentPos++;
                    variantOK = false;
                    break;
                }
            }
            if(!variantOK)
                continue;
            
            for (i=0; i < noGen; i++)
            {
                fscanf(vcf, "%s", data);
                
                p = strtok( data, ":" );
                for(int m=0; m < gt_pos; m++)
                    p = strtok( NULL, ":" );
                
                //get first/left chromosome
                r = strtok(p, "|" );
                if(r == NULL)
                    r = p;
                if (strcmp(r,".") == 0)
                    r[0] = '0';
                gt[0] = atoi(r);
                if (gt[0] > noAlt)
                {
                    printf("*** WARNING***\n");
                    printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                    gt[0] = 0;
                }
                //get second/right chromosome
                r = strtok(NULL, "|" );
                if(r == NULL)
                    r = p;
                if (strcmp(r,".") == 0)
                    r[0] = '0';
                gt[1] = atoi(r);
                if (gt[1] > noAlt)
                {
                    printf("*** WARNING***\n");
                    printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                    gt[1] = 0;
                }
                
                
                for(int no = 0; no < 2; no++)
                    if(gt[no] > 0)
                    {
                        if(deleteSeq[no][i] > 0)
                        {
                            deleteSeq[no][i]--;
                        }
                        else
                        {
                            writeChar(vt[gt[no]-1].alt[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                        }
                    }
                    else
                    {
                        if(deleteSeq[no][i] > 0)
                        {
                            deleteSeq[no][i]--;
                        }
                        else
                        {
                            writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                        }
                    }
            }
            
        }
        else if (strstr(info, "VT=INDEL") != NULL)  //VT is an INDEL
        {
            for (a = 0; a < noAlt; a++)
            {
                if ((255&vt[a].alt[0]) != (255&refC[0]))
                {
                    {
                        printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                        printf("\tPosition: %d. It is an INDEL and first characters of ALT and REF do not match.\n", pos);
                        
                        while (getc(vcf) != '\n')
                            ;
                        fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                        p = strtok(format, ":");
                        gt_pos = 0;
                        while (strcmp(p,"GT") != 0) {
                            p = strtok(NULL, ":");
                            gt_pos++;
                        }
                        samePos = 0;
                        if(pos == currentPos)
                            samePos = 1;
                        else
                            currentPos++;
                        variantOK = false;
                        break;
                    }
                }
            }
            if(!variantOK)
                continue;
            
            if(strlen(refC)==1) //insertion
            {
                //                   strlen(refC) < strlen(vt[0].alt) &&
                
                if(samePos && prevVariant == INS)
                {
                    doubleInsertion=true;
                }
                else
                    doubleInsertion=false;
                
                
                prevVariant = INS;
                
                for (i=0; i < noGen; i++)
                {
                    
                    if(doubleInsertion)
                    {
                        prevIns[0] = insertSeq[0][i];
                        prevIns[1] = insertSeq[1][i];
                    }
                    else
                    {
                        prevIns[0] = 0;
                        prevIns[1] = 0;
                    }
                    
                    
                    
                    fscanf(vcf, "%s", data);
                    
                    p = strtok( data, ":" );
                    for(int m=0; m < gt_pos; m++)
                        p = strtok( NULL, ":" );
                    
                    //get first/left chromosome
                    r = strtok(p, "|" );
                    if(r == NULL)
                        r = p;
                    if (strcmp(r,".") == 0)
                        r[0] = '0';
                    gt[0] = atoi(r);
                    if (gt[0] > noAlt)
                    {
                        printf("*** WARNING***\n");
                        printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                        gt[0] = 0;
                        
                    }
                    //get second/right chromosome
                    r = strtok(NULL, "|" );
                    if(r == NULL)
                        r = p;
                    if (strcmp(r,".") == 0)
                        r[0] = '0';
                    gt[1] = atoi(r);
                    if (gt[1] > noAlt)
                    {
                        printf("*** WARNING***\n");
                        printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                        gt[1] = 0;
                    }
                    
                    for(int no = 0; no < 2; no++)
                        if(gt[no] > 0)
                        {
                            if(!samePos)
                            {
                                if(deleteSeq[no][i] > 0)
                                {
                                    deleteSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                            
                            
                            for (j = 1; j < strlen(vt[gt[no]-1].alt); j++)
                            {
                                
                                if(doubleInsertion && insertSeq[i] > 0)
                                {
                                    insertSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(vt[gt[no]-1].alt[j], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                            insertSeq[no][i]=(strlen(vt[gt[no]-1].alt) - 1) > prevIns[no] ? (strlen(vt[gt[no]-1].alt) - 1) : prevIns[no];
                        }
                        else
                        {
                            
                            if(!doubleInsertion)
                            {
                                insertSeq[no][i]=0;
                            }
                            if(!samePos)
                            {
                                if(deleteSeq[no][i] > 0)
                                {
                                    deleteSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                        }
                }
                
            }
            else //if (strlen(refC)>1) //deletion
            {
                //strlen(vt[0].alt) && strlen(vt[0].alt)==
                prevVariant = DEL;
                
                for (i=0; i < noGen; i++)
                {
                    fscanf(vcf, "%s", data);
                    p = strtok( data, ":" );
                    for(int m=0; m < gt_pos; m++)
                        p = strtok( NULL, ":" );
                    
                    //get first/left chromosome
                    r = strtok(p, "|" );
                    if(r == NULL)
                        r = p;
                    if (strcmp(r,".") == 0)
                        r[0] = '0';
                    gt[0] = atoi(r);
                    if (gt[0] > noAlt)
                    {
                        printf("*** WARNING***\n");
                        printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                        gt[0] = 0;
                    }
                    //get second/right chromosome
                    r = strtok(NULL, "|" );
                    if(r == NULL)
                        r = p;
                    if (strcmp(r,".") == 0)
                        r[0] = '0';
                    gt[1] = atoi(r);
                    if (gt[1] > noAlt)
                    {
                        printf("*** WARNING***\n");
                        printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                        gt[1] = 0;
                    }
                    
                    for(int no = 0; no < 2; no++)
                        if (gt[no] == 0)
                        {
                            if(!samePos)
                            {
                                if(deleteSeq[no][i] > 0)
                                {
                                    deleteSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                        }
                        else
                        {
                            if(!samePos)
                            {
                                if(deleteSeq[no][i] > 0)
                                {
                                    deleteSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                            
                            deleteSeq[no][i] = (deleteSeq[no][i]>(strlen(refC)-1)) ? deleteSeq[no][i] :(strlen(refC)-1) ;
                        }
                    
                }
            }
            
        }
        else if (strstr(info, "VT=SV") != NULL)  //VT is a SV
        {
            if (strstr(info, "SVTYPE=DEL") != NULL)
            {
                
                endPos[0] = '\0';
                if (strstr(info, ";END=") != NULL)
                {
                    ptr=strstr(info, ";END=");
                    j=5;
                    while (ptr[j] != ';' &&  ptr[j] != '\0')
                    {
                        endPos[j-5] =  ptr[j];
                        j++;
                    }
                    endPos[j-5]='\0';
                    
                }
                else if (strncmp(info, "END=", 4) == 0)
                {
                    j=4;
                    while (info[j] != ';' &&  info[j] != '\0')
                    {
                        endPos[j-4] =  info[j];
                        j++;
                    }
                    endPos[j-4]='\0';
                }
                else
                {
                    printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                    printf("\tPosition: %d. SVTYPE=DEL and no END in info: %s\n", pos, info);
                    
                    while (getc(vcf) != '\n')
                        ;
                    fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                    p = strtok(format, ":");
                    gt_pos = 0;
                    while (strcmp(p,"GT") != 0) {
                        p = strtok(NULL, ":");
                        gt_pos++;
                    }
                    samePos = 0;
                    if(pos == currentPos)
                        samePos = 1;
                    else
                        currentPos++;
                    continue;
                    
                }
                
                delLen = atoi(endPos)-pos;
                
                
                for (a = 0; a < noAlt; a++)
                {
                    if (strncmp(vt[a].alt, "<DEL>", 5) == 0) {
                        //isInsertion=false;
                        //prevVariant = DEL;
                        vt[a].alt[0] = refC[0];
                        vt[a].alt[1] = '\0';
                        vt[a].isInsertion=false;
                        
                    }
                    else if(strlen(vt[a].alt) == 1 && vt[a].alt[0]==refC[0])
                    {
                        //isInsertion=false;
                        //prevVariant = DEL;
                        vt[a].alt[0] = refC[0];
                        vt[a].alt[1] = '\0';
                        
                        vt[a].isInsertion=false;
                        
                    }
                    else if(strlen(vt[a].alt) > 1 && vt[a].alt[0]==refC[0])
                    {
                        if(samePos && prevVariant == INS)
                        {
                            vt[a].doubleInsertion = true;
                        }
                        else
                        {
                            vt[a].doubleInsertion = false;
                        }
                        vt[a].isInsertion=true;
                        
                        //prevVariant = INS;
                    }
                    else {
                        printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                        printf("\tPosition: %d. SVTYPE=DEL and for given REF, wrong ALT data: %s\n", pos, altList);
                        
                        while (getc(vcf) != '\n')
                            ;
                        fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                        p = strtok(format, ":");
                        gt_pos = 0;
                        while (strcmp(p,"GT") != 0) {
                            p = strtok(NULL, ":");
                            gt_pos++;
                        }
                        samePos = 0;
                        if(pos == currentPos)
                            samePos = 1;
                        else
                            currentPos++;
                        variantOK = false;
                        break;
                    }
                }
                if(!variantOK)
                    continue;
                
                prevVariant = DEL;
                for (a = 0; a < noAlt; a++)
                {
                    if(strlen(vt[a].alt) > 1)
                        prevVariant = INS;
                }
                
                for (i=0; i < noGen; i++)
                {
                    fscanf(vcf, "%s", data);
                    
                    p = strtok( data, ":" );
                    for(int m=0; m < gt_pos; m++)
                        p = strtok( NULL, ":" );
                    
                    //get first/left chromosome
                    r = strtok(p, "|" );
                    if(r == NULL)
                        r = p;
                    if (strcmp(r,".") == 0)
                        r[0] = '0';
                    gt[0] = atoi(r);
                    if (gt[0] > noAlt)
                    {
                        printf("*** WARNING***\n");
                        printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                        gt[0] = 0;
                        
                    }
                    //get second/right chromosome
                    r = strtok(NULL, "|" );
                    if(r == NULL)
                        r = p;
                    if (strcmp(r,".") == 0)
                        r[0] = '0';
                    gt[1] = atoi(r);
                    if (gt[1] > noAlt)
                    {
                        printf("*** WARNING***\n");
                        printf("\tPosition: %d, GT data for sample no %d is wrong: %s, no variant will be placed in the consensus sequence here\n", currentPos, i+1, data);
                        gt[1] = 0;
                    }
                    
                    for(int no = 0; no < 2; no++)
                        if (gt[no] == 0)
                        {
                            if(!samePos)
                            {
                                if(deleteSeq[no][i] > 0)
                                {
                                    deleteSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                            
                            
                            for (a = 0; a < noAlt; a++)
                            {
                                if(vt[a].isInsertion && !vt[a].doubleInsertion)
                                {
                                    insertSeq[no][i]=0;
                                }
                            }
                        }
                        else //gt>0
                        {
                            if(!samePos)
                            {
                                if(deleteSeq[no][i] > 0)
                                {
                                    deleteSeq[no][i]--;
                                }
                                else
                                {
                                    writeChar(refC[0], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                }
                            }
                            
                            deleteSeq[no][i] = (deleteSeq[no][i] > delLen) ? deleteSeq[no][i] : delLen ;
                            
                            if(vt[gt[no]-1].isInsertion)
                            {
                                if(vt[gt[no]-1].doubleInsertion)
                                {
                                    prevIns[no] = insertSeq[no][i];
                                }
                                else
                                {
                                    prevIns[no] = 0;
                                }
                                for (j = 1; j < strlen(vt[gt[no]-1].alt); j++)
                                {
                                    
                                    if(vt[gt[no]-1].doubleInsertion && insertSeq[no][i] > 0)
                                    {
                                        insertSeq[no][i]--;
                                    }
                                    else
                                    {
                                        writeChar(vt[gt[no]-1].alt[j], buf[no][i], bufPos[no][i], linePos[no][i], resName[no][i]);
                                    }
                                }
                                insertSeq[no][i]=(strlen(vt[gt[no]-1].alt) - 1) > prevIns[no] ? (strlen(vt[gt[no]-1].alt) - 1) : prevIns[no];
                            }
                        }
                    
                    
                }
            }
            else
            {
                printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
                printf("\tUnsupported variant VT=SV at position: %d  and not SVTYPE=DEL, INFO field:  %s\n",  pos, info);
                
                while (getc(vcf) != '\n')
                    ;
                fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
                p = strtok(format, ":");
                gt_pos = 0;
                while (strcmp(p,"GT") != 0) {
                    p = strtok(NULL, ":");
                    gt_pos++;
                }
                samePos = 0;
                if(pos == currentPos)
                    samePos = 1;
                else
                    currentPos++;
                continue;
                
            }
        }
        else
        {
            printf("*** UNSUPPORTED CASE, variant will be omitted. ***\n");
            printf("\tUnsupported variant at position: %d  (neither SNP, INDEL nor SV), INFO field: %s \n", pos, info);
            
            while (getc(vcf) != '\n')
                ;
            fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
            p = strtok(format, ":");
            gt_pos = 0;
            while (strcmp(p,"GT") != 0) {
                p = strtok(NULL, ":");
                gt_pos++;
            }
            samePos = 0;
            if(pos == currentPos)
                samePos = 1;
            else
                currentPos++;
            continue;
            
        }
		
		
		
		fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, altList, qual, filter, info, format);
		
        p = strtok(format, ":");
        gt_pos = 0;
        while (strcmp(p,"GT") != 0) {
            p = strtok(NULL, ":");
            gt_pos++;
        }
        
        
		samePos = 0;
		if(pos == currentPos)
			samePos = 1;
		else
			currentPos++;
        
        free(vt);
	}
	
	chRef = getNextSymbol(ref);
	
	while (!feof(ref))
	{
		
		for (i=0; i < noGen; i++)
		{
			if(deleteSeq[0][i] > 0)
			{
				deleteSeq[0][i]--;
			}
			else
			{
				writeChar(chRef, buf[0][i], bufPos[0][i], linePos[0][i], resName[0][i]);
			}
			
			if(deleteSeq[1][i] > 0)
			{
				deleteSeq[1][i]--;
			}
			else
            {
				writeChar(chRef, buf[1][i], bufPos[1][i], linePos[1][i], resName[1][i]);
			}
		}
		currentPos++;
		chRef = getNextSymbol(ref);
	}
	
	
	for (i=0; i < noGen; i++)
	{
        
        file = fopen(resName[0][i], "a");
        fwrite(buf[0][i], 1, bufPos[0][i], file);
        fclose(file);
        file = fopen(resName[1][i], "a");
        fwrite(buf[1][i], 1, bufPos[1][i], file);
        fclose(file);

        free(buf[0][i]);
        free(buf[1][i]);
	}
    free(buffer);
	
	
    
    printf("Processed 100%% of reference sequence\n");
	fclose(ref);
	fclose(vcf);
    return 0;
}
