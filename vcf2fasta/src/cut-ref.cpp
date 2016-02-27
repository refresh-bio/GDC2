//
//  main.cpp
//  cutChrXY
//
//  Created by Agnieszka Danek on 21/1/13.
//  Copyright (c) 2013 Agnieszka Danek. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <string.h>
#include <stdio.h>

#define LINE_WIDTH 70


char getNextSymbol(FILE *file) {
	char ch;
	while(true) {
		ch = fgetc(file);
		if (ch != '\n')
			break;
	}
	return ch;
	
}

int main(int argc, const char * argv[])
{

    FILE *ref, *out;
    
    char ch;
    int start, end, cur = 0, lnCnt = 0;
    
    
	start = atoi(argv[3]);
    end = atoi(argv[4]);

    
    printf("input reference: %s\n", argv[1]);
	ref = fopen(argv[1], "r");
	if (ref == NULL) {
		printf("Cannot open %s\n", argv[1]);
		exit(8);
	}
    
    printf("output: %s\n", argv[2]);
	out = fopen(argv[2], "w");
	if (out == NULL) {
		printf("Cannot open %s\n", argv[2]);
		exit(8);
	}
    
    
    ch=fgetc(ref);
	//go through first line of ref
	while (ch != '\n') {
		
		fprintf(out, "%c", ch);
        ch=fgetc(ref);
	}
    fprintf(out, " %d-%d\n", start, end);
    
       
    cur = 0;
	while (cur < start){
		
		ch = getNextSymbol(ref);
        cur++;
    }
    
    while (cur >= start && cur <=end)
    {
        fwrite (&ch, 1, 1, out);
        lnCnt++;
        if (lnCnt == LINE_WIDTH) {
            lnCnt=0;
            ch = '\n';
            fwrite(&ch, 1, 1, out);
        }
        
        ch = getNextSymbol(ref);
        cur++;
    }
    
        
    
    fclose(ref);
    fclose(out);
    
    return 0;
}

