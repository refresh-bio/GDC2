#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <iostream>

#define NO_GENOMES 1092
#define MAL_START 2699521
#define MAL_END 154931043


int main (int argc, char * const argv[]) {
    
    
    
	FILE * vcf;
	FILE *Xfem, *Xmal, *Xmal1, *Xmal2;
	
	
    bool isMale[1092] = {true, false, false, false, true, false, true, false, false, true, true, false, false, true, true, true, true, true, false, true, false, false, false, false, false, false, true, false, false, true, false, true, false, false, false, true, false, true, true, true, true, true, true, false, true, true, false, true, true, false, true, true, false, true, true, false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true, false, false, false, true, false, false, false, false, false, false, true, true, true, false, true, false, false, false, true, true, false, false, false, true, false, false, false, true, false, false, false, true, true, false, true, false, false, false, true, false, true, false, false, false, true, true, true, false, false, true, false, false, false, true, true, true, false, false, false, false, false, true, false, false, true, false, false, false, true, false, false, false, false, true, true, false, true, false, true, true, false, false, true, false, false, false, true, false, false, false, false, true, false, true, false, false, false, true, false, true, true, false, true, false, false, false, false, true, false, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, true, false, false, true, true, true, false, false, true, false, true, false, true, false, true, false, true, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, true, false, false, true, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, true, false, true, false, true, false, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, true, true, false, true, true, false, true, false, false, true, true, false, true, true, false, true, true, false, true, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, true, false, true, false, false, true, true, true, false, true, false, false, true, false, true, false, true, false, true, true, false, false, true, false, true, false, true, true, false, false, true, false, true, false, false, false, true, true, false, true, true, false, true, true, false, true, false, true, false, true, true, false, true, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, false, false, false, true, false, true, false, true, false, false, false, false, false, true, true, false, true, false, true, true, false, false, false, false, true, true, true, false, true, true, true, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, true, false, true, true, true, false, false, false, false, false, false, true, true, true, true, true, false, false, false, false, false, true, true, false, true, true, true, true, true, false, false, false, true, true, true, true, true, true, true, true, true, true, true, false, false, false, true, false, true, false, true, false, true, false, true, false, false, true, true, false, false, true, false, true, false, false, true, true, true, false, false, true, false, false, false, true, true, false, false, false, true, true, true, true, false, false, true, true, false, true, false, true, false, false, true, false, false, false, true, true, true, true, true, false, true, true, true, false, true, true, false, false, true, false, false, true, true, true, true, false, true, true, true, true, false, false, true, true, true, false, true, true, false, true, false, true, true, true, false, false, true, true, true, true, true, false, true, true, false, false, true, false, false, true, true, false, true, false, true, false, false, true, true, false, false, true, false, true, false, false, true, false, true, true, false, true, false, false, true, true, false, false, true, false, true, true, false, true, false, true, false, false, true, true, false, true, false, true, false, true, false, false, true, false, true, true, false, true, true, true, false, true, true, false, false, false, true, true, true, false, false, false, false, true, false, true, false, true, true, true, false, true, false, true, true, true, true, true, true, true, true, false, false, true, false, true, true, true, true, false, false, true, true, false, false, true, false, false, false, false, false, true, true, true, false, false, false, false, false, false, false, false, true, true, false, false, true, false, true, true, true, true, false, false, true, false, false, true, false, false, false, false, false, false, false, false, false, false, false, true, false, true, false, false, true, false, true, false, false, true, false, true, false, true, false, true, true, false, true, false, false, true, false, false, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, false, true, false, true, false, true, false, true, false, true, false, true, false, true, false, true, true, false, true, false, true, false, true, true, false, false, true, false, true, false, true, false, true, true, false, true, false, false, true, true, false, false, false, true, false, false, true, false, false, false, false, false, false, false, false, true, false, true, false, true, true, true, true, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, true, true, false, true, true, true, true, false, true, true, true, true, false, false, false, true, false, true, false, true, true, true, true, false, false, false, true, true, true, false, false, true, true, false, true, false, true, true, false, false, true, true, false, false, true, false, false, false, true, false, false, false, false, false, true, true, true, false, true, false, true, false, true, false, true, false, false, true, false, true, false, true, true, false, false, true, true, true, true, false, true, true, true, false, false, false, false};
    
    
	int  i;
	char chrom[15];
    
	
	int pos;
	char * resultFileName;
	
	
    char  name[1000], refC[1000000], alt[1000000], filter[20], info[1000], format[100], data[200], qual[20];
    
	printf("input vcf: %s\n", argv[1]);
	vcf = fopen(argv[1], "r");
	if (vcf == NULL) {
		printf("Cannot open %s\n", argv[1]);
		exit(8);
	}
	

	resultFileName = new char[strlen(argv[1])+ 20];
    
    resultFileName[0]='\0';
    strcat(resultFileName, argv[1]);
    strcat(resultFileName, "-mal");
    strcat(resultFileName, "\0");
    Xmal = fopen(resultFileName, "w");
    if (Xmal == NULL) {
        printf("Cannot open %s\n", resultFileName);
        exit(8);
    }

    resultFileName[0]='\0';
    strcat(resultFileName, argv[1]);
    strcat(resultFileName, "-mal1");
    strcat(resultFileName, "\0");
    Xmal1 = fopen(resultFileName, "w");
    if (Xmal1 == NULL) {
        printf("Cannot open %s\n", resultFileName);
        exit(8);
    }

    resultFileName[0]='\0';
    strcat(resultFileName, argv[1]);
    strcat(resultFileName, "-mal2");
    strcat(resultFileName, "\0");
    Xmal2 = fopen(resultFileName, "w");
    if (Xmal2 == NULL) {
        printf("Cannot open %s\n", resultFileName);
        exit(8);
    }

    resultFileName[0]='\0';
    strcat(resultFileName, argv[1]);
    strcat(resultFileName, "-fem");
    strcat(resultFileName, "\0");
    Xfem = fopen(resultFileName, "w");
    if (Xfem == NULL) {
        printf("Cannot open %s\n", resultFileName);
        exit(8);
    }

    
    char c;

    fscanf(vcf, "%s", data);
    while (strncmp(data, "#CHROM", 6) != 0)
    {
        fprintf(Xmal, "%s", data);
        fprintf(Xmal1, "%s", data);
        fprintf(Xmal2, "%s", data);
        fprintf(Xfem, "%s", data);
        
        c = fgetc(vcf);
        while (c != '\n') {
            fputc(c, Xmal);
            fputc(c, Xmal1);
            fputc(c, Xmal2);
            fputc(c, Xfem);
            c = fgetc(vcf);
        }
        fputc(c, Xmal);
        fputc(c, Xmal1);
        fputc(c, Xmal2);
        fputc(c, Xfem);
        
        fscanf(vcf, "%s", data);
    }
    

    
	fscanf(vcf, "%s %s %s %s %s %s %s %s", data, name, refC, alt, qual, filter, info, format);
    fprintf(Xmal, "#CHROM\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", data, name, refC, alt, qual, filter, info, format);
    fprintf(Xmal1, "#CHROM\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", data, name, refC, alt, qual, filter, info, format);
    fprintf(Xmal2, "#CHROM\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", data, name, refC, alt, qual, filter, info, format);
    fprintf(Xfem, "#CHROM\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", data, name, refC, alt, qual, filter, info, format);
    
    for (i=0; i < NO_GENOMES; i++)
	{
		fscanf(vcf, "%s", data);
        if (isMale[i])
        {
            fprintf(Xmal, "\t%s", data);
            fprintf(Xmal1, "\t%s", data);
            fprintf(Xmal2, "\t%s", data);
        }
        else
            fprintf(Xfem, "\t%s", data);
    }
    
    fprintf(Xmal, "\n");
    fprintf(Xmal1, "\n");
    fprintf(Xmal2, "\n");
    fprintf(Xfem, "\n");
    

	
	fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, alt, qual, filter, info, format);
	    
    while (!feof(vcf))
	{
        
        fprintf(Xfem, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chrom, pos, name, refC, alt, qual, filter, info, format);
        
        if (pos < MAL_START)
            fprintf(Xmal1, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chrom, pos, name, refC, alt, qual, filter, info, format);
        else if (pos <= MAL_END)
            fprintf(Xmal, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chrom, pos, name, refC, alt, qual, filter, info, format);
        else
            fprintf(Xmal2, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", chrom, pos, name, refC, alt, qual, filter, info, format);
        

        
		for (i=0; i < NO_GENOMES; i++)
        {
            fscanf(vcf, "%s", data);
            
            if (isMale[i])
            {
                if (pos < MAL_START)
                    fprintf(Xmal1, "\t%s", data);
                else if (pos <= MAL_END)
                    fprintf(Xmal, "\t%s", data);
                else
                    fprintf(Xmal2, "\t%s", data);
            }
            else
                fprintf(Xfem, "\t%s", data);
        }
        
        if (pos < MAL_START)
            fprintf(Xmal1, "\n");
        else if (pos <= MAL_END)
            fprintf(Xmal, "\n");
        else
            fprintf(Xmal2, "\n");
       
            fprintf(Xfem, "\n");
        
		fscanf(vcf, "%s %d %s %s %s %s %s %s %s", chrom, &pos, name, refC, alt, qual, filter, info, format);
		
	}
	
    fclose(vcf);
	fclose(Xmal);
	fclose(Xmal1);
	fclose(Xmal2);
	fclose(Xfem);
    return 0;
}
