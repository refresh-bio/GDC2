/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#include <vector>
#include <iostream>
#include <string>
#include "defs_gdc2.h"
#include "defs.h"
#include "c1stage.h"
#include "timer.h"
#include "../libs/zlib.h"

using namespace std;

vector<string> file_names;
enum t_working_mode working_mode;

// input parameters
int32
extr_col_no,
extr_from,
extr_to,
ref_col_no,
min_match_1st_len,
min_after_match,
line_width,
threads_count,
compression_level;
bool null_option;

uint64_t max_memory;

uchar compress_mode;
bool indel2;
string archive_name;

CMemoryMonitor* memon;

FILE *rc_file; // for arithmetic coding

void show_usage();
uint32 list_archive_files(bool print);


void save_byte(int x)
{
	putc(x, rc_file);
}

// ***************************************************************
int read_byte()
{
	int a = getc(rc_file);
	return a;
}

void compress()
{
	Compressor1StageParams params;
	params.inputFilenameVector = &file_names;
	params.min_match_1st_len = min_match_1st_len;
	params.indel2 = indel2;
	params.min_after_match = min_after_match;
	params.producer_thread_count = threads_count - 1;
	params.input_line_width = line_width;
	params.memon = memon;
	params.archive_name = archive_name;
	params.rc_file = &rc_file;
	params.compression_level = compression_level;
    
	Compressor1Stage compressor(params);
	
	compressor.compress();
}

void decompress()
{
	Compressor1StageParams params;
	params.min_match_1st_len = min_match_1st_len;
	params.indel2 = indel2;
	params.min_after_match = min_after_match;
    params.producer_thread_count = 1;
    params.consumer_thread_count = threads_count- 1;
	params.input_line_width = line_width;
	params.memon = memon;
	params.archive_name = archive_name;
	params.rc_file = &rc_file;
	params.compression_level = compression_level;
    params.inputFilenameVector = &file_names;
    params.files_in_archive = list_archive_files(false);
    params.null_option = null_option;
    
	Compressor1Stage compressor(params);
	
	compressor.decompress();
}

uint32 list_archive_files(bool print)
{
    FILE *ref;
    
#ifdef WIN32
	errno_t err;
	err = fopen_s(&ref, (archive_name+".gdc2_ref").c_str(), "rb");
	if (err != 0)
    {
        cout << "Cannot open " << archive_name << ".gdc2_ref" << endl;
        exit(1);
    }
#else
	ref = fopen((archive_name + ".gdc2_ref").c_str(), "rb");
	if (ref == NULL)
	{
        cout << "Cannot open " << archive_name << ".gdc2_ref" << endl;
        exit(1);
    }
#endif
    
    int c;
    string filename;
    
    while((c = getc(ref)) != '\0')
	{
		filename.push_back(c);
	}
    
    fclose(ref);
    
#ifdef COMPRESS_DESC
	gzFile_s *in_desc;
#else
	FILE *in_desc;
#endif

	// Read archive description
#ifdef COMPRESS_DESC
	in_desc = gzopen((archive_name + ".gdc2_desc").c_str(), "rb");
#else
#ifdef WIN32
	fopen_s(&in_desc, (archive_name + ".gdc2_desc").c_str(), "rb");
#else
	in_desc = fopen((archive_name + ".gdc2_desc").c_str(), "rb");
#endif
#endif   
	if(!in_desc)
    {
        cout << "Cannot open " << archive_name << ".gdc2_desc" << endl;
        exit(1);
    }
    
    string s;
    uint32 counter = 1;

	int size, _line_width;
    uchar eol_type;
    char new_line_at_end_of_file;
    
    if(print)
    {
        cout << "List of files in archive ``" << archive_name << "'':" << endl;
        cout << "0. (reference) " << filename << endl;
	}
     my_fread(&compression_level, 1, 1, in_desc);
    
	while (1)
	{
		while ((c = my_getc(in_desc)) != EOF)
		{
			if (c != 0)
				s.push_back(c);
			else
				break;
		}
		if (s.length() == 0)
			break;

        my_fread(&size, 1, 4, in_desc);
        my_fread(&eol_type, 1, 1, in_desc);
        my_fread(&new_line_at_end_of_file, 1, 1, in_desc);
        my_fread(&_line_width, 1, 4, in_desc);

        while ((c = my_getc(in_desc)) != 1)
		{
			// Skip sequence name
			while ((c = my_getc(in_desc)) != 0) //read info about sequences in file
                ;
			;
			// Skip size
			for (int i = 0; i < 4; ++i)
				my_getc(in_desc);
		}

        //list
        if(print)
            cout << counter << ". " << s.c_str() << endl;
        counter++;
        s.clear();
	}
    
	my_fclose(in_desc);
    return counter-1;
}

bool parse_parameters(int argc, char *argv[])
{
	if(argc < 3)
    {
        return false;
    }
    
	if(strcmp(argv[1], "c") == 0)
		working_mode = wm_compress;
	else if(strcmp(argv[1], "d") == 0)
		working_mode = wm_decompress;
	else if(strcmp(argv[1], "l") == 0)
		working_mode = wm_list;
	else
		return false;
    
	extr_col_no = -1;
	extr_from   = -1;
	extr_to		= -1;
	ref_col_no  = -1;
    
	min_match_1st_len = 15;
    min_after_match = 4;
	compress_mode = 1;
    line_width = 70;
    indel2 = false;
    null_option = false;
	threads_count = 4; // default 4
	max_memory = 1024; // default 1024MB
	compression_level = 10; // default level: 10
    
    int32 i;
    for(i = 2; i < argc && argv[i][0] == '-'; ++i)
	{
		if(strncmp(argv[i], "-mp", 3) == 0)
		{
			char *p1 = strchr(argv[i], ',');
			if(!p1)
				continue;
			
			*p1 = 0;
			min_match_1st_len = atoi(&argv[i][3]);
            min_after_match = atoi(p1+1);
		}
        else if(strncmp(argv[i], "-i2", 3) == 0)
			indel2 = true;
        else if(strncmp(argv[i], "-lw", 3) == 0)
			line_width = atoi(&argv[i][3]);
		else if (strncmp(argv[i], "-t", 2) == 0)
			threads_count = atoi(&argv[i][2]);
		else if (strncmp(argv[i], "-mm", 3) == 0)
			max_memory = atoi(&argv[i][3]);
		else if (strncmp(argv[i], "-l", 2) == 0)
			compression_level = atoi(&argv[i][2]);
        else if (strncmp(argv[i], "-null", 5) == 0)
			null_option = true;
        else
        {
                cout << "Unknown option \"" << argv[i] << "\" will ge ignored\n";
        }
    }
    
	if (threads_count < 2)
	{
		cerr << "Wrong numer of threads. Minimum is 2" << endl;
		return false;
	}
    
	if (max_memory < 1)
	{
		cerr << "Wrong amount of memory limit" << endl;
		return false;
	}
	else
	{
		max_memory *= 1024 * 1024; // MB=1024KB=1024*1024B
	}
    
	if (compression_level < 1)
		compression_level = 1;
	if (compression_level > 10)
		compression_level = 10;
    
    // No list file name
	if(argc < i )
    {
        return false;
    }
    
    if (i >= argc)
    {
        return 0;
    }
    archive_name = argv[i++];
    cout<< "Archive name: " << archive_name  << endl;;
	
	if(argc == i)
		file_names.clear();
	else
	{
		if(argv[i][0] == '@')
		{
			FILE *f_names = fopen(argv[i]+1, "rt");
			if(!f_names)
			{
				cout << "No file: " << string(argv[i]+1) << "\n";
				return false;
			}
			file_names.clear();
			char name[1024];
			while(fscanf(f_names, "%s", name) != EOF)
				file_names.push_back(string(name));
			fclose(f_names);
		}
		else
		{
			file_names.assign(argv+i, argv+argc);
		}
	}
    
    if(working_mode == wm_compress)
    {
        cout << "Number of files to compress (including reference): " << file_names.size() << endl;
        
        if(file_names.size() == 0)
        {
            cout << "No files to compress!" << endl << endl;
            show_usage();
            return false;
        }
    }
 
	return true;
}

void show_usage()
{
    cout << "Genome Differential Compressor v. 2.0\n";
	cout << "Usage: gdc2 <mode> [options] <archive_name> [@list_of_files | file1_name file2_name ...]\n";
	cout << "  mode - c (compress) archive_name [@list_of_files | file1_name file2_name ...], d (decompress) archive_name [@list_of_files | file1_name file2_name ...], l (list files in archive) archive_name\n";
    cout << "  archive_name - archive name, required\n";
	cout << "  list_of_files - file containing a list of files (1 file per line) to compress\n";
    cout << "  file1_name file2_name file3_name ... - list of files (separated by space) to compress\n";
    cout << " (if some list of files is not provided for decompression, all files will be decompressed)\n";
	cout << "  options:\n";
	cout << "    -mpX,Y - match parameters (default: 15,4):\n";
	cout << "      X - min. len. of 1st part of a match (hashing)\n";
	cout << "      Y - min. len. of 2nd and next parts of a match (after SNP/INS/DEL)\n";
	cout << "    -i2 - additionally look for insertions/deletions (indels) of length 2 (default 1 only)\n";
	//cout << "    -lwX - Set line width in the output file (decompression mode, default 70)\n";
	cout << "    -tX - set number of working threads to X (minimum is 2, default is 4)\n";
	cout << "    -mmX - set memory limit in MB that program can allocate to X (default is 1024MB)\n";
	cout << "    -lX - set compression degree (defining percentage of sequences used in the second level compression) to X, where X is an integer number in range [1-10] (default is 10); ; X*10 percent of sequences will be used in the second level compression. \n";
    cout << "Examples:\n";
	cout << "  gdc2 c my_archive @files_list\n";
	cout << "  gdc2 c -mp15,4 -i2 my_archive @files_list\n";
	cout << "  gdc2 d my_archive @files_list\n";
    cout << "  gdc2 d my_archive\n";
    cout << "  gdc2 l my_archive\n";
}

void count_memory()
{
	uint64_t size = 0;
	for (auto& s : file_names)
	{
		size += s.capacity() * sizeof(char);
	}
	size += sizeof(file_names) + sizeof(string) * file_names.capacity();
	memon->Increase(size);
    
	size = 0;
}

int main(int argc, char * argv[])
{
	SetMemcpyCacheLimit(8*sizeof(char));
    
	if(!parse_parameters(argc, argv))
	{
        show_usage();
        return 0;
	}
    
    
	memon = new CMemoryMonitor(max_memory);
	count_memory();
    
	CStopWatch clock;
	clock.StartTimer();
    
	if(working_mode == wm_compress)
		compress();
	else if(working_mode == wm_decompress)
        decompress();
    else if(working_mode == wm_list)
		list_archive_files(true);
    
	clock.StopTimer();
	cout << clock.GetElapsedTime() << endl;
	
	delete memon;
    
	return 0;
}