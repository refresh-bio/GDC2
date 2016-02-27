/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef COMPRESSOR_1STAGE_H
#define COMPRESSOR_1STAGE_H

#include <vector>
#include <thread>
#include <algorithm>
#include <assert.h>
#include "fasta.h"
#include "hasher.h"
#include "dmutex.h"
#include "p1stage.h"
#include "c2stage.h"
#include "queue.h"
#include "qelem.h"
#include "../libs/asmlib.h"
#include "mmgr.h"
#include "defs_gdc2.h"
#include "bqueue.h"

using namespace std;

struct Compressor1StageParams
{
	vector<string>* inputFilenameVector;
	uint32 min_match_1st_len;
	bool indel2;
    bool null_option;
	uint32 min_after_match;
	uint32 producer_thread_count;
    uint32 consumer_thread_count;
	CMemoryMonitor* memon;
	uint32 input_line_width;
	string archive_name;
	FILE** rc_file;
	uint32 compression_level;
    uint32 files_in_archive;
};

class Compressor1Stage
{
	DiskMutex diskMtx; // mutex to forward to producents (1st stage)
	CRegisteringQueue<string> filenameInputQueue; //producents queue (1st stage with names of files to read)
	CFastaFile *fasta;
	CRegisteringQueue<QElem>* elementQueue;
	CBufferedRegisteringQueue<FileQElem>* fileQueue;
    
	CHasher *hasher;
	uchar *seq_ref;
	uint32 seq_ref_total_size;
    
	uint32 producer_thread_count;
    uint32 consumer_thread_count;
	vector<string>* inputFilenameVector;
	uint32 min_match_1st_len;
	bool indel2; // parameter Producer1Stage
    bool null_option;
	uint32 min_after_match; // parameter Producer1Stage
	CMemoryMonitor* memon;
	uint32 input_line_width;
    uint32 no_seq;
    int32 ref_id_in_list;
    uint32 seq_to_decompress;
    int no_decom_seq = 0;
    std::mutex no_decom_seq_mutex;
    
	uint32 hash_step; // parameter CHasher
	uint32 hash_colisions; // parameter CHasher
	uint32 min_match_ext_len; // parameter CHasher
    
	vector<thread> threads;
    vector<CFastaSequence> data;
	// for compressor2stage
	string archive_name;
	FILE** rc_file;
	uint32 compression_level;
    uint32 files_in_archive;
    
	const int code_N_run		= 1000;
	const int code_n_run		= 1001;
	const int code_other_symbol = 1002;
	const int code_EOF			= 1003;

	void ProcessReference();
    void decompressReference();
	bool addCollectionRef1st(string &col_name, vector<CFastaSequence> &data);
	void runThreads();
	void fillInputFilenameQueue();
	void decompressData();
	void runThreadForDecompress();
    
	void runProducentThread(Producer1StageParams params);
	void runConsumerThread(Compressor2StageParams params);
	void run2StageDecompressor(Compressor2StageParams params);
    
    void decompressStage1Consumer(uint32 tot_len);
    int32 NRunLen(uchar* seq, uint32 pos);
    
    int16 getCode(uchar letter);
    bool getCode3(uchar * data, int16 * code);
    void decode3(int16 code, uchar * dest);
    uchar decode(int16 code);
    
    void print_decompress_info(char * resultFileName);
public:
	Compressor1Stage(Compressor1StageParams params) : filenameInputQueue(1)
	{
		hash_step = 1;
		hash_colisions = 0;
		min_match_ext_len = 0;
        ref_id_in_list = -1;
        
		producer_thread_count = params.producer_thread_count;
		consumer_thread_count = params.consumer_thread_count;
        memon = params.memon;
		inputFilenameVector = params.inputFilenameVector;
		min_match_1st_len = params.min_match_1st_len;
		min_after_match = params.min_after_match;
		indel2 = params.indel2;
        null_option = params.null_option;
		input_line_width = params.input_line_width;
        files_in_archive = params.files_in_archive;
		
		compression_level = params.compression_level;
		archive_name = params.archive_name;
		rc_file = params.rc_file;
        
		fasta = new CFastaFile();
		elementQueue = new CRegisteringQueue<QElem>(producer_thread_count);
		fileQueue = new CBufferedRegisteringQueue<FileQElem>(10, 0.5);
        
		seq_ref = NULL;
		hasher = NULL;
        
		memon->Increase(sizeof(fasta));
		memon->Increase(sizeof(elementQueue));
	}
    
	~Compressor1Stage()
	{
		delete fasta;
		delete elementQueue;
        
		if (seq_ref != NULL)
		{
			delete seq_ref;
			memon->Decrease(sizeof(seq_ref) + sizeof(uchar) * seq_ref_total_size);
		}
        
		if (hasher != NULL)
		{
			delete hasher;
		}
        
		memon->Decrease(sizeof(fasta));
		memon->Decrease(sizeof(elementQueue));
		memon->Decrease(sizeof(fasta));
		memon->Decrease(sizeof(elementQueue));
	}
    
	void compress();
	void decompress();
};

#endif