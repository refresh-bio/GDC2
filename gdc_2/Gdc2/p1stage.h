/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef PRODUCER_1STAGE_H
#define PRODUCER_1STAGE_H

#include <string>
#include <numeric>
#include <algorithm>
#include <assert.h>
#include "dmutex.h"
#include "queue.h"
#include "fasta.h"
#include "hasher.h"
#include "cchasher.h"
#include "qelem.h"
#include "../libs/asmlib.h"
#include "mmgr.h"
#include "bqueue.h"

using namespace std;

struct Producer1StageParams
{
	DiskMutex* diskMtx;
	CBufferedRegisteringQueue<FileQElem>* fileQueue;
	CRegisteringQueue<QElem>* elementQueue;
	CHasher* hasher;
	uchar* seq_ref;
	uint32 min_match_1st_len;
	uint32 min_after_match;
	bool indel2;
	CMemoryMonitor* memon;
	CRegisteringQueue<string>* file_input_queue;
    uint32 no_seq;
};

class Producer1Stage
{
	DiskMutex* _diskMtx;
	CBufferedRegisteringQueue<FileQElem>* fileQueue;
	CRegisteringQueue<QElem>* _elementQueue;
	CRegisteringQueue<string>* file_input_queue;
	vector<CFastaSequence> data;
	CFastaFile* fasta;
	ConcurrentCHasher* hasher;
	CMemoryMonitor* memon;
    uint32 no_seq;
    
	bool readSeq(string filename);
    
	bool addCollectionDif(string &col_name, vector<CFastaSequence> &data, uchar* &outData, int &outDataSize);
	int32 NRunLen(uchar* seq, uint32 pos);
	void writeBE(uint32 tmp, vector<uchar> &rv);
	int writeBElen(uint32 tmp, vector<uchar> &rv);
	uint32 getMatchLen(uchar * seq1, uchar * seq2);
    
	void addOutputDataDoQueue(uchar* outData, int outDataSize, string filename, uchar eol_type, uint32 line_width, uint32 *_seq_sizes, vector<string> &_seq_names, char _new_line_at_end_of_file);
public:
	// start params
	uint32 min_match_1st_len;
	uchar *seq_ref;
	bool indel2;
	uint32 min_after_match;
    
	Producer1Stage(Producer1StageParams p)
	{
		_diskMtx = p.diskMtx;
		fileQueue = p.fileQueue;
		hasher = new ConcurrentCHasher(p.hasher);
		seq_ref = p.seq_ref;
		min_match_1st_len = p.min_match_1st_len;
		min_after_match = p.min_after_match;
		indel2 = p.indel2;
		_elementQueue = p.elementQueue;
		memon = p.memon;
		file_input_queue = p.file_input_queue;
        no_seq = p.no_seq;
		fasta = new CFastaFile();
		fasta->diskMtx = _diskMtx;
	}
    
	void start();
};

#endif