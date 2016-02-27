/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */


#ifndef QELEM_H
#define QELEM_H

#include <string>
#include "defs.h"
#include "fasta.h"

using namespace std;

struct QElem
{
	string filename;
	uchar* data;
    uchar eol_type;
	int dataSize;
    uint32 * seq_sizes;
    uint32 line_width;
    vector<string> seq_names;
    char new_line_at_end_of_file;
};

struct FileQElem
{
	string filename;
	CFastaFile* fFile;
	vector<CFastaSequence>* fData;
};

#endif