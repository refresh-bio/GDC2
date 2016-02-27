/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#include "c1stage.h"

void Compressor1Stage::compress()
{
	ProcessReference(); // sequence read and creating a hash array
    
	runThreads(); // forwarding parameters and start producent threads
    
	fillInputFilenameQueue(); // filling queue with file names
    
	// waiting for producent threads
	for(auto& thread : threads)
	{
        thread.join();
    }
}

/*
 filling queue with file names*/
void Compressor1Stage::fillInputFilenameQueue()
{
	for (unsigned i=1; i<inputFilenameVector->size(); ++i)
	{
		string name = inputFilenameVector->at(i);
		filenameInputQueue.Push(name);
	}
	filenameInputQueue.MarkCompleted();
}

void Compressor1Stage::runProducentThread(Producer1StageParams params)
{
	Producer1Stage p(params);
	p.start();
}

void Compressor1Stage::runConsumerThread(Compressor2StageParams params)
{
	Compressor2Stage c(params);
	c.start_compression();
}

void Compressor1Stage::run2StageDecompressor(Compressor2StageParams params)
{
	Compressor2Stage c(params);
	c.start_decompression();
}

/*
 forwarding parameters and start producent threads
 */
void Compressor1Stage::runThreads()
{
	for (uint32 i=0; i < producer_thread_count; ++i)
	{
		Producer1StageParams params;
		params.diskMtx = &diskMtx;
		params.fileQueue = fileQueue;
		params.hasher = hasher;
		params.seq_ref = seq_ref;
		params.min_match_1st_len = min_match_1st_len;
		params.min_after_match = min_after_match;
		params.indel2 = indel2;
		params.elementQueue = elementQueue;
		params.memon = memon;
        params.no_seq = no_seq;
		params.file_input_queue = &filenameInputQueue;
        
		threads.push_back(thread(&Compressor1Stage::runProducentThread, this, params));
	}
    
	Compressor2StageParams params;
	params.elementQueue = elementQueue;
	params.memon = memon;
	params.archive_name = archive_name;
	params.input_size = inputFilenameVector->size()-1;
	params.rc_file = rc_file;
	params.compression_level_threshold = (uint32)(inputFilenameVector->size()-1)*compression_level*0.1;
    params.compression_level = compression_level;
	params.no_seq = no_seq;
    threads.push_back(thread(&Compressor1Stage::runConsumerThread, this, params));
}


int32 Compressor1Stage::NRunLen(uchar* seq, uint32 pos)
{
	uint32 r = 0;
    
	uchar cN = seq[pos];
    
	if(cN != 'N' && cN != 'n')
		return 0;
    
	for(; seq[pos+r] == cN; ++r)
		;
    
	if(r < (uint32) min_match_1st_len)
		return 0;
	else if(cN == 'N')
		return (int32) r;
	else
		return - ((int32) r);
}

int16 Compressor1Stage::getCode(uchar letter)
{
    switch (letter)
    {
       
        case 'A':
            return 0;
            break;
        case 'C':
            return 1;
            break;
        case 'G':
            return 2;
            break;
        case 'T':
            return 3;
            break;
        case 'a':
            return 4;
            break;
        case 'c':
            return 5;
            break;
        case 'g':
            return 6;
            break;
        case 't':
            return 7;
            break;
        case 'N':
            return 8;
            break;
        case 'n':
            return 9;
            break;
        default:
            return -1;
            break;
    }
}

uchar Compressor1Stage::decode(int16 code1)
{
    switch (code1)
    {
            
        case 0:
            return 'A';
            break;
        case 1:
            return 'C';
            break;
        case 2:
            return 'G';
            break;
        case 3:
            return 'T';
            break;
        case 4:
            return 'a';
            break;
        case 5:
            return 'c';
            break;
        case 6:
            return 'g';
            break;
        case 7:
            return 't';
            break;
        case 8:
            return 'N';
            break;
        case 9:
            return 'n';
            break;
        default:
            return -1;
            break;
    }
}

bool Compressor1Stage::getCode3(uchar * data, int16 * code)
{
    int16 a = getCode(data[0]), b = getCode(data[1]), c = getCode(data[2]);
    if((a & 32768) || (b & 32768) || (c & 32768))
        return false;
    
    *code = a*100 + b*10 + c*1;
    return true;
}

void Compressor1Stage::decode3(int16 code, uchar * dest)
{
    int16 code1 = code/(int16)100;
    dest[0] = decode(code1);
    code = code - code1*100;
    code1 = code/(int16)10;
    dest[1] = decode(code1);
    code = code - code1*10;
    dest[2] = decode(code);
}

/* count number of occurances of each byte */
void countblock(unsigned char *buffer, freq length, freq *counters)
{   unsigned int i;
    /* first zero the counters */
    for (i=0; i<257; i++)
        counters[i] = 0;
    /* then count the number of uccurances of each byte */
    for (i=0; i<length; i++)
        counters[buffer[i]]++;
}

void readcounts(rangecoder *rc, freq *counters)
{   int i;
    for (i=0; i<256; i++)
        counters[i] = decode_short(rc);
}

void Compressor1Stage::ProcessReference()
{
	vector<CFastaSequence> data;
    
	if(!fasta->Open(inputFilenameVector->front()))
	{
		cerr << "Cannot read reference sequence" << endl;
		return;
	}
    
	if (!fasta->Read(data))
	{
		cerr << "Reading reference sequence ended with error" << endl;
		return;
	}
    
	if (!addCollectionRef1st(inputFilenameVector->front(), data))
	{
		cerr << "Cannot create hash array for reference sequence" << endl;
	}
    
#ifdef WIN32
	errno_t err;
	err = fopen_s(rc_file, (archive_name+".gdc2_ref").c_str(), "wb");
	if (err != 0)
    {
        cout << "Cannot open " << archive_name << ".gdc2_ref" << endl;
        exit(1);
    }
#else
	*rc_file = fopen((archive_name+".gdc2_ref").c_str(), "wb");
	if (*rc_file == NULL)
	{
        cout << "Cannot open " << archive_name << ".gdc2_ref" << endl;
        exit(1);
    }
#endif
    
    no_seq = data.size();
    uint32 tot_size = 0;
    
    uint32 max_len = 0;
    
    for(size_t i = 0; i < data.size(); i++)
    {
        tot_size += data[i].size;
        max_len = max_len > data[i].line_len ? max_len : data[i].line_len;
    }
    
    uint32 buf_pos = 0;
    uint64 buf_size = inputFilenameVector->front().length() + 10; //exactly + 7; // + 1 + 1 + 4 + 1; // '\0', eol, line_len, last '\0'
    for(size_t i = 0; i < data.size(); i++)
    {
        buf_size += data[i].seq_name.length() + 5;
	}
    uchar * buf = (uchar*)malloc(buf_size);
    
    buf[0] = '\0';
    strcpy((char*)buf, inputFilenameVector->front().c_str());
    buf_pos += inputFilenameVector->front().length();
    buf[buf_pos++] = '\0';
    
    cout << "Processing reference " << buf << endl;
    //info about eol and line_len
    buf[buf_pos++]= data[0].eol_type;
    memcpy(buf + buf_pos, &max_len, 4);
    buf_pos += 4;
    
    for(size_t i = 0; i < data.size(); i++)
    {
        strcpy((char*)buf + buf_pos, data[i].seq_name.c_str());
        buf_pos += strlen(data[i].seq_name.c_str());
        buf[buf_pos++]='\0';
        memcpy(buf + buf_pos, &data[i].size, 4);
        buf_pos += 4;

    }
    buf[buf_pos++]='\0';
   
    fwrite(buf, 1, buf_pos, *rc_file);
    free(buf);
    
    int16 code;
    int32 len;
    
    rangecoder rc_ref;
    qsmodel qsm;
    int32 syfreq, ltfreq;

	initqsmodel(&qsm, 1004, 12, 2048, NULL, 1);

	start_encoding(&rc_ref, 0, 0);

	for(uint32 i = 0; i < tot_size;)
    {
        if((len = NRunLen(data[0].raw_data, i)) != 0)			// run of Ns or ns
        {
            if(len > 0)			// N run
			{
		        qsgetfreq(&qsm, code_N_run, &syfreq, &ltfreq);
				encode_shift(&rc_ref, syfreq, ltfreq, 12);
				qsupdate(&qsm, code_N_run);
			}
            else				// n run
            {
		        qsgetfreq(&qsm, code_n_run, &syfreq, &ltfreq);
				encode_shift(&rc_ref, syfreq, ltfreq, 12);
				qsupdate(&qsm, code_n_run);

				len = -len;
            }
			encode_byte(&rc_ref, len >> 24);
			encode_byte(&rc_ref, (len >> 16) & 0xff);
			encode_byte(&rc_ref, (len >> 8) & 0xff);
			encode_byte(&rc_ref, len & 0xff);
			i += len;
        }
        else
        {
            if(getCode3(data[0].raw_data+i, &code))		// next three letters are valid
            {
		        qsgetfreq(&qsm, code, &syfreq, &ltfreq);
				encode_shift(&rc_ref, syfreq, ltfreq, 12);
				qsupdate(&qsm, code);

				i += 3;
            }
            else										// some special symbol is found in the next three letters
            {
		        qsgetfreq(&qsm, code_other_symbol, &syfreq, &ltfreq);
				encode_shift(&rc_ref, syfreq, ltfreq, 12);
				qsupdate(&qsm, code_other_symbol);

				encode_byte(&rc_ref, data[0].raw_data[i]);
				i++;
            }
        }
    }
   
	qsgetfreq(&qsm, code_EOF, &syfreq, &ltfreq);
    encode_shift(&rc_ref, syfreq, ltfreq, 12);
    
    done_encoding(&rc_ref);
    deleteqsmodel(&qsm);
	fclose(*rc_file);
	fasta->Close();
}

void Compressor1Stage::decompressReference()
{
#ifdef WIN32
	errno_t err;
	err = fopen_s(rc_file, (archive_name+".gdc2_ref").c_str(), "rb");
	if (err != 0)
    {
        cout << "Cannot open " << archive_name << ".gdc2_ref" << endl;
        exit(1);
    }
#else
	*rc_file = fopen((archive_name + ".gdc2_ref").c_str(), "rb");
	if (*rc_file == NULL)
	{
        cout << "Cannot open " << archive_name << ".gdc2_ref" << endl;
        exit(1);
    }
#endif
    data.push_back(CFastaSequence());
    int c;
    string filename;
    
    while((c = getc(*rc_file)) != '\0')
	{
		filename.push_back(c);
	}
    fread(&data[0].eol_type, 1, 1, *rc_file);
    fread(&data[0].line_len, 1, 4, *rc_file);
    
    uint32 tot_size = 0;
    data[0].seq_name = "";
    while((c = getc(*rc_file)) != '\0')
	{
		data[0].seq_name.push_back(c);
	}
    fread(&data[0].size, 1, 4, *rc_file);
    tot_size += data[0].size;
    
    int j = 1;
    no_seq = 1;
    while ((c = getc(*rc_file)) != '\0')
    {
        data.push_back(CFastaSequence());
        no_seq++;
        data[j].line_len = data[0].line_len;
        data[j].eol_type = data[0].eol_type;
        do
        {
            data[j].seq_name.push_back(c);
        }
        while((c = getc(*rc_file)) != '\0');
        fread(&data[j].size, 1, 4, *rc_file);
        data[j].start_pos = tot_size;
        tot_size += data[j].size;
        j++;
    }

    rangecoder rc_ref;
    qsmodel qsm;
    int32 code, syfreq, ltfreq;
    data[0].raw_data = new uchar[tot_size];
    uint32 buf_pos = 0;
	uchar *buf = data[0].raw_data;
    /* init the model the same as in the compressor */
    initqsmodel(&qsm, 1004, 12, 2048, NULL, 0);
    start_decoding(&rc_ref);
    
    while (1)
    {
        ltfreq = decode_culshift(&rc_ref, 12);
        code = qsgetsym(&qsm, ltfreq);
        qsgetfreq(&qsm, code, &syfreq, &ltfreq);
		decode_update(&rc_ref, syfreq, ltfreq, 1<<12);
		qsupdate(&qsm, code);
		if (code == code_EOF)  /* check for end-of-file */
            break;

		if(code == code_N_run || code == code_n_run)
		{
			int32 len;
			len = decode_byte(&rc_ref);
			len = (len << 8) + decode_byte(&rc_ref);
			len = (len << 8) + decode_byte(&rc_ref);
			len = (len << 8) + decode_byte(&rc_ref);
			uchar c = (code == code_N_run ? 'N' : 'n');
			for(int i = 0; i < len; ++i)
				*buf++ = c;
			buf_pos += len;
		}
		else if(code == code_other_symbol)
		{
			*buf++ = decode_byte(&rc_ref);
			buf_pos++;
		}
		else // valid symbols
		{
			decode3(code, buf);
			buf += 3;
			buf_pos += 3;
		}
	}

    done_decoding(&rc_ref);
    deleteqsmodel(&qsm);
   
    if(inputFilenameVector->size()  != 0)
    {
        for (unsigned i=0; i<inputFilenameVector->size(); ++i)
        {
            if(strcmp(filename.c_str(), inputFilenameVector->at(i).c_str()) == 0)
            {
                ref_id_in_list = i;
                if(!null_option)
                    fasta->Create(filename+".ori");
                else
                    fasta->Create("/dev/null");
                fasta->Write(data);
                cout << "Reference " << filename  <<  ".ori decompressed\n";
                break;
            }
        }
    }
    else
    {
        if(!null_option)
            fasta->Create(filename+".ori");
        else
            fasta->Create("/dev/null");
        
        fasta->Write(data);
        cout << "Reference " << filename  <<  ".ori decompressed\n";
    }
    fclose(*rc_file);
    fasta->Close();
}

/*
 Utworzenie hash array dla sekwencji referencyjnej
 */
bool Compressor1Stage::addCollectionRef1st(string &col_name, vector<CFastaSequence> &data)
{
    vector<uint32> s_pos;
    
	// Copy main ref. seq. data
	seq_ref_total_size = 0;
	for(size_t i = 0; i < data.size(); ++i)
		seq_ref_total_size += data[i].size;
    
    seq_ref = new uchar[seq_ref_total_size+hash_step+(uint32)4];
	A_memcpy(seq_ref+1, data[0].raw_data, seq_ref_total_size);
    
	seq_ref[seq_ref_total_size+1] = '\0';
	fill_n(seq_ref+seq_ref_total_size+2, hash_step+2, '~');
	hasher = new CHasher();
	hasher->SetParams(hash_step, hash_colisions, min_match_1st_len, min_match_ext_len);
	hasher->SetRefSeq(seq_ref, seq_ref_total_size);
	memon->Increase(sizeof(seq_ref) + sizeof(uchar) * seq_ref_total_size);
	return true;
}

void Compressor1Stage::decompress()
{
    decompressReference();
    runThreadForDecompress();
    Compressor1Stage::decompressData();
}

void Compressor1Stage::runThreadForDecompress()
{
	Compressor2StageParams params;
	params.elementQueue = elementQueue;
	params.memon = memon;
	params.archive_name = archive_name;
	params.rc_file = rc_file;
    params.no_seq = no_seq;
    params.inputFilenameVector = inputFilenameVector;
    params.ref_id_in_list = ref_id_in_list;
	threads.push_back(thread(&Compressor1Stage::run2StageDecompressor, this, params));
}

void Compressor1Stage::decompressData()
{
    if(inputFilenameVector->size() > 0)
        seq_to_decompress = inputFilenameVector->size() - (ref_id_in_list >= 0);
    else
        seq_to_decompress = files_in_archive;
    
	uint32 tot_len = 0;
	for (size_t i = 0; i < data.size(); ++i)
		tot_len += data[i].size;
    
    for (uint32 i=0; i< consumer_thread_count ; ++i)
    {
        threads.push_back(thread(&Compressor1Stage::decompressStage1Consumer, this, tot_len));
    }
    // waiting for consumers
	for(auto& thread : threads)
	{
        thread.join();
    }
}


void Compressor1Stage::decompressStage1Consumer(uint32 tot_len)
{
    QElem e;
    FILE *fp;
    char * resultFileName;
	uint32  ref_pos, len;
    int j;
	char ch;
    int pos_dif = 5;
    uint32 buf_pos = 0;
    char *buffer = (char *)malloc (BUF_SIZE+1);
    uint32 file_input_line_width;
    vector<CFastaSequence>  dif_data;
	uchar *raw_data_dif;
	CFastaFile *dif_fasta = new CFastaFile();
	raw_data_dif = new uchar[tot_len + tot_len];
    
    while (!elementQueue->IsCompleted())
    {
        if (elementQueue->Pop(e))
        {
            pos_dif = 0;
            resultFileName = new char[e.filename.size() + 8];
            resultFileName[0] = '\0';
            strcat(resultFileName, e.filename.c_str());
            strcat(resultFileName, ".ori");
            strcat(resultFileName, "\0");
            
            if(!null_option)
                fp = my_fopen(resultFileName, "wb");
            else
                fp = my_fopen("/dev/null", "wb");            
            
            if (fp == NULL)
            {
                cout << "Cannot open " << resultFileName  << endl;
                exit(6);
            }
            setvbuf(fp, NULL, _IOFBF, BUF_SIZE);
            
            for(int i = 0; i < e.dataSize; )
            {
                ch = e.data[i];
                i++;
                if ((ch & 128))
                {
                    ch = ch & 127;
                    raw_data_dif[pos_dif++] = ch;
                }
                else
                {
                    A_memcpy(&ref_pos, &e.data[i], 3); i = i+3;
                    ref_pos = (ref_pos & 0x0000FF00) | ((ref_pos & 0x000000FF) << 16) | ((ref_pos & 0x00FF0000) >> 16) | ((ch << 24) & 0xFF000000);
                    A_memcpy(&len, &e.data[i], 3); i = i+3;
                    len = ((len & 0x000000FF) << 16) | ((len & 0x00FF0000) >> 16) | (len & 0x00FF00);
                    if (ref_pos > POS_MAX)
                    {
                        if (ref_pos == POS_MAX + (uint32)2)
                        {
                            ch = 'N';
                            for (uint32 j = 0; j < len; j++)
                            {
                                
                                raw_data_dif[pos_dif++] = ch;
                            }
                        }
                        else if (ref_pos == POS_MAX + (uint32)1)
                        {
                            ch = 'n';
                            for (uint32 j = 0; j < len; j++)
                            {
                                
                                raw_data_dif[pos_dif++] = ch;
                            }
                        }
                    }
                    else
                    {
                        copy_n(data[0].raw_data + ref_pos - 1, len, raw_data_dif + pos_dif);
                        pos_dif = pos_dif + len;
                    }
                }
            }
            int end = 0, lastEnd = 0;
            file_input_line_width = e.line_width;
            unsigned int limit = BUF_SIZE - file_input_line_width - 1;
            j = 0;
            for(uint32 seq = 0; seq < no_seq; seq++)
            {
                if(seq)
                {
                    lastEnd += e.seq_sizes[seq-1];
                    buffer[buf_pos++] =  e.eol_type;
                }
                buffer[buf_pos++] = '>';
                buffer[buf_pos] = '\0';
                strcat(buffer+buf_pos,  e.seq_names.at(seq).c_str());
                buf_pos = buf_pos +  e.seq_names.at(seq).length();
                buffer[buf_pos++] = '\n';
                end = ( e.seq_sizes[seq] < file_input_line_width)? 0 : e.seq_sizes[seq] - file_input_line_width;
                if( e.eol_type != 255)
                {
                    for (j = lastEnd; j < end+lastEnd;)
                    {
                        if(buf_pos >= limit)
                        {
                            fwrite(buffer, 1, buf_pos, fp);
                            buf_pos = 0;
                        }
                        A_memcpy(buffer+buf_pos, raw_data_dif + j, file_input_line_width);
                        buf_pos = buf_pos + file_input_line_width;
                        buffer[buf_pos++] =  e.eol_type;
                        j = j + file_input_line_width;
                    }
                    if(buf_pos >= limit)
                    {
                        fwrite(buffer,  1, buf_pos, fp);
                        buf_pos = 0;
                    }
                    
                    A_memcpy(buffer+buf_pos, raw_data_dif + j, e.seq_sizes[seq] + lastEnd - j);
                    buf_pos = buf_pos + e.seq_sizes[seq]+lastEnd - j;
                }
                else
                {
                    limit = BUF_SIZE - file_input_line_width - 2;
                    for (j = lastEnd; j < end+lastEnd;)
                    {
                        if(buf_pos >= limit)
                        {
                            fwrite(buffer,  1, buf_pos, fp);
                            buf_pos = 0;
                        }
                        A_memcpy(buffer+buf_pos, raw_data_dif + j, file_input_line_width);
                        buf_pos = buf_pos + file_input_line_width;
                        buffer[buf_pos++] =  0x0d;
                        buffer[buf_pos++] =  0x0a;
                        j = j + file_input_line_width;
                    }
                    if(buf_pos >= limit)
                    {
                        fwrite(buffer,  1, buf_pos, fp);
                        buf_pos = 0;
                    }
                    A_memcpy(buffer+buf_pos, raw_data_dif + j, e.seq_sizes[seq] + lastEnd - j);
                    buf_pos = buf_pos + e.seq_sizes[seq]+lastEnd - j;
                }
            }
            fwrite(buffer, 1, buf_pos, fp);
    
            //add new line ad the end of file
            if(e.new_line_at_end_of_file )//&& e.seq_sizes[no_seq-1] + lastEnd - j != file_input_line_width)
            {
                if(e.eol_type == 0xff)
                {
                    putc(0x0d, fp);
                    putc(0x0a, fp);
                }
                else
                    putc(e.eol_type, fp);
            }
            
            buf_pos = 0;
            
            fclose(fp);
            print_decompress_info(resultFileName);
            delete[] resultFileName;
            delete [] e.seq_sizes;
        }
    }
    delete dif_fasta;
	delete raw_data_dif;
}

void Compressor1Stage::print_decompress_info(char * resultFileName)
{
    std::lock_guard<std::mutex> lock(no_decom_seq_mutex);
    ++no_decom_seq;
    cout << no_decom_seq << " of " << seq_to_decompress << " sequences decompressed [last: " << resultFileName <<  "]" << endl;
}

