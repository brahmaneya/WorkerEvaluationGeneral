#ifndef TEMPORALDATA_H
#define TEMPORALDATA_H

#include<iostream>
#include<stdlib.h>
#include"worker.h"
#include"pool.h"
#include<string>
#include<math.h>
#include"test.h"
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<cassert>

void get_temporal_data(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans);

void get_temporal_data_regular(int& n, int& wn, bool** &ans, bool* &true_ans);

#endif
