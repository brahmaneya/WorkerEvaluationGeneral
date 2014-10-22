#ifndef RTEDATA_H
#define RTEDATA_H

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

void get_rte_data(int& n, int& wn, bool** &attempted, bool** &ans, bool* &true_ans);

void get_rte_data_regular(int& n, int& wn, bool** &ans, bool* &true_ans);

#endif
