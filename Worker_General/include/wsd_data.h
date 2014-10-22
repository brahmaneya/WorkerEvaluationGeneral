#ifndef WSDDATA_H
#define WSDDATA_H

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

void get_wsd_data(int &arity, int& n, int& wn, bool** &attempted, int** &ans, int* &true_ans);

#endif
