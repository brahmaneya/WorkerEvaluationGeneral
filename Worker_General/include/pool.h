#ifndef POOL_H
#define POOL_H
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "worker.h"
#include<Eigen/Dense>
using Eigen::MatrixXd;

class pool {
public:
	virtual worker gen_worker()=0;		
};

class step_pool : public pool{
public: 
	step_pool(double errs[], double probs[], int size);
	worker gen_worker();
	bool sanity_check();//check that arrsize is ok, sum of probabilities is 1, etc. 
private: 
	double *err, *prob;
	int arrsize;
};


class Pool { //For k-ary workers
public:
	virtual Worker gen_Worker()=0;		
};

class Step_Pool : public Pool{
public: 
	Step_Pool(MatrixXd errs[], double probs[], int size);
	Worker gen_Worker();
	bool sanity_check();//check that arrsize is ok, sum of probabilities is 1, etc. 
private: 
	double *prob;
	MatrixXd * err;
	int arrsize;
};

#endif
