#include"pool.h"

step_pool::step_pool(double errs[], double probs[], int size){
	arrsize = size;
	err = new double[size];
	prob = new double[size];
	for(int i=0;i<arrsize;i++){
		err[i] = errs[i];
		prob[i] = probs[i];
	}
}

worker step_pool::gen_worker(){
	int a = rand()%RAND_MAX;
	double cprob = 0;
	double err_rate=-1;
	for(int i=0;i<arrsize;i++){
		cprob += prob[i];
		if(a <= cprob*RAND_MAX){
			err_rate = err[i];
			break;
		}	
	}
	worker w(err_rate);
	return w;
}

bool step_pool::sanity_check(){
	if(arrsize<=0){
		return false;
	}
	double cprob=0;
	for(int i=0;i<arrsize;i++){
		cprob += prob[i];
	}
	if (cprob < 1){
		return false;
	}
	return true;
}


Step_Pool::Step_Pool(MatrixXd errs[], double probs[], int size){
	arrsize = size;
	err = new MatrixXd[size];
	prob = new double[size];
	for(int i=0;i<arrsize;i++){
		err[i] = errs[i];
		prob[i] = probs[i];
	}
}

Worker Step_Pool::gen_Worker(){
	int a = rand()%RAND_MAX;
	double cprob = 0;
	MatrixXd err_rate;
	for(int i=0;i<arrsize;i++){
		cprob += prob[i];
		if(a <= cprob*RAND_MAX){
			err_rate = err[i];
			break;
		}	
	}
	Worker w(err_rate);
	return w;
}

bool Step_Pool::sanity_check(){
	if(arrsize<=0){
		return false;
	}
	double cprob=0;
	for(int i=0;i<arrsize;i++){
		cprob += prob[i];
	}
	if (cprob < 1){
		return false;
	}
	return true;
}
