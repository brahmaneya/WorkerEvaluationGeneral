#ifndef WORKER_H
#define WORKER_H
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<string>
#include<sstream>
#include<cassert>
#include<Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::JacobiSVD;

class worker {
public:
	bool sanity_check(); 	//sanity check on worker. Check that error rate is between 0 and 1
	worker(double err);	
	worker();
	double get_rate();
	void set_rate(double err);
	bool answer();	//generate 1 with probability err_rate and 0 otherwise
	bool answer(bool true_answer);	//give answer not matching true_answer with probability err_rate, and matching true_answer with probablity 1-err_rate 
	friend worker gen_worker(double err);	//create worker with error rate err
	friend worker gen_worker(worker w);	//create duplicate of w
	std::string describe();	//output worker error rate
private:
	double err_rate;
};

worker gen_worker(double err);
worker gen_worker(worker w);

class Worker{
public:
	bool sanity_check(); 	//sanity check on worker. Check that error rate is between 0 and 1
	Worker(MatrixXd err);	
	Worker(int k);
	Worker();
	MatrixXd get_rate();
	int get_arity();
	void set_rate(MatrixXd err);
	void set_arity(int k);
	int answer(int true_answer);	//give answer given true answer 
	friend Worker gen_Worker(MatrixXd err);	//create worker with error rate err
	friend Worker gen_Worker(Worker w);	//create duplicate of w
	
private:
	int arity;
	MatrixXd err_rate;
};

Worker gen_Worker(MatrixXd err);
Worker gen_Worker(Worker w);


#endif
