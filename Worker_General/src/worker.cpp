#include "worker.h"

worker::worker(double err){
	err_rate = err;
	return;
}

worker::worker(){
	err_rate = -1;
	return;
}

bool worker::answer(){
	int a = rand()%RAND_MAX;
	if(double(a) < err_rate*RAND_MAX){
		return true; //worker commits error
	}
	else{
		return false;//worker answers correctly
	}
}

bool worker::answer(bool true_answer){
	int a = rand()%RAND_MAX;
	if(double(a) < err_rate*RAND_MAX){
		return !true_answer; //worker commits error
	}
	else{
		return true_answer;//worker answers correctly
	}
}

double worker::get_rate(){
	return err_rate;
}

void worker::set_rate(double err){
	assert(err >= 0);
	err_rate = err;
}

bool worker::sanity_check(){
	bool chk1 = (0 <= err_rate) && (1 > err_rate);
	return chk1;
}

std::string worker::describe(){
	std::string s = "";
	s += "Error Rate:\t";
	std::ostringstream ss;
	ss << err_rate;
	s += ss.str(); 
	s += "\n";
	return s;
}

worker gen_worker(double err){
	worker w(err);
	return w;
}

worker gen_worker(worker w){
	worker w2(w.get_rate());
	return w2;
}


//k-ary Worker implementation begins here
Worker::Worker(){
	arity = -1;
	return;
}

Worker::Worker(int k){
	arity = k;
	err_rate.resize(k,k);
	for(int i=0;i<k;i++){
		for(int j=0;j<k;j++){
			err_rate(i,j)=-1;
		}
	}
	return;
}

Worker::Worker(MatrixXd err){
	arity = err.outerSize();
	err_rate = err;
	return;
}

MatrixXd Worker::get_rate(){
	return err_rate;
}

int Worker::get_arity(){
	return arity;
}

void Worker::set_rate(MatrixXd err){
	err_rate = err;
	return;
}

void Worker::set_arity(int k){
	arity = k;
	return;
}

int Worker::answer(int true_answer){
	double r = rand()%RAND_MAX;
	r /= RAND_MAX;
	int ans=0;
	while(r > err_rate(true_answer,ans) && ans < arity-1){
		r -= err_rate(true_answer,ans);
		ans++;
	}
	return ans;	
}

bool Worker::sanity_check(){
	bool chk = arity > 0;
	for(int i=0;i<arity;i++){
		double sum = 0;
		for(int j=0;j<arity;j++){
			chk &= err_rate(i,j)>0;
			sum += err_rate(i,j);
		}
		chk &= fabs(sum-1.0)<0.001;
	}
	return chk;
}

Worker gen_Worker(MatrixXd err){
	Worker W(err);
	return W;
}

Worker gen_Worker(Worker w){
	Worker W(w.get_rate());
	return W;
}

