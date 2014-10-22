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

#define f(i,n) for(int i=0;i<n;i++)

using namespace std;

int main(){
	srand(time(NULL));
	int n=100;
	int wn = 7;
	double probs[3] = {0.34,0.33,0.33};
	double errs[3] = {0.1,0.2,0.3};
	step_pool pl(errs,probs,3);
	worker W[wn];
	f(j,wn){
		W[j] = pl.gen_worker();
	} 
	double est[wn],epsilon[wn];
	bool ** attempted;
	bool ** ans;
	attempted = new bool*[wn];
	ans = new bool*[wn];
	bool correct[n];
	f(i,wn){
		attempted[i] = new bool[n];
		ans[i] = new bool[n];
	}
	
	f(i,wn){
		f(j,n){
			attempted[i][j] = true;
			ans[i][j] = W[i].answer(true);
		}
	}
	
	double score = 0;
	double uncertainty = 0;
	double uncertainty_conservative = 0;
	int numIter = 10000;
	double conf = 0.8;
	double density = 1;
	double densities[wn];
	double unc_gold;
	string s = static_cast<ostringstream*>( &(ostringstream() << wn) )->str() + "_" + static_cast<ostringstream*>( &(ostringstream() << n) )->str();
	//ofstream conf_acc(("plotdata/conf_acc_"+s).c_str());
	//ofstream dens_unc(("plotdata/dens_unc_"+s).c_str());
	for(conf=0.1;conf<=1.0;conf+=0.05){
		f(i,wn){
			densities[i] = (0.5*i + (wn-i))/wn;
		}
		score = 0;
		uncertainty = 0;
		uncertainty_conservative = 0;
		unc_gold = 0;
		f(iter,numIter){
			f(j,n){
				correct[j] = true;
			}
			f(i,wn){
				f(j,n){
					attempted[i][j] = (random()%100)<(densities[i]*100);
					ans[i][j] = W[i].answer(correct[j]);
				}
			}	
			N_Worker_nonreg(ans,attempted,wn,n,est,epsilon,conf);//the last (optional) argument determines whether the linear optimization is applied or not	
			f(i,wn){
				if(fabs(est[i]-W[i].get_rate()) < epsilon[i]){
					score+=1;
				}
				uncertainty += epsilon[i];
			}
			
			/* 
			//comparing to interval size of old technique. density needs to be 1 for this.
			n_worker_conservative(ans,wn,n,est,epsilon,conf);
			f(i,wn){
				uncertainty_conservative += epsilon[i];
			}
			*/
			
			f(i,wn){
				Gold_Standard(ans[i],correct,n,est[i],epsilon[i],conf);
				//unc_gold += epsilon[i];
			}
		}
		//conf_acc << conf << '\t' << score/(wn*numIter)<<endl;
		//cout << conf << '\t' << uncertainty/(wn*numIter)<< '\t' << uncertainty_conservative/(wn*numIter)<<endl;
		//dens_unc << density << '\t' << uncertainty/(wn*numIter) << endl;
		cout << conf << '\t' << score/(wn*numIter)<< '\t' << uncertainty/(wn*numIter) << endl;
	}
	
	
	/*
	for(int wwn=3;wwn <= wn;wwn+=40){
		//nonregbin_inter(wwn,n,attempted,ans,est,epsilon,0.9);
		double toteps1 = 0;
		double acc1 = 0;
		f(j,wwn){
			toteps1 += epsilon[j];
			if(fabs(W[j].get_rate()-est[j]) < epsilon[j] ){
				acc1+= 1.0;
			}
		}
		acc1 /= wwn;
		toteps1 /= wwn;
		
		n_worker_solution(wwn,n,ans,est,epsilon,0.9);
		double toteps2 = 0;
		double acc2 = 0;
		f(j,wwn){
			toteps2 += epsilon[j];
			if(fabs(W[j].get_rate()-est[j]) < epsilon[j] ){
				acc2+= 1.0;
			}
		}
		acc2 /= wwn;
		toteps2 /= wwn;
		cout << wwn << '\t' << toteps1 << " "  << acc1 << '\t' << toteps2 << " " << acc2 << endl;
	
	
	}
	*/
	return 0;
}
