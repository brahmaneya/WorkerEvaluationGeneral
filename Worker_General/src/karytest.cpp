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
#include<Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::JacobiSVD;

#define f(i,n) for(int i=0;i<n;i++)

using namespace std;

int main(){
	srand(time(0));
	int arity = 4;
	int poolSize = 3;
	int n = 60;
	int wn = 3;	
	double probs[poolSize];
	probs[0] = 1.0/3; probs[1] = 1.0/3; probs[2] = 1.0/3;
	MatrixXd errs[poolSize];
	f(i,poolSize){
		errs[i].resize(arity,arity);
	}
	
	if(arity == 3){
	errs[0] << 0.6,0.3,0.1,0.1,0.6,0.3,0.3,0.1,0.6;
	errs[1] << 0.8,0.1,0.1,0.2,0.8,0.0,0.0,0.2,0.8;
	errs[2] << 0.9,0.0,0.1,0.1,0.9,0.0,0.0,0.2,0.8;
	}
	else if(arity == 2){	
	errs[0] << 0.9,0.1,0.2,0.8;
	errs[1] << 0.8,0.2,0.1,0.9;
	errs[2] << 0.9,0.1,0.1,0.9;
	}
	else if(arity == 4){
	errs[0] << 0.7,0.1,0.1,0.1,0.1,0.6,0.2,0.1,0.0,0.1,0.8,0.1,0.2,0.1,0.0,0.7;
	errs[1] << 0.8,0.1,0.0,0.1,0.1,0.8,0.0,0.2,0.1,0.1,0.7,0.1,0.0,0.1,0.2,0.7;
	errs[2] << 0.6,0.1,0.2,0.1,0.0,0.7,0.1,0.2,0.1,0.0,0.9,0.0,0.2,0.0,0.0,0.8;
	}
	
	Step_Pool PL(errs,probs,poolSize);
	Worker W[wn];
	f(j,wn){
		W[j] = PL.gen_Worker();
	}
	MatrixXd s(arity,1);
	f(i,arity){
		s(i,0) = 1.0/arity;
	}
	int true_ans[n];
	f(i,n){
		int r = rand()%RAND_MAX;
		int ans = 0;
		while(r > s(ans,0)*RAND_MAX && ans < arity-1){
			r -= s(ans,0)*RAND_MAX;
			ans++;
		}
		true_ans[i] = ans;
	}
	
	bool ** attempted = new bool*[wn];
	int ** answers = new int*[wn];
	f(i,wn){
		answers[i] = new int[n];
		attempted[i] = new bool[n];
	}
	
	f(i,wn){
		f(j,n){
			answers[i][j] = W[i].answer(true_ans[j]);
			attempted[i][j] = rand()%10<8;
		}
	}
	MatrixXd est[wn],epsilon[wn];
	double conf = 0.8;
	double density = 1.0;
	int numIters = 1000/(arity*arity);
	double tot;
	
	MatrixXd P1,P2,P3;
	
	for(conf = 0.1; conf < 1; conf += 0.1){
		tot = 0;
		double uncertainty = 0.0;
		double corrs = 0;
		double corr[wn][arity][arity];
		double total[wn][arity][arity];
		f(j,wn){
			f(i1,arity){
				f(i2,arity){
					corr[j][i1][i2] = 0;
					total[j][i1][i2] = 0;
				}
			}
		}
		
		f(iter,numIters){
			f(j,wn){
				W[j] = PL.gen_Worker();
			}
		
			f(i,n){
				int r = rand()%RAND_MAX;
				int ans = 0;
				while(r > s(ans,0)*RAND_MAX && ans < arity-1){
					r -= s(ans,0)*RAND_MAX;
					ans++;
				}
				true_ans[i] = ans;
			}	
	
			f(i,wn){
				f(j,n){
					if(rand()%100 < density*100){
						answers[i][j] = W[i].answer(true_ans[j]);
						attempted[i][j] = true;
					}
					else{
						attempted[i][j] = false;
					}
				}
			}
	
			Three_Worker_nonreg_Kary(arity, attempted, answers, n, est, epsilon, conf);
			f(j,wn){
				MatrixXd temp = s.cwiseSqrt().asDiagonal() * W[j].get_rate();
				f(i1,arity){
					f(i2,arity){
						if(est[j]((i1,i2) >= 0 || est[j](i1,i2) <= 0) && fabs(est[j](i1,i2)+epsilon[j](i1,i2)) < 4){
							if(epsilon[j](i1,i2)>1000){
								//cout << est[j](i1,i2) << '\t' << epsilon[j](i1,i2) << '\t' << temp(i1,i2) << endl;
							}
							total[j][i1][i2] += 1;
							if(fabs(est[j](i1,i2) - temp(i1,i2)) < epsilon[j](i1,i2)){
								corr[j][i1][i2] += 1;
								corrs += 1;
							}
							else{
								//cout << epsilon[j](i1,i2) << endl;
							}
							uncertainty += epsilon[j](i1,i2);
						}
						else{
							//cout << "error"<<endl;
						}
						tot += 1.0;
					}
				}
				
			}
			
				
		}
		
		cout << conf << '\t' << corrs/tot << '\t' << uncertainty/tot << endl;
	}
	
	
	/*
	f(w,3){
		f(i,arity){
			double tempsum = 0;
			f(j,arity){
				tempsum += est[w](i,j);
			}
			//cout << tempsum*sqr(i,i) << endl;
			f(j,arity){
				est[w](i,j) /= tempsum;
			}
		}
		//cout << est[w] << endl << endl;
	}
	*/
}
