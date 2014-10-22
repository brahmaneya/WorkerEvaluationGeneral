#include<iostream>
#include<stdlib.h>
#include"worker.h"
#include"pool.h"
#include<string>
#include<math.h>
#include"test.h"
#include"datasets.h"
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

/******
We are testing multiple models here: Each model must have a domain for p,d and a function that maps those domains into (1/2,1) (or (0,1)? ), and a model name (for optmization function)

Want to test a parametrized family of models later. But for now, its these four. WIll write one function for each, that takes error frequencies for groups of workers/tasks with equal p/d
and finds best fitting p's d's, and consequent best likelihood.  
******/

//For all functions below: pos is number of correct answers, neg is wrong answers.


double exp_prod(int ps, int ds, int ** pos, int ** neg, double p[], double d[]){//e^(-pd). p,d \in (0,inf) 
	double logprob; 
	f(i,ps){
		p[i] = 50.0/(1+rand()%100);
	}
	f(j,ds){
		d[j] = 50.0/(1+rand()%100);
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log(exp(-p[i]*d[j])) + neg[i][j]*log(1-exp(-p[i]*d[j]));
		}
	}
	
	double learn_rate = 0.00001; //arbitrary. needs to be adjusted based on dataset. this is for rte 11 categories each way
	double improv = 0 - logprob; 
	const double threshold = 0.1;
	while(improv > threshold){
		double pderivs[ps],dderivs[ds];
		f(i,ps){
			pderivs[i] = 0;
		}
		f(j,ds){
			dderivs[j] = 0;
		}
		f(i,ps){
			f(j,ds){
				pderivs[i] += pos[i][j]*(-d[j]) + neg[i][j]*(d[j]*exp(-p[i]*d[j])/(1-exp(-p[i]*d[j])));
				dderivs[j] += pos[i][j]*(-p[i]) + neg[i][j]*(p[i]*exp(-p[i]*d[j])/(1-exp(-p[i]*d[j])));
			}
		}
		
		f(i,ps){
			p[i] += learn_rate*pderivs[i];//cout << i << '\t' << p[i] << '\t' << pderivs[i] << endl; 
			if(p[i] < 0){
				p[i] = 0.00001;
			}
		}
		f(j,ds){
			d[j] += learn_rate*dderivs[j];//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
			if(d[j] < 0){
				d[j] = 0.00001;
			}
		}

		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log(exp(-p[i]*d[j])) + neg[i][j]*log(1-exp(-p[i]*d[j]));
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	
	}
	
	return logprob;
}

double prod(int ps, int ds, int ** pos, int ** neg, double p[], double d[]){ //pd, p,d \in (0,1) 
	double logprob; 
	f(i,ps){
		p[i] = (50+rand()%50)/100.0;
	}
	f(j,ds){
		d[j] = (50+rand()%50)/100.0;
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log(p[i]*d[j]) + neg[i][j]*log(1-p[i]*d[j]);
		}
	}
	
	double learn_rate = 0.00001; //arbitrary. needs to be adjusted based on dataset. this is for rte 11 categories each way
	double improv = 0 - logprob; 
	const double threshold = 0.1;
	while(improv > threshold){
		double pderivs[ps],dderivs[ds];
		f(i,ps){
			pderivs[i] = 0;
		}
		f(j,ds){
			dderivs[j] = 0;
		}
		f(i,ps){
			f(j,ds){
				pderivs[i] += pos[i][j]*(1/p[i]) + neg[i][j]*(-d[j]/(1-p[i]*d[j]));
				dderivs[j] += pos[i][j]*(1/d[j]) + neg[i][j]*(-p[i]/(1-p[i]*d[j]));
			}
		}
		
		f(i,ps){
			p[i] += learn_rate*pderivs[i];//cout << i << '\t' << p[i] << '\t' << pderivs[i] << endl; 
			if(p[i] < 0.01){
				p[i] = 0.01;
			}
			else if (p[i] > 0.999){
				p[i] = 0.999;
			}
		}
		f(j,ds){
			d[j] += learn_rate*dderivs[j];//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
			if(d[j] < 0.01){
				d[j] = 0.01;
			}
			else if (d[j] > 0.999){
				d[j] = 0.999;
			}
		}

		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log(p[i]*d[j]) + neg[i][j]*log(1-p[i]*d[j]);
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	}
	
	return logprob;
}


double pos_prod(int ps, int ds, int ** pos, int ** neg, double p[], double d[]){ //pd+(1-p)(1-d). p,d, \in 0.5,1
	double logprob; 
	f(i,ps){
		p[i] = (50+rand()%50)/99.0;
	}
	f(j,ds){
		d[j] = (50+rand()%50)/99.0;
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log(p[i]*d[j] + (1-p[i])*(1-d[j])) + neg[i][j]*log(p[i]*(1-d[j])+(1-p[i])*d[j]);
		}
	}
	
	double learn_rate = 0.00001; //arbitrary. needs to be adjusted based on dataset. this is for rte 11 categories each way
	double improv = 0 - logprob; 
	const double threshold = 0.1;
	while(improv > threshold){
		double pderivs[ps],dderivs[ds];
		f(i,ps){
			pderivs[i] = 0;
		}
		f(j,ds){
			dderivs[j] = 0;
		}
		f(i,ps){
			f(j,ds){
				pderivs[i] += pos[i][j]*((2*d[j]-1)/(p[i]*d[j] + (1-p[i])*(1-d[j]))) + neg[i][j]*((1-2*d[j])/(p[i]*(1-d[j])+(1-p[i])*d[j]));
				dderivs[j] += pos[i][j]*((2*p[i]-1)/(p[i]*d[j] + (1-p[i])*(1-d[j]))) + neg[i][j]*((1-2*p[i])/(p[i]*(1-d[j])+(1-p[i])*d[j]));
			}
		}
		
		f(i,ps){
			p[i] += learn_rate*pderivs[i];//cout << i << '\t' << p[i] << '\t' << pderivs[i] << endl; 
			if(p[i] < 0.51){
				p[i] = 0.51;
			}
			if(p[i] > 0.99){
				p[i] = 0.99;
			}
		}
		f(j,ds){
			d[j] += learn_rate*dderivs[j];//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
			if(d[j] < 0.51){
				d[j] = 0.51;
			}
			if(d[j] > 0.99){
				d[j] = 0.99;
			}
		}

		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log(p[i]*d[j] + (1-p[i])*(1-d[j])) + neg[i][j]*log(p[i]*(1-d[j])+(1-p[i])*d[j]);
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	
	}
	
	return logprob;
}


double sigmoid(int ps, int ds, int ** pos, int ** neg, double p[], double d[]){ //1/(1+e^(-pd)). p \in (0,inf) d \in (-inf,inf)
	double logprob; 
	f(i,ps){
		p[i] = 50.0/(1+rand()%100);
	}
	f(j,ds){
		d[j] = 50.0/(1+rand()%100);
		if(rand()%2 == 0){
			d[j] = -d[j];
		}
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log(1/(1+exp(-p[i]*d[j]))) + neg[i][j]*log(exp(-p[i]*d[j])/(1+exp(-p[i]*d[j])));
		}
	}
	
	double learn_rate = 0.00001; //arbitrary. needs to be adjusted based on dataset. this is for rte 11 categories each way
	double improv = 0 - logprob; 
	const double threshold = 0.1;
	while(improv > threshold){
		double pderivs[ps],dderivs[ds];
		f(i,ps){
			pderivs[i] = 0;
		}
		f(j,ds){
			dderivs[j] = 0;
		}
		f(i,ps){
			f(j,ds){
				pderivs[i] += pos[i][j]*(d[j]*exp(-p[i]*d[j])/(1+exp(-p[i]*d[j]))) + neg[i][j]*(-d[j]/(1+exp(-p[i]*d[j])));
				dderivs[j] += pos[i][j]*(p[i]*exp(-p[i]*d[j])/(1+exp(-p[i]*d[j]))) + neg[i][j]*(-p[i]/(1+exp(-p[i]*d[j])));
			}
		}
		
		f(i,ps){
			p[i] += learn_rate*pderivs[i];//cout << i << '\t' << p[i] << '\t' << pderivs[i] << endl; 
			if(p[i] < 0){
				p[i] = 0.00001;
			}
		}
		f(j,ds){
			d[j] += learn_rate*dderivs[j];//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
		}

		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log(1/(1+exp(-p[i]*d[j]))) + neg[i][j]*log(exp(-p[i]*d[j])/(1+exp(-p[i]*d[j])));
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	
	}
	
	return logprob;
}


double pos_sigmoid(int ps, int ds, int ** pos, int ** neg, double p[], double d[]){ //1/(1+e^(-pd)). p,d \in (0,inf)
	double logprob; 
	f(i,ps){
		p[i] = 50.0/(1+rand()%100);
	}
	f(j,ds){
		d[j] = 50.0/(1+rand()%100);
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log(1/(1+exp(-p[i]*d[j]))) + neg[i][j]*log(exp(-p[i]*d[j])/(1+exp(-p[i]*d[j])));
		}
	}
	
	double learn_rate = 0.00001; //arbitrary. needs to be adjusted based on dataset. this is for rte 11 categories each way
	double improv = 0 - logprob; 
	const double threshold = 0.1;
	while(improv > threshold){
		double pderivs[ps],dderivs[ds];
		f(i,ps){
			pderivs[i] = 0;
		}
		f(j,ds){
			dderivs[j] = 0;
		}
		f(i,ps){
			f(j,ds){
				pderivs[i] += pos[i][j]*(d[j]*exp(-p[i]*d[j])/(1+exp(-p[i]*d[j]))) + neg[i][j]*(-d[j]/(1+exp(-p[i]*d[j])));
				dderivs[j] += pos[i][j]*(p[i]*exp(-p[i]*d[j])/(1+exp(-p[i]*d[j]))) + neg[i][j]*(-p[i]/(1+exp(-p[i]*d[j])));
			}
		}
		
		f(i,ps){
			p[i] += learn_rate*pderivs[i];//cout << i << '\t' << p[i] << '\t' << pderivs[i] << endl; 
			if(p[i] < 0.001){
				p[i] = 0.001;
			}
		}
		f(j,ds){
			d[j] += learn_rate*dderivs[j];//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
			if(d[j] < 0.001){
				d[j] = 0.001;
			}
		}

		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log(1/(1+exp(-p[i]*d[j]))) + neg[i][j]*log(exp(-p[i]*d[j])/(1+exp(-p[i]*d[j])));
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	
	}
	
	return logprob;
}

double pos_exp_prod(int ps, int ds, int ** pos, int ** neg, double p[], double d[]){//0.5+0.5*e^(-pd). p,d \in (0,inf) 
	double logprob; 
	f(i,ps){
		p[i] = 50.0/(1+rand()%100);
	}
	f(j,ds){
		d[j] = 50.0/(1+rand()%100);
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log(0.5+0.5*exp(-p[i]*d[j])) + neg[i][j]*log(0.5-0.5*exp(-p[i]*d[j]));
		}
	}
	
	double learn_rate = 0.00001; //arbitrary. needs to be adjusted based on dataset. this is for rte 11 categories each way
	double improv = 0 - logprob; 
	const double threshold = 0.1;
	while(improv > threshold){
		double pderivs[ps],dderivs[ds];
		f(i,ps){
			pderivs[i] = 0;
		}
		f(j,ds){
			dderivs[j] = 0;
		}
		f(i,ps){
			f(j,ds){
				pderivs[i] += pos[i][j]*(-0.5*d[j]*exp(-p[i]*d[j])/(0.5+0.5*exp(-p[i]*d[j]))) + neg[i][j]*(0.5*d[j]*exp(-p[i]*d[j])/(0.5-0.5*exp(-p[i]*d[j])));
				dderivs[j] += pos[i][j]*(-0.5*p[i]*exp(-p[i]*d[j])/(0.5+0.5*exp(-p[i]*d[j]))) + neg[i][j]*(0.5*p[i]*exp(-p[i]*d[j])/(0.5-0.5*exp(-p[i]*d[j])));
			}
		}
		
		f(i,ps){
			p[i] += learn_rate*pderivs[i];//cout << i << '\t' << p[i] << '\t' << pderivs[i] << endl; 
			if(p[i] < 0){
				p[i] = 0.00001;
			}
		}
		f(j,ds){
			d[j] += learn_rate*dderivs[j];//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
			if(d[j] < 0){
				d[j] = 0.00001;
			}
		}

		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log(0.5+0.5*exp(-p[i]*d[j])) + neg[i][j]*log(0.5-0.5*exp(-p[i]*d[j]));
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	
	}
	
	return logprob;
}

double exp(double p, double d){
	return exp(-log(p)*log(d));
}

double pos_exp(double p, double d){
	return 0.5+0.5*exp(-log(p)*log(d));
}

double sigmoid(double p, double d){
	return (d>0.5 ? 1:-1) * 1/(1+exp(-log(p)*log(fabs(2*d-1))));
}

double pos_sigmoid(double p, double d){
	return 1/(1+exp(-log(p)*log(d)));
}

double prod(double p, double d){
	return p*d;
}

double pos_prod(double p, double d){
	return 0.5+0.5*p*d;
}

double Optimize(int ps, int ds, int** pos, int** neg, double p[], double d[], double (*func)(double, double) ){//assumes that p,d have to be in (0,1), and function domain is (0,1) or(0.5,1)
	double logprob; 
	f(i,ps){
		p[i] = (1+rand()%100)/101.0;
	}
	f(j,ds){
		d[j] = (1+rand()%100)/101.0;
	}
	
	logprob = 0;
	f(i,ps){
		f(j,ds){
			logprob += pos[i][j]*log((*func)(p[i],d[j])) + neg[i][j]*log(1-(*func)(p[i],d[j]));
		}
	}
	
	double improv = 0 - logprob;
	double threshold = 0.1;
	
	int numIters = 10;
	
	f(iterNo,numIters){
	
		//optimize ds 
		f(j,ds){
			double logprob_local = 0;
			f(i,ps){
				logprob_local += pos[i][j]*log((*func)(p[i],d[j])) + neg[i][j]*log(1-(*func)(p[i],d[j]));
			}cout << j << '\t' << logprob_local << endl;
			double learn_rate = 0.0001; //arbitrary. needs to be adjusted based on dataset. 
			double improv_local = 0 - logprob_local; 
			const double threshold_local = 0.1;
			while(improv_local > threshold_local){cout << j << '\t' << logprob_local << endl;
				double dderiv = 0;
				f(i,ps){
					double lderiv = ((*func)(p[i],d[j]+learn_rate) - (*func)(p[i],d[j]-learn_rate))/(2*learn_rate);
					dderiv += pos[i][j]*lderiv/(*func)(p[i],d[j]) - neg[i][j]*lderiv/(1-(*func)(p[i],d[j]));
				}
				d[j] += learn_rate*dderiv;//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
				if(d[j] < 0.001){
					d[j] = 0.001;
				}
				if(d[j] > 0.999){
					d[j] = 0.999;
				}
				double newlogprob_local = 0;
				f(i,ps){
					newlogprob_local += pos[i][j]*log((*func)(p[i],d[j])) + neg[i][j]*log(1-(*func)(p[i],d[j]));
				}
				improv_local = newlogprob_local - logprob_local;
				logprob_local = newlogprob_local;
			}
		}
		
		//optimize ps
		f(i,ps){
			double logprob_local = 0;
			f(j,ds){
				logprob_local += pos[i][j]*log((*func)(p[i],d[j])) + neg[i][j]*log(1-(*func)(p[i],d[j]));
			}
			double learn_rate = 0.0001; //arbitrary. needs to be adjusted based on dataset. 
			double improv_local = 0 - logprob_local; 
			const double threshold_local = 0.1;
			while(improv_local > threshold_local){
				double pderiv = 0;
				f(j,ds){
					double lderiv = ((*func)(p[i]+learn_rate,d[j]) - (*func)(p[i]-learn_rate,d[j]))/(2*learn_rate);
					pderiv += pos[i][j]*lderiv/(*func)(p[i],d[j]) - neg[i][j]*lderiv/(1-(*func)(p[i],d[j]));
				}
				p[i] += learn_rate*pderiv;//cout << j << '\t' << d[j] << '\t' << dderivs[j] << endl;
				if(p[i] < 0.001){
					p[i] = 0.001;
				}
				if(p[i] > 0.999){
					p[i] = 0.999;
				}
				double newlogprob_local = 0;
				f(j,ds){
					newlogprob_local += pos[i][j]*log((*func)(p[i],d[j])) + neg[i][j]*log(1-(*func)(p[i],d[j]));
				}
				improv_local = newlogprob_local - logprob_local;
				logprob_local = newlogprob_local;
			}
		}
		
	
	
		double newlogprob = 0;
		f(i,ps){
			f(j,ds){
				newlogprob += pos[i][j]*log((*func)(p[i],d[j])) + neg[i][j]*log(1-(*func)(p[i],d[j]));
			}
		}
		improv = newlogprob - logprob;
		logprob = newlogprob;
	}
	
	return logprob;
}


double best_score(int ps, int ds, int** pos, int** neg){ //best logprob possible for any model
	double logprob = 0;
	f(i,ps){
		f(j,ds){
			if(pos[i][j] != 0 && neg[i][j] != 0){
				double p = ((double)pos[i][j])/(pos[i][j]+neg[i][j]);
				logprob += pos[i][j]*log(p) + neg[i][j]*log(1-p);
			}
		}
	}
	return logprob;
}

int main(){
	int n,wn;
	bool* true_ans;
	bool** ans;
	bool** attempted;
	get_imana_data(n, wn, attempted, ans, true_ans);
	
	///May delete parts after this, depending on whether i want to optimize on raw data, or on the aggregates which i compute below. 
	
	f(i,n){
		int corr = 0;
		f(j,wn){
			if(ans[j][i] == true_ans[i]){
				corr++;
			}
		}	
		//cout << i << '\t' << q1[i] << ',' << q2[i] << '\t' << ((double)corr)/wn << endl;
		double acc = ((double)corr)/wn;
		if(acc < 0.5){
			//cout << i << '\t' << q1[i] << ',' << q2[i] << '\t' << ((double)corr)/wn << endl;	
		}
	}


	srand(time(0));
	int wtots[wn];
	double wacc[wn];
	int ttots[n];
	double tacc[n];
	
	f(i,wn){
		wtots[i] = 0;
		wacc[i] = 0;
	}
	f(j,n){
		ttots[j] = 0;
		tacc[j] = 0;
	}
	
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				wtots[i]++;
				ttots[j]++;
				if(ans[i][j] == true_ans[j]){
					wacc[i] += 1;
					tacc[j] += 1.0;
				}
			}
		}
	}
	
	int numCats = 11;
	double t_hist[numCats];
	double w_hist[numCats];
	
	double histacc[numCats][numCats];
	int histtot[numCats][numCats];
	
	f(Windex,numCats){
		w_hist[Windex] = 0;
	}
	f(Tindex,numCats){
		t_hist[Tindex] = 0;	
	}
	
	f(Windex,numCats){
		f(Tindex,numCats){
			histtot[Windex][Tindex] = 0;
			histacc[Windex][Tindex] = 0;
		}
	}
	
	f(i,wn){
		wacc[i] /= wtots[i];
		w_hist[(int)(wacc[i]*(numCats-1))] += 1;
	}
	f(j,n){
		tacc[j] /= ttots[j];
		t_hist[(int)(tacc[j]*(numCats-1))] += 1;
	}
	
	f(i,wn){
		//cout << wacc[i] << '\t' << wtots[i] << endl;
	}
	f(j,n){
		//cout << tacc[j] << '\t' << ttots[j] << endl;
	}
	
	f(i,numCats){
		//cout << i << '\t' << w_hist[i] << '\t' << t_hist[i] << endl;
	}
	//cout << endl;
	
	
	f(i,wn){
		f(j,n){
			if(attempted[i][j]){
				int Windex = (int)(wacc[i]*(numCats-1));
				int Tindex = (int)(tacc[j]*(numCats-1));
				histtot[Windex][Tindex]++;
				if(ans[i][j] == true_ans[j]){
					histacc[Windex][Tindex] += 1;
				}
			}
		}
	}
	
	/*
	for(int i=numCats/2; i < numCats; i++){
		for(int j=numCats/2; j < numCats; j++){
			cout << histacc[i][j] << "," << histtot[i][j] << '\t';
		}
		cout << endl;
	}
	*/
	
	int ps = numCats, ds = numCats;
	double p[ps],d[ds];
	int ** pos = new int*[ps];
	int ** neg = new int*[ps];
	f(i,ps){
		pos[i] = new int[ds];
		neg[i] = new int[ds];
	}
	
	f(i,ps){
		f(j,ds){
			pos[i][j] = histacc[i][j];
			neg[i][j] = histtot[i][j] - histacc[i][j];
		}
	}
	
	double bestp[ps],bestd[ds];
	double maxscore = sigmoid(ps,ds,pos,neg,p,d);
	/*
	f(i,ps){
		bestp[i] = p[i];
	}
	f(j,ds){
		bestd[j] = d[j];
	}
	f(i,1){
		double score = exp_prod(ps,ds,pos,neg,p,d);
		if(maxscore < score){
			maxscore = score;
			f(i,ps){
				bestp[i] = p[i];
			}
			f(j,ds){
				bestd[j] = d[j];
			}
		}
	}
	
	f(i,ps){
		cout << bestp[i] << " ";
	}
	cout << endl;
	f(j,ds){
		cout << bestd[j] << " ";
	}
	cout << endl;
	cout << maxscore << endl;
	*/
	double optscore = Optimize(ps,ds,pos,neg,p,d,exp);
	
	f(i,ps){
		cout << bestp[i] << " ";
	}
	cout << endl;
	f(j,ds){
		cout << bestd[j] << " ";
	}
	cout << endl;
	cout << optscore << endl;
	
	//cout << best_score(ps,ds,pos,neg) << endl;
	
	return 0;
}
