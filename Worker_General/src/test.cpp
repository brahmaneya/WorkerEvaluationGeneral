#include "test.h"

#define f(i,n) for(int i=0;i<n;i++)

using std::cout;
using std::endl;
using std::set;

void generate_questions(bool questions[], int n, double s){
	int r;
	for(int i=0;i<n;i++){
		r = rand()%RAND_MAX;
		if(r < s*RAND_MAX){
			questions[i] = true;
		}
		else{
			questions[i] = false;
		}
	}
}

void generate_answers(worker w[], int wn, int  n, bool questions[], bool** answers){
	for(int i=0;i<n;i++){
		for(int j=0;j<wn;j++){
			answers[i][j] = w[j].answer(questions[i]);	
		}
	}
}

bool weighted_majority(worker w[], int wn, bool ans[], double s){
	double score = log(s/(1-s));
	f(j,wn){
		if(ans[j]){
			score += log((1-w[j].get_rate())/(w[j].get_rate()));
		}
		else{
			score -= log((1-w[j].get_rate())/(w[j].get_rate()));
		}
	}
	return (score > 0);
	
}

bool majority(int wn, bool ans[]){
	double score = 0;
	f(j,wn){
		if(ans[j]){
			score += 1;
		}
		else{
			score -= 1;
		}
	}
	return (score > 0);
}


double majority_error(double error_rates[], int n){
	double poly[n+1];
	for(int i=0;i<n+1;i++){
		poly[i]=0.0;
	}
	poly[0] = 1;
	for(int i=0;i<n;i++){
		for(int j=n;j>0;j--){
			poly[j] = error_rates[i]*poly[j] + (1-error_rates[i])*poly[j-1];
		}
		poly[0] = poly[0] * error_rates[i];
	}
	double err=0.0;
	for(int i=0;2*i<n;i++){
		err += poly[i];
	}
	return err;
}

double majority_error_derivative(double error_rates[], int n, int k){
	double poly[n+1];
	for(int i=0;i<n+1;i++){
		poly[i]=0.0;
	}
	poly[0] = 1;
	for(int i=0;i<n;i++){
		if(i==k){
			continue;
		}
		for(int j=n;j>0;j--){
			poly[j] = error_rates[i]*poly[j] + (1-error_rates[i])*poly[j-1];
		}
		poly[0] = poly[0] * error_rates[i];
	}
	return poly[n/2];
}



void expectation_maximisation(int wn, int n, bool ** answers, double est[]){
	double ans[n];
	double threshold = 0.00001; //will keep iterating till worker ability starts changing slowing than this per worker
	double diff = 1; //will record difference between iterations, and stop when this falls below threshold
	double s = 0.55; //estimated selectivity
	//because of symmetry, there are polar local optima. so if selectivity is on the wrong side of 0.5, their algo converges to 1-error_rates
	f(i,n){
		ans[i] = s; //initialize answers as per selectivity
	}
	f(j,wn){
		est[j]=1; //arbitrary initialization for error rates
	}
	double err;
	double tr,fs; 
	double sel;
	while(diff > threshold){
		diff=0;
		f(j,wn){
			err=0;
			f(i,n){
				if(answers[i][j]){
					err += 1-ans[i];
				}
				else{
					err += ans[i];
				}
			}
			err = err/n;
			diff += fabs(est[j]-err);
			est[j] = err;
		}
		diff = diff/wn;
		//now the other step
		sel = 0;
		f(i,n){
			tr = s;
			fs = 1 - s;
			f(j,wn){
				if(answers[i][j]){
					tr *= 1-est[j];
					fs *= est[j]; 
				}
				else{
					fs *= 1-est[j];
					tr *= est[j];
				}
			}
			ans[i] = tr/(tr+fs);
			sel += ans[i];
		}
		s = sel/n;	
	}
	
}

void n_worker_conservative(bool** answers, int wn, int n, double est[], double epsilon[], double conf){ // only in this n_worker, the order of variables in answers is correct (worker, then task)
	int r;
	double z = erfinv(conf); //because for every worker, we will estimate 3 variables and use them to get worker error rate 
	
	int sums[n][wn];
	for(int i=0;i<n;i++){
		sums[i][0] = (int)answers[0][i]; 
		for(int j=1;j<wn;j++){
			sums[i][j] = sums[i][j-1] + (int)answers[j][i];
		}
	}
	
	double corr[wn][3];
	for(int j=0;j<wn;j++){
		bool b[3];
		int majc1,majc2;
		double q12=0, q23=0, q13=0;
		double est12, est13, est23;
		double epsilon12, epsilon13, epsilon23;
		double upp,low;
		for(int i=0;i<n;i++){
			b[0] = answers[j][i];
			if(j < (wn-1)/2){ 
				majc1 = sums[i][(wn-1)/2] - (int)answers[j][i];
				b[1] = (majc1 >= (wn+1)/4);
				majc2 = sums[i][wn-1]-sums[i][(wn-1)/2];
				b[2] = (majc2 > (wn-1-(wn-1)/2)/2);
			}
			else{
				majc1 = sums[i][(wn-3)/2];
				b[1] = (majc1 >= (wn+1)/4);
				majc2 = sums[i][wn-1]-sums[i][(wn-3)/2]-(int)answers[j][i];
				b[2] = (majc2 > (wn-1-(wn-1)/2)/2);
			}
			q12 += (int)(b[0]==b[1]);
			q13 += (int)(b[0]==b[2]);
			q23 += (int)(b[1]==b[2]);
		}	
		est23 = (q23 + z*z/2)/(n + z*z);
		epsilon23 = (z*sqrt(q23*(n-q23)/(n) + (z*z)/(4) ))/(n + z*z);
		est13 = (q13 + z*z/2)/(n + z*z);
		epsilon13 = (z*sqrt(q13*(n-q13)/(n) + (z*z)/(4) ))/(n + z*z);
		est12 = (q12 + z*z/2)/(n + z*z);
		epsilon12 = (z*sqrt(q12*(n-q12)/(n) + (z*z)/(4) ))/(n + z*z);
	
		double l12 = est12-epsilon12-0.5;
		if(l12 < 0.01){
			l12 = 0.01;
		}
		double l23 = est23-epsilon23-0.5;
		if(l23 < 0.01){
			l23 = 0.01;
		}
		double l13 = est13-epsilon13-0.5;
		if(l13 < 0.01){
			l13 = 0.01;
		}
		double h12 = est12+epsilon12-0.5;
		if(h12 < 0.01){
			h12 = 0.01;
		}
		double h23 = est23+epsilon23-0.5;
		if(h23 < 0.01){
			h23 = 0.01;
		}
		double h13 = est13+epsilon13-0.5;
		if(h13 < 0.01){
			h13 = 0.01;
		}
		
		upp = 0.5 - sqrt(l12*l13/(2*h23)); //danger, of negative square root because of broad limits
		low = 0.5 - sqrt(h12*h13/(2*l23)); 
		est[j] = (upp+low)/2;
		epsilon[j] = (upp-low)/2;
	}
}

void Gold_Standard(bool ans[], bool correct[], int n, double &est, double &epsilon, double conf){ 
	double corr = 0; //number of correct answers
	int r;
	bool true_ans;
	for (int i=0;i<n;i++){
		true_ans = correct[i];
		if(ans[i] == true_ans){
			corr += 1;
		}
	}
	double z = erfinv(conf); 
	est = 1 - (corr + z*z/2)/(n + z*z);
	epsilon = (z*sqrt(corr*(n-corr)/(n) + (z*z)/(4) ))/(n + z*z);
}

void Three_Worker(bool** answers, int n, double est[], double epsilon[], double conf){ 
	double z = erfinv(conf); 
	double q12=0, q23=0, q13=0;
	
	for(int i=0;i<n;i++){
		if(answers[0][i]==answers[1][i]){
			q12 += 1;
		}
		if(answers[0][i]==answers[2][i]){
			q13 += 1;
		}
		if(answers[1][i]==answers[2][i]){
			q23 += 1;
		}
	}
	
	
	double est23 = (q23 + z*z/2)/(n + z*z);
	double epsilon23 = (z*sqrt(q23*(n-q23)/(n) + (z*z)/(4) ))/(n + z*z);
	double devq23 = epsilon23/z;
	double est13 = (q13 + z*z/2)/(n + z*z);
	double epsilon13 = (z*sqrt(q13*(n-q13)/(n) + (z*z)/(4) ))/(n + z*z);
	double devq13 = epsilon13/z;
	double est12 = (q12 + z*z/2)/(n + z*z);
	double epsilon12 = (z*sqrt(q12*(n-q12)/(n) + (z*z)/(4) ))/(n + z*z);
	double devq12 = epsilon12/z;
	
	double meanp1 = 0.5 - sqrt((est12-0.5)*(est13-0.5)/(2*(est23-0.5)));
	double meanp2 = 0.5 - sqrt((est12-0.5)*(est23-0.5)/(2*(est13-0.5)));
	double meanp3 = 0.5 - sqrt((est23-0.5)*(est13-0.5)/(2*(est12-0.5)));
	
	double corr11 = (meanp1*(1-meanp1)*(2*est23-1))/n;
	double corr22 = (meanp2*(1-meanp2)*(2*est13-1))/n;
	double corr33 = (meanp3*(1-meanp3)*(2*est12-1))/n;
	
	double deriv112 = sqrt((est13-0.5)/(8*(est12-0.5)*(est23-0.5)));//deriv of p1 formula wrt q12
	double deriv113 = sqrt((est12-0.5)/(8*(est13-0.5)*(est23-0.5)));//deriv of p1 formula wrt q13
	double deriv123 = -sqrt((est12-0.5)*(est13-0.5)/(8*(est23-0.5)*(est23-0.5)*(est23-0.5)));//deriv of p1 formula wrt q23
	double deriv212 = sqrt((est23-0.5)/(8*(est12-0.5)*(est13-0.5)));//deriv of p2 formula wrt q12
	double deriv213 = -sqrt((est12-0.5)*(est23-0.5)/(8*(est13-0.5)*(est13-0.5)*(est13-0.5)));//deriv of p2 formula wrt q13
	double deriv223 = sqrt((est12-0.5)/(8*(est13-0.5)*(est23-0.5)));//deriv of p2 formula wrt q23
	double deriv312 = -sqrt((est23-0.5)*(est13-0.5)/(8*(est12-0.5)*(est12-0.5)*(est12-0.5)));//deriv of p3 formula wrt q12
	double deriv313 = sqrt((est23-0.5)/(8*(est12-0.5)*(est13-0.5)));//deriv of p3 formula wrt q13
	double deriv323 = sqrt((est13-0.5)/(8*(est12-0.5)*(est23-0.5)));//deriv of p3 formula wrt q23
	
	double devp1 = sqrt(deriv112*deriv112*devq12*devq12 + deriv113*deriv113*devq13*devq13 + deriv123*deriv123*devq23*devq23 + deriv112*deriv113*corr11 + deriv112*deriv123*corr22 + deriv113*deriv123*corr33);
	double devp2 = sqrt(deriv212*deriv212*devq12*devq12 + deriv213*deriv213*devq13*devq13 + deriv223*deriv223*devq23*devq23 + deriv212*deriv213*corr11 + deriv212*deriv223*corr22 + deriv213*deriv223*corr33);
	double devp3 = sqrt(deriv312*deriv312*devq12*devq12 + deriv313*deriv313*devq13*devq13 + deriv323*deriv323*devq23*devq23 + deriv312*deriv313*corr11 + deriv312*deriv323*corr22 + deriv313*deriv323*corr33);
	
	est[0] = meanp1;
	epsilon[0] = z*devp1;
	est[1] = meanp2;
	epsilon[1] = z*devp2;
	est[2] = meanp3;
	epsilon[2] = z*devp3;
	return;	
	
}

void Three_Worker_nonreg(bool** answers, bool** attempted, int n, double est[], double epsilon[], double conf){ 
	double z = erfinv(conf);
	double q12=0, q23=0, q13=0;
	int count12=0,count13=0,count23=0,count123=0;
	
	for(int i=0;i<n;i++){
		if(attempted[0][i] && attempted[1][i]){
			count12++;
			if(answers[0][i]==answers[1][i]){
				q12 += 1;
			}
		}
		if(attempted[0][i] && attempted[2][i]){
			count13++;
			if(answers[0][i]==answers[2][i]){
				q13 += 1;
			}
		}
		if(attempted[1][i] && attempted[2][i]){
			count23++;
			if(answers[1][i]==answers[2][i]){
				q23 += 1;
			}
		}
		
		if(attempted[0][i] && attempted[1][i] && attempted[2][i]){
			count123++;
		}
	}
	
	
	double est23 = q23/count23;//(q23 + z*z/2)/(count23 + z*z);
	double epsilon23 = (z*sqrt(q23*(count23-q23)/(count23) + (z*z)/(4) ))/(count23 + z*z);
	double devq23 = epsilon23/z;
	double est13 = q13/count13;//(q13 + z*z/2)/(count13 + z*z);
	double epsilon13 = (z*sqrt(q13*(count13-q13)/(count13) + (z*z)/(4) ))/(count13 + z*z);
	double devq13 = epsilon13/z;
	double est12 = q12/count12;//(q12 + z*z/2)/(count12 + z*z);
	double epsilon12 = (z*sqrt(q12*(count12-q12)/(count12) + (z*z)/(4) ))/(count12 + z*z);
	double devq12 = epsilon12/z;
	
	if(est12 < 0.51){
		est12 = 0.51;
	}
	if(est23 < 0.51){
		est23 = 0.51;
	}
	if(est13 < 0.51){
		est13 = 0.51;
	}
	
	double meanp1 = 0.5 - sqrt((est12-0.5)*(est13-0.5)/(2*(est23-0.5)));
	double meanp2 = 0.5 - sqrt((est12-0.5)*(est23-0.5)/(2*(est13-0.5)));
	double meanp3 = 0.5 - sqrt((est23-0.5)*(est13-0.5)/(2*(est12-0.5)));
	
	double corr11 = count123*(meanp1*(1-meanp1)*(2*est23-1))/(count12*count13);
	double corr22 = count123*(meanp2*(1-meanp2)*(2*est13-1))/(count12*count23);
	double corr33 = count123*(meanp3*(1-meanp3)*(2*est12-1))/(count13*count23);
	
	double deriv112 = sqrt((est13-0.5)/(8*(est12-0.5)*(est23-0.5)));//deriv of p1 formula wrt q12
	double deriv113 = sqrt((est12-0.5)/(8*(est13-0.5)*(est23-0.5)));//deriv of p1 formula wrt q13
	double deriv123 = -sqrt((est12-0.5)*(est13-0.5)/(8*(est23-0.5)*(est23-0.5)*(est23-0.5)));//deriv of p1 formula wrt q23
	double deriv212 = sqrt((est23-0.5)/(8*(est12-0.5)*(est13-0.5)));//deriv of p2 formula wrt q12
	double deriv213 = -sqrt((est12-0.5)*(est23-0.5)/(8*(est13-0.5)*(est13-0.5)*(est13-0.5)));//deriv of p2 formula wrt q13
	double deriv223 = sqrt((est12-0.5)/(8*(est13-0.5)*(est23-0.5)));//deriv of p2 formula wrt q23
	double deriv312 = -sqrt((est23-0.5)*(est13-0.5)/(8*(est12-0.5)*(est12-0.5)*(est12-0.5)));//deriv of p3 formula wrt q12
	double deriv313 = sqrt((est23-0.5)/(8*(est12-0.5)*(est13-0.5)));//deriv of p3 formula wrt q13
	double deriv323 = sqrt((est13-0.5)/(8*(est12-0.5)*(est23-0.5)));//deriv of p3 formula wrt q23
	
	double devp1 = sqrt(deriv112*deriv112*devq12*devq12 + deriv113*deriv113*devq13*devq13 + deriv123*deriv123*devq23*devq23 + deriv112*deriv113*corr11 + deriv112*deriv123*corr22 + deriv113*deriv123*corr33);
	double devp2 = sqrt(deriv212*deriv212*devq12*devq12 + deriv213*deriv213*devq13*devq13 + deriv223*deriv223*devq23*devq23 + deriv212*deriv213*corr11 + deriv212*deriv223*corr22 + deriv213*deriv223*corr33);
	double devp3 = sqrt(deriv312*deriv312*devq12*devq12 + deriv313*deriv313*devq13*devq13 + deriv323*deriv323*devq23*devq23 + deriv312*deriv313*corr11 + deriv312*deriv323*corr22 + deriv313*deriv323*corr33);
	
	est[0] = meanp1;
	epsilon[0] = z*devp1;
	est[1] = meanp2;
	epsilon[1] = z*devp2;
	est[2] = meanp3;
	epsilon[2] = z*devp3;
	
	//if(est[0] < -0.05){cout << est12 << '\t' << est13 << '\t' << est23 << endl;}
	
	return;		
}


void N_Worker(bool ** answers, int wn, int n,  double est[], double epsilon[], double conf){
	double z = erfinv(conf); 
	
	bool ** ans = new bool*[3];
	f(i,3){
		ans[i] = new bool[n];
	}
	double tempest[3],tempeps[3];
	
	int sums[wn][n];
	for(int i=0;i<n;i++){
		sums[0][i] = (int)answers[0][i]; 
		for(int j=1;j<wn;j++){
			sums[j][i] = sums[j-1][i] + (int)answers[j][i];
		}
	}
	
	for(int j=0;j<wn;j++){
		bool b[3];
		int majc1,majc2;
		for(int i=0;i<n;i++){
			b[0] = answers[j][i];
			if(j < (wn-1)/2){ //need to correct indices
				majc1 = sums[(wn-1)/2][i] - (int)answers[j][i];
				b[1] = (majc1 >= (wn+1)/4);
				majc2 = sums[wn-1][i]-sums[(wn-1)/2][i];
				b[2] = (majc2 > (wn-1-(wn-1)/2)/2);
			}
			else{
				majc1 = sums[(wn-3)/2][i];
				b[1] = (majc1 >= (wn+1)/4);
				majc2 = sums[wn-1][i]-sums[(wn-3)/2][i]-(int)answers[j][i];
				b[2] = (majc2 > (wn-1-(wn-1)/2)/2);
			}
			f(t,3){
				ans[t][i] = b[t];
			}
		}	
		Three_Worker(ans,n,tempest,tempeps,conf);
		est[j] = tempest[0];
		epsilon[j] = tempeps[0];		
	}
}

void N_Worker_nonreg(bool ** answers, bool** attempted, int wn, int n,  double est[], double epsilon[], double conf, bool noopt){
	double z = erfinv(conf); 
	
	int count[wn][wn];
	double q[wn][wn];
	double estq[wn][wn];
	
	f(j1,wn){
		f(j2,wn){
			count[j1][j2] = 0;
			q[j1][j2] = 0.0;
		}
	}
	
	f(i,n){
		f(j1,wn){
			if(!attempted[j1][i]){
				continue;
			}
			f(j2,wn){
				if(attempted[j2][i]){
					count[j1][j2]++;
					if(answers[j1][i] == answers[j2][i]){
						q[j1][j2] += 1;
					}
				}
			}
		}
	}
	
	f(j1,wn){
		f(j2,wn){
			estq[j1][j2] = q[j1][j2]/count[j1][j2];//(q[j1][j2]+z*z/2)/(count[j1][j2]+z*z);
			if(estq[j1][j2] < 0.51){
				estq[j1][j2] = 0.51;
			}
		}
	}
	
	double tempest[3],tempeps[3];
	bool **tempattempted = new bool*[3];
	bool **tempans = new bool*[3];
	
	double tempmean;
	
	f(jt,wn){
		double countt[wn][wn];
		f(j1,wn){
			f(j2,wn){
				countt[j1][j2] = 0;
			}
		}
		f(i,n){
			if(!attempted[jt][i]){
				continue;
			}
			f(j1,wn){
				if(!attempted[j1][i]){
					continue;
				}
				f(j2,wn){
					if(attempted[j2][i]){
						countt[j1][j2]++;
					}
				}
			}
		}
		
		set<int> Remaining;
		f(i,wn){
			if(i != jt && count[jt][i] >= 5){ //arbitrary lower bound (5 instead of 1 because Gaussian assumption works only for big nos). (< 5) common tasks => data not enough
				Remaining.insert(i);
			}
		}
		
		
		int maxPairs = (wn-1)/2;
		int numPairs = 0;
		int WorkerFirst[maxPairs];
		int WorkerSecond[maxPairs];
		
		double Covar[maxPairs][maxPairs];
		double means[maxPairs];
		
		tempmean=0;
		f(pairNo,maxPairs){
			
			bool done = false;
			while(!done){
				if(Remaining.size() < 2){
					break;
				}
				int nextWorker = rand()%Remaining.size();
				set<int>::const_iterator it(Remaining.begin());
	  			advance(it,nextWorker);
	  			nextWorker = *it;
	  			Remaining.erase(nextWorker);
	  			
	  			for(set<int>::iterator rit = Remaining.begin(); rit != Remaining.end(); rit++){
	  				int otherWorker = *rit;
	  				if(count[nextWorker][otherWorker] < 5){
	  					continue;
	  				}
	  				else{
	  					done = true;
	  					WorkerFirst[pairNo] = nextWorker;
	  					WorkerSecond[pairNo] = otherWorker;
	  					Remaining.erase(otherWorker);
	  					break;
	  				}
	  			}
	  			
			}
			if(!done){
				break;
			}
			
			/*
			bool done = false;
			while(!done){
				int workerNo = rand()%wn;
				if(taken[workerNo] ==  false){
					done = true;
					taken[workerNo] = true;
					WorkerFirst[pairNo] = workerNo;
				}
			}
			done = false;
			while(!done){
				int workerNo = rand()%wn;
				if(taken[workerNo] ==  false){
					done = true;
					taken[workerNo] = true;
					WorkerSecond[pairNo] = workerNo;
				}
			}
			*/
			
			
			tempattempted[0] = attempted[jt];
			tempans[0] = answers[jt];
			tempattempted[1] = attempted[WorkerFirst[pairNo]];
			tempans[1] = answers[WorkerFirst[pairNo]];
			tempattempted[2] = attempted[WorkerSecond[pairNo]];
			tempans[2] = answers[WorkerSecond[pairNo]];
			
			Three_Worker_nonreg(tempans,tempattempted,n,tempest,tempeps,conf);
			Covar[pairNo][pairNo] = (tempeps[0]/z)*(tempeps[0]/z);
			means[pairNo] = tempest[0];
			tempmean += means[pairNo];
			//cout << means[pairNo]<<"  ";
			numPairs = pairNo+1;
		}
		tempmean /= numPairs;//cout<<jt << '\t' << numPairs << endl;
		//cout << endl<<tempmean << endl << endl;
		f(pair1,numPairs){
			f(pair2,numPairs){
				if(pair1==pair2){
					continue;
				}
				double deriv11 = sqrt((estq[1][WorkerFirst[pair1]]-0.5)/(8*(estq[1][WorkerSecond[pair1]]-0.5)*(estq[WorkerFirst[pair1]][WorkerSecond[pair1]]-0.5)));
				double deriv12 = sqrt((estq[1][WorkerSecond[pair1]]-0.5)/(8*(estq[1][WorkerFirst[pair1]]-0.5)*(estq[WorkerFirst[pair1]][WorkerSecond[pair1]]-0.5)));
				double deriv21 = sqrt((estq[1][WorkerFirst[pair2]]-0.5)/(8*(estq[1][WorkerSecond[pair2]]-0.5)*(estq[WorkerFirst[pair2]][WorkerSecond[pair2]]-0.5)));
				double deriv22 = sqrt((estq[1][WorkerSecond[pair2]]-0.5)/(8*(estq[1][WorkerFirst[pair2]]-0.5)*(estq[WorkerFirst[pair2]][WorkerSecond[pair2]]-0.5)));
				double correl = 0;
				correl += deriv11*deriv21*(tempmean*(1-tempmean)*countt[WorkerFirst[pair1]][WorkerFirst[pair2]]*(2*estq[WorkerFirst[pair1]][WorkerFirst[pair2]] -1)/(count[jt][WorkerFirst[pair1]]*count[jt][WorkerFirst[pair2]])); 
				correl += deriv12*deriv21*(tempmean*(1-tempmean)*countt[WorkerSecond[pair1]][WorkerFirst[pair2]]*(2*estq[WorkerSecond[pair1]][WorkerFirst[pair2]] -1)/(count[jt][WorkerSecond[pair1]]*count[jt][WorkerFirst[pair2]])); 
				correl += deriv11*deriv22*(tempmean*(1-tempmean)*countt[WorkerFirst[pair1]][WorkerSecond[pair2]]*(2*estq[WorkerFirst[pair1]][WorkerSecond[pair2]] -1)/(count[jt][WorkerFirst[pair1]]*count[jt][WorkerSecond[pair2]])); 
				correl += deriv12*deriv22*(tempmean*(1-tempmean)*countt[WorkerSecond[pair1]][WorkerSecond[pair2]]*(2*estq[WorkerSecond[pair1]][WorkerSecond[pair2]] -1)/(count[jt][WorkerSecond[pair1]]*count[jt][WorkerSecond[pair2]])); 
				
				Covar[pair1][pair2] = correl;
			}
		}
		//cout << jt << '\t' << numPairs << endl;
		
		double a[numPairs];
		f(i,numPairs){
			a[i] = 1.0/numPairs;
		}
		
		if(!noopt){
			MatrixXd Covars(numPairs,numPairs);
			MatrixXd as(numPairs,1);
			f(i1,numPairs){
				as(i1,0) = 1.0;
				f(i2,numPairs){
					Covars(i1,i2) = Covar[i1][i2];
				}
			}
			as = Covars.inverse()*as;
			f(i,numPairs){
				a[i] = as(i,0);
			}
		}
		
		double mean = 0;
		double asum = 0.0;
		f(i,numPairs){
			mean += a[i]*means[i];
			asum += a[i];
		}
		mean /= asum;
		double dev = 0;
		f(i1,numPairs){
			f(i2,numPairs){
				dev += a[i1]*a[i2]*Covar[i1][i2]/(asum*asum);
			}
		} 
		if(dev<0.0001){
			dev = 0.0001;
		}
		dev = sqrt(dev);
		est[jt] = mean;
		epsilon[jt] = z*dev;
		//hack
		if(est[jt] < epsilon[jt]){
			est[jt] = epsilon[jt];
		}
	}
	
}


double erfinv(double conf, double threshold){
	double x=0, y=8; //erf is sufficiently close to 1 even for y >=5
	double z=(x+y)/2;
	while(fabs(erf(z)-conf) >= threshold){
		if (erf(z) > conf){
			y = z;
			z = (x+y)/2;
		}
		else{
			x = z;
			z = (x+y)/2;
		}
	}
	return z*sqrt(2);
}


void Gold_Standard_nonreg_Kary(int arity, bool attempted[], int ans[], int correct[], int n, MatrixXd &est, MatrixXd &epsilon, double conf){
	double qcounts[arity];
	double z = erfinv(conf);
	MatrixXd acounts(arity,arity);
	f(i,arity){
		qcounts[i] = 0;
		f(j,arity){
			acounts(i,j)=0;
		}
	}
	f(i,n){
		if(attempted[i]){
			qcounts[correct[i]] += 1;
			acounts(correct[i],ans[i]) += 1;
		}
	}
	est.resize(arity,arity);
	epsilon.resize(arity,arity);
	f(i,arity){
		f(j,arity){
			est(i,j) = (acounts(i,j) + z*z/2)/(qcounts[i]+z*z);
			epsilon(i,j) = (z*sqrt(acounts(i,j)*(qcounts[i]-acounts(i,j))/(qcounts[i]) + (z*z)/(4) ))/(qcounts[i] + z*z);
		}
	}
	return;
}

MatrixXd MatSqrt(MatrixXd input){
	JacobiSVD<MatrixXd> J = input.jacobiSvd(1<<2 | 1<<4);
	MatrixXd U = J.matrixU();
	MatrixXd V = J.matrixV();
	MatrixXd Sing = J.singularValues();
	Sing = Sing.cwiseSqrt();
	MatrixXd sol = U*Sing.asDiagonal()*V.transpose();
	return sol; 
}	
	
void Mat3Solve(int arity, MatrixXd est[], double *** est_counts, int tot[]){
	MatrixXd Q[3];
	f(i,3){
		Q[i] = MatrixXd::Zero(arity,arity);
		est[i].resize(arity,arity);
	}	
	f(i,arity){
		f(j,arity){
			f(k,arity){
				Q[0](i,j) += tot[3]*est_counts[i][j][k];
				Q[1](j,k) += tot[3]*est_counts[i][j][k];
				Q[2](k,i) += tot[3]*est_counts[i][j][k];
			}
			
			Q[0](i,j) += tot[0]*est_counts[i][j][arity];
			Q[1](i,j) += tot[1]*est_counts[arity][i][j];
			Q[2](i,j) += tot[2]*est_counts[j][arity][i];
		}
	}
	
	Q[0] /= tot[0]+tot[3];
	Q[1] /= tot[1]+tot[3];
	Q[2] /= tot[2]+tot[3];
	
	est[0] = MatSqrt(Q[2].transpose() * Q[1].inverse() * Q[0].transpose());
	est[1] = est[0].transpose().inverse() * Q[0];
	est[2] = est[0].transpose().inverse() * Q[2].transpose();
	
	MatrixXd P = MatrixXd::Zero(arity,arity);
/*
	MatrixXd Pold;
	MatrixXd temp2old;
	MatrixXd temp1old;
	MatrixXd triplesold;
	MatrixXd Uold;
	double old;
*/
	f(i,1){//can do f(i,arity). had some nan problems there
		MatrixXd Triples(arity,arity);
		f(j,arity){
			f(k,arity){
				Triples(j,k) = est_counts[j][k][i];
			}
		}
		const MatrixXd temp1 = est[0].transpose().inverse() * Triples * est[1].inverse();
		JacobiSVD<MatrixXd> J = temp1.jacobiSvd(1<<2 | 1<<4);
		MatrixXd U = J.matrixU();
		MatrixXd temp2 = U.inverse() * est[0];
		if(!(temp2(0,0)>0 || temp2(0,0) <= 0)){
			//cout << "U, est\n" << U << "\n\n" << est[0] << endl;
		}
		
		if(!(est[0](0,0)>0 || est[0](0,0) <= 0)){
			//cout  << "Qs\n" << Q[1] << "\n\n" << Q[1].inverse() <<endl<< endl;
			f(i1,arity){
				f(i2,arity){
			//		cout << i1 << ", " << i2 << '\t' << tot[1]*est_counts[arity][i1][i2];
					f(k,arity){
						//cout << " + " << tot[3]*est_counts[k][i1][i2];
					}
			//		cout << endl;
				}
			}
		}
		
		f(j,arity){
			double jsum = 0;
			f(k,arity){
				jsum += temp2(j,k);
			}
			if(jsum < 0){
				f(k,arity){
					temp2(j,k) = -temp2(j,k);
				}
			}
		}
		bool max_done[arity];
		f(j,arity){
			max_done[j] = false;
		}
		
		f(j,arity){
			int max_index=-1;
			double max=temp2(j,0)-100;
			
			f(k,arity){
				if(max_done[k] == false){
					max_index = k;
					max = temp2(j,k);
					break;	
				}
			}
			
			f(k,arity){
				if(max_done[k]){
					continue;
				}
				if(temp2(j,k) > max){
					max_index = k;
					max = temp2(j,k);
				}
			}
			
			max_done[max_index] = true;
			f(k,arity){
				P(max_index,k) += temp2(j,k)/1;//if we have f(i,arity), then divide by arity instead of 1 here
			}
		}
		/*
		MatrixXd temp = P.transpose().inverse();
		MatrixXd tem = (temp2/arity).transpose().inverse();
		MatrixXd temold = (temp2old/arity).transpose().inverse();
		
		if(i==1 && !(temp(0,0)>=0 || temp(0,0) < 0)){
			MatrixXd PTemp = temp2+temp2old;
			MatrixXd P2Temp = U.inverse() + Uold.inverse();
			cout << temp2 << "\n\n" << U.inverse()*est[0] << "\n\n" << temp2old << "\n\n" << Uold.inverse()*est[0] << endl << endl; 
			cout << PTemp << "\n\n" << P2Temp*est[0] << "\n\n" << P2Temp.transpose().inverse() << endl << endl;
			exit(0);
		}
		Pold = P;
		temp1old = temp1;
		triplesold = Triples;
		temp2old = temp2;
		Uold = U;
		*/
	}
	est[0] = P;//if(!(P(0,0) > 0 || P(0,0) <= 0)){cout << Q[0] << "\n" << Q[1] << "\n" << Q[2] << endl;exit(0);}
	est[1] = est[0].transpose().inverse() * Q[0];
	est[2] = est[0].transpose().inverse() * Q[2].transpose();
}	

double Covariance_for_3WnK(int arity, int i1, int i2, int i3, int j1, int j2, int j3, double *** est_counts, int tot[], double z = 0.0){ //only handles case where 2 or more have answered each set. else returns 0
	if((i1 == arity && i2 == arity) || (i2 == arity && i3 == arity) || (i3 == arity && i1 == arity)){
		return 0;
	}
	if((j1 == arity && j2 == arity) || (j2 == arity && j3 == arity) || (j3 == arity && j1 == arity)){
		return 0;
	}
	if((i1-arity)*(j1-arity) == 0 and i1 != j1){
		return 0;
	}
	if((i2-arity)*(j2-arity) == 0 and i2 != j2){
		return 0;
	}
	if((i3-arity)*(j3-arity) == 0 and i3 != j3){
		return 0;
	}
	int total = 0;
	if(i3 == arity){
		total = tot[0];
	}
	else if(i1 == arity){
		total = tot[1];
	}
	else if(i2 == arity){
		total = tot[2];
	}
	else{
		total = tot[3];
	}
	if(total == 0){
		return 0;
	}
	else if(i1==j1 && i2==j2 && i3==j3){//changing this to make it similar to wilson score, rather than direct normal approx.
		//return est_counts[i1][i2][i3]*(1-est_counts[i1][i2][i3])/total;
		return ( (est_counts[i1][i2][i3]*(1-est_counts[i1][i2][i3])/total) + z*z/(4*total*total) ) / ((1+z*z/total)*(1+z*z/total));
	}				
	else{
		return -est_counts[i1][i2][i3]*est_counts[j1][j2][j3]/total;
	}
}
	

void Three_Worker_nonreg_Kary(int arity, bool ** attempted, int ** answers, int n, MatrixXd est[], MatrixXd epsilon[], double conf){
	double z = erfinv(conf);
	int tot[4];
	f(i,4){//0,1,2 are for 12,23,31. 3 is for 123
		tot[i] = 0;
	}
	double dev_counts[arity+1][arity+1][arity+1];// last index value (=arity) is if task not attempted by that worker 
	double *** est_counts = new double**[arity+1]; //we're not bothering with questions attempted by less than 2 workers. 
	f(i,arity+1){
		est_counts[i] = new double*[arity+1];
		f(j,arity+1){
			est_counts[i][j] = new double[arity+1];
			f(k,arity+1){
				est_counts[i][j][k] = 0;
				dev_counts[i][j][k] = 0;
			}
		}
	}
	
	f(i,n){
		if(attempted[0][i] && attempted[1][i] && !attempted[2][i]){
			tot[0]++;
			est_counts[answers[0][i]][answers[1][i]][arity] += 1;
		}
		if(!attempted[0][i] && attempted[1][i] && attempted[2][i]){
			tot[1]++;
			est_counts[arity][answers[1][i]][answers[2][i]] += 1;
		}
		if(attempted[0][i] && !attempted[1][i] && attempted[2][i]){
			tot[2]++;
			est_counts[answers[0][i]][arity][answers[2][i]] += 1;
		}
		if(attempted[0][i] && attempted[1][i] && attempted[2][i]){
			tot[3]++;
			est_counts[answers[0][i]][answers[1][i]][answers[2][i]] += 1;
		}
	}
	
	double derivatives[3][arity][arity][arity+1][arity+1][arity+1]; //derivatives[i1][i2][i3][j1][j2][j3] is deriv of est[i1](i2,i3) wrt. term est_counts[j1](j2,j3)
	f(i1,3){
		f(i2,arity){
			f(i3,arity){
				f(j1,arity+1){
					f(j2,arity+1){
						f(j3,arity+1){
							derivatives[i1][i2][i3][j1][j2][j3] = 0.0;
						}
					}
				}
			}
		}
	}
	
	f(i,arity){
		f(j,arity){//don't have Wilson Score approximation here, so won't work well for edge values, and small tot values. 
			   //The other assignments to dev do use the wilson score interval half-size. the half size in wilson score is bigger for border values of p, as should occur here. 
			if(tot[3] != 0){
				f(k,arity){
					//dev_counts[i][j][k] = sqrt(est_counts[i][j][k]*(tot[3]-est_counts[i][j][k])/(tot[3]))/tot[3];
					dev_counts[i][j][k] = sqrt(est_counts[i][j][k]*(tot[3]-est_counts[i][j][k])/(tot[3]*tot[3]*tot[3]) + z*z/(4*tot[3]*tot[3]))/(1+z*z/tot[3]);
					est_counts[i][j][k] = est_counts[i][j][k]/tot[3];
				}
			}
			if(tot[0] != 0){
				//dev_counts[i][j][arity] = sqrt(est_counts[i][j][arity]*(tot[0]-est_counts[i][j][arity])/(tot[0]))/tot[0];
				dev_counts[i][j][arity] = sqrt(est_counts[i][j][arity]*(tot[0]-est_counts[i][j][arity])/(tot[0]*tot[0]*tot[0]) + z*z/(4*tot[0]*tot[0]))/(1+z*z/tot[0]);
				est_counts[i][j][arity] = est_counts[i][j][arity]/tot[0];
			}
			if(tot[1] != 0){
				//dev_counts[arity][i][j] = sqrt(est_counts[arity][i][j]*(tot[1]-est_counts[arity][i][j])/(tot[1]))/tot[1];
				dev_counts[arity][i][j] = sqrt(est_counts[arity][i][j]*(tot[1]-est_counts[arity][i][j])/(tot[1]*tot[1]*tot[1]) + z*z/(4*tot[1]*tot[1]))/(1+z*z/tot[1]);
				est_counts[arity][i][j] = est_counts[arity][i][j]/tot[1];
			}
			if(tot[2] != 0){
				//dev_counts[j][arity][i] = sqrt(est_counts[j][arity][i]*(tot[2]-est_counts[j][arity][i])/(tot[2]))/tot[2];
				dev_counts[j][arity][i] = sqrt(est_counts[j][arity][i]*(tot[2]-est_counts[j][arity][i])/(tot[2]*tot[2]*tot[2]) + z*z/(4*tot[2]*tot[2]))/(1+z*z/tot[2]);
				est_counts[j][arity][i] = est_counts[j][arity][i]/tot[2];
			}
			
		}
	}
	
	Mat3Solve(arity, est, est_counts, tot);
	
	f(d1,arity+1){
		f(d2,arity+1){
			f(d3,arity+1){
				if((d1 == arity && d2 == arity) || (d2 == arity && d3 == arity) || (d3 == arity && d1 == arity)){
					f(i1,3){
						f(i2,arity){
							f(i3,arity){
								derivatives[i1][i2][i3][d1][d2][d3] = 0;
							}
						}
					}
					continue;
				}
				int total = 0;
				if(d3 == arity){
					total = tot[0];
				}
				else if(d1 == arity){
					total = tot[1];
				}
				else if(d2 == arity){
					total = tot[2];
				}
				else{
					total = tot[3];
				}
				if(total == 0){
					f(i1,3){
						f(i2,arity){
							f(i3,arity){
								derivatives[i1][i2][i3][d1][d2][d3] = 0;
							}
						}
					}
					continue;
				}
				
				MatrixXd temp_upper[3],temp_lower[3],temp_orig[3];
				
				Mat3Solve(arity, temp_orig, est_counts, tot);
				
				double delta = dev_counts[d1][d2][d3];
				if(delta == 0){
					delta = 0.001;
				}
				
				est_counts[d1][d2][d3] += delta;
				Mat3Solve(arity, temp_upper, est_counts, tot);
				
				est_counts[d1][d2][d3] -= 2*delta;
				Mat3Solve(arity, temp_lower, est_counts, tot);
				
				est_counts[d1][d2][d3] += delta;
				f(i1,3){
					f(i2,arity){
						f(i3,arity){
							derivatives[i1][i2][i3][d1][d2][d3] = (temp_upper[i1](i2,i3) - temp_lower[i1](i2,i3))/(2*delta);
							//if(!(derivatives[i1][i2][i3][d1][d2][d3] > 0 || derivatives[i1][i2][i3][d1][d2][d3] <= 0)){cout << est[i1] << endl << temp_upper[i1] << endl << temp_lower[i1] << endl << temp_upper[i1]-temp_lower[i1]<<'\t'<<est_counts[d1][d2][d3]<<'\t'<<d1<<","<<d2<<","<<d3<<'\t' << delta << endl;exit(0);}
						}
					}
				}
			}
		}
	}
	
	f(i1,3){
		epsilon[i1].resize(arity,arity);
		f(i2,arity){
			f(i3,arity){
				int total;
				epsilon[i1](i2,i3) = 0.0;
				double A,B;
				total = tot[3];
				f(j1,arity){
					f(j2,arity){
						f(j3,arity){
							if(total==0){
								j1=j2=j3=arity;
								continue;
							}
							A += derivatives[i1][i2][i3][j1][j2][j3] * derivatives[i1][i2][i3][j1][j2][j3] * ( (est_counts[j1][j2][j3]*(1-est_counts[j1][j2][j3])/total) + z*z/(4*tot[3]*total) ) / ((1+z*z/total)*(1+z*z/total));
							A += derivatives[i1][i2][i3][j1][j2][j3] * derivatives[i1][i2][i3][j1][j2][j3] * est_counts[j1][j2][j3]/sqrt(total) * est_counts[j1][j2][j3]/sqrt(total);
							B += derivatives[i1][i2][i3][j1][j2][j3] * est_counts[j1][j2][j3]/sqrt(total);
						}
					}
				}
				if(A - B*B < 0){
					epsilon[i1](i2,i3) += 0;
				}
				else{
					epsilon[i1](i2,i3) += A - B*B;
				}
				
				A=B=0;
				total = tot[0];
				f(j1,arity){
					f(j2,arity){
						if(total==0){
							j1=j2=arity;
							continue;
						}
						A += derivatives[i1][i2][i3][j1][j2][arity] * derivatives[i1][i2][i3][j1][j2][arity] * ( (est_counts[j1][j2][arity]*(1-est_counts[j1][j2][arity])/total) + z*z/(4*tot[3]*total) ) / ((1+z*z/total)*(1+z*z/total));
						A += derivatives[i1][i2][i3][j1][j2][arity] * derivatives[i1][i2][i3][j1][j2][arity] * est_counts[j1][j2][arity]/sqrt(total) * est_counts[j1][j2][arity]/sqrt(total);
						B += derivatives[i1][i2][i3][j1][j2][arity] * est_counts[j1][j2][arity]/sqrt(total);

					}
				}
				if(A - B*B < 0){
					epsilon[i1](i2,i3) += 0;
				}
				else{
					epsilon[i1](i2,i3) += A - B*B;
				}
								
				A=B=0;
				total = tot[1];
				f(j2,arity){
					f(j3,arity){
						if(total==0){
							j2=j3=arity;
							continue;
						}
						A += derivatives[i1][i2][i3][arity][j2][j3] * derivatives[i1][i2][i3][arity][j2][j3] * ( (est_counts[arity][j2][j3]*(1-est_counts[arity][j2][j3])/total) + z*z/(4*tot[3]*total) ) / ((1+z*z/total)*(1+z*z/total));
						A += derivatives[i1][i2][i3][arity][j2][j3] * derivatives[i1][i2][i3][arity][j2][j3] * est_counts[arity][j2][j3]/sqrt(total) * est_counts[arity][j2][j3]/sqrt(total);
						B += derivatives[i1][i2][i3][arity][j2][j3] * est_counts[arity][j2][j3]/sqrt(total);
					}
				}
				if(A - B*B < 0){
					epsilon[i1](i2,i3) += 0;
				}
				else{
					epsilon[i1](i2,i3) += A - B*B;
				}
								
				A=B=0;
				total = tot[2];
				f(j1,arity){
					f(j3,arity){
						if(total==0){
							j1=j3=arity;
							continue;
						}
						A += derivatives[i1][i2][i3][j1][arity][j3] * derivatives[i1][i2][i3][j1][arity][j3] * ( (est_counts[j1][arity][j3]*(1-est_counts[j1][arity][j3])/total) + z*z/(4*tot[3]*total) ) / ((1+z*z/total)*(1+z*z/total));
						A += derivatives[i1][i2][i3][j1][arity][j3] * derivatives[i1][i2][i3][j1][arity][j3] * est_counts[j1][arity][j3]/sqrt(total) * est_counts[j1][arity][j3]/sqrt(total);
						B += derivatives[i1][i2][i3][j1][arity][j3] * est_counts[j1][arity][j3]/sqrt(total);
					}
				}
				if(A - B*B < 0){
					epsilon[i1](i2,i3) += 0;
				}
				else{
					epsilon[i1](i2,i3) += A - B*B;
				}
				
				epsilon[i1](i2,i3) = z*sqrt(epsilon[i1](i2,i3));
				if(!(est[i1](i2,i3) > 0 || est[i1](i2,i3) <= 0)){
					//cout << i1 << '\t' << i2 << '\t' << i3 << endl;
				}
			}
		}
	}
	
	/*
	//Change of plan, to make this more efficient. Going to express the whole sum as A - B*B. 
	f(i1,3){
		epsilon[i1].resize(arity,arity);
		f(i2,arity){
			f(i3,arity){
				epsilon[i1](i2,i3) = 0.0;
				f(j1,arity+1){
					f(j2,arity+1){
						f(j3,arity+1){
							f(k1,arity+1){
								f(k2,arity+1){
									f(k3,arity+1){
										double covar = Covariance_for_3WnK(arity, j1, j2, j3, k1, k2, k3, est_counts, tot, z);
										epsilon[i1](i2,i3) += derivatives[i1][i2][i3][j1][j2][j3] * derivatives[i1][i2][i3][k1][k2][k3] * covar;
									}
								}
							}
						}
					}
				}
				epsilon[i1](i2,i3) = z*sqrt(epsilon[i1](i2,i3));
			}
		}
	}
	*/

	return;
}


void NEWWWW_Three_Worker_nonreg_Kary(int arity, bool ** attempted, int ** answers, int n, MatrixXd est[], MatrixXd epsilon[], double conf){
	int tot123=0;
	int tot[3];
	MatrixXd counts[3],est_counts[3],dev_counts[3]; //0 is 12, 1 is 23, 2 is 31
	
	double rands[arity];
	f(i,arity){
		rands[i] = ((double) rand() / (RAND_MAX)) - 0.5;
	}
	MatrixXd special(arity,arity);
	MatrixXd Triples[arity];
	f(i,arity){
		Triples[i].resize(arity,arity);
	}
	
	f(i,3){
		tot[i] = 0;
		counts[i].resize(arity,arity);
		est_counts[i].resize(arity,arity);
		dev_counts[i].resize(arity,arity);
	}
	double est_counts123[arity][arity][arity], dev_counts123[arity][arity][arity];
	double z = erfinv(conf);
	
	f(i,arity){
		f(j,arity){
			special(i,j) = 0;
			f(k,3){
				counts[k](i,j) = 0;
			}
			f(k,arity){
				est_counts123[i][j][k] = 0;
			}
		}
	}
	
	f(i,n){
		if(attempted[0][i] && attempted[1][i]){
			tot[0]++;
			counts[0](answers[0][i],answers[1][i]) += 1;
		}
		if(attempted[1][i] && attempted[2][i]){
			tot[1]++;
			counts[1](answers[1][i],answers[2][i]) += 1;
		}
		if(attempted[2][i] && attempted[0][i]){
			tot[2]++;
			counts[2](answers[2][i],answers[0][i]) += 1;
		}
		if(attempted[0][i] && attempted[1][i] && attempted[2][i]){
			tot123++;
			est_counts123[answers[0][i]][answers[1][i]][answers[2][i]] += 1;
			
			special(answers[0][i],answers[1][i]) += rands[answers[2][i]];
			Triples[answers[2][i]](answers[0][i],answers[1][i]) += 1;
		}
	}
	
	double derivatives[3][arity][arity][3][arity][arity]; //derivatives[i1][i2][i3][j1][j2][j3] is deriv of est[i1](i2,i3) wrt. term est_counts[j1](j2,j3)
	f(i1,3){
		f(i2,arity){
			f(i3,arity){
				f(j1,3){
					f(j2,arity){
						f(j3,arity){
							derivatives[i1][i2][i3][j1][j2][j3] = 0.0;
						}
					}
				}
			}
		}
	}
	
	
	f(i,arity){
		f(j,arity){//don't have Wilson Score approximation here, so won't work well for edge values, and small tot values. 
			f(k,3){
				est_counts[k](i,j) = counts[k](i,j)/tot[k];
				dev_counts[k](i,j) = sqrt(counts[k](i,j)*(tot[k]-counts[k](i,j))/(tot[k]))/tot[k];
			}
			
			f(k,arity){
				dev_counts123[i][j][k] = sqrt(est_counts123[i][j][k]*(tot123-est_counts123[i][j][k])/(tot123))/tot123;
				est_counts123[i][j][k] = est_counts123[i][j][k]/tot123;
			}
		}
	}
	
	f(i,3){
		est[i] = MatSqrt(est_counts[(i+2)%3].transpose() * est_counts[(i+1)%3].inverse() * est_counts[i].transpose());
	}
	
	//testing new idea
	est[1] = est[0].transpose().inverse() * est_counts[0];
	special =  est[0].transpose().inverse() * special * est[1].inverse();
	MatrixXd toAvg = MatrixXd::Zero(arity,arity);
	
	f(i,arity){
		MatrixXd temp1 = est[0].transpose().inverse() * Triples[i] * est[1].inverse();
		JacobiSVD<MatrixXd> J = temp1.jacobiSvd(1<<2 | 1<<4);
		MatrixXd U = J.matrixU();
		MatrixXd temp2 = U.inverse() * est[0];
		f(j,arity){
			int max_index = -1;
			double max = -100;
			f(k,arity){
				if(temp2(j,k) > max){
					max_index = k;
					max = temp2(j,k);
				}
			}
			
			f(k,arity){
				toAvg(max_index,k) += temp2(j,k)/arity;
			}
		}
	
	}
	JacobiSVD<MatrixXd> J = special.jacobiSvd(1<<2 | 1<<4);
	MatrixXd U = J.matrixU();
	est[1] = toAvg;
	
	f(i,3){
		epsilon[i] = est[i];
	}
	return;
	//testing ends
	
	MatrixXd temp[3][2];
	
	f(j1,3){
		f(j2,arity){
			f(j3,arity){
				est_counts[j1](j2,j3) += z*dev_counts[j1](j2,j3);
				f(i1,3){
					temp[i1][0] = MatSqrt(est_counts[(i1+2)%3].transpose() * est_counts[(i1+1)%3].inverse() * est_counts[i1].transpose());
				}
				est_counts[j1](j2,j3) -= 2*z*dev_counts[j1](j2,j3);
				f(i1,3){
					temp[i1][1] = MatSqrt(est_counts[(i1+2)%3].transpose() * est_counts[(i1+1)%3].inverse() * est_counts[i1].transpose());
				}
				est_counts[j1](j2,j3) -= 2*z*dev_counts[j1](j2,j3);
				
				f(i1,3){
					f(i2,arity){
						f(i3,arity){
							derivatives[i1][i2][i3][j1][j2][j3] = (temp[i1][0](i2,i3) - temp[i1][1](i2,i3))/(2*z*dev_counts[j1](j2,j3));
						}
					}
				}
			}
		}
	}
	
	f(i1,3){
		epsilon[i1].resize(arity,arity);
		f(i2,arity){
			f(i3,arity){
				epsilon[i1](i2,i3) = 0.0;
				f(j1,3){
					f(j2,arity){
						f(j3,arity){
							f(k1,3){
								f(k2,arity){
									f(k3,arity){
										double covar = 0.0;
										if(j1==k1){
											if(j2==k2 && j3==k3){
												covar = dev_counts[j1](j2,j3)*dev_counts[j1](j2,j3);
											}
											else{
												(-est_counts[j1](j2,j3)*est_counts[k1](k2,k3))/tot[j1];
											}
										}
										else if(k1 == (j1+1)%3){
											if(j3==k2){
												double temp;
												if(j1==0){
													temp = est_counts123[j2][j3][k3];
												}
												else if(j1==1){
													temp = est_counts123[k3][j2][j3];
												}
												else{
													temp = est_counts123[j3][k3][j2];
												}
												tot123*(temp-est_counts[j1](j2,j3)*est_counts[k1](k2,k3))/(tot[j1]*tot[k1]);
											}
											else{
												tot123*(-est_counts[j1](j2,j3)*est_counts[k1](k2,k3))/(tot[j1]*tot[k1]);
											}
										}
										else{//j1 == (k1+1)%3
											if(k3==j2){
												double temp;
												if(k1==0){
													temp = est_counts123[k2][k3][j3];
												}
												else if(k1==1){
													temp = est_counts123[j3][k2][k3];
												}
												else{
													temp = est_counts123[k3][j3][k2];
												}
												tot123*(temp-est_counts[j1](j2,j3)*est_counts[k1](k2,k3))/(tot[j1]*tot[k1]);
											}
											else{
												tot123*(-est_counts[j1](j2,j3)*est_counts[k1](k2,k3))/(tot[j1]*tot[k1]);
											}
										}
										epsilon[i1](i2,i3) += derivatives[i1][i2][i3][j1][j2][j3] * derivatives[i1][i2][i3][k1][k2][k3] * covar;
									}
								}
							}
						}
					}
				}
				
				epsilon[i1](i2,i3) = z*sqrt(epsilon[i1](i2,i3));
			}
		}
	}

	return;
}



