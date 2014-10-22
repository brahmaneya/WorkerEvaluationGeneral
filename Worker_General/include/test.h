#ifndef TEST_H
#define TEST_H

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
#include<set>
#include<cassert>
#include<Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::JacobiSVD;
//Tests to estimate workers error rate

void generate_questions(bool questions[], int n, double s = 0.5);//populate questions array with random true answers as per selectivity s

void generate_answers(worker w[], int wn, int  n, bool questions[], bool** answers);//populate answers array with answers from workers

void expectation_maximisation(int wn, int n, bool ** answers, double est[]); //expectation maximisation estimate based on worker answers. could also return estimate for true answers, but easy to do once you have worker reliability estimates

double majority_error(double error_rates[], int n); //error rate of majority of n workers, with error rates given by error_rates[]

double majority_error_derivative(double error_rates[], int n, int k); //derivative of majority error rate of n workers, with respect to error rate of k^th worker

double erfinv(double conf, double threshold=0.000001);

bool weighted_majority(worker w[], int wn, bool ans[], double s = 0.5); // find answer using weighted majority

bool majority(int wn, bool ans[]); //find majority answer

//Old methods

void n_worker_conservative(bool** answers, int wn, int n, double est[], double epsilon[], double conf=0.9);

/////////New Efficient Functions and Generalizations begin here

void Gold_Standard(bool ans[], bool correct[], int n, double& est, double& epsilon, double conf=0.9); //n is number of questions, ans are workers answers, correct are correct answers, conf is confidence, est will be assigned 

void Three_Worker(bool** answers, int n, double est[], double epsilon[], double conf = 0.9);//answers is the double array of worker answers. answer[j][i] is j^th workers answer to i^th question

void N_Worker(bool ** answers, int wn, int n, double est[], double epsilon[], double conf = 0.9); //questions is the array of correct answers

void Three_Worker_nonreg(bool** answers, bool ** attempted, int n, double est[], double epsilon[], double conf = 0.9);//answers is the double array of worker answers. answer[j][i] is j^th workers answer to i^th question

void N_Worker_nonreg(bool ** answers, bool ** attempted, int wn, int n, double est[], double epsilon[], double conf = 0.9, bool noopt = false); //questions is the array of correct answers
//noopt flag means dont optimize the linear weights (a's)

void Gold_Standard_nonreg_Kary(int arity, bool attempted[], int ans[], int correct[], int n, MatrixXd &est, MatrixXd &epsilon, double conf = 0.9); 

void Three_Worker_nonreg_Kary(int arity, bool ** attempted, int ** answers, int n, MatrixXd est[], MatrixXd epsilon[], double conf = 0.9);

MatrixXd MatSqrt(MatrixXd input); // Symmetric Square root of a positive Matrix

#endif
