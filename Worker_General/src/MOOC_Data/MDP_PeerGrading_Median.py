import csv
import sys
import random
import time
import numpy as np
import operator
from collections import defaultdict
#import matplotlib.pyplot as plt

goldSet ={} # submission_id: grade
assignmentNum={} # submission_id: assignment_no
graderType = {} # grader_id: grader_type
grade={} # submission_id, grader_id, grade 
finalGrade={} # submission_id: peer_grade, staff_grade
allGrades = defaultdict(list) # submission_id: [all (graderId, grade) assigned]
totalEval = defaultdict(int) #Number of evaluations done by the grader
avg_groundTruth = {} # Single Class Ground Set
ind_idClass = {} # Individual Error Distribution IdClass
ind_e = {} # Individual Error Distribution
ind_costClass = {}
ind_probClass = {}    
groundTruth = {} # Individual Error Distribution Ground Set
Prior = {}
majoritySet = {}


pen1={}
for i in range(1,7):
    for j in range(1,7):
        if j==i or j==i-1 or j==i+1:
            pen1[(i,j)] =  0
        else: pen1[(i,j)] = 1
        
pen2={}
for i in range(1,7):
    for j in range(1,7):
        if j==i :
            pen2[(i,j)] =  0
        else: pen2[(i,j)] = 1

pen3={}
for i in range(1,7):
    for j in range(1,7):
        pen3[(i,j)]=abs(j-i)
    

m = 5 #5
upper = 100.0
c_gran = 1
p_gran = 10
p_size = int(upper/p_gran)+1
c_size = int((m+1)/c_gran)
beta = 1
factor = float(p_gran)/float(upper)

# Optimizing penalty function
pen = pen3
# Reporting penalty function
penReport = pen3

# Per Grade threshold to consider Average Distribution
# instead of the Individual error distribution
confidence_threshold = 10 #20
# Worker threshold to remove low-participating graders 
worker_threshold = 90 #90

def mapG(grade):
    if (grade=='A*'): return 1
    elif (grade=='A'): return 2
    elif (grade=='B'): return 3
    elif (grade=='C'): return 4
    elif (grade=='D'): return 5
    elif (grade=='F'): return 6
    else:
        # Some median grades are absent
        print 'NULL Grade' 
        return random.choice(range(1,7))

# Required since median is taken as average and 
# hence fractional scores in the peer grade
def mgrade(s):
    if s=="-1":
        return "NULL"
    score = float(s)
    if (score>=4.5): return 'A*'
    elif(score>=3.5): return 'A'
    elif(score>=2.5): return 'B'
    elif(score>=1.5): return 'C'
    elif(score>=0.5): return 'D'
    return 'F'

# Calculating max probability grade for a particular
# submission while forming the ground set
def maxGrade(gradeList, idClass, e):
    global Prior
    if gradeList==[]: return 'NULL'
    maxV = 0.0
    maxG = ""
    total = 0.0
    
    for grade in ['A*','A','B','C','D','F']:
        p_evidence_grade = 1000000.0
        count = 0;
        for evidence in gradeList:
            if (totalEval[evidence[0]]<worker_threshold): continue
            if (count==m): break
            p_evidence_grade = p_evidence_grade * e[idClass[evidence[0]]][(evidence[1],grade)]
            count +=1
        numerator = p_evidence_grade * Prior[grade]
        total += numerator
        
        if(numerator >= maxV): 
            maxV = numerator
            maxG = grade
    
    if total==0: confidence=0
    else: confidence = float(maxV)/float(total)
    return maxG, confidence


# For plotting distribution (average disagreement per grader)
# from dictionary containing disagreement of each grader
def confidenceVectors(confidence):
    confidence_data = []
    for key in confidence.keys():
        confidence_data.append(confidence[key])
        
    factor = 10
    tempData = []
    for i in confidence_data:
        tempData.append(int(i*factor))
    
    frequencies = defaultdict(int)
    total = 0
    for num in tempData:
            frequencies[num] += 1
            total += 1
    x_iter = []
    y_iter = []
    for key in frequencies.keys():
        x_iter.append(key/float(factor))
        y_iter.append(frequencies[key]/float(total))
    
    return (x_iter,y_iter)

def loadData(filename):
    global groundTruth
    global totalEval
    global ind_e
    global ind_idClass
    global avg_e
    original = file(filename, 'rU')
    reader = csv.reader(original)
    
    # Format of the CSV file
    # "submission_id",assignment_num,"grader_id","grader_group",grade,peer_grade,staff_grade
    temp1 = 0
    temp2 = 0
    for row in reader:
        #print row
        if(int(row[1])==1): continue
        temp1+=1
        assignmentNum[row[0]] = int(row[1])
        graderType[row[2]] = row[3]
        if(row[3]=='staff'):
            temp2 += 1
            goldSet[row[0]] = mgrade(row[4])
            finalGrade[row[0]] = (mgrade(row[4]),mgrade(row[6]))  
        else:
            if (mgrade(row[5]) != "NULL"):
                finalGrade[row[0]] = (mgrade(row[5]),mgrade(row[6]))  
        grade[(row[0],row[2])] = mgrade(row[4])
    print "Total Submissions: ", len(finalGrade.keys())
    print "Total Evaluations: ", temp1
    print "Total Staff Evaluations: ", temp2
    
    ### Computing Prior
    totalPrior = float(0)
    for i in ['A*','A','B','C','D','F']: Prior[i]=float(0)
    
    for i in finalGrade.keys():
        if (finalGrade[i][0]=="NULL"): continue
        Prior[finalGrade[i][0]] += 1
        totalPrior += 1
#    for i in goldSet.keys():
#        Prior[goldSet[i]] += 1
#        totalPrior += 1
    for i in ['A*','A','B','C','D','F']: 
        Prior[i]=Prior[i]/float(totalPrior)
        #Prior[i] = 1.0/6.0

    for (subId,graderId) in grade.keys():
        allGrades[subId].append((graderId,grade[(subId,graderId)]))  

    
    idClass = {}
    countGrader = 0
    for graderId in graderType.keys():
        idClass[graderId] = countGrader
        countGrader += 1
        
    # Calculating Individual Error Distribution
    e = {}
    count = {}
    for grader in graderType.keys():
        e[idClass[grader]] = {}
        count[idClass[grader]] = {}
      
    for grader in graderType.keys():
        for i in ['A*','A','B','C','D','F']:
            count[idClass[grader]][i] = 0.0
            for j in ['A*','A','B','C','D','F']:
                e[idClass[grader]][(i,j)]=0.0
                
    for i in finalGrade.keys():
        for (k,x) in allGrades[i]:
            count[idClass[k]][finalGrade[i][0]]+=1.0
            e[idClass[k]][ (x,finalGrade[i][0]) ] +=1.0
#    for i in goldSet.keys():
#        for (j,k) in grade.keys():
#            if (i==j):
#                count[k][goldSet[i]] += 1.0 
#                e[k][ (grade[(j,k)],goldSet[i]) ] +=1.0
                
    for grader in graderType.keys():
        for i in ['A*','A','B','C','D','F']:
            totalEval[grader] += count[idClass[grader]][i]
    
    avg_idClass = singleClass()
    avg_e = makeClass(avg_idClass,1)[2]
    avg_confidence = {}
    for subId in allGrades.keys():
        (avg_groundTruth[subId],avg_confidence[subId]) = maxGrade(allGrades[subId], avg_idClass,avg_e)

    for grader in graderType.keys():    
        for i in ['A*','A','B','C','D','F']:
            for j in ['A*','A','B','C','D','F']:
                if (count[idClass[grader]][i]<confidence_threshold):
                    e[idClass[grader]][(j,i)] = avg_e[0][(j,i)]
                else:
                    e[idClass[grader]][(j,i)] = e[idClass[grader]][(j,i)]/float(count[idClass[grader]][i])
                e[idClass[grader]][(mapG(j),mapG(i))] = e[idClass[grader]][(j,i)]
        
#    print "\nAvg Error Distribution\n"
#    for i in ['A*','A','B','C','D','F']:
#        for j in ['A*','A','B','C','D','F']: 
#            print (j,i),'---->',avg_e[0][(j,i)]
#        print ''
#                        
#    print "\nIndividual Error Dist\n"
#    for grader in graderType.keys():
#        print 'Grader id ', grader
#        for i in ['A*','A','B','C','D','F']:
#            for j in ['A*','A','B','C','D','F']: 
#                print (j,i),'---->',e[grader][(j,i)]
#            print '' 
        
    ind_e = e
    ind_idClass = idClass
    numGraders= len(graderType.keys())
    for i in graderType.keys():
        ind_costClass[ind_idClass[i]]=1.0
        ind_probClass[ind_idClass[i]] =1.0/float(numGraders)
        
    confidence = {}
    for subId in allGrades.keys():
        (groundTruth[subId],confidence[subId]) = maxGrade(allGrades[subId], idClass,e)
    
# Calculating the MDP strategy for given alpha and Grader Classes
def mdp(pen,m,alpha,numClass, probClass, costClass, e, memo_cumE, memo_nextState):
    
    ## Defining States
    statePts=[]
    for p1 in np.arange(0,p_size):
        for p2 in np.arange(0,p_size):
            for p3 in np.arange(0,p_size):
                for p4 in np.arange(0,p_size):
                    for p5 in np.arange(0,p_size):
                        if p1+p2+p3+p4+p5>p_size-1: continue
                        for c in np.arange(0,m+1,c_gran):
                                statePts.append((p1,p2,p3,p4,p5,c))
    
    ### Defining Actions
    action = {}
    action[1]='stop'
    action[2]='continue'
    
    #### Computing Reward: Reward(s) when s is terminating
    reward = np.zeros((p_size,p_size,p_size,p_size,p_size,c_size))
    minProb = np.zeros((p_size,p_size,p_size,p_size,p_size,c_size))
    termAnswer = np.zeros((p_size,p_size,p_size,p_size,p_size,c_size))
    
#    for p1 in np.arange(0,p_size):
#        for p2 in np.arange(0,p_size):
#            for p3 in np.arange(0,p_size):
#                for p4 in np.arange(0,p_size):
#                    for p5 in np.arange(0,p_size):
#                        if p1+p2+p3+p4+p5>p_size-1: continue
#                        p6 = p_size-1-p1-p2-p3-p4-p5
#                        (temp1,temp2) = min( (p2*pen[(1,2)] + p3*pen[(1,3)] + p4*pen[(1,4)] + p5*pen[(1,5)] + p6*pen[(1,6)],1), 
#                                             (p1*pen[(2,1)] + p3*pen[(2,3)] + p4*pen[(2,4)] + p5*pen[(2,5)] + p6*pen[(2,6)],2),
#                                             (p1*pen[(3,1)] + p2*pen[(3,2)] + p4*pen[(3,4)] + p5*pen[(3,5)] + p6*pen[(3,6)],3),
#                                             (p1*pen[(4,1)] + p2*pen[(4,2)] + p3*pen[(4,3)] + p5*pen[(4,5)] + p6*pen[(4,6)],4),
#                                             (p1*pen[(5,1)] + p2*pen[(5,2)] + p3*pen[(5,3)] + p4*pen[(5,4)] + p6*pen[(5,6)],5),
#                                             (p1*pen[(6,1)] + p2*pen[(6,2)] + p3*pen[(6,3)] + p4*pen[(6,4)] + p5*pen[(6,5)],6)                
#                                             )
#                        for c in np.arange(0,m+1,c_gran):
#                                (minProb[p1][p2][p3][p4][p5][c],termAnswer[p1][p2][p3][p4][p5][c]) = (temp1,temp2)    
#                                reward[p1][p2][p3][p4][p5][c] = alpha*c + beta*(float(minProb[p1][p2][p3][p4][p5][c])*factor)
#                                
#                                
    for (p1,p2,p3,p4,p5,c) in statePts:
        p6 = p_size-1-p1-p2-p3-p4-p5
        (minProb[p1][p2][p3][p4][p5][c],termAnswer[p1][p2][p3][p4][p5][c]) = min( (p2*pen[(1,2)] + p3*pen[(1,3)] + p4*pen[(1,4)] + p5*pen[(1,5)] + p6*pen[(1,6)],1), 
        (p1*pen[(2,1)] + p3*pen[(2,3)] + p4*pen[(2,4)] + p5*pen[(2,5)] + p6*pen[(2,6)],2),
        (p1*pen[(3,1)] + p2*pen[(3,2)] + p4*pen[(3,4)] + p5*pen[(3,5)] + p6*pen[(3,6)],3),
        (p1*pen[(4,1)] + p2*pen[(4,2)] + p3*pen[(4,3)] + p5*pen[(4,5)] + p6*pen[(4,6)],4),
        (p1*pen[(5,1)] + p2*pen[(5,2)] + p3*pen[(5,3)] + p4*pen[(5,4)] + p6*pen[(5,6)],5),
        (p1*pen[(6,1)] + p2*pen[(6,2)] + p3*pen[(6,3)] + p4*pen[(6,4)] + p5*pen[(6,5)],6)                
        )
        # TODO: We do not need to calculate for every cost since only dependent on posterior
        reward[p1][p2][p3][p4][p5][c] = alpha*c + beta*(float(minProb[p1][p2][p3][p4][p5][c])*factor)
    
    ### Computing Value Function
    def computeValue(state,graderClass,valueFunc):
        value = 0.0
        p =[0]*7
        (p[1],p[2],p[3],p[4],p[5],c) = state
        c = c+costClass[graderClass]
        for g in range(1,7):
            (p1,p2,p3,p4,p5) = memo_nextState[p[1]][p[2]][p[3]][p[4]][p[5]][graderClass][g]
            #print (p[1],p[2],p[3],p[4],p[5]), graderClass, g, (p1, p2, p3, p4, p5, c)
            #print "CumE:", cum_e[graderClass][g], "New Temp: ", cumE(state,graderClass,g)
            value += memo_cumE[p[1]][p[2]][p[3]][p[4]][p[5]][graderClass][g]*valueFunc[p1][p2][p3][p4][p5][c]
        return value

    #### Computing Strategy    
    pi = np.zeros((p_size,p_size,p_size,p_size,p_size,c_size))
    value = np.zeros((p_size,p_size,p_size,p_size,p_size,c_size))
    
    c = m
    for p1 in np.arange(0,p_size):
        for p2 in np.arange(0,p_size):
            for p3 in np.arange(0,p_size):
                for p4 in np.arange(0,p_size):
                    for p5 in np.arange(0,p_size):
                        if p1+p2+p3+p4+p5>p_size-1: continue
                        pi[p1][p2][p3][p4][p5][c]=1
                        value[p1][p2][p3][p4][p5][c]=reward[p1][p2][p3][p4][p5][c]

            
    #sample=0.01                
    for c in range(m-c_gran,-1,-c_gran):
        for p1 in np.arange(0,p_size):
            for p2 in np.arange(0,p_size):
                for p3 in np.arange(0,p_size):
                    for p4 in np.arange(0,p_size):
                        for p5 in np.arange(0,p_size):
                            if p1+p2+p3+p4+p5>p_size-1: continue
                            
#                            temp = 0
#                            for graderClass in range(numClass):
#                                if random.random() < sample:
#                                    temp += (1/sample)*probClass[graderClass]*computeValue((p1,p2,p3,p4,p5,c),graderClass,value)

                            temp = 0
                            for graderClass in range(numClass):
                                temp += probClass[graderClass]*computeValue((p1,p2,p3,p4,p5,c),graderClass,value)
                            
                            if (reward[p1][p2][p3][p4][p5][c]<=temp):
                                pi[p1][p2][p3][p4][p5][c]=1 # STOP
                                value[p1][p2][p3][p4][p5][c]=reward[p1][p2][p3][p4][p5][c]
                            else:
                                pi[p1][p2][p3][p4][p5][c]=2 # CONTINUE
                                value[p1][p2][p3][p4][p5][c]=temp
    
    return (pi,termAnswer, minProb)


# Simulating Strategy on Real Data
def runStrategy(pi,termAnswer, penReport, idClass, probClass, costClass, e):
    stratAns = {}
    totalQuestions = 0
    stratQuestions = 0
    
    # Calculating next state depending on graders answer g
    # belonging to graderClass at a given state
    def nextState(state,graderClass,g):
        p =[0]*7
        p_new =[0]*7
        (p[1],p[2],p[3],p[4],p[5],c) = state
        p[6] = p_size-1-(p[1]+p[2]+p[3]+p[4]+p[5])
        if (p[6]<0):
            print 'p6 negative and state', state
            sys.exit(0)
        
        c_new = c+costClass[graderClass]
        norm = 0.0
        for grade in range(1,7):
            p_new[grade]= e[graderClass][(g,grade)]*p[grade]
            norm += p_new[grade]
            
        for grade in range(1,7):
            if norm==0: 
                p_new[grade]=1.0/6.0
            else:
                p_new[grade]=p_new[grade]/float(norm) 
        return (int(int(p_new[1]*upper)/p_gran),int(int(p_new[2]*upper)/p_gran),int(int(p_new[3]*upper)/p_gran),int(int(p_new[4]*upper)/p_gran),int(int(p_new[5]*upper)/p_gran),c_new)
    
    
    for sub_id in allGrades.keys():
        grades_old = allGrades[sub_id]
        
        # Eliminating Low Confidence Workers while simulating
        grades = []
        for (grader,gradeGiven) in grades_old:
            if totalEval[grader]>=worker_threshold:
                grades.append((grader,gradeGiven))
        #AGP
        if sub_id in goldSet.keys(): continue

        # Considering only the first m grades and
        # eliminating all items with <m grades
        #AGP 
        temp = min(m, len(grades))
        size = temp-1
        

        answer =''
        (a,b,c,d,g,f) = (int(int(Prior['A*']*upper)/p_gran),int(int(Prior['A']*upper)/p_gran),int(int(Prior['B']*upper)/p_gran),int(int(Prior['C']*upper)/p_gran),int(int(Prior['D']*upper)/p_gran),0)

        answerOver = False
        while (answer ==''):
            #print (a,b,c,d,g,f), size
            if (pi[a][b][c][d][g][f]==2):
                #print 'continuing'
                if(size== -1):
                    # Will never reach here since we are eliminating <m items
                    #print 'answers over', sub_id
                    #answerOver = True
                    #break;
                    answer = termAnswer[a][b][c][d][g][f]
                else: 
                    if grades[size][1]=='A*': (a,b,c,d,g,f)=nextState((a,b,c,d,g,f),idClass[grades[size][0]],1)
                    elif grades[size][1]=='A':(a,b,c,d,g,f)=nextState((a,b,c,d,g,f),idClass[grades[size][0]],2)
                    elif grades[size][1]=='B':(a,b,c,d,g,f)=nextState((a,b,c,d,g,f),idClass[grades[size][0]],3)
                    elif grades[size][1]=='C':(a,b,c,d,g,f)=nextState((a,b,c,d,g,f),idClass[grades[size][0]],4)
                    elif grades[size][1]=='D':(a,b,c,d,g,f)=nextState((a,b,c,d,g,f),idClass[grades[size][0]],5)
                    else: (a,b,c,d,g,f)=nextState((a,b,c,d,g,f),idClass[grades[size][0]],6)
                    size-=1
            else: 
                #print 'terminating'
                answer = termAnswer[a][b][c][d][g][f]
        
        if(answerOver==False):
            totalQuestions += temp
            stratQuestions += temp-(size+1)
            stratAns[sub_id] = int(answer)
        #print 'Answer for', sub_id,':', stratAns[sub_id], '\n' 

    # Calculating reporting Cost and Error using penReport
    count = 0
    total = 0
    for sub_id in stratAns.keys():
        j = int(stratAns[sub_id])
        i = int(mapG(groundTruth[sub_id]))
        count += penReport[(i,j)]
        total += 1
    
    print 'Running Strategy Stats', stratQuestions, totalQuestions, count, total
    return ( (float(stratQuestions)/ float(totalQuestions))*100.0 ,(float(count)/total))


# Forming idClass belonging to single class
# all graders considered identical
def singleClass():
    idClass = {}
    for graderId in graderType.keys():
        ## Comment the IF part if you don't want to eliminate 
        ## low-confidence workers from Error matrix calculating
        if (totalEval[graderId]<worker_threshold):
            idClass[graderId]=-1
        else:
            idClass[graderId] = 0
    return idClass


# Forming idClass based on variance
# Dividing into equal numClasses
def varianceClass(numClass):
    idClass = {}
    score = defaultdict(int)
    total = defaultdict(int)
    for grader in graderType.keys():
        score[grader] = 0
        total[grader] = 0
        
    # Calculating Variances of each grader
    for subId in finalGrade.keys():
        for (grader,grade_given) in allGrades[subId]:
            i = mapG(finalGrade[subId][0])
            j = mapG(grade_given)
            total[grader]+=1
            score[grader]+=abs(i-j)
            
    for key in score.keys():
        ## Comment the below IF part if you don't want to eliminate 
        ## low-confidence workers from Error matrix calculating
        if (total[key]<worker_threshold):
            idClass[key]=-1
            del score[key]
            continue
        
        if total[key]!=0: 
            score[key] = score[key]/float(total[key])
    
    length = len(score.keys())
    print 'dividing', length, 'into ', numClass, 'classes'
    size = length/numClass   
    sorted_score = sorted(score.iteritems(), key=operator.itemgetter(1))
    for classNo in range(0,numClass):
        if(classNo==numClass-1):
            for grader, scoreGrader in sorted_score[classNo*size:]:
                idClass[grader]=classNo
        else:
            for grader, scoreGrader in sorted_score[classNo*size:(classNo+1)*size]:
                idClass[grader]=classNo
        
    return idClass


# Forming idClass based on bias
# Dividing into equal numClasses
def biasClass(numClass):
    idClass = {}
    score = defaultdict(int)
    total = defaultdict(int)
    for grader in graderType.keys():
        score[grader] = 0
        total[grader] = 0
        
    # Calculating Bias of each grader
    for subId in finalGrade.keys():
        for (grader,grade_given) in allGrades[subId]:
            if graderType[grader]=='staff': continue
            i = mapG(finalGrade[subId][0])
            j = mapG(grade_given)
            total[grader]+=1
            score[grader]+=i-j
    
    for key in score.keys():
        ## Comment the below IF part if you don't want to eliminate 
        ## low-confidence workers from Error matrix calculating
        if (total[key]<worker_threshold):
            idClass[key]=-1
            del score[key]
            continue
        
        if total[key]!=0: 
            score[key] = score[key]/float(total[key])
    
    length = len(score.keys())
    print 'dividing', length, 'into ', numClass, 'classes'
    size = length/numClass   
    sorted_score = sorted(score.iteritems(), key=operator.itemgetter(1))
    for classNo in range(0,numClass):
        if(classNo==numClass-1):
            for grader, scoreGrader in sorted_score[classNo*size:]:
                idClass[grader]=classNo
        else:
            for grader, scoreGrader in sorted_score[classNo*size:(classNo+1)*size]:
                idClass[grader]=classNo
        
    return idClass


# Getting error distribution, probability according to idClass
def makeClass(idClass, numClass):
    global avg_e
    numGrader=0.0
    count = [0.0]*numClass
    for grader in idClass.keys():
        if (idClass[grader]==-1):continue
        numGrader += 1
        count[idClass[grader]] +=1
    costClass=[]
    probClass=[]
    for it in range(numClass):
        print 'Size of Class', it,':',count[it]
        costClass.append(1)
        probClass.append(count[it]/float(numGrader))
        
    # calculating Error distribution of Grader Classes
    # specified by the idClass
    e = []
    for it in range(numClass):
        e.append({})
    
    for it in range(numClass):
        
        for i in ['A*','A','B','C','D','F']:
            for j in ['A*','A','B','C','D','F']:
                e[it][(i,j)]=0.0
        count ={}
        for i in ['A*','A','B','C','D','F']: count[i]=0
        # Fill-up e[it] using gold_set data
#        for i in goldSet.keys():
#            for (j,k) in grade.keys():
#                if (i==j and idClass[k]==it and graderType[k]!='staff'):
#                    count[goldSet[i]] += 1.0 
#                    e[it][ (grade[(j,k)],goldSet[i]) ] +=1.0  
                    
        for i in finalGrade.keys():
            for (k,x) in allGrades[i]:
                if( idClass[k]==it and graderType[k]!='staff'):
                    count[finalGrade[i][0]]+=1.0
                    e[it][ (x,finalGrade[i][0]) ] +=1.0            
                          
        for i in ['A*','A','B','C','D','F']:
            for j in ['A*','A','B','C','D','F']:
                if (count[i]==0):
                    print 'No grades seen for class', it, 'for grade', i
                    e[it][(j,i)] = avg_e[0][(j,i)]
                    #e[it][(j,i)] = 1.0/6.0
                else:
                    e[it][(j,i)] = e[it][(j,i)]/float(count[i])
                e[it][(mapG(j),mapG(i))] = e[it][(j,i)]
        
    return (probClass, costClass, e)


# Plotting Random Strategy cost-error curve
def plotRandomCostError():
    plotV = {}       
    for cost in [0.3,0.35,0.4,0.45,0.5,0.55,0.6]:
        randomGroundTruth = {}
        randomConfidence = {}
        randomAllGrades = defaultdict(list)
            
        for subId in allGrades.keys():
            #AGP
            if subId in goldSet.keys(): 
                continue 
            for (grader,gradeGiven) in allGrades[subId]:
                if (random.random()<cost):
                    randomAllGrades[subId].append((grader,gradeGiven))
            
        for subId in randomAllGrades.keys():
            (randomGroundTruth[subId],randomConfidence[subId]) = maxGrade(randomAllGrades[subId],ind_idClass,ind_e)
        
        count = 0
        total = 0
        for sub_id in randomGroundTruth.keys():
            j = int(mapG(groundTruth[sub_id]))
            i = int(mapG(randomGroundTruth[sub_id]))
            count += penReport[(i,j)]
            total += 1
            
        print cost, (cost*100,count/float(total))
        plotV[count/float(total)] = cost*100
    
    errorV = []
    costV = []
    for key in sorted(plotV.iterkeys()):
        errorV.append(key)
        costV.append(plotV[key])
    return (costV, errorV)

def median(grades):
    grades_new = []
    
    grade_to_int = {}
    grade_to_int["A*"] = 0
    grade_to_int["A"] = 1
    grade_to_int["B"] = 2
    grade_to_int["C"] = 3
    grade_to_int["D"] = 4
    grade_to_int["F"] = 5
    int_to_grade = {}
    for x in grade_to_int.keys():
        int_to_grade[grade_to_int[x]] = x
        
    for g in grades:
        grades_new.append(grade_to_int[g])
    grades_new.sort()
    
    return int_to_grade[grades_new[len(grades_new)/2]]

# Plotting Random Strategy cost-error curve
def plotRandomMedianCostError():
    plotV = {}       
    for cost in [0.3,0.35,0.4,0.45,0.5,0.55,0.6]:
        randomMedianGroundTruth = {}
        randomAllGrades = defaultdict(list)
            
        for subId in allGrades.keys():
            #AGP
            if subId in goldSet.keys(): 
                continue 
  
            for (grader,gradeGiven) in allGrades[subId]:
                if (random.random()<cost):
                    randomAllGrades[subId].append(gradeGiven)
            
        for subId in randomAllGrades.keys():
            randomMedianGroundTruth[subId] = median(randomAllGrades[subId])
        
        count = 0
        total = 0
        for sub_id in randomMedianGroundTruth.keys():
            j = int(mapG(groundTruth[sub_id]))
            i = int(mapG(randomMedianGroundTruth[sub_id]))
            count += penReport[(i,j)]
            total += 1
            
        print cost, (cost*100,count/float(total))
        plotV[count/float(total)] = cost*100
    
    errorV = []
    costV = []
    for key in sorted(plotV.iterkeys()):
        errorV.append(key)
        costV.append(plotV[key])
    return (costV, errorV)

# Plotting cost-error curve for different alpha, different algorithms
def plotCostError(numClass,idClass, filename):
    (probClass, costClass, e) = makeClass(idClass,numClass)       
#    
#    probClass = ind_probClass
#    costClass = ind_costClass
#    e = ind_e
    
    # Computing Grader Probability of giving a particular grade at a given state
    def cumE(state,graderClass,g):
        p =[0]*7
        (p[1],p[2],p[3],p[4],p[5]) = state
        p[6] = p_size-1 -(p[1]+p[2]+p[3]+p[4]+p[5])
        
        temp=0
        for j in range(1,7):
            temp += e[graderClass][(g,j)]*p[j]
        return (float(temp)*float(p_gran)/float(upper))
    
     #### Computing Next State
    def nextState(state,graderClass,g):
        p =[0]*7
        p_new =[0]*7
        (p[1],p[2],p[3],p[4],p[5]) = state
        p[6] = p_size-1 -(p[1]+p[2]+p[3]+p[4]+p[5])
        if (p[6]<0):
            print 'p6 negative and state', state
            sys.exit(0)
        
        norm = 0.0
        for grade in range(1,7):
            p_new[grade]= e[graderClass][(g,grade)]*p[grade]
            norm += p_new[grade]
            
        for grade in range(1,7):
            if norm==0: 
                p_new[grade]=1.0/6.0
            else:
                p_new[grade]=p_new[grade]/float(norm)  
        return (int(int(p_new[1]*upper)/p_gran),int(int(p_new[2]*upper)/p_gran),int(int(p_new[3]*upper)/p_gran),int(int(p_new[4]*upper)/p_gran),int(int(p_new[5]*upper)/p_gran))


    memo_cumE = np.zeros((p_size,p_size,p_size,p_size,p_size,numClass,7))
    dt = np.dtype("i4,i4,i4,i4,i4")
    memo_nextState = np.zeros((p_size,p_size,p_size,p_size,p_size,numClass,7),dt)
    
    print "Computation start", time.time()
    for p1 in np.arange(0,p_size):
        for p2 in np.arange(0,p_size):
            #print "p3 done", time.time()
            for p3 in np.arange(0,p_size):
                for p4 in np.arange(0,p_size):
                    for p5 in np.arange(0,p_size):
                        if p1+p2+p3+p4+p5>p_size-1: continue
                        for graderClass in range(numClass):
                            for gradeGiven in range(1,7):
                                memo_cumE[p1][p2][p3][p4][p5][graderClass][gradeGiven] = cumE((p1,p2,p3,p4,p5),graderClass,gradeGiven)
                                memo_nextState[p1][p2][p3][p4][p5][graderClass][gradeGiven] = nextState((p1,p2,p3,p4,p5),graderClass,gradeGiven)
    print "Computation ends", time.time()
    
    plotV = {}
    #for alpha in range(0,40):
    #    alpha = float(alpha)/100.0
    for alpha in [0.1,0.2,0.3,0.4,0.5]:
        temp = mdp(pen, m, alpha, numClass, probClass, costClass, e, memo_cumE, memo_nextState)
        #temp = mdp(pen, m, alpha, numClass, ind_probClass, ind_costClass, ind_e, memo_cumE,memo_nextState)
        
        pi = temp[0]
        termAnswer = temp[1]
        if(filename=='single'):
            # Comment/Uncomment if want to simulate for individual error dist for single class too
            #(cost,error)=runStrategy(pi,termAnswer, penReport, idClass, numClass, probClass, costClass, e)
            (cost,error)=runStrategy(pi,termAnswer, penReport, ind_idClass, ind_probClass, ind_costClass, ind_e)
        else:
            #(cost,error)=runStrategy(pi,termAnswer, penReport, idClass, numClass, probClass, costClass, e)
            (cost,error)=runStrategy(pi,termAnswer, penReport, ind_idClass, ind_probClass, ind_costClass, ind_e)
        print alpha, (cost,error)
        plotV[error] = cost
    
    errorV = []
    costV = []
    for key in sorted(plotV.iterkeys()):
        errorV.append(key)
        costV.append(plotV[key])

    return (costV, errorV)


# Reporting the disagreement of the median algorithm with our ground set
def reportMedianCostError(penReport):
    count = 0
    total = 0
    for sub_id in finalGrade.keys():
        j = int(mapG(finalGrade[sub_id][0]))
        i = int(mapG(groundTruth[sub_id]))
        count += penReport[(i,j)]
        total += 1
    print '\nRunning Median Strategy Stats', count, total, float(count)/total   


# Plotting the Theoretical expected cost-error curve for different alpha's, different algorithms
def plotTheoreticalCostError(numClass, idClass, filename):
    (probClass, costClass, e) = makeClass(idClass,numClass)
    def nextState(state,graderClass,g):
        p =[0]*7
        p_new =[0]*7
        (p[1],p[2],p[3],p[4],p[5],c) = state
        p[6] = p_size-1-(p[1]+p[2]+p[3]+p[4]+p[5])
        if (p[6]<0):
            print 'p6 negative and state', state
            sys.exit(0)
        
        c_new = c+costClass[graderClass]
        norm = 0.0
        for grade in range(1,7):
            p_new[grade]= e[graderClass][(g,grade)]*p[grade]
            norm += p_new[grade]
            
        for grade in range(1,7):
            if norm==0: 
                p_new[grade]=1.0/6.0
            else:
                p_new[grade]=p_new[grade]/float(norm) 
        return (int(int(p_new[1]*upper)/p_gran),int(int(p_new[2]*upper)/p_gran),int(int(p_new[3]*upper)/p_gran),int(int(p_new[4]*upper)/p_gran),int(int(p_new[5]*upper)/p_gran),c_new)
    def cumE(state,graderClass,g):
        p =[0]*7
        (p[1],p[2],p[3],p[4],p[5],c) = state
        p[6] = p_size-1 -(p[1]+p[2]+p[3]+p[4]+p[5]) +c-c
        
        temp=0
        for j in range(1,7):
            temp += e[graderClass][(g,j)]*p[j]
        return (float(temp)*float(p_gran)/float(upper))
    def expectedCostError(pi, minProb):
        prob = np.zeros((p_size,p_size,p_size,p_size,p_size,c_size))
        (a,b,c,d,g,f) = (int(int(Prior['A*']*upper)/p_gran),int(int(Prior['A']*upper)/p_gran),int(int(Prior['B']*upper)/p_gran),int(int(Prior['C']*upper)/p_gran),int(int(Prior['D']*upper)/p_gran),0)
        #print 'Prior State', (a,b,c,d,g,f)
        prob[a][b][c][d][g][f] = 1.0
        for c in np.arange(0,m,c_gran):
            for p1 in np.arange(0,p_size):
                for p2 in np.arange(0,p_size):
                    for p3 in np.arange(0,p_size):
                        for p4 in np.arange(0,p_size):
                            for p5 in np.arange(0,p_size):
                                if p1+p2+p3+p4+p5>p_size-1: continue
                                if (pi[p1][p2][p3][p4][p5][c]==1): continue
                                if (prob[p1][p2][p3][p4][p5][c]==0): continue
                                
                                #print  (p1,p2,p3,p4,p5,c), ':', prob[p1][p2][p3][p4][p5][c]
                                for graderClass in range(numClass):
                                    for grade in range(1,7):
                                        (u,n,o,p,q,r) = nextState((p1,p2,p3,p4,p5,c),graderClass,grade)
                                        prob[u][n][o][p][q][r] += prob[p1][p2][p3][p4][p5][c]*probClass[graderClass]*cumE((p1,p2,p3,p4,p5,c),graderClass,grade)
                                        #print 'Grader Class, Grade, Prob, Next State: ',graderClass,grade,(u,n,o,p,q,r) , prob[u][n][o][p][q][r]
    
        expect_cost = 0.0 
        expect_error= 0.0
        factor = 1.0
        for c in np.arange(0,c_size):
                for p1 in np.arange(0,p_size):
                    for p2 in np.arange(0,p_size):
                        for p3 in np.arange(0,p_size):
                            for p4 in np.arange(0,p_size):
                                for p5 in np.arange(0,p_size):
                                    if p1+p2+p3+p4+p5>p_size-1: continue
                                    if (pi[p1][p2][p3][p4][p5][c]==2): continue
                                    expect_cost += prob[p1][p2][p3][p4][p5][c]*c*factor
                                    expect_error += prob[p1][p2][p3][p4][p5][c]*(float(minProb[p1][p2][p3][p4][p5][c])*float(p_gran)/float(upper))*factor

        return (float(expect_cost)/factor, float(expect_error)/factor)
        
    plotV = {}
    #for alpha in range(0,102,2):
    #   alpha = float(alpha)/100.0
    for alpha in [0.005,0.01,0.1,0.2,0.3]:
        temp = mdp(pen, m, alpha, numClass, probClass, costClass, e)
        pi = temp[0]
        minProb=temp[2]
        
        (cost,error) = expectedCostError(pi,minProb)
        print alpha, (cost,error)
        plotV[error] = cost
    
    errorV = []
    costV = []
    for key in sorted(plotV.iterkeys()):
        errorV.append(key)
        costV.append(plotV[key])
    return (costV, errorV)


# Plot Graph 1, Random vs Single VS Variance vs Bias
def main1():
    loadData('evaluation_data.csv')
    reportMedianCostError(penReport)
    
    print '\nRandom Strategy Analysis'
    (rand_cost,rand_error) = plotRandomCostError()
    
    print '\nRandom Median Strategy Analysis'
    (randM_cost,randM_error) = plotRandomMedianCostError()
    
    time0 = time.time()
    print '\nVariance Analysis'
    numClass=2
    idClass = varianceClass(numClass)
    #(var_Tcost, var_Terror) = plotTheoreticalCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    #print 'Switching'
    (var_cost, var_error) = plotCostError(numClass,idClass,'var_2')
    
    time1 = time.time()
    print time1-time0
    print '\nBias Analysis'
    #Bias Plots
    numClass=3
    idClass = biasClass(numClass)
    #(bias_Tcost, bias_Terror) = plotTheoreticalCostError(numClass,idClass, 'bias_3_pen3pen2_varyGround.pdf')
    #print 'Switching'
    (bias_cost, bias_error) = plotCostError(numClass,idClass, 'bias_3')
    
    time2 = time.time()
    print time2-time1
    print '\nSingle Analysis'
    numClass=1
    idClass = singleClass()
    #(sing_Tcost, sing_Terror) = plotTheoreticalCostError(numClass,idClass,'single_pen3pen2_varyGround.pdf')
    #print 'Switching'
    (sing_cost, sing_error) = plotCostError(numClass,idClass,'single')
    
    time3 = time.time()
    print time3-time2
    
    plt.figure()
    plt.plot(rand_cost, rand_error, '^-',randM_cost, randM_error, '^-', sing_cost, sing_error, 'o-', var_cost, var_error,'*-', bias_cost, bias_error,'s-')
    plt.ylabel('Distance Disagreement with Ground (%)')
    plt.xlabel('Fractional Cost Used (%)')
    plt.legend(('Random Algo','Random Median','Single Class','Variance Class_2','Bias Class_9'))
    plt.title(' Graph1- Random vs Median vs Single vs Variance vs Bias')
    #plt.savefig('Graph1')
    plt.show()
    plt.close()


# Plot Graph 2, varying Variance Classes
def main2():
    loadData('evaluation_data.csv')
    reportMedianCostError(penReport)
    
    time0 = time.time()
    print '\nVariance Analysis- Class 2'
    numClass=2
    idClass = varianceClass(numClass)
    (var_cost2, var_error2) = plotCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    
    time1 = time.time()
    print time1-time0
    print '\nVariance Analysis - Class 4'
    numClass=4
    idClass = varianceClass(numClass)
    (var_cost4, var_error4) = plotCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    
    time2 = time.time()
    print time2-time1
    print '\nVariance Analysis - Class 8'
    numClass=8
    idClass = varianceClass(numClass)
    (var_cost8, var_error8) = plotCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    
    time3 = time.time()
    print time3-time2
    
    plt.figure()
    plt.plot(var_cost2, var_error2, 'o-', var_cost4, var_error4,'*-', var_cost8, var_error8,'s-')
    plt.ylabel('Distance Disagreement with Ground (%)')
    plt.xlabel('Fractional Cost Used (%)')
    plt.legend(('2 Class','4 Class','8 Class'))
    plt.title('Graph2 - Varying Variance Classes')
    plt.savefig('Graph2')
    plt.show()
    plt.close()
    
    
# Plot Graph 3, varying Bias Classes
def main3():
    loadData('evaluation_data.csv')
    reportMedianCostError(penReport)
    
    time0 = time.time()
    print '\Bias Analysis- Class 3'
    numClass=3
    idClass = biasClass(numClass)
    (var_cost2, var_error2) = plotCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    
    time1 = time.time()
    print time1-time0
    print '\nBias Analysis - Class 6'
    numClass=6
    idClass = biasClass(numClass)
    (var_cost4, var_error4) = plotCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    
    time2 = time.time()
    print time2-time1
    print '\nBias Analysis - Class 9'
    numClass=9
    idClass = biasClass(numClass)
    (var_cost8, var_error8) = plotCostError(numClass,idClass,'var_9_pen3pen2_varyGround.pdf')
    
    time3 = time.time()
    print time3-time2
    
    plt.figure()
    plt.plot(var_cost2, var_error2, 'o-', var_cost4, var_error4,'*-', var_cost8, var_error8,'s-')
    plt.ylabel('Distance Disagreement with Ground (%)')
    plt.xlabel('Fractional Cost Used (%)')
    plt.legend(('3 Class','6 Class','9 Class'))
    plt.title('Graph 3 - Varying Bias Classes')
    plt.savefig('Graph3')
    plt.show()
    plt.close()

main1()
#main2()
#main3()