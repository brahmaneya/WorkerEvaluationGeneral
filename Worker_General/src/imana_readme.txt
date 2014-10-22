Paper: Learning from Crowds in the Presence of Schools of Thought, ACM SIGKDD 2012
Author: Yuandong Tian (yuandong@andrew.cmu.edu) and Jun Zhu (dcszj@mail.tsinghua.edu.cn)

The dataset is in dataset.txt. It is a 402-by-64 matrix. Each row is a response for a unique worker. You can use Matlab to read it. 

The meaning of each dimension of the response is as follows:

• Dimension 1 − 12: The response of task sky.
• Dimension 13 − 24: The response of task building.
• Dimension 25 − 36: The response of task computer.
• Dimension 37 − 40: The responses of object counting missions. 
• Dimension 41 − 52: The responses of task beauty1.
• Dimension 53 − 64: The responses of task beauty2.

For each task except for the counting task, the ordering of images in Fig. 1-5 in *images.pdf* (from left to right then from top to bottom) corresponds to the 1st-12th bits of 12 dimensional responses. 

There is no strict ground truth for these tasks. Plausible ground truths judged by authors are (1 = "is")

Task sky:        1 0 0 1 1 0 1 1 0 0 0 1
Task building:   0 1 0 0 0 0 1 0 1 1 1 0
Task computer:   0 1 1 0 1 1 0 0 0 1 0 0
Task beauty1:    0 0 1 1 1 1 1 0 0 0 0 1
Task beauty2:    1 1 1 1 1 0 0 1 0 0 0 0

The ground truth of the counting task is 65 5 8 27, as listed in the paper. 

For more details, please reference our paper. 

##################
Ground truths according to Manas
Task sky:	 1 1 0 0 1 0 1 1 0 0 0 1
Task building:	 0 1 0 0 0 0 1 1 1 1 1 0
Task computer:	 0 1 1 0 1 1 1 1 1 0 0 1
Task beauty1:	 



