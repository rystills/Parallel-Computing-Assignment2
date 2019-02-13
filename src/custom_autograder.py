from __future__ import print_function
import os, cmd, sys

numTests = 1

def main():
    ### Usage: python custom_autograder.py #ranks ###
    
    os.system("mpicc \-o cla.out cla.c")
    #this is the path to the folder containing all of the provided test cases and correct output
    path = "/home/parallel/spring-2019/stillr/project2/test_cases"
    #run all the test cases, sticking output in O0-O8
    for i in range(numTests):
        os.system('mpirun -np {0} ./cla.out {1}/t{2}.txt {3}/O{4}.txt'.format(sys.argv[1], path,i,path,i))
    #read all the test results to two lists
    userOutput = []
    expectedOutput = []
    for i in range(numTests):
        with open('{0}/O{1}.txt'.format(path,i),"r") as f:
            userOutput.append(f.readlines())
        with open('{0}/A{1}.txt'.format(path,i),"r") as f:
            expectedOutput.append(f.readlines())
    #compare the test results
    print("Test1: Pass")
    testsFailed = 0
    for i in range(numTests):
        if (userOutput[i] == expectedOutput[i]):
            print("Test{0}: Pass".format(i+2))
        else:
            print("Test{0}: Fail\n{1}{2}".format(i+2,userOutput[i][0],expectedOutput[i][0]),end='')
            testsFailed +=1
    print("Total: {0}/{1} Tests Passed".format(numTests+1-testsFailed,numTests+1))
        
if __name__ == "__main__":
    main()