import os, cmd

def main():
    os.system("gcc cla.c")
    #this is the path to the folder containing all of the provided test cases and correct output
    path = "C:\\Users\\Ryan\\Desktop\\school files\\RPI\\spring 2019\\Parallel Computing\\programming\\Parallel-Computing-Assignment1\\test_cases"
    #run all the test cases, sticking output in O0-O8
    for i in range(8):
        os.system('a.exe < "{0}\\t{1}.txt" > "{2}\\O{3}.txt"'.format(path,i,path,i))
    #read all the test results to two lists
    userOutput = []
    expectedOutput = []
    for i in range(8):
        with open('{0}\\O{1}.txt'.format(path,i),"r") as f:
            userOutput.append(f.readlines())
        with open('{0}\\A{1}.txt'.format(path,i),"r") as f:
            expectedOutput.append(f.readlines())
    #compare the test results
    print("Test1: Pass")
    testsFailed = 0
    for i in range(8):
        if (userOutput[i] == expectedOutput[i]):
            print("Test{0}: Pass".format(i+2))
        else:
            print("Test{0}: Fail\n{1}{2}".format(i+2,userOutput[i][0],expectedOutput[i][0]),end='')
            testsFailed +=1
    print("Total: {0}/{1} Tests Passed".format(9-testsFailed,9))
        
if __name__ == "__main__":
    main()