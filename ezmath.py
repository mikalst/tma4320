import numpy 
import math

def bisectionMethod(mathFunction, start, end, tolerance):
    print("<--- Bisection Method --->")
    if (mathFunction(start)*mathFunction(end)<0):
        countIter = 0
        while (end-start)/2 > tolerance:
            middle = (start+end)/2
            if (mathFunction(middle)==0):
                print("iter {}: x = {}".format(countIter,(start+end)/2))
                return middle
            elif (mathFunction(start)*mathFunction(middle)<0):
                end = middle
            else:
                start = middle
            print("iter {}: x = {}".format(countIter,(start+end)/2))
            countIter += 1
        print("iter {}: x = {}".format(countIter,(start+end)/2))
        return (start+end)/2
    else:
        print("Your interval does not appear to contain a root, please choose a smaller one")
        return False

def fixedPointIteration(mathFunction, start, tolerance, solution=0, method=2, stepSize = 0.001):
    a = solution
    if (method == 0):
        print("<--- fpi DUMB --->")
        def derivatedMathFunction(x): return (mathFunction(x+stepSize)-mathFunction(x))/stepSize
    elif (method == 1):
        print("<--- fpi Richardsson's --->")
        def derivatedMathFunction(x): return (mathFunction(x+stepSize)-mathFunction(x))/stepSize
   elif (method == 2):
        print("<--- fpi Newton's --->")
        def derivatedMathFunction(x): return (mathFunction(x+stepSize)-mathFunction(x))/stepSize
    elif (method == 3):
        print("<--- fpi Chebyshev's --->")
        def derivatedMathFunction(x): return (mathFunction(x+stepSize)-mathFunction(x))/stepSize
    countIter = 0




def matrixDeterminant(x):
    matrixA = [[1,4,7,x],[4,5,x,6],[7,x,8,9],[x,10,11,12]]
    myDet = np.linalg.det(matrixA)-1000
    return myDet

# <--- Math Functions --->
def a(x):
    return 2*pow(x,3)-6*x-1

def b(x):
    return math.exp(x-2)+pow(x,3)-x

def c(x):
    return 1+5*x-6*pow(x,3)-math.exp(2*x)

def harmonic(x):
    return 1/x

# <--- End of math --->

def main():
    bisectionMethod(matrixDeterminant,0, 20, 0.001)
    return 0
    
if __name__ == "__main__":
    main()
