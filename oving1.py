from bokeh.plotting import figure, output_file, show
import numpy as np

def bisectionMethod(mathFunction, start, end, tolerance):
    if (mathFunction(start)*mathFunction(end)<0):
        while (end-start)/2 > tolerance:
            middle = (start+end)/2
            if (mathFunction(middle)==0):
                return middle
            elif (mathFunction(start)*mathFunction(middle)<0):
                end = middle
            else:
                start = middle
        print("The approximated root is x:",(start+end)/2)
        return (start+end)/2
    else:
        print("Your interval does not appear to contain a root, please choose a smaller one")
        return False

def fixedPointIteration(mathFunction, start, tolerance):
    return 0

def matrixDeterminant(x):
    matrixA = [[1,4,7,x],[4,5,x,6],[7,x,8,9],[x,10,11,12]]
    myDet = np.linalg.det(matrixA)-1000
    return myDet

def harmonic(x):
    return 1/x

def plotFunction(mathFunction, x_start, x_end, dataPoints):
    stepSize = (x_end - x_start)/dataPoints
    x = np.arange(x_start, x_end, stepSize)
    y = []
    for i in x:
        y.append(mathFunction(i))
    y_start = min(y)
    y_end = max(y)
    output_file("output.html") 
    p = figure(x_range=(x_start,x_end), y_range=(y_start, y_end))
    p.line(x,y)
    show(p)

def main():
    plotFunction(matrixDeterminant,-20, 10, 1000)
    print(bisectionMethod(matrixDeterminant,0, 20, 0.00001))
    return 0
    
if __name__ == "__main__":
    main()
