from bokeh.plotting import figure, output_file, show
import numpy as np

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
   return 0
    
if __name__ == "__main__":
    main()
