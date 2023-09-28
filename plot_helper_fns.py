# libraries
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


# drawline is used to draw a line on the current plot
def draw_line(p1, p2, color, ax):
    # an axis can be provided. If it isn't, it will be created
    if ax==None:
        ax = plt.gca()

    l = mlines.Line2D([p1[0],p2[0]], [p1[1],p2[1]], color=color)
    ax.add_line(l)
    return l

# drawline is used to draw a line on the current plot
def draw_lines(poly, color, ax=None):
    if ax==None:
        ax = plt.gca()
        
    for i in range(0,len(poly)-1):
        draw_line(poly[i], poly[i+1], color, ax)