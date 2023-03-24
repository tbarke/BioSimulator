import math
import log
from matplotlib.patches import Rectangle
l = log.log()

def getColorCode(color):
    color_code =  [(0, 0, 0), (255, 255, 255), (255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 255, 0), (0, 255, 255), (255, 0, 255), (192, 192, 192), (128, 128, 128), (128, 0, 0), (128, 128, 0), (0, 128, 0), (128, 0, 128), (0, 128, 128), (0, 0, 128), (244, 230, 30), (32, 147, 140), (68, 2, 86)]
    colors = ['Black', 'White', 'Red', 'Lime', 'Blue', 'Yellow', 'Cyan', 'Magenta', 'Silver', 'Gray','Maroon', 'Olive', 'Green', 'Purple', 'Teal', 'Navy', 'y1', 'g1', 'p1']
    for i in range(len(colors)):
        if colors[i] == color:
            return color_code[i]
    return color

def findcolor(max, min, color_scheme, num, absolute = True):
    portion = (num - min) / (max - min)
    index = portion * (len(color_scheme) - 1)
    index_int_bottom = math.floor(index)
    if index_int_bottom == len(color_scheme) - 1:
        col = getColorCode(color_scheme[index_int_bottom])
        if absolute:
            return col[0]/255, col[1]/255, col[2]/255
        else:
            return col[0], col[1], col[2]
    index_int_top = index_int_bottom + 1
    index = index - index_int_bottom
    smallColor = getColorCode(color_scheme[index_int_bottom])
    largetColor = getColorCode(color_scheme[index_int_top])
    new_x = 0
    new_y = 0
    new_z = 0
    if smallColor[0] > largetColor[0]:
        new_x = (largetColor[0] + (smallColor[0] - largetColor[0]) * (1-index)) / 255
    else:
        new_x = (smallColor[0] + (math.fabs(largetColor[0] - smallColor[0]) * index)) / 255

    if smallColor[1] > largetColor[1]:
        new_y = (largetColor[1] + (smallColor[1] - largetColor[1]) * (1-index)) / 255
    else:
         new_y = (smallColor[1] + (math.fabs(largetColor[1] - smallColor[1]) * index)) / 255
    if smallColor[2] > largetColor[2]:
        new_z = (largetColor[2] + (smallColor[2] - largetColor[2]) * (1-index)) / 255
    else:
        new_z = (smallColor[2] + (math.fabs(largetColor[2] - smallColor[2]) * index)) / 255
    if absolute:
        return [new_x, new_y, new_z]
    else:
        return [new_x*255, new_y*255, new_z*255]

def findColor4D(colors, x1, y1):

    def dist(x1, y1, x2, y2):
        return math.sqrt(math.fabs(math.pow(y2-y1, 2)) +math.fabs(math.pow(x2-x1, 2)))

    #x2 = x1
    xplane1_dist = dist(x1, y1, x1, 0)
    #y2 = y1
    yplane1_dist = dist(x1, y1, 0, y1)
    #x2 = x1
    xplane2_dist = 1 - xplane1_dist
    #y2 = y1
    yplane2_dist = 1 - yplane1_dist

    red_infl = (1 - xplane1_dist) * (1-yplane1_dist)
    blue_infl = (1 - xplane2_dist) * (1-yplane1_dist)
    green_infl =(1 - xplane1_dist) * (1-yplane2_dist)
    yellow_infl = (1 - xplane2_dist) * (1-yplane2_dist)


    def multColor(infl, color):
        new_color = [0,0,0]
        new_color[0] = color[0] * infl
        new_color[1] = color[1] * infl
        new_color[2] = color[2] * infl
        return new_color

    new_red  = multColor(red_infl, colors[0])
    new_blue = multColor(blue_infl, colors[1])
    green_blue = multColor(green_infl, colors[2])
    yellow_blue = multColor(yellow_infl, colors[3])

    def color_combine(colors):
        length = len(colors)
        new_color = [0,0,0]
        for i in range(length):
            new_color[0] += colors[i][0]
            new_color[1] += colors[i][1]
            new_color[2] += colors[i][2]

        return new_color

    return color_combine([new_red, new_blue, green_blue, yellow_blue])

def drawRectangles(boxes, colors_4D):
    index = 1/boxes
    rects = []
    for i in range(boxes):
        i_index = i * index
        for j in range(boxes):
            j_index = j*index
            rects.append(Rectangle((i_index, j_index), index, index, facecolor=findColor4D(colors_4D, i_index, j_index )))
    return rects

def d3color(emp_red, emp_blue):
    color_Scheme1 = ['Black', 'Red']
    color_Scheme2 = ['Black', 'Blue']
    color_red = findcolor(50, 0,color_Scheme1, emp_red, absolute=False)
    color_blue = findcolor(50,0,color_Scheme2, emp_blue, absolute=False)
    color_Scheme3 = [color_red, color_blue]
    dist = 0
    if color_red[0] + color_blue[2] >0:
        dist = color_blue[2] / (color_red[0] + color_blue[2])
    color_purp = findcolor(1,0,color_Scheme3, dist)
    return color_purp