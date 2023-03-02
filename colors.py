import math
import log
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