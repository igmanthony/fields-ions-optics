import numpy as np
# import matplotlib.pyplot as plt

# start70 = np.array([
#     [  0,-70,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
# ])

# start100 = np.array([
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,  0,  0],
#     [  0,  0,  0,100,  0],
# ])

# # the value of each cell is its number of neighbors
# boundary_conditions = np.array([
#     [  2,  3,  3,  3,  2],
#     [  3,  4,  4,  4,  3],
#     [  3,  4,  4,  4,  3],
#     [  3,  4,  4,  4,  3],
#     [  2,  3,  3,  3,  2],
# ])

# # this needs to be 2x2 bigger to allow pasting the start
# # four times - left, right, up, and down
# temporary_canvas = np.zeros((7, 7))

# for i in range(10_000):
#     temporary_canvas[:-2,1:-1] += start70 # paste left
#     temporary_canvas[2:,1:-1] += start70 # paste right
#     temporary_canvas[1:-1,:-2] += start70 # paste up
#     temporary_canvas[1:-1,2:] += start70 # paste down
#     temporary_canvas[1:-1,1:-1] /= boundary_conditions
#     start70 = np.copy(temporary_canvas[1:-1,1:-1])
#     start70[0,1] = -70
#     start70[4, 3] = 0
#     temporary_canvas.fill(0)

# for i in range(10_000):
#     temporary_canvas[:-2,1:-1] += start100 # paste left
#     temporary_canvas[2:,1:-1] += start100 # paste right
#     temporary_canvas[1:-1,:-2] += start100 # paste up
#     temporary_canvas[1:-1,2:] += start100 # paste down
#     temporary_canvas[1:-1,1:-1] /= boundary_conditions
#     start100 = np.copy(temporary_canvas[1:-1,1:-1])
#     start100[0,1] = 0
#     start100[4, 3] = 100
#     temporary_canvas.fill(0)

# print(start70)
# print(start100)
# plt.imshow(start70, cmap="turbo")
# plt.show()
# plt.imshow(start100, cmap="turbo")
# plt.show()


data = np.array([[ 13.23529412,   0.0       ,   25.0      ,    38.23529412,  44.11764706],
                 [ 26.47058824,  26.47058824,  36.76470588,  45.58823529  ,50.0        ],
                 [ 39.70588235,  42.64705882,  50.0       ,   57.35294118 , 60.29411765],
                 [ 50.        ,  54.41176471,  63.23529412,  73.52941176  ,73.52941176],
                 [ 55.88235294,  61.76470588,  75.        , 100.0         , 86.76470588]] )

print(data/2)