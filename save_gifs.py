'''
Created on Jun 25, 2021

@author: jiadongc
'''
import imageio
writer = []
filenames = []
for i in range(0,8):
    filenames.append(str(i))
print(filenames)
for filename in filenames:
    image = imageio.imread("C://Users//jiadongc//classes//myResearch//prelim//gif_figs2//"+filename+".png")
    writer.append(image)
exportname = "output11.gif"
kargs = { 'duration': 0.8 }
imageio.mimsave(exportname, writer, 'GIF', **kargs)