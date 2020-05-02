import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from numba import njit


def readfile(fname):
    with open(fname, 'r') as f:
        print('Working on '+fname+'\nAnalyzing mesh...')
        Nnode = np.nan
        Nelem = np.nan
        Ncomm = 0
        for line in f:
            if line.startswith('#') or line.startswith(':'):
                Ncomm += 1
            if ':NodeCount' in line:
                Nnode = np.int64(line.split()[-1])
            if ':ElementCount' in line:
                Nelem = np.int64(line.split()[-1])
            if np.isnan(Nnode) is False and np.isnan(Nelem) is False:
                break
    node = np.loadtxt(fname,comments=[':','#'],max_rows=Nnode)
    elem = np.loadtxt(fname,comments=[':','#'],skiprows=Nnode+Ncomm,dtype=np.int64)
    print('NodeCount %d' % len(node))
    print('ElementCount %d' % len(elem)) 
    return node, elem


@njit
def erosion(node, elem, rdm, risk, current_iteration_num, coupling_period):
    if np.mod(current_iteration_num, coupling_period) == 0:
        num_of_element_of_interest = 0
        element_of_interest = np.zeros(len(elem),dtype=np.int64)
        for i_element in np.random.permutation(np.arange(len(elem))):
            i1 = elem[i_element,0] - 1
            i2 = elem[i_element,1] - 1
            i3 = elem[i_element,2] - 1
            if (node[i1,2] == 0 or node[i2,2] == 0 or node[i3,2] == 0) and (node[i1,2] == 1 or node[i2,2] == 1 or node[i3,2] == 1):
                element_of_interest[num_of_element_of_interest] = i_element
                num_of_element_of_interest += 1
        if num_of_element_of_interest != 0:
            for i_element in element_of_interest[:num_of_element_of_interest]:
                i1 = elem[i_element,0] - 1
                i2 = elem[i_element,1] - 1
                i3 = elem[i_element,2] - 1
                bernoulli = 0
                if rdm[i_element] < risk[i_element]:
                    bernoulli = 1
                for i in [i1, i2, i3]:
                    if node[i,2] == 1 and bernoulli == 1:
                        node[i,2] = 0
    return node[:,2]


def evolution(fname, Niterations, Nrealizations, coupling_period):
    if Niterations < 1:
        Niterations = 1
    if Nrealizations < 1:
        Nrealizations = 1
    node, elem = readfile(fname)
    z = np.copy(node[:, 2])
    ave = np.zeros(len(node))
    for realization in range(Nrealizations):
        print('Realization '+str(realization+1))
        node[:, 2] = np.copy(z)
        for current_iteration_num in range(Niterations):
            rdm = np.random.random(len(elem))
            r = np.ones(len(elem))/2
            risk = 0.5*(1 + erf(4*(r - 1)/np.sqrt(r**2 + 1)/np.sqrt(2)))
            node[:,2] = erosion(node, elem, rdm, risk, current_iteration_num, coupling_period)
            sys.stdout.write('\rIteration '+str(current_iteration_num))
            sys.stdout.flush()
        print()
        ave += node[:,2]
    ave /= Nrealizations
    print()
    return node, ave


def main(Niterations, Nrealizations, onlyIC=False):
    mycmap = np.ones((21,4))
    mycmap[:,0] = np.linspace(to_rgba('navy')[0], 244/256, 21)
    mycmap[:,1] = np.linspace(to_rgba('navy')[1], 242/256, 21)
    mycmap[:,2] = np.linspace(to_rgba('navy')[2], 205/256, 21)
    mycmap = ListedColormap(mycmap)
    
    fList = ['BOTTOM25.t3s', 'BOTTOM50.t3s', 'BOTTOM100.t3s']
    if onlyIC is False:
        for couplingList in [[1, 1, 1], [1, 2, 4]]:
            plt.figure(figsize=(9,12))
            for i, fname in enumerate(fList):
                node, ave = evolution(fname, Niterations, Nrealizations, couplingList[i])
                x, y, z = node[:,0], node[:,1], node[:,2]
                for j, zplot, ttl in zip([1, 2], [z, ave], ['Single Realization', 'Ensemble Average of '+str(Nrealizations)+' Realizations']):
                    plt.subplot(int(str(len(fList))+'2'+str(j+i*2)))
                    plt.tripcolor(x,y,zplot,cmap=mycmap)
                    plt.colorbar(fraction=0.03,ticks=[0,1])
                    plt.xlim(np.min(x),np.max(x))
                    plt.ylim(np.min(y),np.max(y))
                    plt.xticks(np.arange(np.min(x),np.max(x)+1,10))
                    plt.yticks(np.arange(np.min(y),np.max(y)+1,10))
                    plt.title(ttl+'\nMesh size: '+str(float(fname.split('BOTTOM')[1].split('.t3s')[0])/100)+'    CP: '+str(couplingList[i]))
                    plt.tight_layout()
            plt.savefig(''.join(str(couplingList)[1:-1].split(', '))+'.png',dpi=300)
            plt.close()
    else:
        plt.figure(figsize=(9,12))
        for i, fname in enumerate(fList):
            node, elem = readfile(fname)
            print()
            x, y, z = node[:,0], node[:,1], node[:,2]
            
            plt.subplot(int(str(len(fList))+'2'+str(2+i*2)))
            plt.tripcolor(x,y,z,cmap=mycmap)
            plt.colorbar(fraction=0.03,ticks=[0,1])
            plt.xlim(np.min(x),np.max(x))
            plt.ylim(np.min(y),np.max(y))
            plt.xticks(np.arange(np.min(x),np.max(x)+1,10))
            plt.yticks(np.arange(np.min(y),np.max(y)+1,10))
            plt.title('Initial Condition')
            plt.tight_layout()

            plt.subplot(int(str(len(fList))+'2'+str(1+i*2)))
            patches = []
            for i_element in range(len(elem)):
                tri = np.zeros((3,2))
                for j in range(3):
                    tri[j] = node[elem[i_element,j]-1][:2]
                patches.append(Polygon(tri, closed=True))
            p = PatchCollection(patches, edgecolors='k', linewidths=.1, facecolors='none')
            plt.gca().add_collection(p)
            plt.xlim(np.min(x),np.max(x))
            plt.ylim(np.min(y),np.max(y))
            plt.xticks(np.arange(np.min(x),np.max(x)+1,10))
            plt.yticks(np.arange(np.min(y),np.max(y)+1,10))
            plt.title('Mesh size: '+str(float(fname.split('BOTTOM')[1].split('.t3s')[0])/100))
            plt.tight_layout()
##            inset_axes(plt.gca(), height='35%', width='35%', loc=3)
##            patches = []
##            for i_element in range(len(elem)):
##                tri = np.zeros((3,2))
##                for j in range(3):
##                    tri[j] = node[elem[i_element,j]-1][:2]
##                patches.append(Polygon(tri, closed=True))
##            p = PatchCollection(patches, edgecolors='k', linewidths=.1, facecolors='none')
##            plt.gca().add_collection(p)
##            plt.xlim(22,28)
##            plt.ylim(22,28)
##            plt.xticks([])
##            plt.yticks([])
        plt.savefig('ic.png',dpi=300)
        plt.show()

  
if __name__ == '__main__':
    main(Niterations=101, Nrealizations=100, onlyIC=False)

