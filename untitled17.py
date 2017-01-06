
import matplotlib.pylab as plt
import numpy as np
class Ising(object):
    def __init__(self, H, T, length=20):
        self.length =length
        self.H = H
        self.T = T
        self.system = []
        for i in range(self.length):
            self.system.append([1]*self.length)
    def surrounding(self,_posi):  
        if 1 <= _posi[0] <= self.length-2 and 1<=_posi[1]<=self.length-2:
            return [_posi[0]-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],_posi[1]+1]
        if _posi[0]==0 and (1<=_posi[1]<=self.length-2):
            return [self.length-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],_posi[1]+1]
        if _posi[0]==self.length-1 and  (1<=_posi[1]<=self.length-2):
            return [_posi[0]-1,_posi[1]],[0,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],_posi[1]+1]
        if 1<=_posi[0]<=self.length-2 and _posi[1]==0:
            return [_posi[0]-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],self.length-1],[_posi[0],_posi[1]+1]
        if 1<=_posi[0]<=self.length-2 and _posi[1]==self.length-1:
            return [_posi[0]-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],0]
        if _posi==[0,0]:
            return [self.length-1,0],[1,0],[0,self.length-1],[0,1]
        if _posi==[self.length-1,0]:
            return [self.length-2,0],[0,0],[self.length-1,self.length-1],[self.length-1,1]
        if _posi==[0,self.length-1]:
            return [self.length-1,self.length-1],[1,self.length-1],[0,self.length-2],[0,0]
        if _posi==[self.length-1,self.length-1]:
            return [self.length-2,self.length-1],[0,self.length-1],[self.length-1,self.length-2],[self.length-1,0]           
    def MCstep(self,_system, _H=0., _T=0.5):  # sweep all spins -- Monte Carlo step
        for i in range(self.length):
            for j in range(self.length):
                self.temp1,self.temp2,self.temp3,self.temp4=self.surrounding([i,j])
                self.temp=np.array([_system[self.temp1[0]][self.temp1[1]],_system[self.temp2[0]][self.temp2[1]],_system[self.temp3[0]][self.temp3[1]],_system[self.temp4[0]][self.temp4[1]]])
                self.delta_e = sum(self.temp)*_system[i][j]*2.+2.*_system[i][j]*_H
                if self.delta_e <=0:
                    _system[i][j]=-_system[i][j]
                else :
                    self.temp=np.exp(-self.delta_e/_T)
                    self.temp1=np.random.rand()
                    if self.temp1<=self.temp:
                        _system[i][j]=-_system[i][j]
    def energy_ave(self,_system,_H=0.):
        self.energy=0.
        for i in range(self.length):
            for j in range(self.length):
                self.temp1,self.temp2,self.temp3,self.temp4=self.surrounding([i,j])
                self.temp=np.array([_system[self.temp1[0]][self.temp1[1]],_system[self.temp2[0]][self.temp2[1]],_system[self.temp3[0]][self.temp3[1]],_system[self.temp4[0]][self.temp4[1]]])
                self.energy=self.energy+(-sum(self.temp)*_system[i][j])/2.- _system[i][j]*_H
        return self.energy/self.length**2
    def plot_system(self,_sys):
        pass
    def plot_spinvstime(self,_times):
        self.time=np.linspace(0,_times,_times+1)
        self.mag=[np.mean(self.system)]
        self.ene=[self.energy_ave(self.system,self.H)]
        for i in range(_times):
            self.MCstep(self.system, self.H, self.T)
            self.mag.append(np.mean(self.system))
            self.ene.append(self.energy_ave(self.system,self.H))
        np.fig=plt.figure(figsize=(13,12))
        ax1=plt.axes([0.1,0.2,0.35,0.7])
        ax2=plt.axes([0.55,0.2,0.35,0.7])
        ax1.plot(self.time,self.mag,'-b')
        ax2.plot(self.time,self.ene,'-r')
        ax1.set_xlabel("Time / M-C step",fontsize=15)
        ax1.set_ylabel("Magnetization",fontsize=15)
        ax2.set_xlabel("Time / M-C step",fontsize=15)
        ax2.set_ylabel("Energy",fontsize=15)
        ax1.set_title("Fluctuation: Magnetization vs. time",fontsize=18)
        ax2.set_title("Fluctuation: Energy vs. time",fontsize=18)
        ax1.set_ylim(-1.5,1.5)
        ax2.set_ylim(-3,3)
        plt.show()
a=Ising(0,0.5)
b=Ising(0,1.5)
a=Ising(0,2.0)
b=Ising(0,2.25)
a=Ising(0,3.0)
b=Ising(0,4.0)
a.plot_spinvstime(200)
