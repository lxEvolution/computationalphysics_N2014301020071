
import matplotlib.pylab as plt
import numpy as np
class Ising():
    def __init__(self, H=0, T=0.5, length=20):
        self.length = length
        self.H =H
        self.T = T
        self.system =np.ones((self.length,self.length))
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
    def MCstep(self,_system, _H=0., _T=0.5):  
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
    def first_order(self): 
        self.H=list(np.linspace(10,-10,60))
        self.M=[]
        print (self.H)
        for i in self.H:
            self.temp_time=20
            self.temp_h=0
            for j in range(self.temp_time):
                self.MCstep(self.system, i, 2.1)
                self.temp_h=self.temp_h+np.mean(self.system)
            self.mfile=open(r"d:\somewhat%.2f.txt"%i,'w')
            for i in range(self.length):
                for j in range(self.length):
                    print ( self.mfile, self.system[i][j])
            self.mfile.close()
            self.M.append(self.temp_h/self.temp_time)
        plt.figure(figsize=(5,5))
        ax=plt.subplot(111)
        ax.plot(self.H[0:60],self.M[0:60],'o-b',markersize=3)
        ax.plot(self.H[60:119],self.M[60:119],'o-r',markersize=3)
        ax.set_xlabel("Magnetic Field",fontsize=15)
        ax.set_ylabel("Magnetization",fontsize=15)
        ax.set_title("First-order phase transition",fontsize=18)
        ax.set_xlim(-10.2,10.2)
        ax.set_ylim(-1.1,1.1)
        plt.show()

cmp=Ising()
cmp.first_order()
