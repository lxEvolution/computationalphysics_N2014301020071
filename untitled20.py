
import matplotlib.pylab as plt
import numpy as np

class Ising():
   def __init__(self, H=0, T=2.5, length=20):
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
                self.temp=array([_system[self.temp1[0]][self.temp1[1]],_system[self.temp2[0]][self.temp2[1]],_system[self.temp3[0]][self.temp3[1]],_system[self.temp4[0]][self.temp4[1]]])
                self.delta_e = sum(self.temp)*_system[i][j]*2.+2.*_system[i][j]*_H
                if self.delta_e <=0:
                    _system[i][j]=-_system[i][j]
                else :
                    self.temp=exp(-self.delta_e/_T)
                    self.temp1=random.rand()
                    if self.temp1<=self.temp:
                        _system[i][j]=-_system[i][j]
    def energy_ave(self,_system,_H=0.):   
        self.energy=0.
        for i in range(self.length):
            for j in range(self.length):
                self.temp1,self.temp2,self.temp3,self.temp4=self.surrounding([i,j])
                self.temp=array([_system[self.temp1[0]][self.temp1[1]],_system[self.temp2[0]][self.temp2[1]],_system[self.temp3[0]][self.temp3[1]],_system[self.temp4[0]][self.temp4[1]]])
                self.energy=self.energy+(-sum(self.temp)*_system[i][j])/2.- _system[i][j]*_H
        return self.energy/self.length**2
    def corrrelation(self):   
        self.i=range(1,self.length/2+1)
        self.fi=[0.]*len(self.i)
        self.n=2000
        for j in range(self.n):
            self.MCstep(self.system,0.,self.T)
            if j >20 :
                for i in self.i:
                    for k in range(self.length):
                        for r in range(self.length/2):
                            self.fi[i-1]=self.fi[i-1]+self.system[k][r]*self.system[k][r+i]
        for i in range(len(self.fi)):
            self.fi[i]=self.fi[i]/((self.n-21)*(self.length**2/2))
i_dis=range(1,11)
fi_dis=[]
for i in [1.0,1.8,2.25,3.,6.]:
    cmp=ISING(0.,i,20)
    cmp.corrrelation()
    fi_dis.append(cmp.fi)
type=['r','m','g','y','c']
plt.figure(figsize=(6,6))
ax=plt.subplot(111)
print len(fi_dis)
for i in range(5):
    ax.plot(i_dis,fi_dis[i],'o--'+type[i],markersize=5)
ax.set_xlabel("Distance",fontsize=15)
ax.set_ylabel("Correlation",fontsize=15)
ax.set_title("Correlation function",fontsize=18)
ax.set_xlim(0,11)
ax.set_ylim(-0.1,1.1)
plt.show()