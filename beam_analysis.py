import numpy as np
from numpy import sin, cos, tan, pi
from sympy import Symbol, symbols, solveset, linsolve
from matplotlib import pyplot as plt

class beam():
    unit='m' # kN,N-m,N/m
    length=1
    ds=0.01 #m
    LoadTyp={'point','moment','GDL','UDL','rTDLrmax','rTDLlmax'}
    # Load specification load,dirn,type,pos  dirn:0 to 360 degree, {cw,ccw} pos:x,[xs,xe],
    constraintDef={'hinge':{'Rx':True,'Ry':True,'Mz':False},'fixed':{'Rx':True,'Ry':True,'Mz':True},'roller':{'Rx':False,'Ry':True,'Mz':False},'free':{'Rx':False,'Ry':False,'Mz':False}}
    dirn={'up':0,'down':1,'ccw':0,'cw':1}
    def __init__(self,length=1,units='m',ds=0.01):
        self.unit=units
        self.length=length
        self.ds=ds
        self.sections={} #dict({0:{'load':None,'dirn':None,'typ':None,'pos':None}})
        self.secCount=0
        self.constriants=None
        
    def add_constraints(self,s,e):
        self.constraints=[s,e]
        
    def add_section(self,s,e,load=None,dirn=None,typ=None,pos=None, label=None):
        self.sections.update({self.secCount:{'load':load,'dirn':dirn,'typ':typ,'pos':pos,'start':s,'endd':e,'label':label}})
        self.secCount=len(self.sections)
        
    def calc_reactions(self):
        #inds=[val for i,val in enumerate(self.constraintDef if val==self.constraints[0])
        #inde=[val for i,val in enumerate(self.constraintDef if val==self.constraints[1])

        ## Find which reactions exists
        se=[0,1] # start and end reaction notations
        a=[]
        for i in se:
            x=self.constraintDef[self.constraints[i]]
            #print(x)
            for j in x:
                if x[j]:
                    a.append(Symbol(j+str(i)))
        #print(a)
        # Append reactions in global force/mom vector
        Fx=[]
        Fy=[]
        Maz=[]
        for x in a:
            cha=None
            cha=[y for y in str(x)]
            if 'R' in cha:
                if 'x' in cha:
                    Fx.append([x,None,float(cha[-1])*self.length]) # (symbol,value,position)
                elif 'y' in cha:
                    Fy.append([x,None,float(cha[-1])*self.length])
            elif 'M' in cha:
                Maz.append([x,None,float(cha[-1])*self.length])
        print(Fx)
        print(Fy)
        print(Maz)

        # Append the loads defined in the sections in the global force/mom vector
        for i in self.sections:
            if self.sections[i]['typ']=='point':
                #,'genericDistributed','uniformDistributed','moment','rTriangleRmax','rTraingleLmax'}
                if self.isanumber(self.sections[i]['dirn']):
                    ang=self.sections[i]['dirn']
                    if not (ang<0 & ang>360):
                        cfact=self.conformSignConvention(ang)
                        fx=cos(ang*pi/180)*cfact[0]*self.sections[i]['load']
                        fy=sin(ang*pi/180)*cfact[1]*self.sections[i]['load']
                        Fx.append([Symbol('Px'+str(i)),fx,self.sections[i]['pos']])
                        Fy.append([Symbol('Py'+str(i)),fy,self.sections[i]['pos']])
                    else:
                        print('Improper Angle for point load')
                else:
                    print('Point load should be provided with an angle of the vector as the direction.')
                    print('Point load is omitted.')
            elif self.sections[i]['typ']=='moment':
                if isinstance(self.sections[i]['dirn'],str):
                   msense=self.sections[i]['dirn']
                   if msense=='cw':
                        msign=float(-1)
                   elif msense=='ccw':
                        msign=1.0
                   else:
                        print('Unknown moment sense. Moment omitted.')
                        break;
                   Maz.append([Symbol('Mz'+str(i)),msign*self.sections[i]['load'],self.sections[i]['pos']])

                else:
                    print('Moment should be provided with cw or ccw as string for the direction.')
                    print('Moment load is omitted.')
                    
            elif self.sections[i]['typ']=='UDL':
                if self.isanumber(self.sections[i]['dirn']):
                    ang=self.sections[i]['dirn']
                    if (ang==90 or ang==270):
                        if ang==90:
                            Lsign=-1.0
                        elif ang==270:
                            Lsign=1.0
                        else:
                            print('Unknown angle for UDLoad vector. UDL omitted.')
                            break;
                        fx=0.0
                        fy=Lsign*self.sections[i]['load']*(self.sections[i]['pos'][1]-self.sections[i]['pos'][0])
                        loc=(self.sections[i]['pos'][1]+self.sections[i]['pos'][0])/2
                        Fx.append([Symbol('UDLx'+str(i)),fx,loc,self.sections[i]['load'],self.sections[i]['pos']])
                        Fy.append([Symbol('UDLy'+str(i)),fy,loc,Lsign*self.sections[i]['load'],self.sections[i]['pos']])
                    else:
                        print('Improper Angle for UDL')
                else:
                    print('UDL should be provided with an angle of the load vector (90 or 270) as the direction.')
                    print('UDL is omitted.')
                    
            elif self.sections[i]['typ']=='rTDLrmax': # load in N/m, symbol rTDLa, where a denotes rmax.
                if self.isanumber(self.sections[i]['dirn']):
                    ang=self.sections[i]['dirn']
                    if (ang==90 or ang==270):
                        if ang==90:
                            Lsign=-1.0
                        elif ang==270:
                            Lsign=1.0
                        else:
                            print('Unknown angle for UDLoad vector. UDL omitted.')
                            break;
                        fx=0
                        fy=Lsign*0.5*self.sections[i]['load']*(self.sections[i]['pos'][1]-self.sections[i]['pos'][0])
                        loc=self.sections[i]['pos'][1]-(self.sections[i]['pos'][1]-self.sections[i]['pos'][0])*(1/3)
                        Fx.append([Symbol('rTDLax'+str(i)),fx,loc,self.sections[i]['load'],self.sections[i]['pos']]) # just for completeness
                        Fy.append([Symbol('rTDLay'+str(i)),fy,loc,Lsign*self.sections[i]['load'],self.sections[i]['pos']])
                    else:
                        print('Improper Angle for TDL')
                else:
                    print('TDL should be provided with an angle of the load vector (90 or 270) as the direction.')
                    print('TDL is omitted.')
                    
            elif self.sections[i]['typ']=='XYZ':
                if self.isanumber(self.sections[i]['dirn']):
                    ang=self.sections[i]['dirn']
                    if not (ang<0 & ang>360):
                        cfact=self.conformSignConvention(ang)
                        fx=cos(ang*pi/180)*cfact(0)*self.sections[i]['load']
                        fy=sin(ang*pi/180)*cfact(1)*self.sections[i]['load']
                        Fx.append([Symbol('Px'+str(i)),fx,self.sections[i]['pos']])
                        Fy.append([Symbol('Py'+str(i)),fy,self.sections[i]['pos']])
                    else:
                        print('Improper Angle for point load')
                else:
                    print('XYZ load should be provided with an angle of the vector as the direction.')
                    print('XYZ load is omitted.')
                    
            else:
                print('Unknown load type')

        print(Fx)
        print(Fy)
        print(Maz)

        # Sum Fx=0
        v=[]
        eqnFx=0
        for f in Fx:
            eqnFx += f[0]
            if not (f[1]==None):
                v.append((str(f[0]),f[1]))
        print('Sum Fx = 0')
        print(eqnFx)
        print(v)
        eqnFx_kn=eqnFx.subs(v)
        eqnFxsol=solveset(eqnFx_kn,Fx[0][0])
        # (x0, y0) = next(iter(sol)) # get an iterator based on object and pull one by one using next command
        if (type(eqnFxsol)=='FiniteSet'):
            print(str(Fx[0][0])+'= %d' % (eqnFxsol.args[0]))
            Fx[0][1]=eqnFxsol.args[0]
            
        # sum Fy=0
        vy=[]
        eqnFy=0
        for f in Fy:
            eqnFy += f[0]
            if not (f[1]==None):
                vy.append((str(f[0]),f[1]))

        print('Sum Fy = 0')
        print(eqnFy)
        print(vy)
        eqnFy_kn=eqnFy.subs(vy)
        #eqnFysol=solveset(eqnFy_kn)
        print(eqnFy_kn)
        
        # sum Mz=0
        vm=[]
        eqnM0=0
        for f in Fy:
            if not (f[1]==None):
                if not (str(f[0])=='Ry0'):
                    vm.append((str(f[0]),f[1]))
                    eqnM0 += f[1]*f[2]
            else:
                eqnM0 += f[0]*f[2]
                
        if not Maz==[]:
            for mz in Maz:
                    if not (mz[1]==None):
                        vm.append((str(mz[0]),mz[1]))
                    eqnM0 += mz[0]
        print('Sum Mz = 0')
        print(eqnM0)
        print(vm)
        eqnM0_kn=eqnM0.subs(vm) # substitute knonw moments
        print(eqnM0_kn)
        #eqnFysol=solveset(eqnFy_kn)
        
        # solve equations and find unknowns.
#       forces=linsolve([eqnFy_kn,eqnM0],(Fy[0][0],Fy[1][0]))
        unkn=[]
        for ff in Fy:
            if ff[1] is None:
                unkn.append(ff[0])
        for mm in Maz:
            if mm[1] is None:
                unkn.append(mm[0])
        print(unkn)
        
        forces=linsolve([eqnFy_kn,eqnM0_kn],unkn)
# TBD print forces thru loop check if forces exist based on the prioir list.
        #print(str(Fy[0][0])+','+str(Fy[1][0])+'= %f, %f' %(forces.args[0][0],forces.args[0][1]))
        print(forces)
        cnt=0
        for u in unkn:
            for k in Fy:
                if u==k[0]:
                    if k[1]==None: # yet another check
                       k[1]=forces.args[0][cnt]
                       # add some measure toavoid multiple variables being assigned.
                       # break the loop if assigning operation done.
            for kk in Maz:
                if u == kk[0]:
                    if kk[1]==None: # yet another check
                       kk[1]=forces.args[0][cnt]
            cnt+=1
            
        print(Fy)
        print(Maz)
        
    # calculate shear force
        dx=0.001
        LB=np.linspace(0,self.length,num=1000) #0:dx:self.length
        V=[]
        for xx in LB:
            if xx==0:
                V.append(Fy[0][1])
            #elif xx==LB[-1]:
            #    V.append(Fy[1][1])
            else:
                Vxx=0
                for ff in Fy:
                    if len(ff)==3: # point load
                        if (ff[2] < xx):
                            Vxx += ff[1]
                    elif len(ff)==5: # UDL,rTDL
                        if ff[4][0]<xx: # if starting point less han current section
                            if not (ff[4][1]<xx): #if ending point not less than the current section
                                if ('U' in str(ff[0])): # UDL N/m * m
                                    Vxx +=ff[3]*(xx-ff[4][0])
                                elif (('T' in str(ff[0])) & ('a' in str(ff[0]))): #rTDLrmax
                                    Vxx +=0.5*(xx-ff[4][0])*ff[3]*((xx-ff[4][0])/(ff[4][1]-ff[4][0]))
                            else: #if we are past the UDL or TDL use the equivalent load.
                                if ('U' in str(ff[0])): # UDL N/m * m
                                    Vxx +=ff[1]
                                elif (('T' in str(ff[0])) & ('a' in str(ff[0]))): #rTDLrmax
                                    Vxx +=ff[1] 
                                
                V.append(Vxx)

        plt.plot(LB,V)
        
            
    # calculate bending moment
        dx=0.001
        LB=np.linspace(0,self.length,num=1000) #0:dx:self.length
        L=self.length
        M=[]
        for xx in LB:
            if xx==0:
                if Maz==[]:
                    M.append(0)
                else:
                    M.append(-1.0*Maz[0][1]) # Mzz = -(sum of all moments until each section)
            else:
                Mzz=0 # Mzz at section is assumed to be positive
                for ff in Fy:
                    if len(ff)==3: # point load
                        if (ff[2] < xx):
                            Mzz += ff[1]*(-1.0)*(xx-ff[2]) #Moment about any section to the right of the loads =>(-1.0) for negation effect.
                    elif len(ff)==5: # UDL,rTDL
                        if ff[4][0]<xx: # if starting point less han current section
                            if not (ff[4][1]<xx): # if ending point not less han current section
                                if ('U' in str(ff[0])): # UDL N/m * m
                                    Mzz +=(-1.0)*ff[3]*(xx-ff[4][0])*((xx-ff[4][0])/2)
                                elif (('T' in str(ff[0])) & ('a' in str(ff[0]))): #rTDLrmax
                                    totalht=ff[4][1]-ff[4][0]
                                    Mzz +=(-1.0)*0.5*(xx-ff[4][0])*ff[3]*((xx-ff[4][0])/(ff[4][1]-ff[4][0]))*((1/3)*totalht)
                            else: # if we are past the udl or trdl use equivalent load to calculate the moment.
                                Mzz+=(-1.0)*ff[1]*(xx-ff[2]) # same for both UDL or TDL since the storage format is same.
                            
                
                for mm in Maz:
                    if (mm[2] < xx):
                        Mzz +=mm[1]
                        
                M.append(-Mzz) # Mzz = -(sum of all moments until each section)
                
        plt.plot(LB,M)
        plt.show()
                
    @staticmethod
    def isanumber(a):
        try:
            float(repr(a))
            outBool= True
        except:
            outBool= False
        return outBool
    
    @staticmethod
    def conformSignConvention(ang):
        # update to ASTC or use difference II:180-ang +costh for x and -sin th for y.
        conv=[(-1,0),(-1,-1),(0,-1),(1,-1),(1,0),(1,1),(0,1),(-1,1)]
        if ang==0:
            return conv[0]
        elif ang==90:
            return conv[2]
        elif ang==180:
            return conv[4]
        elif ang==270:
            return conv[6]
        elif (ang>0 & ang<90):
            return conv[1]
        elif (ang>90 & ang<180):
            return conv[3]
        elif (ang>180 & ang<270):
            return conv[5]
        elif (ang>270 & ang<360):
            return conv[7]
            
        
if __name__=="__main__":
    b=beam() # m
# test 1
    #b.add_section(0,0.3,10,90,'point',0.15) # pointing down

# test 2 (verified Hibbeler Example 7.10)
    if False:
        b.length=2
        print('Beam Length = '+str(b.length)+' '+b.unit)
        b.add_constraints('fixed','free')
        print('constraints = '+ str(b.constraints))
        b.add_section(0,1.2,400,90,'UDL',[0,1.2])
        b.add_section(1.2,2,600,90,'point',2)
        b.add_section(1.2,2,100,'cw','moment',2)
        print(b.sections)
        b.calc_reactions()
        
# test 3 (verified Hibbeler Example 7.11)
    if 1:
        b.length=8
        print('Beam Length = '+str(b.length)+' '+b.unit)
        b.add_constraints('roller','roller')
        print('constraints = '+ str(b.constraints))
        b.add_section(0,2,2,90,'point',2)
        b.add_section(2,4,3,90,'point',4)
        b.add_section(4,6,2,90,'point',6)
        print(b.sections)
        b.calc_reactions()

# test 4 (Hibbeler Example 7.12) (check error: moment dosent go to zero in the right end.
    if 0:
        b.length=20
        print('Beam Length = '+str(b.length)+' '+b.unit)
        b.add_constraints('hinge','hinge')
        print('constraints = '+ str(b.constraints))
        b.add_section(0,10,600,90,'point',10)
        b.add_section(10,15,4000,'cw','moment',15)
        print(b.sections)
        b.calc_reactions()

# test 5 (Hibbeler Example 7.8) check error . maximum of moment dosent match
    if 0:
        b.length=9
        print('Beam Length = '+str(b.length)+' '+b.unit)
        b.add_constraints('hinge','roller')
        print('constraints = '+ str(b.constraints))
        b.add_section(0,9,6,90,'rTDLrmax',[0,9])
        print(b.sections)
        b.calc_reactions()       

# symplify, subs, lambdify f=lambdify(x,expr,"numpy")=> f(x=a)
        
