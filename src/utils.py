#-*-coding:Utf-8-*-
__author__ ="Virginie Lollier"
__version__ = "1.0.1"
__license__ = "BSD"

import datetime


class Atom:
   
    dicomass={'Ag':106.905095,'Al':26.981541,'Ar':39.962383,'As':74.921596,'Au':196.96656,'B':11.009305,'Ba':137.905236,'Be':9.012183,'Bi':208.980388,'Br':78.918336,'C':12,'=C':12,'Ca':39.962591,'Cd':113.903361,'Ce':139.905442,'Cl':34.968853,'Co':58.933198,'Cr':51.94051,'Cs':132.905433,'Cu':62.929599,'Dy':163.929183,'Er':165.930305,'Eu':152.921243,'F':18.998403,'Fe':55.934939,'Ga':68.925581,'Gd':157.924111,'Ge':73.921179,'H':1.007825,'He':4.002603,'Hf':179.946561,'Hg':201.970632,'Ho':164.930332,'I':126.904477,'In':114.903875,'Ir':192.962942,'K':38.963708,'Kr':83.911506,'La':138.906355,'Li':7.016005,'Lu':174.940785,'Mg':23.985045,'Mn':54.938046,'Mo':97.905405,'N':14.003074,'=N':14.003074,'Na':22.98977,'Nb':92.906378,'Nd':141.907731,'Ne':19.992439,'Ni':57.935347,'O':15.994915,'=O':15.994915,'Os':191.961487,'P':30.973763,'Pb':207.976641,'Pd':105.903475,'Pr':140.907657,'Pt':194.964785,'Rb':84.9118,'Re':186.955765,'Rh':102.905503,'Ru':101.904348,'S':31.972072,'Sb':120.903824,'Sc':44.955914,'Se':79.916521,'Si':27.976928,'Sm':151.919741,'Sn':119.902199,'Sr':87.905625,'Ta':180.948014,'Tb':158.92535,'Te':129.906229,'Th':232.038054,'Ti':47.947947,'Tl':204.97441,'Tm':168.934225,'U':238.050786,'V':50.943963,'W':183.950953,'X':125.904281,'Xe':131.904148,'Y':88.905856,'Yb':173.938873,'Zn':63.929145,'Zr':89.904708}
    dicovalence={'C':4,'O':2,'N':3,'S':6,'P':5,'Cl':1,'I':1,'Se':1,'Te':1,'Br':1,'F':1,'':0,'Li':1,'K':1,'Cs':1}
    dicocharge={'H':1,'Li':1,'Cl':-1,'Na':1,'K':1,'Cs':1}
    
    @staticmethod
    def mass(atom):
        if atom in Atom.dicomass.keys():
            return Atom.dicomass[atom]
        else:
            return -1
    @staticmethod
    def valence(atom):
        if atom in Atom.dicovalence.keys():
            return Atom.dicovalence[atom]
        else:
            return -1  
    
    @staticmethod
    def charge(atom):
        if atom in Atom.dicocharge:
            return Atom.dicocharge[atom]
        else:
            return 0
            
        

class Logger:
    LEVEL=0
    CODE_LEVEL={0:"INFO ",1:"WARNING ",2:"ERROR "}
    
    @staticmethod
    def debug(msg,lvl=0):        
        if lvl>=Logger.LEVEL:
            if lvl in Logger.CODE_LEVEL.keys():
                print(Logger.CODE_LEVEL[lvl]+" "+str(datetime.datetime.now())+" " +str(msg))
            else:
                print(Logger.CODE_LEVEL[0]+" "+str(datetime.datetime.now())+" " +str(msg))
    
    
            

