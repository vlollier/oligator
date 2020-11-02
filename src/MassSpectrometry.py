#-*-coding:Utf-8-*-


__author__ ="Virginie Lollier"
__version__ = "1.0.0"
__license__ = "BSD"


import re,os
# for array intersection (pb sur set() & set() qui ne conserve pas l'ordre des éléments)
import numpy as np
import random

import platform
import networkx as nx
import matplotlib.pyplot as plt
#from .Graphes import *
from utils import Logger,Atom

def plotGose(G,ax=None):
    """
    """        
    coloromap={'pentose':'salmon','hexose':'lightblue'}
    colorose=[]
    for io in G.nodes():
        nct=G.nodes[io]['nbc']
        
        if nct in coloromap.keys():
            colorose.append(coloromap[nct])
        else:
            colorose.append('grey')

    pos = nx.spring_layout(G)
    nx.draw(G,pos,ax=ax)
    nx.draw_networkx_nodes(G,pos,node_color=colorose)
    nx.draw_networkx_labels(G,pos,font_size=14)
    edge_labels = nx.get_edge_attributes(G,'bound')        
    nx.draw_networkx_edge_labels(G, pos,edge_labels)                       


    if ax==None:
        plt.show()   



class Spectrum:
    ## from Domon & Costello
    positive={"B":0,"C":2,"X":[1],"Y":2,"Z":0}
    negative={"B":-2,"C":0,"A":[-1],"X":[-1],"Y":0,"Z":-2}    
    default={"B":0,"C":0,"A":0,"X":0,"Y":0,"Z":0}    
    
    def __init__(self,Gatom,Gose,start,dicoMass=None):
        glabels=Gose.copy()
        self.gmass=Gatom.copy()
        
        self.randomint={"[A,X]":[3000.0, 10000.0],"[B,C,Y,Z]":[10.0, 4000.0]}
        
        
        if dicoMass!=None:
            for n,m in dicoMass.items():
                Mass.dicomass[n]=m
                Mass.valence[n]=1
        else:
            Mass.dicomass=Atom.dicomass.copy()
            Mass.valence=Atom.dicovalence.copy()
        fragment=Mass( self.gmass,glabels,start)   
        
        self.ions=[]
        fragment.fragose(self.ions)
        fragment.fraglink(self.ions)
        fragment.fragzero(self.ions)
        self.add_rdm_intensity("A")
        
        Logger.debug("neutral mass of entire molecule %.2f"%fragment.precurmass())
        nx.draw(glabels)
        
    @staticmethod
    def __adduct__(add,dctions):
        ions=[]
        sign=re.sub("[^+-]","",add)
        atom=re.sub("[+-]","",add)
        
        if atom in Atom.dicocharge.keys():            
            mass=Atom.mass(atom)
            if sign=="-":
                mass*=-1
            for item in dctions:
                if isinstance(item,dict):
                    ion=item.copy()
                else:
                    ion={}
                    
                    elts=item.split(Mass.field_sep)
                    if len(elts)==3:
                        ion["mz"]=float(elts[0])
                        ion["intensity"]=float(elts[1])
                        ion["name"]=elts[2]                
                
                ion["mz"]+=mass                
                if isinstance(item,dict):
                    ions.append(ion)
                else:
                    ions.append("%f%s%f%s%s"%(ion["mz"],Mass.field_sep,ion["intensity"],Mass.field_sep,ion["name"]))                
        return ions
        
    
         
    
    def __dicomode__(dicomode,dicoion,atomcharge):
        """
        """
        ions=[]
        
        for item in dicoion:  
            ion={}
            if isinstance(item,dict):
                ion=item.copy()
            else:
                ion={}
                
                elts=item.split(Mass.field_sep)
                if len(elts)==3:
                    ion["mz"]=float(elts[0])
                    ion["intensity"]=float(elts[1])
                    ion["name"]=elts[2]
                
            if re.search("^C",ion["name"]):
                ion["mz"]+=dicomode["C"]*Atom.dicomass[atomcharge]
            elif re.search("^B",ion["name"]):
                ion["mz"]+=dicomode["B"]*Atom.dicomass[atomcharge]
            elif re.search("^Z",ion["name"]):
                ion["mz"]+=dicomode["Z"]*Atom.dicomass[atomcharge]
            elif re.search("^Y",ion["name"]):
                ion["mz"]+=dicomode["Y"]*Atom.dicomass[atomcharge]                                   
                
            elif re.search("}A_{",ion["name"]):
                ion["mz"]+=dicomode["A"]*Atom.dicomass[atomcharge]
               
            elif re.search("}X_{",ion["name"]):
                ion["mz"]+=dicomode["X"]*Atom.dicomass[atomcharge]
               
            
            if isinstance(item,dict):
                ions.append(ion)
            else:
                ions.append("%f%s%f%s%s"%(ion["mz"],Mass.field_sep,ion["intensity"],Mass.field_sep,ion["name"]))
            
    
        return ions
    
    #randomize intensities by default
    def get_ions(self):
        """
        """
        lsions=[]
        intensite=0
        for ion in self.ions:
            for regexpion,minmax in self.randomint.items():
                if re.search(regexpion,ion["name"])!=None:
                    intensite=random.uniform(minmax[0], minmax[1]) 
            if intensite>0:
                lsions.append("%0.3f;%0.3f;%s"%(ion["mz"],intensite,ion["name"]))
                
        return lsions
    
    def get_peaks(self,intensity_rule,isrange):
        """
        """
        lsions=[]
        intensite=0       
        
        
        for ion in self.ions:
            for regexpion,irule in intensity_rule.items():
                if re.search(regexpion,ion["name"])!=None:
                    if isrange:
                        intensite=random.uniform(irule[0], irule[1]) 
                    else:
                        intensite=irule
                       
                        
            if intensite>0:  
                lsions.append("%0.3f;%0.3f;%s"%(ion["mz"],intensite,ion["name"]))            


        return lsions
    
    def print(self):
        """
        """
        for i in self.lsions:            
            Logger.debug(i)
     
    def getseparator(self):
        """
        """
        return Mass.field_sep

    def add_rdm_intensity(self,typion,cfrom=None,cto=None,mini=0.0,maxi=0.0):
        """
        """
        if typion in ["A","X"]:
            if cfrom==None and cto==None:
                self.randomint[typion]=[mini,maxi]
            else:
                regexpion="^\\^{%i,%i}%s"%(cfrom,cto,typion)
                self.randomint[regexpion]=[mini,maxi]
    
            
    def ratio_intensity(self,dicoratio,maxi=100.0):
        """
        """
        lsions=[]
        intensite=maxi
        for ion in self.ions:
            for regexpion,ratio in dicoratio.items():
                if re.search(regexpion,ion["name"])!=None:
                    intensite*=ratio
            if intensite>0:
                lsions.append("%0.3f;%0.3f;%s"%(ion["mz"],intensite,ion["name"]))
                
        return lsions        
   

class Mass:    
    
    """
     Mass: partitions atom graph and computes masses from subgraphs
    """
    
    dicomass={}
    valence={}
   
        
    field_sep=';'
    
    def __init__(self,graphatom,graphose,reducing_end=None):
        self.Gatom=graphatom
        desoxynodes=[]   
        
        for n in self.Gatom.nodes():
            if self.Gatom.nodes[n]["atom"]=="":
                desoxynodes.append(n)
        for n in desoxynodes:
            self.Gatom.remove_node(n)       
        
        
        self.Gose=graphose
        self.dicoanhydro={}
        self.cyclatoms=self.__cycles__()
        self.greek=945        
        self.diconame={}
        self.re_osenum={}
        self.nre_osenum={}
        self.lsObound=self.__getObound__()
        if reducing_end!=None:
            self.set_redend(reducing_end)
        else:
            self.get_firstose()
            
    
       
    def get_firstose(self):
        """
        used when ose graph has not reducing_end attribute (get 1st ancestor of tree graph)
        """
        if nx.is_tree(self.Gose):
            [self.set_redend(n) for n,d in self.Gose.in_degree() if d==0] 
    
    def __getObound__(self):
        """
        """
        ls=[]
        for n in self.Gatom.nodes():
            descr=self.Gatom.nodes[n]
            if descr['atom']=='O':
                neigh=list(self.Gatom[n])
                if len(neigh)==2:                    
                    descr1=self.Gatom.nodes[neigh[0]]
                    descr2=self.Gatom.nodes[neigh[1]]
                    if descr1['ose']!=descr2['ose']:
                        ls.append(n) 
                        edge=(descr1['ose'],descr2['ose'])
                        if edge in self.Gose.edges():
                            self.Gose.edges[edge]["o"]=n   
                        else:
                            edge=(descr2['ose'],descr1['ose'])
                            try:
                                self.Gose.edges[edge]["o"]=n   
                            except:
                                Logger.debug("problem MS/MS : unable to find osidic binding")
                                
        return ls
    
    def __cycles__(self):
        """
        """
        anhy_edges=[]
        self.dicoanhydro={}
        cyclodico={}
        for edge in self.Gatom.edges():
            if "link" in self.Gatom.edges[edge]:                
                if self.Gatom.edges[edge]["link"]=="anhydro":
                    anhy_edges.append(edge)
                
        cycle_nodes= nx.cycle_basis(self.Gatom, 1)      
        
        
        
        for i in range(0,len(cycle_nodes)-1):
            for j in range(i+1,len(cycle_nodes)):                
                s1=set(cycle_nodes[i])
                s2=set(cycle_nodes[j])
                inter=list(s1.intersection(s2))
                if len(inter)>0:
                    if len(s1)<len(s2):
                        self.dicoanhydro[self.Gatom.nodes[inter[0]]["ose"]]=cycle_nodes[i]
                    else:
                        self.dicoanhydro[self.Gatom.nodes[inter[0]]["ose"]]=cycle_nodes[j]
                        
       
        if len(self.dicoanhydro.keys())>0:
            if len(anhy_edges)>0:
                G=self.Gatom.copy()
                for e in anhy_edges:
                    G.remove_edge(*e)
                cycle_nodes= nx.cycle_basis(G, 1)  
            else:
                Logger.debug("Error anhydro links are missing in edge attribute")
                
        for o in self.Gose.nodes():
            cyclodico[o]=[]
            for cycle in cycle_nodes:
                for n in cycle:
                    descr=self.Gatom.nodes[n]
                    if descr['ose']==o:
                        cyclodico[o].append(n) 
                    if descr['atom']=='O':
                        # for intra-cyclic ion labels
                        self.Gatom.nodes[n]['cnum']=0 
                        
        return cyclodico
    
    def set_redend(self,reducing_end) :
        """
        """        
        self.r=reducing_end
        self.assign_branch(reducing_end,'') 
        self.__assign_osenum__(reducing_end,1)
        self.__parse_nrepath__()
    
    # branch naming from ose start
    def assign_branch(self,start,name):
        """
        branch naming from ose start
        """
        desc=self.tri(self.Gose[start])        
        d=len(desc)
        self.diconame[start]=name
        
        if d>0:
            if d==1:
                self.assign_branch(desc[0],name)
            else:
                # alpha beta gamma...
                if name=="":
                    for i in range(0,len(desc)):
                        name=chr(self.greek+i)      
                        self.assign_branch(desc[i],name)
                        
                # ' '' '''...
                else:
                    for i in range(0,len(desc)):
                        name+="'"                        
                        self.assign_branch(desc[i],name)                                    
                
       
    
    def get_nbdesc(self,ose,num):
        """
        """
        for n in self.Gose[ose]:
            if n!=ose:
                self.get_nbdesc(n,num+1)
        return num
    
    def get_weight(self,ose,w):
        """
        """
        
        for n in self.Gose[ose]:
            if n!=ose:
                w=self.get_weight(n,w)+self.get_omass(ose)
        return w        
    
    
    def tri(self,lso):
        """
        branch naming according the number of descendants {nbdesc:[oses]} (should be the computed masses)
        """
        lstrie=[]
        lsd={}
        for n in lso:
            #k=self.get_nbdesc(n,0)
            k=self.get_weight(n,0)
            if k in lsd.keys():            
                lsd[k].append(n)
            else:
                lsd[k]=[n]
        ktri=list(lsd.keys())
        ktri.sort(reverse=True)
        for i in ktri:
            for elt in lsd[i]:
                lstrie.append(elt)
        return lstrie
    
    def get_omass(self,ose):
        """
        """        
        ls=self.get_oseatoms_degree(ose,0,0)  
        m=self.compute_mass(ls)
        m=1
        return m

    

    def compute_mass(self,lsatom_bond):
        """
        """
        m=0
        for ab in lsatom_bond:
            atom=re.sub("[0-9]","",ab[0])
            links=ab[1]
            if atom.find("=")!=-1:
                atom=re.sub("=","",atom)                
                m-=self.dicomass["H"]*2
               
            nh=self.valence[atom]-links
            
            if nh<0:
                Logger.debug("valence probleme! "+ atom+" "+str(links),1)
            else:
                m+=self.dicomass[atom]+self.dicomass['H']*nh
        return m
        
    def __assign_osenum__(self,node,num):
        """
        """
        self.re_osenum[node]=num
        num+=1                 
        desc=self.Gose[node]
        if len(desc)>0:
            for n in desc:
                self.__assign_osenum__(n,num)
    
    
    
    def __parse_nrepath__(self):
        """
        numeroter les oses en partant des extremites (non-reducing ends)
        pour les troncs communs, numerote à partir de la + longue branche        
        """
        ions={}
        paths=[]
        upath=[]
        g=self.Gose
        for i in g.nodes():
            if len(g[i])==0:
                paths.append(nx.shortest_path(g,source=self.r,target=i))
        if len(paths)>1:
            upath=unique_path(paths)            
        else:
            upath=[paths[0]]
        
        for p in upath:             
            p.reverse()
            self.__assign_nre_osenum(p,0)    
    
    def __assign_nre_osenum(self,path,num):                 
        """
        """
        if num<len(path):
            node=path[num]            
            self.nre_osenum[node]=num+1
            num+=1                
            self.__assign_nre_osenum(path,num)        
    
    def fragzero(self,outlist):
        zlab="Z_{0}"
        ylab="Y_{0}"
        
        blab="B_{"+str(self.nre_osenum[self.r])+self.diconame[self.r]+"}"
        clab="C_{"+str(self.nre_osenum[self.r])+self.diconame[self.r]+"}"        
        
        c1node=None
        for n in self.Gatom.nodes():
            if self.Gatom.nodes[n]["ose"]==self.r and self.Gatom.nodes[n]["atom"]=="C" and self.Gatom.nodes[n]["cnum"]==1 and not self.Gatom.nodes[n]["substit"]:
                c1node=n
                break
         
        neigh=list(self.Gatom[c1node])
        right=[]
        for n in neigh:
            if self.Gatom.nodes[n]["substit"] and self.Gatom.nodes[n]["atom"]=="O":                
                right.append(n)
                            
        if len(right)==1:                 
            morceaux=self.partitionne([(c1node,right[0])])
            node_left=morceaux[c1node]
            node_right=morceaux[right[0]]     
            if len(node_right)>1:    
                mleft=self.__calcmass__(node_left)
                mright=self.__calcmass__(node_right)  
                
                outlist.append({"mz":mleft+self.dicomass["O"],"name":clab})
                outlist.append({"mz":mleft,"name":blab})                          
                
                outlist.append({"mz":mright-self.dicomass["O"],"name":zlab})
                outlist.append({"mz":mright,"name":ylab})                  
    
    def fraglink(self,outlist):
        """
        """        
        if self.r!=None:
            for l in self.lsObound:
            
                edges=[]
                [edges.append(i) for i in self.Gatom[l]]
                
                if len(edges)==2:                                
                    ose0=self.Gatom.nodes[edges[0]]['ose']
                    ose1=self.Gatom.nodes[edges[1]]['ose']
                    atom0=self.Gatom.nodes[edges[0]]['atom']
                    atom1=self.Gatom.nodes[edges[1]]['atom']
                    left=l
                    right=edges[0]
                    ## quel sens???                
                    if (ose0,ose1) in self.Gose.edges():                    
                        zlab="Z_{"+str(self.re_osenum[ose0])+self.diconame[ose1]+"}"
                        ylab="Y_{"+str(self.re_osenum[ose0])+self.diconame[ose1]+"}"
                        
                        blab="B_{"+str(self.nre_osenum[ose1])+self.diconame[ose1]+"}"
                        clab="C_{"+str(self.nre_osenum[ose1])+self.diconame[ose1]+"}"
                        
                        right=edges[0]                    
                        
                    else:                    
                        zlab="Z_{"+str(self.re_osenum[ose1])+self.diconame[ose0]+"}"
                        ylab="Y_{"+str(self.re_osenum[ose1])+self.diconame[ose0]+"}"
                        blab="B_{"+str(self.nre_osenum[ose0])+self.diconame[ose0]+"}"                    
                        clab="C_{"+str(self.nre_osenum[ose0])+self.diconame[ose0]+"}"                    
                        right=edges[1]
                        
                    morceaux=self.partitionne([(left,right)])
                    node_left=morceaux[left]
                    node_right=morceaux[right]                
                                       
                    mleft=self.__calcmass__(node_left)
                    mright=self.__calcmass__(node_right)                 
                    
                    
                    
                    outlist.append({"mz":mleft,"name":clab})
                    outlist.append({"mz":mleft-self.dicomass["O"],"name":blab})                          
                    
                    outlist.append({"mz":mright,"name":zlab})
                    outlist.append({"mz":mright+self.dicomass["O"],"name":ylab})                        
                else:                    
                    Logger.debug('erreur nombre de connexions: %s'%(len(edges)),2)
                
        else:            
            Logger.debug("erreur: pas d'extremité non réductrice désignée dans la structure",2)
                    
            
            
    def fragose(self,outlist):
        """
        """
        if self.r!=None:
            for o in self.Gose.nodes():                   
                edges=[]   
                for e in self.casse_cycle(self.cyclatoms[o]) :
                    if not self.protect_edge(e,o):
                        edges.append(e)
                
                branch=self.diconame[o]
                nre=str(self.nre_osenum[o])
                re=str(self.re_osenum[o]-1)
                xion="X_{"+re+branch+"}"
                aion="A_{"+nre+branch+"}"
                
               
                
                # trouver le C du cycle le + pres de l'extremite reductrice (startc)
                startc=min(self.cyclatoms[o])  # pour eviter une valeur nulle
                if o!=self.r:
                    ascendant=nx.shortest_path(self.Gose,self.r,o)
                    ascendant.reverse()
                    ascendant=ascendant[1]
                
                    olink=self.Gose.edges[ascendant,o]["o"]
                    for n in self.Gatom[olink]:
                        if self.Gatom.nodes[n]['ose']==o:
                            startc=n 
                else:
                    for i in self.cyclatoms[o]:
                        if self.Gatom.nodes[i]['cnum']==1:
                            startc=i
                            break
                
                
                
                # calculer la masse des ions X et A (les ions X contiennent le startc)
                for casse in edges:                    
                    morceaux=self.partitionne(casse)   
                    massx=0
                    massa=0
                    for g in morceaux.keys():                        
                        if startc in morceaux[g]:
                            massx=self.__calcmass__(morceaux[g])-Mass.dicomass["H"]
                        else:
                            massa=self.__calcmass__(morceaux[g])-Mass.dicomass["H"]
                        
                    
                    lblion=""
                    edge=casse[0]                    
                    edge=[self.Gatom.nodes[edge[0]]['cnum'],self.Gatom.nodes[edge[1]]['cnum']]
                    edge.sort()
                    ncc=len(self.cyclatoms[o])-1
                    if edge!=[0,ncc]:                        
                        lab0=min(edge[0]-startc,edge[1]-startc)+startc                  
                    else:                        
                        lab0=ncc
                    
                    edge=casse[1]                                                 
                    edge=[self.Gatom.nodes[edge[0]]['cnum'],self.Gatom.nodes[edge[1]]['cnum']]
                    edge.sort()    
                    if edge!=[0,ncc]:
                        lab1=min(edge[0]-startc,edge[1]-startc)+startc 
                    else:
                        lab1=ncc
                    lblion=str(min(lab0,lab1))+","+str(max(lab0,lab1))
                    
                    outlist.append({"mz":massx,"name":"^{%s}%s"%(lblion,xion)})
                    outlist.append({"mz":massa,"name":"^{%s}%s"%(lblion,aion)})
        else:            
            Logger.debug("erreur: pas d'extremité non réductrice désignée dans la structure",2)    
        
    
        
    def precurmass(self):
        """
        """
        return self.__calcmass__(list(self.Gatom.nodes()))
    
    def __calcmass__(self,morceau):  
        """
        """
        mass=0
        for n in morceau:
            atom=re.sub("[0-9]","",self.Gatom.nodes[n]['atom'])                        
            if atom.find("=")!=-1:
                atom=re.sub("=","",atom)                
                mass-=self.dicomass["H"]*2            
            nh=self.valence[atom]-nx.degree(self.Gatom,n)
            if nh<0:
                print("valence probleme! "+ atom+" "+str(nx.degree(self.Gatom,n)),1)
            else:                        
                mass+=self.dicomass[atom]+self.dicomass['H']*nh    
        return mass   
  
    def protect_edge(self,edge,ose):
        if ose in self.dicoanhydro.keys():
            cyc=self.dicoanhydro[ose]
            s1=set(edge[0])
            s2=set(edge[1])
            scyc=set(cyc)
            if s1.issubset(scyc) or s2.issubset(scyc):
                Logger.debug(str(edge)+" protected")
                return True
            else:
                return False
        return False
        
   
    def casse_cycle(self,path):
        """
        """
        ncc=len(path)
        points=[]
        edges=[]
        for i in range(0,ncc-1):
            for j in range(i+1,ncc):
                points.append([i,j])
        
        for p in points:
            edge1=(path[p[0]],path[p[0]+1])
            if p[1]==ncc-1:
                edge2=(path[p[1]],path[0])
            else:
                edge2=(path[p[1]],path[p[1]+1])
        
            edges.append([edge1,edge2])
    
        return edges
                            
    def get_oseatoms_degree(self,iose,fromc,toc):
        """
        """
        result=[]
        nodes=[] 
        
        for inode in self.Gatom.nodes():
            description=self.Gatom.nodes[inode]
            if description["ose"]==iose:               
                if fromc==toc:
                    nodes.append(inode)
                else:
                    if description["cnum"] in range(fromc,toc+1):
                        nodes.append(inode)
        
        for n in nodes:    
            binome=(self.Gatom.nodes[n]['atom'],nx.degree(self.Gatom,n))
            result.append(binome)
        return result
    
    
    
    def partitionne(self,edges):
        """
        """
        result={}
        Gprime=nx.Graph()
        Gprime.add_nodes_from(self.Gatom.nodes())
        Gprime.add_edges_from(self.Gatom.edges())
        
        Gprime.remove_edges_from(edges)        
        
        
        for e in edges:            
            deb=e[0]
            fin=e[1] 
                        
            result[deb]=self.__subgroup__(Gprime,deb)
            result[fin]=self.__subgroup__(Gprime,fin)
                 
            
        return result    
    
    def __subgroup__(self,G,start):
        """
        """
        group=[start]
        self.__collect4__(G,start,start,group)
        return group 
    
    def __collect4__(self,G,start,node,ls):         
        """
        """
        for n in G[node]:                
            if n !=start and n not in ls:                          
                ls.append(n)
                self.__collect4__(G,node,n,ls)
        return ls         
  
    
    
def tri_path(paths,pref,tri):          
    """
    internal function of the module
    """
    tri.append(pref)
    paths.remove(pref)
    
    for ip in range(0,len(paths)):
        paths[ip]=diff(paths[ip],pref)  
        
    if len(paths)>1:        
        tri_path(paths,paths[longest(paths)],tri)
    else:
        #tri.append(paths[0])
        if len(paths)>0:
            tri.append(paths[0])
            
def diff(t1,t2):
    result=[]
    for i in t1:
        if i not in t2:
            result.append(i)
    
    return result

def longest(paths):
    l=0
    pref=0
    for i in range(0,len(paths)):
        p=paths[i]
        if len(p)>l:
            pref=i
            l=len(p)
    return pref


def unique_path(paths):
    result=[]
    tri_path(paths,paths[longest(paths)],result)
    return result
