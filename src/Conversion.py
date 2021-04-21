#-*-coding:Utf-8-*-
__author__ ="Virginie Lollier"
__version__ = "1.0.1"
__license__ = "BSD"

import re,os,string,sys
import networkx as nx


from GlycanTopology import OsidicBond,OseModel,Substitution,SubstitutionLibrary,SingletonTopo
import matplotlib.pyplot as plt
import platform
from networkx.drawing.nx_agraph import graphviz_layout
from utils import Logger,Atom


RefForm={4:'furane',5:"pyrane"}
"""
"""
FormRef={'furane':4,"pyrane":5}
"""
"""
RefNbC={10:"decose",9:"nonose",8:"octose",7:"heptose",6:"hexose",5:"pentose",4:"tetrose",3:"triose"}
"""
"""
NbCRef={"decose":10,"nonose":9,"octose":8,"heptose":7,"hexose":6,"pentose":5,"tetrose":4,"triose":3}
"""
"""
# from monosaccharideDB
IsomerCode={"1":"L","2":"D","0":""}
"""
"""
CodeIsomer={"L":"1","D":"2","":"X"}
"""
"""

# non-exhaustive list...
WurcsTerminalCarbs=[("m",""),("h","O"),("C","CO/1O/1."),("o","=O"),("A","O/1=O")]
WurcsNonTermCarbs=[("d",""),("C","CO/1O"),("O","=O")]

   

def get_nodes(G,**kwargs):
    """
    search nodes having the given attributes
    :param G: graph where nodes are searched
    :param kwargs: atom= , in_cycle=
    
    :return: a list of nodes
    :rtype: list
    """
    nodes=[]
    if G!=None:
        for n in G.nodes():
            match=True
            for k,v in kwargs.items():
                if k in G.nodes[n] and v==G.nodes[n][k]:
                    match&=True
                else:
                    match&=False
            if match:
                nodes.append(n)
    return nodes

def tricycle(tab,firstelt,sens):
    """
    
    """
    result=[]    
    cptab=tab.copy()[::sens]        
        
    index=cptab.index(firstelt)
    for i in range(index,len(cptab)):
        result.append(cptab[i])
    for i in range(0,index):
        result.append(cptab[i])
    
    return result

def plotG(G):
    pos=nx.spring_layout(G)
    nx.draw(G,pos)
    nx.draw_networkx_labels(G,pos,nx.get_node_attributes(G,'atom')  )
    plt.show()

def smi2graph(smi,G=None,starti=0):
    """
    """
    if G==None:
        # debug
        #print(smi)
        
        #init
        G=nx.Graph(startnode=starti)
        atoms=re.findall("=?[A-Z][a-z]?[0-9]*",smi)
        debcyc=-1
        fincyc=-1
        
        node_cycles={}
        
        i=starti
        for elt in atoms:
            if elt!="R":
                val=Atom.dicovalence[re.sub("[=0-9]","",elt)]                
                #if re.search("=",elt):
                    #val-=1    
            else:
                val=2
            
            G.add_node(i,atom=elt,valence=val)
            if re.search(elt,smi):
                start,end=re.search(elt,smi).span()    
                smi=smi[0:start]+str(i)+","+smi[end:]       
                
                if re.search("[0-9]",elt):                        
                    # marche pour les substituants où il n'y a plus de 10 cycles (C32 -> cycles 3 et 2)
                    cycnums=re.findall("[0-9]",elt)
                    for cycnum in cycnums:
                        if cycnum in node_cycles.keys():                        
                            node_cycles[cycnum].append(i)
                        else:
                            node_cycles[cycnum]=[i]   
                    
            i+=1
        for cyc in node_cycles.keys():
            edge=node_cycles.get(cyc)
            #print("cycle "+str(edge[0])+" "+str(edge[1]))
            G.add_edge(*edge,close_cyc=True)                        

            
        smi2graph(smi,G)
        
    else:            
        have_branch=re.search("\([0-9,]+\)",smi)
        if have_branch:
            branche=re.search("\([0-9,]+\)",smi).group()
            nodes=re.findall("[0-9]+",branche)            
            edge=re.search("[0-9]+,\("+nodes[0],smi).group()
            edge=re.findall("[0-9]+",edge)                
            G.add_edge(int(edge[0]),int(edge[1]))
            
            for n in range(0,len(nodes)-1):
                G.add_edge(int(nodes[n]),int(nodes[n+1]))
            smi=re.sub(re.escape(branche),"",smi)
            smi2graph(smi,G)
        
        elif len(re.findall("[0-9]+",smi))>0:
            nodes=re.findall("[0-9]+",smi)
            for n in range(0,len(nodes)-1):                    
                G.add_edge(int(nodes[n]),int(nodes[n+1]))        
          
    return G    

#def smi2graph(smi,G=None,starti=0):
    #"""
    #"""
    #if G==None:
        #print(smi)
        ##init
        #G=nx.Graph(startnode=starti)
        #atoms=re.findall("=?[A-Z][0-9a-z]?",smi)
        #debcyc=-1
        #fincyc=-1
        
        #i=starti
        #for elt in atoms:
            #G.add_node(i,atom=elt)
            #if re.search(elt,smi):
                #start,end=re.search(elt,smi).span()    
                #smi=smi[0:start]+str(i)+","+smi[end:]            
                #if re.search("[0-9]",elt):
                    #if debcyc>0:
                        #fincyc=i
                    #else:
                        #debcyc=i                
            #i+=1
        #if debcyc>0 and fincyc>0:
            #print(str(debcyc)+" "+str(fincyc))
            #G.add_edge(int(debcyc),int(fincyc))        
        #smi2graph(smi,G)
        
    #else:
        #have_branch=re.search("\([0-9,]+\)",smi)
        #if have_branch:
            #branche=re.search("\([0-9,]+\)",smi).group()
            #nodes=re.findall("[0-9]+",branche)            
            #edge=re.search("[0-9]+,\("+nodes[0],smi).group()
            #edge=re.findall("[0-9]+",edge)
            #G.add_edge(int(edge[0]),int(edge[1]))
            #for n in range(0,len(nodes)-1):
                #G.add_edge(int(nodes[n]),int(nodes[n+1]))
            #smi=re.sub(re.escape(branche),"",smi)
            #smi2graph(smi,G)
        
        #elif len(re.findall("[0-9]+",smi))>0:
            #nodes=re.findall("[0-9]+",smi)
            #for n in range(0,len(nodes)-1):
                #G.add_edge(int(nodes[n]),int(nodes[n+1]))
            
    #return G


def __collect4__(G,start,node,ls):         
    """
    recursively add in a list all the connected nodes from a start
    
    :return: a list of node numbers
    """
    for n in G[node]:                
        if n !=start and n not in ls:                          
            ls.append(n)
            __collect4__(G,node,n,ls)
    return ls     

class GlycoCTEncoder:
    """
    """
    
    dicomods={"desoxy":"d","oxydation":"a","keto":"keto","double_bond":"en"}
    def __init__(self,gatom,gose):   
        """
        """
        self.gatom=gatom.copy()
        self.gose=gose
        gctose={}
        self.reslist=[]
        self.linlist=[]
        lso=[]
        for n in self.gatom.nodes():
            self.gatom.nodes[n]["atom"]=re.sub("[0-9]*","",self.gatom.nodes[n]["atom"])
        
        for o in gose.nodes():
            ot=ose_template(o)
            lso.append(ot)
            ot.set_txtnbC(gose.nodes[o]["nbc"])
            ot.set_txtform(gose.nodes[o]["form"])
            startc=111
            endc=-1
            for n in get_nodes(self.gatom,atom="C",in_cycle=True,ose=o)    :
                cnum=self.gatom.nodes[n]["cnum"]
                if cnum<startc:
                    startc=cnum
                if cnum>endc:
                    endc=cnum
            ot.set_cycle_start(startc)
            ot.set_cycle_end(endc)
            
            for n in get_nodes(self.gatom,atom="C",substit=False,ose=o):
                cnum=self.gatom.nodes[n]["cnum"]
                if "isomer" in self.gatom.nodes[n]:
                    ciso=self.gatom.nodes[n]["isomer"]
                    ot.set_iso(cnum,ciso)
                    
                else:
                    ot.set_iso(cnum,"")
                    
            
            for e in gose.edges():
                bond=re.findall("[0-9]+",gose.edges[e]["bond"])
                if o==int(bond[0]):
                    ot.mods[int(bond[1])-1]=-1
                elif o==int(bond[3]):
                    ot.mods[int(bond[2])-1]=-1
            ot.set_mod(endc,-1)        
            for cnum in range(1,ot.nbc()+1):
                if ot.get_mod(cnum)!=-1 and cnum!=endc:
                    subformul=GlycoCTEncoder.getsubformule(self.gatom,o,cnum)
                    ot.set_mod(cnum,subformul)
            
            
    
                
        resnum=0
        corres_id={}
        for o in lso:            
            resnum+=1
            self.reslist.append("%i%s"%(resnum,self.fromOseToGlycoCT(o)))
            corres_id[o.identifier]=resnum
            for icarb in range(0,len(o.mods)):
                if o.mods[icarb]!=-1 and o.mods[icarb]!="":
                    subname=SubstitutionLibrary.get_subformul(o.mods[icarb]).name.lower()
                    if subname!="hydroxy" and subname not in GlycoCTEncoder.dicomods.keys():
                        resnum+=1
                        self.reslist.append("%is:%s"%(resnum,subname))
                        corres_id[(o.identifier,icarb+1)]=resnum                        
           
        
        linlist=[]
        
        for e in gose.edges():            
            bond=[]
            for b in re.findall("[0-9]+",gose.edges[e]["bond"]):             
                bond.append(int(b))
            linlist.append("%io(%s+%s)%sd"%(corres_id[bond[0]],bond[1],bond[2],corres_id[bond[3]]))
            
        for ref,val in corres_id.items():
            if isinstance(ref,tuple):                
                linlist.append("%so(%s+1)%sn"%(corres_id[ref[0]],ref[1],val))                 
        
        linlist=GlycoCTEncoder.__sort_asnum__(linlist)
        linnum=0
        for lin in linlist:
            linnum+=1
            self.linlist.append("%i:%s"%(linnum,lin))
                     
            
    
    def getsubformule(g,ose,cnum):
        """
        """
        formul=""
        degat={}        
        nodes=get_nodes(g,ose=ose,cnum=cnum,substit=True)        
        
        if len(nodes)>0:             
            for n in nodes:                
                deg=0
                atom=re.sub("[0-9]+","",g.nodes[n]["atom"])
                if re.search("=",atom):
                    cnode=get_nodes(g,ose=ose,cnum=cnum,substit=False,atom="C")[0]
                    if (cnode,n) in list(g.edges()):
                        deg=1
                    else:
                        deg=2
                    atom=re.sub("=","",atom)
                deg+=nx.degree(g,n)   
                if atom in degat.keys():
                    degat[atom]+=1
                else:
                    degat[atom]=1
                
                nh=Atom.valence(atom)-deg
                if "H" in degat.keys():                                      
                    degat["H"]+=nh
                else:
                    degat["H"]=nh
                    
        for a in sorted(degat.keys()):
            if degat[a]>0:
                formul+=a+str(degat[a])
        return formul
   
    def fromOseToGlycoCT(self, ot):
        """
        """
        name=""
       
        c1=""
        if ot.get_iso(ot.cyclebounds[0])=="D":
            c1="a"
        elif ot.get_iso(ot.cyclebounds[0])=="L":
            c1="b"
        else:
            c1="x"
       
        txtiso=""
        name1="ukn"
        for icarb in range(ot.cyclebounds[0],ot.cyclebounds[1]-1):   
            if ot.get_iso(icarb+1) in ["D","L"]:
                txtiso+=ot.get_iso(icarb+1)
            else:
                txtiso+=" "
        # pentose pyrane
        if ot.cyclebounds[1]==ot.nbc():
            if txtiso[:-1] in IupacName._iupac:
                name1=txtiso[-1:]+IupacName._iupac[txtiso[:-1]]    
        else:
            if txtiso in IupacName._iupac:
                name1=ot.get_iso(ot.cyclebounds[1])+IupacName._iupac[txtiso]
            
        name1=name1.lower()
        txtiso=""
        name2=None
        # exclusivement pour l'acide neuraminic! connais pas la regle à appliquer pour répondre à d'autres cas (demander à un chimiste)...
        for cnum in range(ot.cyclebounds[1]+2,ot.nbc()):
            if ot.get_iso(cnum)!="x":
                txtiso+=ot.get_iso(cnum)
        if txtiso!="":
            name2=txtiso[-1:]+IupacName._iupac[txtiso[:-1]] 
            name2=name2.lower()
       
        if name2:
            name="b:%s-%s-%s-%s-%i:%i"%(c1,name2,name1,RefNbC[ot.nbc()][:3].upper(),ot.cyclebounds[0],ot.cyclebounds[1]) 
        else:        
            name="b:%s-%s-%s-%i:%i"%(c1,name1,RefNbC[ot.nbc()][:3].upper(),ot.cyclebounds[0],ot.cyclebounds[1]) 
        
        
        
        # parse mods[1] for linkage carbon and type 
        for icarb in range(0,len(ot.mods)):
            if ot.mods[icarb]!=-1:
                # desoxy
                if ot.mods[icarb]=="":
                    name+="|%i:d"%(icarb+1)               
                elif SubstitutionLibrary.get_subformul(ot.mods[icarb]):
                    subname=SubstitutionLibrary.get_subformul(ot.mods[icarb]).name.lower()                   
                    if subname in GlycoCTEncoder.dicomods.keys():
                        name+="|%i:%s"%(icarb+1,GlycoCTEncoder.dicomods[subname])
        return name
    
    def __sort_asnum__(ls):
        """
        """
        oo=[]
        nbo=[]
        for i in ls:
            nbo.append(int(re.findall("^[0-9]+",i)[0]))
        
        trio=sorted(nbo)
        for elt in trio:    
            index=nbo.index(elt)
            oo.append(ls[index])
            nbo[index]=""
            
        return oo
        
class GlycoCTDecoder:
    """
    """
    patc1="[abox]"
    patose="[ldx][a-z]{3}"
    patnbc="[A-Z]{3}"
    patring="[1-9]:[1-9]"
    patterm="[1-9]:(d|a|aldi|keto|en|geminal|sp|sp2)"    
    
    def __init__(self,residues,links):
        """
        """
        self.gatom=nx.Graph()
        self.gose=nx.DiGraph()
        gct_ids={}
        c1free=[]
        self.check=None
        
        for line in residues:
            
            gctid=re.sub("[a-z]","",line[:line.find(":")])
            gct_ids[gctid]=GlycoCTDecoder.readResidue(line)
            if isinstance(gct_ids[gctid],ose_template):
                c1free.append(gctid)
        
        #apply ose modifs first
        for line in links:
            ifrom,cfrom,cto,ito=re.findall("[0-9]+",line)
            linkage=re.findall("[odhn]+",line)
            if isinstance(gct_ids[ifrom],ose_template) and isinstance(gct_ids[ito],Substitution):
                gct_ids[ifrom].set_mod(int(cfrom),gct_ids[ito].identifier)                    
                GlycoCTDecoder.applymod2gatom(gct_ids[ifrom],int(cfrom)) 
        
            
            elif isinstance(gct_ids[ifrom],Substitution) and isinstance(gct_ids[ito],ose_template): 
                gct_ids[ito].set_mod(int(cto),gct_ids[ifrom].identifier)                    
                GlycoCTDecoder.applymod2gatom(gct_ids[ito],int(cto))                    
                
            #else:
                #Logger.debug("glycoct format error on line %s"%line,2)  
        
        # apply osidic binding second
        removeO=-1
        for line in links:
            ifrom,cfrom,cto,ito=re.findall("[0-9]+",line)
            linkage=re.findall("[odhn]+",line)
            
            if isinstance(gct_ids[ifrom],ose_template) and isinstance(gct_ids[ito],ose_template):
                ofrom=gct_ids[ifrom]
                oto=gct_ids[ito]
                if int(cfrom)==ofrom.cyclebounds[0]:
                    c1free.remove(ifrom) 
                    
                if int(cto)==oto.cyclebounds[0]:
                    c1free.remove(ito)                                    
                
                #osidic binding
                edge=[-1,-1]
                opos=-1
                if not int(ifrom) in list(self.gose.nodes()):
                    self.gatom=nx.disjoint_union(self.gatom,ofrom.gatom)
                    self.gose.add_node(int(ifrom),form=ofrom.get_txtform(),nbc=ofrom.get_txtnbc())
                if not(int(ito)) in list(self.gose.nodes()):
                    self.gatom=nx.disjoint_union(self.gatom,oto.gatom)
                    self.gose.add_node(int(ito),form=oto.get_txtform(),nbc=oto.get_txtnbc())                        
                                      
                if linkage[0]=="d" and linkage[1]=="o":
                    remove=get_nodes(self.gatom,ose=int(ifrom),atom="O",cnum=int(cfrom),substit=True,in_cycle=False)
                    if len(remove)==1:
                        
                        removeO=remove[0]
                        edge[0]=get_nodes(self.gatom,ose=int(ifrom),atom="C",cnum=int(cfrom),substit=False)[0]
                        edge[1]=get_nodes(self.gatom,ose=int(ito),atom="O",cnum=int(cto),substit=True)[0]
                        opos=edge[1]
                elif linkage[1]=="d" and linkage[0]=="o":
                    remove=get_nodes(self.gatom,ose=int(ito),atom="O",cnum=int(cto),substit=True,in_cycle=False)
                    if len(remove)==1:
                        
                        removeO=remove[0]
                        edge[1]=get_nodes(self.gatom,ose=int(ito),atom="C",cnum=int(cto),substit=False)[0]
                        edge[0]=get_nodes(self.gatom,ose=int(ifrom),atom="O",cnum=int(cfrom),substit=True)[0]
                        opos=edge[0]                        
                else:
                    Logger.debug("GlycoCT: failed to connect oses %s and %s"%(ifrom,ito),1)
                    
                if opos>0 and removeO>0:
                    self.gatom.remove_node(removeO)
                    self.gose.add_edge(int(ifrom),int(ito),bond=re.sub("[odhn]+","",line),opos=opos)
                    self.gatom.add_edge(*edge)
                else:                    
                    Logger.debug("GlycoCT: failed to connect oses %s and %s"%(ifrom,ito),2)
                    self.check="obond"
        
        if len(c1free)==1:
            self.gose.graph["reducing_end"]=int(c1free[0])
            self.gatom.graph["reducing_end"]=get_nodes(self.gatom,ose=int(c1free[0]),cnum=gct_ids[c1free[0]].cyclebounds[0],atom="C",in_cycle=True)[0]
        else:            
            self.gose.graph["reducing_end"]=1
            self.gatom.graph["reducing_end"]=get_nodes(self.gatom,ose=self.gose.graph["reducing_end"],cnum=gct_ids["1"].cyclebounds[0],atom="C",in_cycle=True)[0]
        
       
    
    def checkformat(res,lin):
        """
        """
        
        linpattn="[0-9]+[odhn]\([0-9]\+[0-9]+\)[0-9]+[odhn]"
        idres=[]
        
        if len(res)==0:
            return "no residue found"
        if len(lin)==0:
            return "no link found"
        
        for r in res:            
            
            if re.match("^[0-9]+b:.*",r):
                elts=r[r.find(":")+1:].split("-")
                check=True
                for elt in elts:
                    if re.match(GlycoCTDecoder.patc1,elt):
                        check=True
                    elif re.match(GlycoCTDecoder.patnbc,elt):
                        check=True
                    elif re.match(GlycoCTDecoder.patose,elt):
                        check=True
                    elif re.search(GlycoCTDecoder.patring,elt):
                        check=True
                    else:
                        check=False
                if check:
                
                    descr=r.split(":")
                    idres.append(re.sub("[^0-9]","",descr[0]))
                else:
                    return "malformed RES: %s"%r                
            elif re.search("[0-9]+[s]:[\-a-z]+",r):
                descr=r.split(":")
                idres.append(re.sub("[^0-9]","",descr[0]))                
            else:
                return "unsupported form"
        for l in lin:
            if re.match(linpattn,l):
                o1=l.split("+")[0]
                o2=l.split("+")[1]
                o1=o1.split("(")[0]
                o2=o2.split(")")[1]
                o1=re.sub("[^1-9]","",o1)
                o2=re.sub("[^1-9]","",o2)
                if o1 not in idres:
                    return "LIN id %s  not in residue list"%o1
                if o2 not in idres:
                    return "LIN id %s  not in residue list"%o2
            else:
                return "unsupported LIN: %s"%l
        return None
    
    def readResidue(line,inode=1):
        """
        """
        sep=line.find(":")
        dicomods={"d":"desoxy","a":"oxydation","keto":"keto","en":"double_bond"}
        
        if sep>0:
            residue=line[sep+1:]
            identif=int(line[:sep-1])    
            restype=line[sep-1]
            if restype=="b":
                descr=residue.split("-")  
                ot=ose_template(identif)
                
                if len(descr)>=4: 
                    c1iso=None
                    carbiso=""
                    ring=""
                    nbc=""    
                    closure=""
                    termmods=[]
                    # collect description
                    for elt in descr:
                        if re.match(GlycoCTDecoder.patnbc,elt):
                            nbc=6
                            for k in NbCRef.keys():
                                if elt==k.upper()[:3]:
                                    nbc=NbCRef[k]
                                    break
                        elif re.match(GlycoCTDecoder.patc1,elt):
                            if elt=="a":
                                c1iso="D"
                            elif elt=="b":
                                c1iso="L"
                            else:
                                c1iso=""
                        elif re.match(GlycoCTDecoder.patose,elt):
                            end=re.sub("x","",elt[0]).upper()
                            for k,v in IupacName._iupac.items():                                
                                if elt[1:]==v.lower():  
                                    if carbiso=="":
                                        carbiso=k+end
                                    else:
                                        carbiso=k+end+"x"+carbiso
                                    break
                        elif re.search(GlycoCTDecoder.patring,elt):
                            closure=re.search(GlycoCTDecoder.patring,elt).group(0)
                            termmods=elt.split("|")[1:]
                    
                    # design template
                    if nbc!="":                        
                        ot.set_nbc(int(nbc))
                        if closure!="":
                            cycle=closure.split(":")
                            ot.set_form(int(cycle[1])-int(cycle[0])+1)
                            ot.set_cycle_start(int(cycle[0]))
                            ot.set_cycle_end(int(cycle[1]))                        
                            
                            if c1iso!=None:
                                ot.set_iso(int(cycle[0]),c1iso)
                                
                            if carbiso!="":                                 
                                for i in range(0,len(carbiso)):                                            
                                    ot.set_iso(int(cycle[0])+i+1,carbiso[i]) 
                                    
                                ot.get_graphatom(inode,False)
                                # says terminal or cycle carbons of ose are modified                                                    
                                if len(termmods)>0:
                                    for carbmod in range(0,len(termmods)):
                                        if re.match(GlycoCTDecoder.patterm,termmods[carbmod]):
                                            carb=int(termmods[carbmod].split(":")[0])
                                            mod=termmods[carbmod].split(":")[1]     
                                            ## keto type is already known by cycle bounds
                                            if mod!="keto":
                                                if mod=="d":
                                                    ot.mods[carb-1]=SubstitutionLibrary.get_subfromname(dicomods["d"]).identifier
                                                    ot.set_iso(carb,"")
                                                    GlycoCTDecoder.applymod2gatom(ot,carb)                                                    
                                                elif mod in dicomods.keys():
                                                    ot.mods[carb-1]=SubstitutionLibrary.get_subfromname(dicomods[mod]).identifier
                                                    GlycoCTDecoder.applymod2gatom(ot,carb)                                                    
                                           
                                                else:
                                                    Logger.debug("unknown glycoCT key %s"%mod,2)    
                                         
                       
                        
                        
                        
                        return ot
            elif restype=="s":
                return SubstitutionLibrary.get_subfromname(residue.lower())
            else:
                return None
    
    def applymod2gatom(ot,cnum):    
        """
        """
        if isinstance(ot.mods[cnum-1],int):
            notation=SubstitutionLibrary.getSub(ot.mods[cnum-1]).smiles
            if notation=="":
                removeO=ot.get_nodes(substit=True,cnum=cnum,atom="O")[0]
                ot.gatom.remove_node(removeO)                
            else:
                # oxydation par exemple (=O)O , on force un noeud racine pour avoir C(=O)O dans gatom
                if notation[0]=="(":
                    notation="Q"+notation
                    
                smidecod=SmilesDecoder(notation,False)
                modgraph=smidecod.Gatom.copy()
                remove=-1
                carbnode=-1
                bout=None                
               
                for node in modgraph.nodes():
                    modgraph.nodes[node]["ose"]=ot.identifier
                    modgraph.nodes[node]["cnum"]=cnum
                    modgraph.nodes[node]["substit"]=True
                    modgraph.nodes[node]["in_cycle"]=False
                    if modgraph.nodes[node]["atom"]=="Q":
                        bout=node
                
                if bout:
                    modgraph.remove_node(bout)
                for subg in nx.connected_component_subgraphs(modgraph):                       
                    modgraph.nodes[min(list(subg.nodes()))]["connector"]=True
                                  
                ot.gatom=nx.disjoint_union(ot.gatom,modgraph)
                removeO=ot.get_nodes(substit=True,cnum=cnum,atom="O")[0]
                ot.gatom.remove_node(removeO)
                connectors=ot.get_nodes(substit=True,cnum=cnum,connector=True)
                carbnode=ot.get_nodes(substit=False,atom="C",cnum=cnum)[0]
                for con in connectors:
                    ot.gatom.add_edge(carbnode,con)
                
                            

    
class IupacName:
    """
    """
    
    #dictionnaire des noms IUPAC
    _iupac = {}
    _iupac["DDD"] = "All"
    _iupac["LDD"] = "Alt"
    _iupac["DLD"] = "Glc"
    _iupac["LLD"] = "Man"
    _iupac["DDL"] = "Gul"
    _iupac["LDL"] = "Ido"
    _iupac["DLL"] = "Gal"
    _iupac["LLL"] = "Tal"
    _iupac["DD"] = "Rib"
    _iupac["LD"] = "Ara"
    _iupac["DL"] = "Xyl"
    _iupac["LL"] = "Lyx"
    _iupac["D"] = "erythro"
    _iupac["L"] = "Tro"
    _iupac[""] = "Gro"
    
    _isocolr = {}
    _isocolr["DDD"] = "purple"
    _isocolr["LDD"] = "magenta"
    _isocolr["DLD"] = "blue"
    _isocolr["LLD"] = "green"
    _isocolr["DDL"] = "darkorange"
    _isocolr["LDL"] = "saddlebrown"
    _isocolr["DLL"] = "gold"
    _isocolr["LLL"] = "cyan"
    _isocolr["DD"] = "magenta"
    _isocolr["LD"] = "green"
    _isocolr["DL"] = "darkorange"
    _isocolr["LL"] = "gold"
    _isocolr["D"] = ""
    _isocolr["L"] = ""
    _isocolr[""] = ""    
    
    #dictionnaire des noms des modifications
    _substitution = {}
    _substitution["hydroxy_hydroxy_hydroxy_hydroxy_hydroxy_hydroxy_"]= ""
    _substitution["hydroxy_n-acetyl_hydroxy_hydroxy_hydroxy_hydroxy_"]= "NAc"
    _substitution["hydroxy_amino_hydroxy_hydroxy_hydroxy_hydroxy_"]= "N"
    _substitution["hydroxy_hydroxy_hydroxy_hydroxy_hydroxy_carboxyl_"]= "A"
    _substitution["hydroxy_n-sulfate_hydroxy_hydroxy_hydroxy_hydroxy_"]= "S"
    _substitution["hydroxy_hydroxy_hydroxy_hydroxy_hydroxy_hydrogen_"]= ""
    _substitution["hydroxy_hydroxy_hydroxy_hydroxy_hydroxy_"]= ""
    _substitution["hydroxy_n-acetyl_hydroxy_hydroxy_hydroxy_"]= "NAc"
    _substitution["hydroxy_amino_hydroxy_hydroxy_hydroxy_"]= "N"
    _substitution["hydroxy_hydroxy_hydroxy_hydroxy_carboxyl_"]= "A"
    _substitution["hydroxy_n-sulfate_hydroxy_hydroxy_hydroxy_"]= "S"
     
    
    @staticmethod
    def iupacName(form, carbon, modif) :
        """
        """
        #determine le nom IUPAC 
        
        #nb de carbone et forme du cycle
        nbC = len(modif)
        c1 = carbon[0] #positionnement du groupement OH sur le carbone 1
        c5 = carbon[nbC-2] #positionnement du groupement OH sur le carbone 5
        modif_c6 = modif[nbC-1] #modification sur le carbone 6
        
        #configuration absolue et configuration du carbone anomerique
        carb1 = IupacName.anomericCarbon(c1)
        carb5 = IupacName.absoluteConfiguration(c5)
        
        
        #forme du cycle 
        if (form=="pyrane"):
            formcycle = "p"
        else :
            formcycle = "f"        
        
        
        #recupere la position des groupements sur chaque carbones 
        isomerie = ""
        for i in range(1, nbC-2) :
            isomerie = isomerie + carbon[i]
        #recupere les substitutions subit sur tous les carbones 
        substitution = ""
        for i in range(len(modif)):   
            substitution = substitution + modif[i] + "_"
        
        #notation des substitutions et de l'isomerie
        subs = IupacName.substitutions(substitution)
        iso = IupacName.typeOfOse(isomerie, nbC, modif_c6)
        
        
        #on retourne le nom IUPAC avec les informations retenues plus haut 
        #si les substitutions ne sont pas connues on donne le type de l'ose (hexose ou pentose)
        if subs == "unknown" :
            if(nbC==6):
                return carb1 + carb5 + "HEX" + formcycle
            elif(nbC==5):
                return carb1 + carb5 + "PEN" + formcycle
            else :
                return carb1 + carb5 + "OSE" + formcycle
        else :
            return carb1 + "-" + carb5 + "-" + iso + formcycle + " " + subs
      
    def getIUPACName(isomerie):
        """
        """
        if isomerie in IupacName._iupac:
            return IupacName._iupac[isomerie]
        else :
            return "ose"       
        
    @staticmethod        
    def typeOfOse(isomerie, nbC, modif):
        """
        """
        if isomerie in IupacName._iupac :
            #recherche s'il y a une substitution sur le carbone 6
            if (nbC == 6 and modif=="hydrogen"):
                if IupacName._iupac[isomerie]=="Gal" :
                    return "Fuc"
                elif IupacName._iupac[isomerie]=="Man" :
                    return "Rha"
                elif IupacName._iupac[isomerie]=="Glc" :
                    return "Qui"            
                else :            
                    return IupacName._iupac[isomerie]
            else :
                return IupacName._iupac[isomerie]
        else :
            if(nbC==6):
                return "HEX"
            elif(nbC==5):
                return "PEN"
            else :
                return "OSE"
    
    
    @staticmethod        
    def anomericCarbon(position):
        """        
        :return: the configuration of then anomeric carbon
        
        """
        #retourne la configuration du carbone anomerique
        #position du goupement sur le carbone 1 qui donne la position alpha ou beta
        if position=="D":
            return "alpha"
        elif position == "L":
            return "beta"
        else :
            return "x"  
    
    
    @staticmethod        
    def absoluteConfiguration(position):  
        """
        assign ose series (D or L)
        :param position: carbon number ending the ose cycle
        """
        #position du groupement sur l'avant dernier carbone qui donne la serie de l'ose (D ou L)
        if position=="D":
            return "D"
        elif position=="L" :
            return"L"
        else : 
            return "X"        
    
    
    @staticmethod
    def substitutions(modif):
        """        
        """
        #recherche dans le dictionnaire des substitutions la notation des modifications subit par l'ose 
        #si le dictionnaire ne contient pas cette notation on renvoie "unknown"
        if modif in IupacName._substitution :
            return IupacName._substitution[modif]
        else :
            return "unknown" 

class SmilesEncoder:  
    """
    """
    def __init__(self,graph):        
        self.G=graph.copy()
       
        
        if self.G.size()>0:
            self.nbcy=0
            self.cycles=nx.cycle_basis(self.G)   
            typecyc=SmilesDecoder.generic_type_cycle(self.G,self.cycles)
            anhdro=[]
            id_anhydro=0
            rmcyc=[]
            icyc=0
            for cyc in self.cycles:
                if typecyc[icyc]=="ose":
                    for node in cyc:
                        if graph.nodes[node]["atom"]=="O" and graph.nodes[node]["cnum"]!=0:                        
                            neigh=list(graph[node])
                            anhydro=(node,neigh[0])                        
                            if anhydro not in anhdro:
                            
                                anhdro.append(anhydro)
                                self.G.remove_edge(*anhydro)
                                id_anhydro+=1
                                self.G.nodes[node]["anh"]=id_anhydro
                                self.G.nodes[neigh[0]]["anh"]=id_anhydro                            
                            break
                icyc+=1
                
            self.sub_cycles={}
            subids={}
            for n in self.G.nodes():      
                if self.G.nodes[n]["substit"]:
                    atom=self.G.nodes[n]["atom"]
                    ose=self.G.nodes[n]["ose"]
                    cnum=self.G.nodes[n]["cnum"]
                    subid=str(ose)+"_"+str(cnum)
                    if n in subids.keys():
                        subids[subid].append(n)
                    else:
                        subids[subid]=[n]
                    if re.search("[0-9]",atom):                        
                        # marche pour les substituants où il n'y a plus de 10 cycles (C32 -> cycles 3 et 2)
                        cycnums=re.findall("[0-9]",atom)                        
                        for cycnum in cycnums:   
                            cycnum=subid+"_"+cycnum
                            if cycnum in self.sub_cycles.keys():                                    
                                self.sub_cycles[cycnum].append(n)                                
                            else:
                                self.sub_cycles[cycnum]=[n]
                                
                                
            if len(self.sub_cycles)>0:                
                for cyc in self.sub_cycles.keys():
                    if len(self.sub_cycles.get(cyc))==2:
                        self.G.remove_edge(*self.sub_cycles.get(cyc))
                        for n in self.sub_cycles.get(cyc):
                            if "subcycle" in self.G.nodes[n]:
                                self.G.nodes[n]["subcycle"]+=","+cyc
                            else:
                                self.G.nodes[n]["subcycle"]=cyc
                    else:
                        print("rror")           
                
                c2=nx.cycle_basis(self.G)   
                for c in self.cycles:
                    if c in c2:
                        for elt in c:
                            self.G.nodes[elt]["in_cycle"]=True                        
                    else:
                        for elt in c:
                            self.G.nodes[elt]["in_cycle"]=False                        
                    #if c not in c2:
                        #for elt in c:
                            #self.G.nodes[elt]["in_cycle"]=False
                    
                self.cycles=c2.copy()
            
            Logger.debug("number of nodes: "+str(len(self.G.nodes()))+", number of cycles: "+str(len(self.cycles)))      
            
            
            # supprimer les numeros de cycle qui seront réaffectés
            for n in self.G.nodes():
                self.G.nodes[n]['atom']=re.sub("[0-9]","",self.G.nodes[n]['atom'])
                self.G.nodes[n]['cend']=False
          
            #ajouter l'attribut cycle end
            for cyc in self.cycles:
                tricycle=sorted(cyc,key=lambda n:self.G.nodes[n]['cnum'],reverse=True)
                self.G.nodes[tricycle[0]]['cend']=True
            
            #nx.write_graphml(self.G,"/tmp/smienc_gatom.graphml")
            
    def parse_graphcycles(self):
        start=0
        g=self.G.copy()
        
        if "reducing_end" in self.G.graph.keys():
            start=self.G.graph['reducing_end']
        else:
            # find reducing end        
            start=self.find_start()                
         
            
            
    def in_cycle(self,inode): 
        """
        """
        for cycle in self.cycles:         
            if inode in cycle:
                return True
        return False
    
    
    
    # quand l'info n'est pas un attribut du graphe
    def find_start(self):
        """
        when this information is not speficied by graph attribute, search any C1 (or C2 if ketose) without osidic binding
        """
        lsC1=[]
        oco=['O','C','O']
        oco.sort()  
        
        start=list(self.G.nodes())[0]
        # tous les noeuds
        for n in self.G.nodes():
            # qui sont dans 1 cycle et dont l'atome est C
            if self.in_cycle(n) and self.G.nodes[n]['atom']=='C':
                neigh=[]
                atoms=[]
                # ont pour voisin(s)...
                for v in self.G[n]:
                    neigh.append(v)
                    atoms.append(self.__getatom__(v))
                    
                atoms.sort()
                # ... 2 oxygenes
                if atoms==oco:
                    # dont l'un des 2 n'est pas lie a 1 cycle (peut-etre ajouter patron de cycle osidique)
                    for a in range(0,len(atoms)):
                        if atoms[a]=='O' and not self.in_cycle(neigh[a]):                            
                            aubout=False
                            nextnode=neigh[a] 
                            prevnode=n
                            while not aubout:
                                voisinage=list(self.G[nextnode])
                                voisinage.remove(prevnode)
                                if len(voisinage)==0:
                                    start=n
                                    aubout=True
                                    
                                else:
                                    for node in voisinage:
                                        if self.in_cycle(node):
                                            aubout=True
                                    prevnode=nextnode
                                    nextnode=voisinage[0]
        return start
                
    def __getatom__(self, inode):
        """
        :return: atom name +/- smiles isomerie if available of the given node number
        """
        if inode in self.G.nodes():
            note=self.G.nodes[inode]['atom']            
            if "isomer" in self.G.nodes[inode].keys():                    
                if self.G.nodes[inode]["isomer"]=="L":  
                    # cf SMILES specification for ring closure when notation writing starts with ose 1, carbon 1
                    if self.G.nodes[inode]['cend']:
                        note="["+note+"@@"+"H]"                
                    else:                        
                        note="["+note+"@"+"H]"
                elif self.G.nodes[inode]["isomer"]=="D":                    
                    if self.G.nodes[inode]['cend']:
                        note="["+note+"@"+"H]"
                    else:
                        note="["+note+"@@"+"H]"                
            
            return note
    
    # cyclic atoms first
    def __sortneigh__(self,neigh):
        """
        sort nodes such that cycle nodes are first elements
        
        :param neigh: list of nodes
        
        :rtype: list
        """
        cyn=[]
        noncyn=[]
        for n in neigh:
            if self.in_cycle(n):
                cyn.append(n)
            else:
                noncyn.append(n)
        cyn.sort()
        noncyn.sort()
        return cyn+noncyn
                 
    def branch(self,inode,lsnode):
        smi=";"+str(inode)+";"
        #smi=self.__getatom__(inode)
        lsnode.append(inode)
        neigh=list(set(self.G[inode])-set(lsnode))        
        neigh=self.__sortneigh__(neigh)  
        is_end=False
        if inode==self.osend: 
            self.osend=None 
            is_end=True   
            
        if len(neigh)>0:
            if self.in_cycle(inode): 
                if self.osend==None and not is_end:
                    cycle_neigh=[]
                    atom_neigh=[]
                    for n in neigh:
                        if self.in_cycle(n):
                            cycle_neigh.append(n)
                            atom_neigh.append(self.__getatom__(n))                    
                    if 'O' in atom_neigh:
                        self.osend=cycle_neigh[atom_neigh.index('O')]                        
                        
                    else:                                                                             
                        self.osend=max(cycle_neigh) 
                    neigh.remove(self.osend)  
                            
            elif "subcycle" in self.G.nodes[inode]:
                neigh=list(set(self.G[inode])-set(lsnode))        
                neigh.sort()                
                            
            for n in neigh[:-1]:                    
                smi+="("+self.branch(n,lsnode)+")"
            smi+=self.branch(neigh[-1:][0],lsnode)
        
        return smi            
            
    
                
    
    # graph networkx avec attribut "atom"
    # build smiles notation from graph starting from given node number (should be the reducing end)  
    def graph2smi(self):
        """
        build smiles notation from atom graph starting (the reducing end)  
        """
        self.osend=None
        self.subends={}
        self.cnum=0   
        self.redend=None
        
         
            
        if "reducing_end" in self.G.graph.keys():
            self.redend=self.G.graph['reducing_end']
        else:
            # find reducing end        
            self.redend=self.find_start()                
        
            
        
        Logger.debug("node number of reducing end: "+str(self.redend))
        ls=[]
        lsrank={}
        
        s0=self.branch(self.redend,ls)
        
        fract=re.findall("\(?;[0-9]+;\)?",s0)
        
        dicoli={}
        #for i in range(0,len(fract)):
            #dicoli[fract[i]]=ls[i]        
        if len(ls)==len(fract):
            for i in range(0,len(ls)):
                dicoli[ls[i]]=fract[i]
                lsrank[ls[i]]=i
        # collect cycle nodes from graph attributes            
        nanh={}
        subcyc={}
        for n in self.G.nodes():
            attr=self.G.nodes[n]
            if "anh" in attr:
                if attr["anh"] in nanh.keys():
                    nanh[attr["anh"]].append(n)
                else:
                    nanh[attr["anh"]]=[n]
            elif "subcycle" in attr:
                nums=re.split(",",attr["subcycle"])
                for s in nums:
                    if s in subcyc.keys():
                        subcyc[s].append(n)
                    else:
                        subcyc[s]=[n]
            
        cycles=nx.cycle_basis(self.G)
        for a in nanh.keys():
            cycles.append(nanh.get(a))
        for s in subcyc.keys():
            cycles.append(subcyc.get(s))
        
        
        
            
        for icyc in range(0,len(cycles)):
            cyc=cycles[icyc]              
            c1=sorted(cyc,key= lambda inode:lsrank.get(inode))
            cycles[icyc]=[c1[0],c1[-1]]
        cycles=sorted(cycles,key=lambda c:lsrank.get(c[0]))    
        numcyc=0
        #cyclenodes=[]
        for cyc in cycles:
            numcyc+=1
            if numcyc==10:
                numcyc=1
            for cnode in cyc:
                value=dicoli.get(cnode) 
                if re.search("\)$",value):
                    dicoli[cnode]=value[0:-1]+str(numcyc)+")"                   
                    
                else:
                    dicoli[cnode]+=str(numcyc)
                
                
        smi=""
        for inode in ls:
            note=dicoli.get(inode)            
            inode=re.search(";[0-9]+;",note).group()
            inode=re.sub(";","",inode)
            atom=self.__getatom__(int(inode))
            smi+=re.sub(";"+inode+";",atom,note)
            

                
        return smi    



   
    
class SmilesDecoder:
    """
    """
    # https://en.wikipedia.org/wiki/CPK_coloring
    coloratomap={"C":"grey","O":"red","N":"darkblue","S":"yellow","P":"orange","Cl":"green","I":"darkviolet"}
    #patatom="=?[A-Z][0-9a-z]?"
    patatom="=?[A-Z][a-z]?[0-9]*"
    RefForm={4:'furane',5:"pyrane"}
    RefNbC={5:"pentose",6:"hexose",9:"non"}
    
    
        
    def __init__(self,notation,ose=True):
        """
        """
        self.smi=notation
        if re.search("@",notation)==None:
            self.smi=notation
            self.isonode=[]
        else:
            isosmi=re.sub("[^a-zA-Z@]","",notation)
            isosmi=re.sub("@@H","2",isosmi)
            isosmi=re.sub("@H","1",isosmi)            
            self.isonode=re.findall(self.patatom,isosmi)
            self.smi=re.sub("[@H\\[\\]]","",notation) 
            
        
        # nodes
        self.nodatoms=re.findall(self.patatom,self.smi)
       
        
        self.edgatoms=[]
        self.Gatom=nx.Graph()
        self.Ose_graphs=[]
        self.ose_list=[]
       
        smi1=re.sub("[=0-9]+","",self.smi)     
        
        self.collect_nodatom()
        self.collect_edgeatom(smi1,list(range(1,len(self.nodatoms)+1)))
        # extraction des atomes constitutifs des cycles trouves dans le graphe    
        self.cycle_nodes= nx.cycle_basis(self.Gatom, 1)  
        self.anhydro_edge=[]
        
        # typer les cycles
        self.typcycl=self.type_cycle()
       
        # check intersection between cycles for internal bond (anhydro function)
        for i in range(0,len(self.cycle_nodes)-1):
            for j in range(i+1,len(self.cycle_nodes)):                    
                if self.typcycl[i]=="ose" or self.typcycl[j]=="ose":
                #if self.typcycl[i]!="benz" or self.typcycl[j]!="benz":
                    s1=set(self.cycle_nodes[i])
                    s2=set(self.cycle_nodes[j])
                    inter=s1.intersection(s2)
                    if len(inter)!=0:
                        bond=[]
                        onode=-1
                        cnode=-1
                        oselsn=[]
                        if len(s1)<len(s2):
                            bond=list(self.cycle_nodes[i])
                            oselsn=list(self.cycle_nodes[j])
                        else:
                            bond=list(self.cycle_nodes[j])
                            oselsn=list(self.cycle_nodes[i])
                        for node in bond:
                            if re.sub("[0-9]+","",self.Gatom.nodes[node]["atom"])=="O":
                                onode=node
                                break
                        if onode>0:
                            neigh=list(self.Gatom[onode])
                            for node in neigh:
                                if re.sub("[0-9]+","",self.Gatom.nodes[node]["atom"])=="C" and node in oselsn:
                                    cnode=node
                                    break
                        
                            if cnode>0:
                                edge=(onode,cnode)
                                self.anhydro_edge.append(edge)
                                self.Gatom.remove_edge(*edge)
        
        if len(self.anhydro_edge)>0:
            self.cycle_nodes= nx.cycle_basis(self.Gatom, 1)  
            self.typcycl=self.type_cycle()        
        if ose:
            try:
                self.parse_cycles()
                for anhy_edge in self.anhydro_edge:
                    self.Gatom.add_edge(*anhy_edge,link="anhydro")
                
            except ValueError:
                
                Logger.debug("unable to decode notation",2)
        
       
    def collect_nodatom(self):
        """
        creates the gatom nodes with the corresponding atom and isomer attributes
        """
        self.coloratoms=[]
        
        ## nodes
        
        for inode in range(0,len(self.nodatoms)):
            atom=re.sub("[=0-9\\(\\)]*","",self.nodatoms[inode])
            isomer=None
            if len(self.isonode)>0:
                isomer=re.sub("[A-Z][a-z]?","",self.isonode[inode])
                
            if atom in self.coloratomap:
                self.coloratoms.append(self.coloratomap[atom])
            else:
                self.coloratoms.append("grey")                
            
            if isomer and isomer in IsomerCode:
                self.Gatom.add_node(inode+1,atom=self.nodatoms[inode],in_cycle=False,isomer=IsomerCode[isomer])        
            else:
                self.Gatom.add_node(inode+1,atom=self.nodatoms[inode],in_cycle=False)        
                
    
    def __debranche__(self,branch,smi1,tabn):    
        """
        
        """
        start=smi1.find(branch)
        end=start+len(branch)
        B=len(re.findall(self.patatom,smi1[0:start]))
        N=len(re.findall(self.patatom,smi1[0:end]))-1    
        
        remove=[]
        for inode in (range(B,N)):
            self.Gatom.add_edge(tabn[inode],tabn[inode+1])
            remove.append(tabn[inode+1])
        
        for i in remove:        
            tabn.remove(i)    
        smi1=smi1[0:start+1]+smi1[end:]       
        return(smi1,tabn)
        
    def collect_edgeatom(self,smi1,tabn):
        """
        """
        patbranch="[A-Za-z]\\([A-Za-z]+\\)"    
        branchages=re.findall(patbranch,smi1)
        
        if len(branchages)>0:
            for branch in branchages:
                (smi1,tabn)=self.__debranche__(branch,smi1,tabn)                
            self.collect_edgeatom(smi1,tabn)
        else:
            if len(tabn)>0:
                for n in range(0,len(tabn)-1):
                    self.Gatom.add_edge(tabn[n],tabn[n+1])        
            lsnum=[]
            cnumycle={}
            for n in re.findall("[0-9]",self.smi):
                if n not in lsnum:
                    lsnum.append(n)
                    
            for n in lsnum:
                p=[]
                for inode in range(0,len(self.nodatoms)):
                    if re.search(n,self.nodatoms[inode]):    
                        p.append(inode+1)
                cnumycle[n]=p
                for i in range(0,len(p),2):
                    ideb=p[i]
                    ifin=p[i+1]
                    
                    Logger.debug("Cycle:%s %i - %s %i"%(self.nodatoms[ideb-1],ideb,self.nodatoms[ifin-1],ifin))
                    self.Gatom.add_edge(ideb,ifin)            
                  
           
    
   
   
            
    def get_atom(self,node):
        """
        :return: the "atom" attribute of th given node number
        """
        return re.sub('[0-9]','',self.Gatom.nodes[node]['atom'])
    
    
             
    
   
   
    def __collect4__(self,G,start,node,ls):         
        """
        recursively add in a list all the connected nodes from a start
        
        :return: a list of node numbers
        """
        for n in G[node]:                
            if n !=start and n not in ls:                          
                ls.append(n)
                self.__collect4__(G,node,n,ls)
        return ls        
   
    def check_obond(self,ls,node):
        """
        """        
        ls.append(node)
        if self.get_atom(node)=="C":
            if self.Gatom.nodes[node]["in_cycle"]==True:
                return node
            else:
                for n in list(self.Gatom[node]):
                    if n not in ls:
                        return self.check_obond(ls,n)
        return None
   
   
    def parse_cycles(self):
        """
        central function building the ose descriptions from the cycles detection in gatom
        """
        self.Gose=nx.DiGraph()
        icycle=0
        lsOB={}
        lsC1={}
        self.ose_list=[]
        cyc=[]
        
        ##cycle as subsituent?
        #subcyc=None
        ## test nb of double bonds
        #for cycle in self.cycle_nodes:
            #db=0            
            #for node in cycle:                
                #atom=self.Gatom.nodes[node]["atom"]        
                #if re.search("=",atom):
                    #db+=1
            #if db>0:
                #subcyc=cycle
                
        #if subcyc!=None:
            #self.cycle_nodes.remove(subcyc)
        
        subcycles=[] 
        for icycle in range(0,len(self.cycle_nodes)):  
            if self.typcycl[icycle]=="ose":
                for node in self.cycle_nodes[icycle]:            
                    cyc.append(node)
            else:
                subcycles.append(self.cycle_nodes[icycle])
        if len(subcycles)>0:
            for sub in subcycles:
                self.cycle_nodes.remove(sub)
    
        #coc
        coc=[]
        for node in self.Gatom.nodes():
            if node in cyc:
                self.Gatom.nodes[node]["in_cycle"]=True
            else:
                self.Gatom.nodes[node]["in_cycle"]=False
                
        for icycle in range(0,len(self.cycle_nodes)-1):
            for jcycle in range(icycle+1,len(self.cycle_nodes)):   
                #if self.typcycl[icycle]=="ose" and self.typcycl[jcycle]=="ose":
                start=self.cycle_nodes[icycle][0]
                end=self.cycle_nodes[jcycle][0]
                path=nx.shortest_path(self.Gatom,start,end)
                if len(path)>=3:
                    for i in range(0,len(path)-2):
                        triplet=path[i:i+3]
                        atoms=[]
                        for i in range(0,3):
                            atoms.append(self.get_atom(triplet[i]))
                        if atoms==["C","O","C"] and triplet[1] not in cyc:                                        
                            exist=False
                            for lob in coc:  
                                if set(triplet).issubset(set(lob)):
                                
                                    exist=True
                                    break
                                
                            if not exist:
                                coc.append(triplet)
                                
        
        gprime=self.Gatom.copy()
        g_template_oses=[]
        ox2add=[]
        # find coc that are not osidic binding
        rm_coc=[]
        for _,o in enumerate(coc):      
            
            right=self.check_obond([o[1]],o[0])
            left=self.check_obond([o[1]],o[2])
            if right==None or left==None:
                
                rm_coc.append(o)
            else:
                
                if o[1] in list(gprime.nodes()):
                    gprime.remove_node(o[1])
                else:
                    print("error in rm ox node: "+str(o[1]))
        
        for o in rm_coc:
            coc.remove(o)
        
        for o in coc:   
        
            right=o[0]
            left=o[2]                       
            node_right=self.__collect4__(gprime,right,right,[right])
            node_right.sort()            
            if node_right not in g_template_oses:
                g_template_oses.append(node_right)                
                ox2add.append([o[1]])
            else:
                i=g_template_oses.index(node_right)                
                ox2add[i].append(o[1])
            
            node_left=self.__collect4__(gprime,left,left,[left])
            node_left.sort()
            if node_left not in g_template_oses:
                g_template_oses.append(node_left)
                ox2add.append([o[1]])
            else:
                i=g_template_oses.index(node_left)
                ox2add[i].append([o[1]])
        
        ose_list=[]
        
        #! redondance pour les oses au milieu
        for i in range(0,len(g_template_oses)):
            ose_id=i+1
            g_template_oses[i]+=ox2add[i]
            
            ot=ose_template(ose_id)            
            self.ose_list.append(ot)
            gprime=self.Gatom.subgraph(g_template_oses[i]).copy()            
            ot.get_graph_cnum(gprime)
            
            lsC1[ose_id]=ot.getcyclebound_nodes()[0]
            
            
            
            for node in gprime.nodes():  
                cnum=gprime.nodes[node]["cnum"]
                
                if ot.get_iso(cnum) in ["D","L"] and not gprime.nodes[node]["substit"]:                    
                    #print(ot.get_iso(cnum))
                         # =================ISOMER THING ================
                    if not self.checkisomer(node,gprime):
                        ot.inverse_iso(cnum)
                        self.Gatom.nodes[node]["isomer"]=ot.get_iso(cnum)
                        #print("iso changed for node %i , cnum %i "%(node,cnum))
                         # ===========================================
                    
                        
                self.Gatom.nodes[node]["cnum"]=cnum
                self.Gatom.nodes[node]["ose"]=gprime.nodes[node]["ose"]
                self.Gatom.nodes[node]["substit"]=gprime.nodes[node]["substit"]
            #nx.write_graphml(self.Gatom,"gatom.graphml")
            self.Gose.add_node(ose_id,form=ot.get_txtform(),nbc=ot.get_txtnbc())     
        
                
                
        # set edges of ose graph and find reducing_end
        redend_candidates=list(self.Gose.nodes())
        for o in coc:
            oleft=None
            oright=None
            carbright=None
            carbleft=None
            for ot in self.ose_list:
                if o[0] in ot.gatom.nodes():                
                    oright=ot                
                    carbright=ot.gatom.nodes[o[0]]["cnum"]
                if o[2] in ot.gatom.nodes():
                    oleft=ot
                    carbleft=ot.gatom.nodes[o[2]]["cnum"]
            if oleft!=None and oright!=None:                  
                lsOB[o[1]]=OsidicBond(oright.identifier,oleft.identifier,carbright,carbleft)
                ### est-ce que le C1 d'un cetose represente une extremité réductrice (bien que non impliqué dans les liaisons O, c'est C2 a la place)??? 
                ## problem avec le C1Free de SingletonTopo
                ##if carbright==1:
                    ##redend_candidates.remove(oright.identifier)
                ##if carbleft==1:
                    ##redend_candidates.remove(oleft.identifier)  
                if carbright==oright.cyclebounds[0]:
                    redend_candidates.remove(oright.identifier)
                if carbleft==oleft.cyclebounds[0]:
                    redend_candidates.remove(oleft.identifier)  
                
                      
        
                            
        if len(redend_candidates)==1:            
            Logger.debug("youpi, SmilesDecoder found reducing end for ose number: "+str(redend_candidates[0]))
            self.Gose.graph["reducing_end"]=redend_candidates[0]
            self.Gatom.graph["reducing_end"]=lsC1[redend_candidates[0]]
                
            # assign edges direction of Gose
            self.__create_directedges__(lsOB,redend_candidates[0])
            
            
           
        else:
            Logger.debug(" could not find reducing end",2)
            Logger.debug("try to start with first cycle detected",2)
            self.Gose.graph["reducing_end"]=0
            
            # assign edges direction of Gose
            self.__create_directedges__(lsOB,0)      
            
        
        
     
    def __create_directedges__(self,lsob,parent):   
        """
        """
        for oposition,bond in lsob.items():
            connect=False
            if bond.parent_ose==parent:
                connect=True
            elif bond.child_ose==parent:
                bond.setParent(parent)           
                connect=True
                
                
            if connect:
                self.Gose.add_edge(bond.parent_ose,bond.child_ose,bond=bond.getAttributString(),opos=oposition)
                
                dd=lsob.copy()
                dd.pop(oposition)                
                self.__create_directedges__(dd,bond.child_ose)        
       
    @staticmethod
    def check_smiles(notation):
        """
        """
        if re.search("[\\(\\)]",notation):
            open_bracket=len(re.findall("\\(",notation))
            close_bracket=len(re.findall("\\)",notation))
            if open_bracket!=close_bracket:                
                Logger.debug("nb brackets inconsistent "+notation,2)
                return False
            #cycles=re.findall("[0-9]+",notation)
            #! pas plus de 10 cycles (9+1->1)
            cycles=re.findall("[0-9]",notation)
            if len(cycles)%2>0:                
                Logger.debug("nb cycles inconsistent "+notation,2)
                return False
            else:
                cycles.sort()
                for i in range(0,len(cycles)-1,2):
                    if cycles[i]!=cycles[i+1]:
                        Logger.debug("open/close cycle inconsistent "+str(cycles),2)
                        return False 
            
        return True
    
    # "3 neighbors that follow a tetrahedral specification" and implicit H
    def checkisomer(self,tetrahnode,G):
        cycle=None
        cnums=[]
        onode=-1
        # check the asumption @=1=L(beta) and @@=2=D(alpha)
        check=True
        for c in nx.cycle_basis(G):
            if tetrahnode in c:
                cycle=c
                break
        if cycle:
            for n in cycle:
                cnums.append(G.nodes[n]["cnum"])
                if re.match("O[0-9]?",self.nodatoms[n-1]):
                    onode=n
                    
        if G.nodes[tetrahnode]["cnum"]==max(cnums):            
            if onode>tetrahnode:
                if onode-tetrahnode==1:
                    check=False
                else:
                    check=True            
            else:
                check=True
            
        
                
        return check
            
    def generic_type_cycle(g,cycles):        
        typecycle={}
        icycle=0
        for cyc in cycles:
            atoms=[]
            compo={}
            for n in cyc:
                atom=re.sub("[0-9]","",g.nodes[n]["atom"])
                if atom in compo.keys():
                    compo[atom]+=1
                else:
                    compo[atom]=1
                    
            typec=""
            for elt in compo.keys():
                typec+=","+elt+str(compo.get(elt))
            typec=re.sub(",","",typec,count=1)            
            
            if "O" in compo.keys():  
                if set(compo.keys())=={"O","C"}:                    
                    typecycle[icycle]="ose"
                    #if compo.get("O")==1 and compo.get("C") in [5,6]:                    
                        #typecycle[icycle]="ose"
                    #else:                        
                        #typecycle[icycle]=typec
                else:
                    typecycle[icycle]=typec
            elif "=C" in compo.keys():
                if compo.get("=C")==3 :
                    if compo.get("C")==3:
                        typecycle[icycle]="benz"
                    elif compo.get("N")==1 :
                        typecycle[icycle]="benz"
                    else:                
                        typecycle[icycle]="nd"                        
                else:
                    typecycle[icycle]=typec
            else:
                typecycle[icycle]=typec         
            icycle+=1
        
        return typecycle        
        
    def type_cycle(self):         
        val=SmilesDecoder.generic_type_cycle(self.Gatom,self.cycle_nodes)
        return val
    
        
class WurcsDecoder:
    """
    """
    CarbonDescriptor_term={"A":"CO/2=O","d":"","h":"CO","m":"C","o":"C=O"}
    CarbonDescriptor_cyc={"d":""}
    def __init__(self,line):
        """
        """
        self.wurcs=line
        # demo
        if line=="":
            self.wurcs="2.0/3,3,2/[a2112h-1b_1-5_4*NCC/3=O][a1221h-1a_1-5][a2112h-1a_1-5]/1-2-3/a2-b1_a3-c1"
        parsing=self.wurcs
        wurcs_types=re.findall("\[[^\[\]]*\]",parsing)
        
        ose_types=[]
        temp=0
        for wtype in wurcs_types:
            parsing=re.sub(re.escape(wtype),"",parsing)
                   
            omod=re.sub("[\[\]]","",wtype)
            omod=omod.split("_")
            # h...h     
            if re.match("^h[12]+h$",omod[0]) :
                iso=omod[0]
                cycle="1-"+str(len(re.findall("[12]",iso))+1)
                modif=[]
                if len(omod)>1:
                    modif=omod[1:]
                omod=[iso,cycle]
                for m in modif:
                    omod.append(m)
                
            otype=self.create_ose(omod[0],omod[1],temp)
            temp+=1
            
            for mod in omod[2:]:
                cmod=mod.split("*")
                otype.set_mod(int(cmod[0]),cmod[1])
                
            ose_types.append(otype)    
            
        elts=parsing.split("/")
        ose_order=elts[3].split("-")
        ose_bonds=elts[4].split("_")
        
        self.Gose=nx.DiGraph()
        self.Gatom=nx.Graph()
        iatom=0
        
        # WURCS notation starts from reducing end
        self.Gatom.graph["reducing_end"]=iatom
        self.Gose.graph["reducing_end"]=1
        
        for o in range(0,len(ose_order)):
            ose=ose_types[int(ose_order[o])-1].copy(o+1)
            
            self.Gose.add_node(o+1,form=RefForm[ose.form],nbc=RefNbC[len(ose.carbons)])
            g=self.modgatom(ose,iatom)
            iatom+=max(g.nodes())+1
            if len(self.Gatom.nodes())==0:
                self.Gatom=g
            else:
                self.Gatom=nx.disjoint_union(self.Gatom,g)
        
        for bond in ose_bonds:
            connex=bond.split("-")
            ofrom=self.__decode_52__(re.sub("[0-9]","",connex[0]))
            cfrom=int(re.sub("[a-zA-Z]+","",connex[0]))
            oto=self.__decode_52__(re.sub("[0-9]","",connex[1]))
            cto=int(re.sub("[a-zA-Z]+","",connex[1]))
            Logger.debug("ose from:%s,carbon from:%s - ose to:%s,carbon to:%s"%(ofrom,cfrom,oto,cto))
            oxfrom=-1
            oxto=-1
            
            
            for node in self.Gatom.nodes():                
                if self.Gatom.nodes[node]["atom"]=="O":
                    if self.Gatom.nodes[node]["ose"]==ofrom and self.Gatom.nodes[node]["cnum"]==cfrom:
                        oxfrom=node
                    if self.Gatom.nodes[node]["ose"]==oto and self.Gatom.nodes[node]["cnum"]==cto:
                        oxto=node
                    
               
                    
            if oxfrom>=0 and oxto>=0:
                neigh=list(self.Gatom[oxfrom])
                if len(neigh)==1 :
                    node=neigh[0]
                    if self.Gatom.nodes[node]["atom"]=="C" and self.Gatom.nodes[node]["ose"]==ofrom and self.Gatom.nodes[node]["cnum"]==cfrom:
                        self.Gatom.remove_node(oxfrom)
                        self.Gatom.add_edge(node,oxto)                         
                        self.Gose.add_edge(ofrom,oto,bond="%i(%i+%i)%i"%(ofrom,cfrom,cto,oto),opos=oxto)
            else:
                Logger.debug("coc problem",2)
        
                
        
    def create_ose(self,iso,cycle,i):
        """
        """
        ose=ose_template(i)        
        iso=re.sub("[\[\]]","",iso).split("-") 
        closure=re.sub("[\[\]]","",cycle).split("-")
        ose.cyclebounds=[int(closure[0]),int(closure[1])]
        ose.set_form(int(closure[1])-int(closure[0])+1)        
        ose.set_nbc(len(iso[0]))
        ls_iso=list(iso[0])
        if len(iso)==1:
            c1=1
            ls_iso[c1-1]="0"
        else:
            c1= int(iso[1][0])
            if re.search("[0-5]b",iso[1]):            
                ls_iso[c1-1]="1"
            elif re.search("[0-5]a",iso[1]):
                ls_iso[c1-1]="2"
            else:
                ls_iso[c1-1]="0"   
        
        Lseries=ls_iso[int(closure[1])-1]=="1"    
        for c in range(0,len(ls_iso)):
            if re.search("[0-2]",ls_iso[c]):  
                if Lseries and c<int(closure[1])-1 and c>0:
                    if ls_iso[c]=="1":
                        ose.set_iso(c+1,IsomerCode["2"])
                    elif ls_iso[c]=="2":
                        ose.set_iso(c+1,IsomerCode["1"])
                    else:
                        ose.set_iso(c+1,IsomerCode[ls_iso[c]])
                else:
                    ose.set_iso(c+1,IsomerCode[ls_iso[c]])
            elif c in range(ose.cyclebounds[0]-1,ose.cyclebounds[1]):
                cyc=list(self.CarbonDescriptor_cyc.keys())
                for mod in cyc:
                    if re.search(mod,ls_iso[c]):
                        ose.set_mod(c+1,self.CarbonDescriptor_term[mod])                    
            elif ls_iso[c]!="h":
                termin=list(self.CarbonDescriptor_term.keys())
                termin.remove("h")                 
                for term in termin:
                    if re.search(term,ls_iso[c]):
                        ose.carbons[c]=self.CarbonDescriptor_term[term][0]
                        modterm=""                        
                        for elt in self.CarbonDescriptor_term[term][1:]:
                            if elt.isdigit():
                                modterm+=str(int(elt)-1)
                            else:
                                modterm+=elt
                                
                        ose.set_mod(c+1,modterm)
            
        return ose
        
    
    def modgatom(self,ose,iatom):
        """
        """
        g=ose.get_graphatom(iatom,False)
        imod=max(g.nodes())+iatom+1
       
        
        #for alin in ose.mods:
        for icarb in range(0,len(ose.mods)):
            cnum=icarb+1
            alin=ose.mods[icarb]
            main_nodes=[]
            if alin!="O": #OH is not a modif
                #cnum=ose.mods.index(alin)+1
                connector=-1
                remove=-1
                branches=alin.split("/")
                main=list(branches[0])
                for node in g.nodes():
                    if len(g.nodes[node])>0:
                        if g.nodes[node]["cnum"]==cnum:
                            if g.nodes[node]["atom"]=="C" and g.nodes[node]["substit"]==False:#ose.carbons[ose.mods.index(alin)]  :
                                connector=node
                                main_nodes.append(node)
                            elif  g.nodes[node]["atom"]==alin :
                                remove=node
                    else:                        
                        print("bizarre 1 noeud sans attribut %i dans l'ose %i, carbon? %i"%(node,ose.identifier,cnum))
                        
                if connector>=0 and alin!="":
                    for atom in main:                        
                        g.add_node(imod,atom=atom,in_cycle=False,cnum=cnum,ose=ose.identifier,substit=True)
                        g.add_edge(connector,imod)            
                        connector=imod
                        main_nodes.append(imod)
                        imod+=1
                        
                    for mod in branches[1:]:
                        b=re.findall("=?[0-9A-Za-z]",mod)
                        # *NCC/3=O => main=[NCC], indides=[0,1,2]                        
                        connector=main_nodes[int(b[0])-1]
                        for atom in b[1:]:                            
                            g.add_node(imod,atom=atom,in_cycle=False,cnum=cnum,ose=ose.identifier,substit=True)
                            g.add_edge(connector,imod)            
                            connector=imod
                            imod+=1
                            
                    
                if remove>0:
                    g.remove_node(remove)
        
        return g
            
   
                    
    # recursive decoding
    def __decode_52__(self,code,it=0,num=0):    
        """
        recursive coding the ascii letter(s) into a number
        """
        if it==len(code):        
            return num      
        else:
            # digit
            d=len(code)-it-1
            # char index
            c=string.ascii_letters.index(code[it])+1        
            return num+self.__decode_52__(code,it+1,(52**d)*c)
    
    
    
class WurcsEncoder:
    """
    """
    VERSION="2.0"
    """
    wurcs version
    """
                  
    def __init__(self,gatom,gose):
        """
        """
        self.dictwurcs_t={}
        self.dictwurcs_nt={}
        
        for elt in WurcsTerminalCarbs:
            self.dictwurcs_t[elt[1]]=elt[0]
        for elt in WurcsNonTermCarbs:
            self.dictwurcs_nt[elt[1]]=elt[0]    
            
        self.gatom=gatom
        self.gose=gose
        self.lsot=self.create_osetemplates()        
        self.notation="WURCS="+WurcsEncoder.VERSION     
                  
        self.encode_wurcs()
         
            
    def create_osetemplates(self):
        """
        """
        lsot=[]
        nb_type=0
        for o in self.gose.nodes():
            ot=ose_template(o)            
            ot.set_form(FormRef[self.gose.nodes[o]["form"]])
            ot.set_nbc(NbCRef[self.gose.nodes[o]["nbc"]])
            ot.gatom=nx.Graph()
            atoms=[]
            opos=[]
            cin=100
            cout=-1
            
            for e in self.gose.edges():
                if o in e:
                    opos.append(self.gose.edges[e]["opos"] )
            
            for n in self.gatom.nodes():
                if self.gatom.nodes[n]["ose"]==o:
                    # to extract subgraph of monosacch
                    atoms.append(n)                    
                    if re.sub("[0-9]","",self.gatom.nodes[n]["atom"])=="C" :
                        cnum=self.gatom.nodes[n]["cnum"]
                        if not self.gatom.nodes[n]["substit"]:
                            if "isomer" in self.gatom.nodes[n]:
                                ot.iso[cnum-1]=self.gatom.nodes[n]["isomer"]  
                            else:
                                ot.iso[cnum-1]=0 
                        if self.gatom.nodes[n]["in_cycle"]:
                            # get cycle bounds                                                    
                            if cnum>cout:
                                cout=cnum
                            if cnum<cin:
                                cin=cnum                        
                        
                        
            for obond in opos:
                if obond not in atoms:
                    # nodes where ose number is wrong (assigned to #ose from reducing end in gatom)
                    atoms.append(obond)
                    
            ot.cyclebounds=[cin,cout]             
            Lseries=ot.iso[cout-1]=="L"
            if Lseries:
                for cnum in range(cin+1,cout):
                    ot.inverse_iso(cnum)
                    
            if len(atoms)>0:
                ot.gatom=self.gatom.subgraph(atoms).copy() 
               
                
                lsot.append(ot)
                
        return lsot
    

    def __get_rankatoms__(self,lsatom,degrees):
        """
        """
        rank=list(range(len(lsatom)))
        for i in range(len(lsatom)):
            for j in range(i+1,len(lsatom)):
                
                mi=Atom.mass(lsatom[i])
                mj=Atom.mass(lsatom[j])
                ci=degrees[i]
                cj=degrees[j]
                if ci<cj:
                    rank[i]=j 
                    rank[j]=i
                    lsatom[j]=lsatom[i]
                    degrees[j]=degrees[i]
                elif ci==cj:
                    if mi>mj:
                        rank[i]=j
                        lsatom[j]=lsatom[i]
                        degrees[j]=degrees[i]
                        rank[j]=i
                                            
                     
        return rank
    
    
    def substit2alin(self,rootnode,ot,cnum,w="",count=0,uplevel="",ls=[]):
        """
        """
        neigh=list(ot.gatom[rootnode])     
                
        #level 0
        for n in ot.gatom[rootnode]:
            if n in ls:
                neigh.remove(n)
            elif not ot.gatom.nodes[n]["substit"]:
                neigh.remove(n)
           
       
        nbvoisin=len(neigh)    
        if nbvoisin==0:            
            w+=uplevel
            return w
        
        else:
            # sort neigh according to connectivity and atomic mass...
            lsatoms=[]
            lsdegrees=[]
            ls+=neigh
            for n in neigh:
                atom=re.sub("[0-9]","",ot.gatom.nodes[n]["atom"])
                degree=nx.degree(ot.gatom,n)
                if re.search("=",atom):
                    atom=re.sub("=","",atom)
                    lsatoms.append(atom)
                    degree+=1
                    
                    degree=Atom.valence(atom)-degree
                    lsdegrees.append(degree)
                else:
                    lsatoms.append(atom)
                    
                    degree=Atom.valence(atom)-degree
                    lsdegrees.append(degree)
            rank=self.__get_rankatoms__(lsatoms,lsdegrees)            
            sorted_neigh=[]
            for i in rank:
                sorted_neigh.append(neigh[i])
                
            #
            count+=1
            w+=self.gatom.nodes[sorted_neigh[0]]["atom"]            
                               
            for n in sorted_neigh[1:]:
                
                lvl2="/%i%s"%(count,ot.gatom.nodes[n]["atom"])
                uplevel+=lvl2
            return self.substit2alin(sorted_neigh[0],ot,cnum,w,count,uplevel,ls)        
        return w+uplevel
        
        
        
    # ose template to wurcs type of ose notation    
    def ot2to(self,ot):
        """
        ose template to wurcs type of ose notation    
        """
        w=""
        isomeri=""
        alins=""
        cnum_chiral=[]
        for node in ot.gatom.nodes():                   
            if re.sub("[0-9]","",ot.gatom.nodes[node]["atom"])=="C" and not ot.gatom.nodes[node]["substit"]:
                cnum=ot.gatom.nodes[node]["cnum"]
                
                # check if this carbon displays chirality
                neigh=list(ot.gatom[node])    
                
                if len(neigh)==3:                    
                    for neighnode in ot.gatom[node]:
                        if ot.gatom.nodes[neighnode]["in_cycle"]:
                            neigh.remove(neighnode)
                        elif re.search("=",ot.gatom.nodes[neighnode]["atom"]):
                            neigh.remove(neighnode)
                        
                            
                       
                    if ot.gatom.nodes[node]["in_cycle"]:
                        if len(neigh)==1 :
                            cnum_chiral.append(cnum)
                    elif nx.degree(self.gatom,node)>2 :
                        d=self.parseatoms(neigh,node)
                        equals=True
                        k=list(d.keys())
                        for i in range(len(k)):                            
                            for j in range(i+1,len(k)):
                                conxi=d[k[i]]
                                conxj=d[k[j]]
                                if len(conxi)==len(conxj):
                                    if set(sorted(conxi))==set(sorted(conxj)):
                                        pass
                                    else:                                        
                                        equals=False
                                        break
                                        
                                else:
                                    equals=False
                                    break
                        if not equals:
                            cnum_chiral.append(cnum)
                        
                        
                        
                # create alin notation
                if cnum!=ot.cyclebounds[1]:
                    substit=self.substit2alin(node,ot,cnum,"",0,"",[])                        
                    
                    if substit!="O":                                
                        ot.mods[cnum-1]=substit 
                        
        # monosaccharide wurcs encoding                    
        for cnum in range(1,len(ot.carbons)+1):
            # substituents linear notation
            alin=ot.mods[cnum-1]
            
           
            #carbon descriptors
            if cnum in range(ot.cyclebounds[0],ot.cyclebounds[1]):
                carbondesc=None
                if alin!="O":
                    for expr,code in self.dictwurcs_nt.items():
                        if re.search("^%s$"%expr,alin):
                            carbondesc=code
                            break
                    if carbondesc==None:
                        alins+="_%i*%s"%(cnum,alin)
                    
                if cnum==ot.cyclebounds[0]:
                    if alin=="O":
                        w+="a"     
                    else:
                        w+="U"
                elif cnum in cnum_chiral :
                    if ot.iso[cnum-1] in CodeIsomer.keys():
                    #if ot.iso[cnum-1]!=0 :
                        w+=CodeIsomer[ot.iso[cnum-1]]
                    else:
                        w+="x"
                else:
                    
                    if carbondesc==None:
                        w+="X"
                    else:
                        w+=carbondesc
               
            else: 
                carbondesc=None
                if alin!="O":
                    for expr,code in self.dictwurcs_t.items():
                        if re.search("^%s$"%expr,alin):
                            carbondesc=code
                            break
                    if carbondesc==None:
                        alins+="_%i*%s"%(cnum,alin)                
                
                
                if cnum in cnum_chiral: 
                    if ot.iso[cnum-1] in CodeIsomer.keys():                
                        w+=CodeIsomer[ot.iso[cnum-1]] 
                    else:
                        w+="x"
                else:
                    if cnum==ot.nbc():
                        if ot.get_mod(cnum)=="O":
                            w+="h"
                        elif ot.get_mod(cnum)=="":
                            w+="m"
                            #print("smkfdh"+str(ot.get_mod(cnum)))
                        else:
                            w+="h"
                    else:
                        if carbondesc==None:
                            w+="X"
                        else:
                            w+=carbondesc
                            
                       
        # cycle stuff    
        if ot.iso[ot.cyclebounds[0]-1]=="L":
            w+="-%ib_%i-%i"%(ot.cyclebounds[0],ot.cyclebounds[0],ot.cyclebounds[1])
        elif ot.iso[ot.cyclebounds[0]-1]=="D":
            w+="-%ia_%i-%i"%(ot.cyclebounds[0],ot.cyclebounds[0],ot.cyclebounds[1])
        else:
            w+="-%ix_%i-%i"%(ot.cyclebounds[0],ot.cyclebounds[0],ot.cyclebounds[1])
            
        return w+alins
        
    def encode_wurcs(self):
        """
        merge informations into a wurcs string
        """
        lstype={}
        monotypes=[]
        #parse gose
        ose_order=[self.gose.graph["reducing_end"]]
        self.parse_topo(self.gose.graph["reducing_end"],ose_order,[])        
        
        
        triot=[]
        for oid in ose_order:
            for ot in self.lsot:
                if ot.identifier==oid:
                    triot.append(ot)
        # ose types
        for ot in triot:
            w=self.ot2to(ot)
            lstype[ot.identifier]=w
            Logger.debug(w)
            if not w in monotypes:
                monotypes.append(w)
        self.notation+="/%i,%i,%i/"%(len(monotypes),len(self.gose.nodes()),len(self.gose.edges()))
        for to in monotypes:
            self.notation+="[%s]"%to
        self.notation+="/"    
        
        
        # list of monosacch types
        for identifier in ose_order:
            self.notation+="%i-"%(monotypes.index(lstype[identifier])+1)
        
        
        self.notation=self.notation[:-1]+"/"
        
        # list of osidic bindings
        for i in range(0,len(ose_order)-1):
            ofrom=ose_order[i]
            for j in range(i+1,len(ose_order)):                
                oto=ose_order[j]
                         
                try:
                    attr=self.gose[ofrom][oto]["bond"]
                    binding=re.split("[+\(\)]",attr)            
                    self.notation+="%s%s-%s%s_"%(self.encode_52(i),binding[1],self.encode_52(j),binding[2])
                except KeyError:
                    Logger.debug("invalid link %i->%i"%(ofrom,oto),0)
                except :                
                    Logger.debug("Unexpected error for link: %i->%i"%(sys.exc_info()[0],ofrom,oto),2)
                
        # final notation without trailing char "_"
        self.notation=self.notation[:-1]        
        Logger.debug(self.notation)
        
        
    
    def parse_topo(self,start,path,waiting):
        """
        """
        neigh=list(self.gose[start])
        
        ##precaution inutile avec DiGraph()
        
        # parcours        
        if len(neigh)>0:
            if len(neigh)==1:    
                path+=neigh
                self.parse_topo(neigh[0],path,waiting)
            else:
                carbs={}
                for o in neigh:
                    if "bond" in self.gose[start][o] and len(re.split("[\(\)+]",self.gose[start][o]["bond"]))==4:
                        carbon=re.split("[\(\)+]",self.gose[start][o]["bond"])                    
                        carbs[int(carbon[1])]=int(carbon[3])
                
                for i in range(1,len(carbs)):
                    waiting.append(carbs[list(carbs.keys())[i]])
                path.append(carbs[list(carbs.keys())[0]])                
                self.parse_topo(carbs[list(carbs.keys())[0]],path,waiting)
        elif len(waiting)==1:
            path+=waiting
            self.parse_topo(waiting[0],path,[])
        elif len(waiting)>1:
            path.append(waiting[0])
            self.parse_topo(waiting[0],path,waiting[1:])
        else:            
            #print(("path",path))
            Logger.debug(str(("path",path)),1)
    
    def parseatoms(self,nodes,exclu):
        """
        """
        dico={}
        for n in nodes:
            subgraph=self.__collect__(self.gatom,exclu,n,[n,exclu])    
            subgraph.remove(exclu)
            for i in range(len(subgraph)):
                subgraph[i]=self.gatom.nodes[subgraph[i]]["atom"]
            dico[n]=subgraph
        return dico
    
    def __collect__(self,G,start,node,ls):     
        """
        """
        for n in G[node]:                
            if n !=start and n not in ls:                          
                ls.append(n)
                self.__collect__(G,node,n,ls)
        return ls      
                    
    # firstly find the corresponding number of digits then decompose the number into letters (starts with 0='a')
    def encode_52(self,num):
        """
        firstly find the corresponding number of digits then decompose the number into letters (starts with 0='a')
        
        :type num: int
        """
        if num<52:
            return string.ascii_letters[num]
        else:
            i=0
            while int(num)/(52**i)>1:    
                i+=1
                
            s="" 
            
            for it in range(i-1,-1,-1):
                d=int(num/52**it)
                m=num%52**it    
                reste=0
                for j in range(0,it,1):
                    reste+=52**j
                    
                if m<reste:
                    d-=1
                    m=num-d*(52**it)
                s+=string.ascii_letters[d-1]
                num=m
            return s
    

class ose_template:    
    """
    Pattern of ose descriptions reading and writing information into an atom graph 
    """
    def __init__(self,ident):
        self.identifier=ident
        self.form=0
        self.carbons=[]
        self.iso=[]
        self.mods=[]
        self.cyclebounds=[-1,-1]
        self.gatom=None
    
    def set_iso(self,cnum,iso):
        """        
        """
        i=cnum-1
        if i in range(0,len(self.iso)):
            if iso in IsomerCode.keys():
                self.iso[i]=IsomerCode[iso]
            elif iso in IsomerCode.values() :
                self.iso[i]=iso
            else:
                self.iso[i]=""
                Logger.debug("isomerie non specified",1)
    
    def inverse_iso(self,cnum):
        iso=self.get_iso(cnum)
        if iso=="D":
            # verifier ok sur smiles
            self.set_iso(cnum,"L")
        elif iso=="L":
            self.set_iso(cnum,"D")
        elif  iso==1:
            self.set_iso(cnum,2)
        elif iso==2:
            self.set_iso(cnum,1)
            
    
                        
                                  
                
    
    def get_iso(self,cnum):
        """
        """
        i=cnum-1
        if i in range(0,len(self.iso)):
            return self.iso[i]
        else:
            return ""
        
    def set_mod(self,cnum,mod):
        """
        """
        i=cnum-1
        if i in range(0,len(self.mods)):
            self.mods[i]=mod
    
    def join_modgraph(self,subgraph,cnum):
        """
        links a substituent graph to the ose graph
        
        :param: subgraph: subgraph of substituant atoms
        """
        if self.gatom !=None:
            cnode=get_nodes(self.gatom,atom="C",cnum=cnum,substit=False)[0]
            Rnode=get_nodes(subgraph,atom="R")
            if len(Rnode)>0:
                neigh=list(subgraph[Rnode[0]])
                subgraph.remove_node(Rnode[0])
                
                for subnode in subgraph.nodes():
                    self.gatom.add_node(subnode,atom=subgraph.nodes[subnode]["atom"],substit=True,cnum=cnum,ose=self.identifier,in_cycle=False)
                    
                   
                for subnode in neigh:                                    
                    self.gatom.add_edge(cnode,subnode)
                
                for subedge in list(subgraph.edges()):
                    self.gatom.add_edge(*subedge)
                                
            else:                
                for subnode in list(subgraph.nodes()):                                    
                    self.gatom.add_node(subnode,atom=subgraph.nodes[subnode]["atom"],substit=True,cnum=cnum,ose=self.identifier,in_cycle=False)
                    
                for subdedge in list(subgraph.edges()):
                    self.gatom.add_edge(*subdedge)
                
                self.gatom.add_edge(cnode,subgraph.graph["startnode"])                
            
        
    def modsmi2graph(self,cnum,smi):
        """
        applies substitution to carbon in ose graph
        :param cnum: carbon number
        :param smi: SMILES notation of the subsituent
        """
        self.set_mod(cnum,smi)
        if self.gatom !=None:    
            
            startinode=max(list(self.gatom.nodes()))+1
            subGraph=None
            if smi=="":
               
                onodes=get_nodes(self.gatom,atom="O",cnum=cnum,substit=True)
                if len(onodes)>0:
                    self.gatom.remove_nodes_from(onodes)
                    Logger.debug("removing OH for ose %i,carbon number %i"%(self.identifier,cnum))
            else:
                if smi[0]=="(":
                    subGraph=smi2graph("R"+smi,subGraph,startinode)                    
                    
                    self.join_modgraph(subGraph,cnum)
                    
                elif re.match("^=?[A-Z][a-z]?$",smi):
                    # if O in osidic binding, node is removed from Gatom
                    cnode=get_nodes(self.gatom,atom="C",cnum=cnum,substit=False)[0]
                    self.gatom.add_node(startinode,atom=smi,substit=True,cnum=cnum,ose=self.identifier,in_cycle=False)                    
                    self.gatom.add_edge(cnode,startinode)
                    
                else:
                    subGraph=smi2graph(smi,subGraph,startinode)                    
                    
                    self.join_modgraph(subGraph,cnum)
                        
                      


    def get_mod(self,cnum):
        """
        """
        i=cnum-1
        if i in range(0,len(self.mods)):
            return self.mods[i]
        else:
            return None
    
    def set_cycle_start(self,cnum):
        """
        cnum=1 for aldose, cnum=2 for ketose, more ?
        """
        if cnum in range(1,len(self.carbons)+1):
            self.cyclebounds[0]=cnum
    
    def set_cycle_end(self,cnum):
        """
        """
        if cnum in range(1,len(self.carbons)+1):
            self.cyclebounds[1]=cnum
    
    def set_txtform(self,name):
        """
        """
        for num,ref in RefForm.items():
            if name==ref:
                self.set_form(num)
            
    def set_txtnbC(self,name):
        """
        :param name: hexose,pentose,etc...
        """
        for num,ref in RefNbC.items():
            if name==ref:
                self.set_nbc(num)
                
        
    def set_form(self,nbC):
        """
        :param nbC: number of carbon in the cycle
        :type nbC: int
        """
        self.form=nbC
        if self.cyclebounds==[-1,-1]:
            self.cyclebounds=[1,nbC]
            
            
    def set_nbc(self,nbc):
        """
        :param nbc: total number of carbon in the ose
        :type nbc: int
        """
        for i in range(0,nbc):
            self.carbons.append("C")
        for c in range(0,nbc):
            self.mods.append("O")
            self.iso.append("0")
    
    def nbc(self):
        """
        get the total number of carbons
        """
        return len(self.carbons)
    
    def get_txtform(self):
        """        
        """
        return RefForm[self.form]

    def get_txtnbc(self):
        """
        """
        return RefNbC[len(self.carbons)]
    
    def copy(self,ident):
        """
        """
        osecp=ose_template(ident)
        osecp.form=self.form
        osecp.carbons=self.carbons
        osecp.iso=self.iso.copy()
        osecp.mods=self.mods.copy()
        osecp.cyclebounds=self.cyclebounds.copy()
        return osecp
    
    # create graph of monosacch atoms (with or without oxygens)
    def get_graphatom(self,inode,wo):
        """
        builds graph of monosacch atoms (with or without oxygens)
        :return: gatom of ose (ose_template member)
        :rtype: nx.graph
        """        
        if self.gatom==None:
            gatom=nx.Graph()
            
            cycstart=-1
            cycend=-1
            carbnodes=[]
            for inum in range(0,len(self.carbons)):  
                cnum=inum+1
                gatom.add_node(inode,atom=self.carbons[inum],cnum=cnum,in_cycle=(cnum in range(self.cyclebounds[0],self.cyclebounds[1]+1)),ose=self.identifier,substit=False)
                
                if cnum==self.cyclebounds[0]:
                    cycstart=inode       
                if cnum==self.cyclebounds[1]:
                    cycend=inode
                else:
                    carbnodes.append(inode)
                    
                if inum in range(0,len(self.iso)):
                    gatom.nodes[inode]["isomer"]=self.iso[inum]
                
                if cnum>1:
                    gatom.add_edge(inode,inode-1) 
                
               
                    
                inode+=1                    
            
            # end with closing the cycle by 1 oxygen
                            
            gatom.add_node(inode,atom="O" ,cnum=cnum,in_cycle=False,ose=self.identifier,substit=True)                    
            gatom.nodes[inode]["in_cycle"]=True
            gatom.nodes[inode]["cnum"]=0
            gatom.nodes[inode]["substit"]=False
            
            gatom.add_edge(inode,cycstart)
            gatom.add_edge(inode,cycend)  
            inode+=1 
            
            # with oxygens
            if not wo:
                for node in carbnodes:  
                    cnum=gatom.nodes[node]["cnum"]                    
                    gatom.add_node(inode,atom=self.mods[cnum-1] ,cnum=cnum,in_cycle=False,ose=self.identifier,substit=True)
                    gatom.add_edge(node,inode)                  
                    inode+=1  
            self.gatom=gatom.copy()
        
        return self.gatom    
    
    def get_graph_cnum(self,gatom):
        #cycle as subsituent?
        subcyc=[]
        ose=None
        cycles=nx.cycle_basis(gatom,min(gatom.nodes()))           
        if len(cycles)>1:
            # test nb of double bonds
            for c in cycles:
                db=0
                for node in c:                
                    atom=gatom.nodes[node]["atom"]        
                    if re.search("=",atom):
                        db+=1
                if db>0:
                    subcyc.append(c)
                else:                    
                    ose=c
        else:
            ose=cycles[0]  
        
        if len(subcyc)>0:        
            gprime=gatom.copy()   
            dico_subcx={}
            
            rmsub=[]
            for isub in range(0,len(subcyc)-1):
                for jsub in range(isub+1,len(subcyc)):
                    inodes=set(subcyc[isub])
                    jnodes=set(subcyc[jsub])
                    inter=inodes.intersection(jnodes)
                    if len(inter)>0:
                        noninter=inodes.union(jnodes)-inter
                        maxpath=[]
                        bindnode=-1
                        for n in noninter:
                            p=nx.shortest_path(gprime,n,ose[0])
                            if len(p)>len(maxpath):
                                maxpath=p
                                bindnode=n
                        if bindnode in jnodes:
                            rmsub.append(subcyc[jsub])
                        elif bindnode in inodes:
                            rmsub.append(subcyc[isub])
                            #rmsub.append(lsubcyc[isub])
                        else:
                            print("error set OT graph ",maxpath)
                        
                   
            if len(rmsub)>0:
                for sub in rmsub:
                    if sub in subcyc:
                        subcyc.remove(sub)
                        
            for sub in subcyc:
                minpath=nx.shortest_path(gprime,min(sub),ose[0])                
                for n in sub:
                    path=nx.shortest_path(gprime,n,ose[0])
                    if len(path)<len(minpath):
                        minpath=path
                link=list(set(minpath) - set(sub))
                link=list(set(link) - set(ose))
                #print(link)
                for i in range(0,len(minpath)-1):  
                    if minpath[i] in sub and minpath[i+1] in link:                          
                        e=(minpath[i],minpath[i+1])
                        #print(e)
                        gprime.remove_edge(*e)
                                           
                        toremove=__collect4__(gprime,minpath[i],minpath[i],[minpath[i]])
                        dico_subcx[path[i+1]]=toremove.copy()                                                
                        gprime.remove_nodes_from(toremove)
                        
                       
                        break
                
        
            self.graph_cnum(gprime,ose)
            
            for node in gatom.nodes():
                if node in gprime.nodes():
                    gatom.nodes[node]["substit"]=gprime.nodes[node]["substit"]
                    gatom.nodes[node]["ose"]=gprime.nodes[node]["ose"]
                    gatom.nodes[node]["cnum"]=gprime.nodes[node]["cnum"]
                    if node in dico_subcx.keys():
                        for subnode in dico_subcx[node]:
                            gatom.nodes[subnode]["substit"]=True
                            gatom.nodes[subnode]["ose"]=gprime.nodes[node]["ose"]
                            gatom.nodes[subnode]["cnum"]=gprime.nodes[node]["cnum"]
        else:
            self.graph_cnum(gatom,ose)
             
    # assign ose description given atom graph (and store annotated graph)        
    def graph_cnum(self,gatom,cycle):        
        """
        assigns ose description given atom graph (and store annotated graph)        
        """
             
        oxnode=None        
        for node in gatom.nodes():
            gatom.nodes[node]["substit"]=False
            gatom.nodes[node]["ose"]=self.identifier
            if node in cycle:
                gatom.nodes[node]["in_cycle"]=True                
                if re.sub("[0-9]+","",gatom.nodes[node]["atom"])=="O":
                    oxnode=node
            else:
                gatom.nodes[node]["in_cycle"]=False
        # search hemi-acetalic node
        carbon1=[]
        
        for c in cycle:
            neigh=list(gatom[c])
            atoms=[]             
            for a in neigh:                                
                atom=re.sub("[0-9]+","",gatom.nodes[a]['atom'])                                
                atoms.append(atom)
                
            
            # OCO when C1 has no substitution on OH, OCOC when neuraminic acid...
            if len(neigh) in [3,4] and atoms.count("O")==2:                
                carbon1.append(c)
                
                    
            # OCN when amidation
            elif len(neigh)==3 and atoms.count("O")==1 and atoms.count("N")==1:
                carbon1.append(c)
            
        
        start_num_carb=1
        carbonEnd=[]
        minc=start_num_carb
        # transform cyclic to linear form for carbon counting
        if len(carbon1)==1 and oxnode!=None:            
            carbonEnd=list(gatom[oxnode])
            carbonEnd.remove(carbon1[0])
            atom_oxnode=gatom.nodes[oxnode]
            gatom.remove_node(oxnode)
            
            # cyclic nodes
            numc=start_num_carb
            for c in nx.shortest_path(gatom,carbon1[0],carbonEnd[0]):                
                self.__parse_noncyclic__(gatom,numc,c,[],0)  
                numc+=1
                
            
            # terminal carbs
            neigh=list(gatom[carbon1[0]])
            for n in neigh:
                if re.sub("[0-9]+","",gatom.nodes[n]["atom"])=="C":
                    if gatom.nodes[n]["in_cycle"]==False:
                        self.__parse_noncyclic__(gatom,gatom.nodes[carbon1[0]]["cnum"],n,[carbon1[0]],-1)                    
                else:
                    gatom.nodes[n]["cnum"]=gatom.nodes[carbon1[0]]["cnum"]
                    
            neigh=list(gatom[carbonEnd[0]])
            for n in neigh:
                if re.sub("[0-9]+","",gatom.nodes[n]["atom"])=="C":
                    if gatom.nodes[n]["in_cycle"]==False:
                        self.__parse_noncyclic__(gatom,gatom.nodes[carbonEnd[0]]["cnum"],n,[carbonEnd[0]],1)                    
                else:
                    gatom.nodes[n]["cnum"]=gatom.nodes[carbonEnd[0]]["cnum"]
                    
            # renumber for "real" backbone start
            for i in gatom.nodes():
                if gatom.nodes[i]["cnum"]<minc:
                    minc=gatom.nodes[i]["cnum"]
            if minc<start_num_carb:
                for i in gatom.nodes():
                    gatom.nodes[i]["cnum"]=gatom.nodes[i]["cnum"]+(start_num_carb-minc)
            # re-cycle
            gatom.add_node(oxnode)
            atom_oxnode["cnum"]=0
            for k,v in atom_oxnode.items():
                gatom.nodes[oxnode][k]=v
                
            gatom.add_edge(*(oxnode,carbon1[0]))  
            gatom.add_edge(*(oxnode,carbonEnd[0]))  
            
        else:               
            Logger.debug("can not find hemi-acetal function")
        
        
        self.gatom=gatom.copy()
        cnums=[]
        nbc=0
        dct_iso={}
        for n in self.gatom.nodes():
            if self.gatom.nodes[n]["cnum"]>nbc:
                nbc=self.gatom.nodes[n]["cnum"]
        self.set_nbc(nbc)        
        
        for c in cycle:
            if re.sub("[0-9]+","",gatom.nodes[c]["atom"])=="C":
                cnums.append(gatom.nodes[c]["cnum"])
            
            if "isomer" in gatom.nodes[c]:      
                cnum=gatom.nodes[c]["cnum"]
                isomer=gatom.nodes[c]["isomer"]
                if isomer in IsomerCode.keys():
                    dct_iso[cnum]=IsomerCode[isomer]
                elif isomer in IsomerCode.values():
                    dct_iso[cnum]=isomer
                    
                else:
                    dct_iso[cnum]=""
                    
        self.cyclebounds=[min(cnums),max(cnums)]
        self.set_form(max(cnums)-min(cnums)+1)
        for k,v in sorted(dct_iso.items()):
            self.iso[k-1]=v
      

        
      
            
    def __parse_noncyclic__(self,g,cnum,node,ls,incr):        
        """
        """
        ls.append(node) 
        
        #if not cyclic should not have number from smiles...!
        if re.sub("[0-9]+","",g.nodes[node]["atom"])=="C":
            cnum+=incr
        else:
            incr=0
           
        g.nodes[node]["cnum"]=cnum          
        # tag node as susbstit when the atom is not a carbon belonging to monosacch backbone
        g.nodes[node]["substit"]=(incr==0 and not g.nodes[node]["in_cycle"])
            
        if len(ls)<100:
            for n in g[node]:            
                if n not in ls and g.nodes[n]["in_cycle"]==False: 
                    self.__parse_noncyclic__(g,cnum,n,ls,incr)   
            return ls      
        else:
            return ls 
        
    # for OsidicBond in graphs
    def get_co_nodes(self,cnum):        
        """
        for OsidicBond in graphs
        """
        cnode,onode=-1,-1
        if self.gatom!=None:
            for n in self.gatom.nodes():
                descr=self.gatom.nodes[n]
                if descr["atom"]=="C" and descr["cnum"]==cnum and descr["substit"]==False:
                    cnode=n
                    neigh=list(self.gatom[n])
                    for nn in neigh:
                        ndescr=self.gatom.nodes[nn]
                        if ndescr["atom"]=="O" and ndescr["cnum"]==cnum and ndescr["substit"]==True:
                            onode=nn
                    break
        
        return cnode,onode
    
    def get_cyclenodes(self):
        """
        """
        if self.gatom!=None:
            return nx.cycle_basis(self.gatom,min(self.gatom.nodes()))
        
    def getcyclebound_nodes(self):
        """
        :return: list of node numbers for start and end cycle
        """
        nodes=[]
        for i in self.gatom.nodes():
            if self.gatom.nodes[i]["in_cycle"]==True:
                if self.gatom.nodes[i]["cnum"]==self.cyclebounds[0]:
                    nodes.append(i)
                elif self.gatom.nodes[i]["cnum"]==self.cyclebounds[1]:
                    nodes.append(i)
        return nodes
    
       
    def get_nodes(self,**kwargs):
        """
        search type of nodes
        
        :param kwargs: dictionnary of node attributes (ose=...)
        
        :return: a list of node numbers
        """
        nodes=[]
        if self.gatom!=None:
            for n in self.gatom.nodes():
                match=True
                for k,v in kwargs.items():
                    if k in self.gatom.nodes[n] and v==self.gatom.nodes[n][k]:
                        match&=True
                    else:
                        match&=False
                if match:
                    nodes.append(n)
        return nodes
    
    def equals(self,ot2):  
        """
        """
        patsmi="(\(?[=A-Z][a-z]?[0-9]?\)?)"
        patform="([A-Z][a-z]?[0-9])"
        patalin="([A-Z]+/[0-9][=A-Z]+[a-z]?)"
        
        if self.nbc!=ot2.nbc:
            return False
        elif self.form!=ot2.form:
            return False
        elif self.cyclebounds!=ot2.cyclebounds:
            return False
        elif self.iso!=ot2.iso:
            return False        
        elif self.mods!=ot2.mods:
            
            for isub in range(len(self.mods)):
                sub1=self.mods[isub]
                if isub>=len(ot2.mods):
                    return False
                else:
                    sub2=ot.mods[isub]
                    if re.search(patsmi,sub1):
                        
                        if sorted(sub1)!=sorted(sub2):
                            return False
                    elif re.search(patform,sub1):
                        if sorted(re.findall(patform,sub1))!=sorted(re.findall(patform,sub2)):
                            return False
                    elif sub1!=sub2:
                        return False            
            return True
        else:
            return True

    def add_anhdydro(self,cnum=(3,6)):
        if self.gatom==None:
            Logger.debug("nothing to do in add_anhydro")
        else:
            # remove H2O between (C3,C6)
            c3=list(get_nodes(self.gatom,atom="C",in_cycle=True,substit=False,cnum=cnum[0]))[0]            
            c6=list(get_nodes(self.gatom,atom="C",substit=False,cnum=cnum[1]))[0]
            
            neigh3=self.gatom[c3]
            neigh6=self.gatom[c6]
            o3=None
            o6=None
            for n in neigh3:
                if self.gatom.nodes[n]["atom"]=="O":
                    if self.gatom.nodes[n]["substit"]==True:
                        if self.gatom.nodes[n]["in_cycle"]==False:
                            if self.gatom.nodes[n]["cnum"]==cnum[0]:
                                o3=n
                  
                         
            for n in neigh6:
                if self.gatom.nodes[n]["atom"]=="O":
                    if self.gatom.nodes[n]["substit"]==True:
                        if self.gatom.nodes[n]["in_cycle"]==False:
                            if self.gatom.nodes[n]["cnum"]==cnum[1]:
                                o6=n            
            
            if o3!=None and o6!=None:
                self.gatom.remove_node(o3)
                self.gatom.add_edge(c3,o6,link="anhydro")
                
                
                
class Topograph:
    prefix_modif="m"
    RefForm={4:'furane',5:"pyrane"}
    RefNbC={5:"pentose",6:"hexose",9:"non"}  
    
    def __init__(self):                
        self.Gatom=nx.Graph()
        # should be directed graph according to reducing end
        self.Gose=nx.DiGraph()  
        self.Gatom.graph['reducing_end']=1
    
    def topo2graph(self):     
        """
        topo2graph
        """
        self.Gatom=nx.Graph()
        # should be directed graph according to reducing end
        self.Gose=nx.DiGraph()  
        self.Gatom.graph['reducing_end']=1
        
        if self.Gatom.size()==0 or self.Gose.size()==0:        
            
            self.OSE_LIST=list(SingletonTopo.OSE_LIST.values())
            root=self.OSE_LIST[0].oseid      
            self.lsB=SingletonTopo.OBOND_LIST
            

            inode=0
            snode=0

            
            templates={}

            for o in self.OSE_LIST:  
                ot=ose_template(o.oseid)
                ot.set_form(o.ncc)
                ot.set_nbc(o.nct)
                ot.set_cycle_start(o.startcycle)
                ot.set_cycle_end(o.startcycle+o.ncc-1)
                self.Gose.add_node(o.oseid,form=ot.get_txtform(),nbc=ot.get_txtnbc(),anhydro=o.anhydro)                
                for c in range(ot.cyclebounds[0],ot.cyclebounds[1]+1):
                    ot.set_iso(c,o.modifs[0][c-1])#?c
                
                # init carbon graph of ose   
                ot.get_graphatom(inode,True) 
                
                # bind substituants to carbons
                for i,mod in enumerate(o.modifs[1]):
                    #ot.set_mod(i+1,SubstitutionLibrary.SUBSTITUTIONS[mod].smiles)
                    if mod>0 :  
                        smi=SubstitutionLibrary.getSub(mod).smiles
                        
                        #if not SubstitutionLibrary.check_smiles(mod):
                            
                            ##compo=Composition(SubstitutionLibrary.getSub(mod).formula,SubstitutionLibrary.getSub(mod).link,o.in_cycle(i+1))
                        
                            
                            ##smi=compo.smiles(compo.start,[])                            
                        
                            #if compo.check_H()!=0:
                                #Logger.debug("Wrong number of H? formula %s gives smiles '%s'"%(SubstitutionLibrary.getSub(mod).formula,smi),2)                                                        
                        
                        # mod should be =-1 for end cycle?
                        if i+1!=ot.cyclebounds[1]:                        
                            ot.modsmi2graph(i+1,smi)
                    
                
                if o.anhydro:
                    ot.add_anhdydro()
                
                self.Gatom.add_nodes_from(ot.gatom.nodes(data=True))
                self.Gatom.add_edges_from(ot.gatom.edges(data=True))

                inode=max(self.Gatom.nodes())+1
                
                templates[o.oseid]=ot
               
            #redend=templates[1].getcyclebound_nodes()[0]  
            redends=list(templates.keys())
           
            # osidic backbone, 1st, numbering oses and carbon in cycles               
            for b in self.lsB:                
                # parent is reducing end side
                osenum_parent=b.parent_ose
                cparent=b.parent_carbon                    
                cchild=b.child_carbon                        
                osenum_child=b.child_ose

                inode_parent=-1
                inode_child=-1
                removeO=-1
                coparent=templates[osenum_parent].get_co_nodes(cparent)
                cochild=templates[osenum_child].get_co_nodes(cchild)
                
                # si la liaison était TJS correctement orientée, il suffirait de mettre redends.remove(osenum_end)
                if cparent==templates[osenum_parent].cyclebounds[0]:
                    redends.remove(osenum_parent)
                if cchild==templates[osenum_child].cyclebounds[0]:
                    redends.remove(osenum_child)
               
                # parent keeps O (co=tuple carbon,oxygen nodes)
                inode_parent=coparent[1]
                removeO=cochild[1]
                inode_child=cochild[0]

                if removeO>0:
                    self.Gatom.remove_node(removeO)
                    
                if inode_parent>=0 and inode_child>=0:
                    # get oxygen nodes
                    
                    self.Gatom.add_edge(inode_parent,inode_child)                                                
                    b.set_opos(inode_parent)
                    
                else:
                    # nodes have not been create when topoOse.modif[cnum]<-1 (indicates a binding in oseModel)
                    if inode_parent<0 :
                        inode_parent=max(list(self.Gatom.nodes()))+1
                        self.Gatom.add_node(inode_parent,atom="O",cnum=cparent,ose=osenum_parent,substit=True,in_cycle=False)
                        self.Gatom.add_edge(inode_parent,coparent[0])
                        self.Gatom.add_edge(inode_parent,inode_child) 
                        b.set_opos(inode_parent)
                        
                            
                    else:                        
                        Logger.debug("problem",2)

            
            if len(redends)==1:
                root=redends[0]
                self.Gose.graph["reducing_end"]=root
                self.Gatom.graph['reducing_end']=templates[root].getcyclebound_nodes()[0]  
            
            else:
            
                Logger.debug("C1 node for reducing end not found",1)
            
            
                
            #Ose graph
            self.__buildGograph__(root,[])    
            
            #nx.write_graphml(self.Gatom,"/tmp/topo2gatom.graphml")
            

       

    

    def plotGose(self,ax=None):
        """
        """        
        coloromap={'pentose':'salmon','hexose':'lightblue'}
        colorose=[]
        for io in self.Gose.nodes():
            nct=self.Gose.nodes[io]['nbc']
            
            if nct in coloromap.keys():
                colorose.append(coloromap[nct])
            else:
                colorose.append('grey')

        pos = graphviz_layout(self.Gose, prog='dot')
        nx.draw(self.Gose,pos,ax=ax)
        nx.draw_networkx_nodes(self.Gose,pos,node_color=colorose)
        nx.draw_networkx_labels(self.Gose,pos,font_size=14)
        edge_labels = nx.get_edge_attributes(self.Gose,'bound')        
        nx.draw_networkx_edge_labels(self.Gose, pos,edge_labels)                       


        if ax==None:
            plt.show()    


    def getNode(self,**kwargs):
        """
        """        
        results=[]
        target=kwargs
        #tocompare={'atom':typeatom,'ose':osenum,'cnum':cnum,'in_cycle':in_cycle}        
        for n in self.Gatom.nodes():            
            descr={}
            for attr in target.keys():
                if attr=="atom":
                    descr[attr]=re.sub("[0-9]*","",self.Gatom.nodes[n][attr])
                else:
                    descr[attr]=self.Gatom.nodes[n][attr]
            if target==descr:
                results.append(n)
        return results


    # without isomeric information
    def __cycle_graph__(self,start,o):  
        """
        """        
        inode=start
        close=False
        isomers=o.modif[0]
        for c in range(1,o.carbon_number+1):
            inode+=1
            self.Gatom.add_node(inode,atom='C',ose=o.ose_number,cnum=c,in_cycle=(c<o.carbon_number))
            if (c-1)<len(isomers):
                self.Gatom.nodes[inode]["isomer"]=isomers[c-1]

            if inode>start+2:
                self.Gatom.add_edge(inode,inode-2)

            inode+=1
            self.Gatom.add_node(inode,atom='O',ose=o.ose_number,cnum=c,in_cycle=False)
            self.Gatom.add_edge(inode-1,inode)
            # ring size
            if o.ose_form=='furane' and c==4:
                close=True
            elif o.ose_form=='pyrane' and c==5:
                close=True
            else:
                close=False

            if close:
                self.Gatom.add_edge(start+1,inode)            
                self.Gatom.nodes[inode]["cnum"]=0
                self.Gatom.nodes[inode]["in_cycle"]=True            

        
        return inode
        
   
   

    
    

    def __buildGograph__(self,root,lsB):
        """
        """        
        reste=list(set(self.lsB)-set(lsB))
        if len(reste):            
            for bond in reste:
                if bond.parent_ose==root:
                    self.Gose.add_edge(bond.parent_ose,bond.child_ose,bond=bond.getAttributString(),opos=bond.onode)
                    lsB.append(bond)
                    self.__buildGograph__(bond.child_ose,lsB)  
                
   
    
    def graph2topo(self):    
        """
        """
        Edges=[]
        Oses=[]
        print(len(self.Gatom))
        if len(self.Gatom)==0 or len(self.Gose)==0:
            # bye bye
            return None,None
        else:
            
            # pff... on a pas le droit de commencer à 0
            if min(self.Gose.nodes())==0:
                start_osecount=1
            else:
                
                start_osecount=0
            
            self.anhydro_edge=[]
            for edge in self.Gatom.edges():
                if "link" in self.Gatom.edges[edge] and self.Gatom.edges[edge]["link"]=="anhydro":
                    self.anhydro_edge.append(edge)
            for edge in self.anhydro_edge:
                self.Gatom.remove_edge(*edge)
    
            # (osenum,cnum):modif name
            dico_modif={}
            dico_isomer={}
            # lesmodifs et l'isomerie ICI!
            for atomnode in self.Gatom.nodes():
                descr=self.Gatom.nodes[atomnode]         
                osenum=descr['ose']
                cnum=descr["cnum"]   
                
                if "isomer" in descr:
                    if osenum in dico_isomer:
            
                        dico_isomer[osenum].append((str(cnum),descr["isomer"]))
                    else:
                        dico_isomer[osenum]=[(str(cnum),descr["isomer"])]
            
    
                # subsituents on carbons in cycle or on last ose carbon            
                elt=re.sub("[0-9]+","",descr["atom"])
                subformul={}
                
                if elt=="C" and not descr["substit"]:
                    # should have hydrogens at least
                    
                    if osenum+start_osecount in dico_modif.keys():
                        dico_modif[osenum+start_osecount].append((cnum,SubstitutionLibrary.NOSUBID))
                    else:
                        dico_modif[osenum+start_osecount]=[(cnum,SubstitutionLibrary.NOSUBID)]
                    
                    # not H
                    clinks=[]
                    substit=self.getSubstitutGraph(cnum,osenum)
                    
                    if substit!=None:
                        subformul=self.formula(substit,-1)
                        for neigh in self.Gatom[atomnode]:                       
                    
                            if neigh in substit.nodes():
                                smi=Topograph.__simple__(substit,neigh,-1)
                                if smi==None:
                                    smi=Topograph.__NOTsimple__(substit,neigh,-1)
                                    clinks.append([subformul,smi]) 
                                else:
                                    clinks.append([subformul,smi]) 
                    else:
                        anh=False
                        # anhydro?
                        for edge in self.anhydro_edge:
                            if atomnode in edge:
                                anh=True
                                break
                            
                        if not anh:
                            #desoxy?
                            subformul["H"]=1                            
                            for neigh in self.Gatom[atomnode]:
                                # modify subformula to say there is no modif in such cases:
                                # end cycle
                                if self.Gatom.nodes[neigh]["cnum"]==0 and cnum>2:
                                    subformul["O"]=1
                                    subformul["H"]=1
                                # carbon in osidic binding
                                if self.Gatom.nodes[neigh]["atom"]=="O" and self.Gatom.nodes[neigh]["cnum"]!=cnum:
                                    subformul["O"]=1                  
                                    subformul["H"]=1
                           
                    
                                
                    if len(clinks)>4:
                        Logger.debug("ose %i carbon atom %i has too much atomic binding: %i"%(osenum,cnum,len(clinks)))
                    
                    smiles=""               
                    for notH in clinks:   
                        if smiles!="":# ? ça risque pas d'arriver...
                            smiles="("+smiles+")"+notH[1]                     
                        else:
                            smiles=notH[1]
                       
                    
                            
                    #print((osenum,cnum,subformul))        
                    ref=SubstitutionLibrary.get_substitution(subformul)
    
                    
                    if ref!=None and ref.name=="keto" : # ask SUBSTITUTIONS if there is constraint of cnum
                        if cnum==1:
                            dico_modif[osenum+start_osecount].append((cnum,ref.identifier))
                    
                        else:
                            ref=None
    
                    if ref!=None:
                        dico_modif[osenum+start_osecount].append((cnum,ref.identifier))
                       
                    elif len(subformul)>0: 
                        link=None
                        if smiles[0]=="(":
                            link=""
                            
                                                 
                        
                        newsub=Substitution(self.__formul2txt__(subformul),link,None,smiles)
                        SubstitutionLibrary.add_substit(newsub)
                        
                        
    
                        dico_modif[osenum+start_osecount].append((cnum,newsub.identifier))
                        
                        #print(dico_modif[osenum+start_osecount])
    
                        
    
            
            for onode in self.Gose.nodes():
                cycc=[]
                for n in self.getNode(ose=onode,in_cycle=True,atom="C"):
                    cycc.append(self.Gatom.nodes[n]["cnum"])
                anh=False
                for edge in self.anhydro_edge:
                    if self.Gatom.nodes[edge[0]]["ose"]==onode:
                        anh=True
                    
                
                ose_model=OseModel(FormRef[self.Gose.nodes[onode]["form"]],NbCRef[self.Gose.nodes[onode]["nbc"]],min(cycc))  
                ose_model.oseid=onode
                ose_model.set_anhydrobond(anh)
                for carbsub in dico_modif[onode+start_osecount]:
                    ose_model.set_modcarb(carbsub[0],carbsub[1])
                    
                
                if onode in dico_isomer:                    
                
                    for carb in dico_isomer[onode]:
                        ose_model.set_isocarb(int(carb[0]),carb[1])
                    
                    
                   
                  
                Oses.append(ose_model)
    
    
            if "reducing_end" in self.Gose.graph.keys():         
                start_ose=self.Gose.graph["reducing_end"]         
                self.get_oseedges(start_ose,Edges,start_osecount)
            else:
                Logger.debug("pas d'ose start",2)
            return Oses,Edges
            
    def __formul2txt__(self,f):  
        """
        """        
        txt=""
        for elt,nb in f.items():
            txt+=elt+str(nb)
        return txt

    def __NOTsimple__(G,node,previous):        
        gprime=G.copy()
        cycnum=0
        #nx.write_graphml(gprime,"/tmp/gprime.graphml") 
        subcycles=nx.cycle_basis(gprime)
        edge_cycles={}
        for n in gprime.nodes():
            atom=gprime.nodes[n]["atom"]
            nums=re.findall("[0-9]",atom)
            if len(nums)>0:
                gprime.nodes[n]["atom"]=re.sub("[0-9]","",atom)    
                for num in nums:                    
                    if num in edge_cycles.keys():
                        edge_cycles[num].append(n)
                    else:
                        edge_cycles[num]=[n]
        for cyc in edge_cycles.keys():
            cycnum+=1
            edge=edge_cycles.get(cyc)
            if len(edge)==2:                
                gprime.remove_edge(*edge)
                gprime.nodes[edge[0]]["atom"]+=str(cycnum)
                gprime.nodes[edge[1]]["atom"]+=str(cycnum)
            else:
                print("uncycling problem in __NOTsimple__")
           
            
        #nx.write_graphml(gprime,"/tmp/gprime.graphml")    
        return Topograph.__simple__(gprime,node,previous)
    
    # convert subsituent graph into SMILES notation (assuming that no cycle is present)
    @staticmethod
    def __simple__(G,node,previous):
        """
        """        
        if len(nx.cycle_basis(G))>0:
            return None
        smi=G.nodes[node]["atom"]
        neigh=list(G[node])
        if previous in neigh:
            neigh.remove(previous)    
        if len(neigh)>0:
            for n in neigh[:-1]:            
                smi+="("+Topograph.__simple__(G,n,node)+")"
            smi+=Topograph.__simple__(G,neigh[-1:][0],node)      

        return smi     


    def get_oseedges(self,start,edges,start_osecount):
        """
        """        
        neigh=self.Gose[start]
        if len(neigh)>0:         
            for ose in neigh:
                bond_attr=self.Gose.edges[(start,ose)]["bond"]
                bond=OsidicBond.getBondFromAttribute(bond_attr)
                o1=int(bond[0])+start_osecount
                c1=int(bond[1])
                c2=int(bond[2])
                o2=int(bond[3])+start_osecount

                edges.append([o1,c1,o2,c2])   
                self.get_oseedges(ose,edges,start_osecount)

    def formula(self,subgraph,start_nh):
        """
        """        
        form=""
        nh=start_nh
        dico_nbatom={}
        paths=list(nx.connected_components(subgraph))
        if len(paths)>1:
            nh=-len(paths)
        #if len(subgraph.nodes())-1>len(subgraph.edges()):
            ## nodes should be inter-connected unless the substit is directly linked to ose carbon
            #nh-=1 
            
            
        for n in subgraph.nodes():
            atom=re.sub("[0-9]","",subgraph.nodes[n]['atom'])
            if atom.find("=")>-1:
                for path in paths:
                    if n in path:
                        if n==sorted(list(path))[0]:
                            nh-=1
                        
                        else:
                            nh-=2
                atom=re.sub("=","",atom)                
                #nh-=2

            if atom in dico_nbatom.keys():
                dico_nbatom[atom]+=1
            else:
                dico_nbatom[atom]=1


           
            nh+=Atom.valence(atom)-nx.degree(subgraph,n) 
        if nh>0:
            dico_nbatom['H']=nh
        return dico_nbatom


    def getSubstitutGraph(self,cnum,osenum):
        """
        """        
        
        subnodes=get_nodes(self.Gatom,ose=osenum,cnum=cnum,substit=True)
        
        
        if len(subnodes)>0:
            return self.Gatom.subgraph(subnodes).copy()
        else:
            return None
        
   

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

    def get_osetype(self,onode):
        descr=self.Gose.nodes[onode]
        otype=""
        if "nbc" in descr:
            otype+=descr["nbc"][0:3]
        if "form" in descr:
            otype+=str(descr["form"][0]).upper()
        return otype

    def get_oseattr(self,ose):
        """
        """        
        descr={}
        descr["form"]=ose.ose_form
        descr["nbc"]=self.RefNbC[ose.carbon_number]
        isomer=""      
        for c in range(ose.carbon_number-1):
            isomer+=ose.modif[0][c]
        descr["isomer"]=isomer
        return descr

    def joinformula(self,form1,form2):
        """
        """        
        formula={}
        for k,v in form1.items():
            if k in form2.keys():
                formula[k]=form1[k]+form2[k]
            else:
                formula[k]=form1[k]

        for k,v in form2.items():
            if k not in form1.keys():
                formula[k]=form2[k]

        return formula
    

    
    def get_ostart(self):
        """
        """        
        if self.Gatom.size()>0 and "reducing_end" in self.Gatom.graph.keys():
            return self.Gatom.nodes[self.Gatom.graph['reducing_end']]["ose"]
        else:
            return 1
    
    def get_nbcycles(self):
        return len(nx.cycle_basis(self.Gatom))

    def get_nboses(self):
        return len(self.Gose.nodes())
        
# tools to convert formula to graph to smiles
class Composition:
    """
    tools to convert formula to graph to smiles
    """
    light=2
    # formula is text
    def __init__(self,formula,link,cin=True):
        self.Gatom=nx.Graph()
        self.groupatom=re.findall("[A-Z][a-z]?[0-9]+",formula)
        self.nH=0
        self.start=None
        self.dicovalence={}
        self.cin=cin
        i=0
        # create nodes
        for elt in self.groupatom:
            atom=re.sub("[0-9]","",elt)
            nb=int(re.sub("[^0-9]","",elt))
            if atom=="H":
                self.nH=nb
            else:     
                val=Atom.valence(atom)
                self.dicovalence[atom]=val
                for n in range(nb):
                    self.Gatom.add_node(i,atom=atom,nh=val,valence=val)                
                    i+=1    
                    
        # create edges
        self.setStart(link)
        
        if self.start!=None:
            self.connectHeavy(self.start)
            self.connectLight()
            self.fill_nh(max(self.dicovalence.values()))
    
    def check_H(self):
        return self.nbh_graph()-self.nH
        
    def smiles(self,start,done): 
        """
        creates recursively the smiles notation from a start node of gatom
        """
        neigh=list(self.Gatom[start])          
        for node in done:
            if node in neigh:
                neigh.remove(node)
        done.append(start)
        smi=self.Gatom.nodes[start]["atom"]
        if len(neigh)>0:
            neigh=sorted(neigh,key=lambda n:self.Gatom.nodes[n]["valence"])  
                     
            for i in range(len(neigh)-1):
                smi+="("+self.smiles(neigh[i],done)+")"
            smi+=self.smiles(neigh[-1],done)
        
        return re.sub("R","",smi)
      
    def setStart(self,atomlink):
        """
        
        """
        #print(atomlink)
        nodename=""
        if len(self.Gatom.nodes())>0:
            if atomlink=="":
                inode=len(self.Gatom.nodes())
                if self.cin:
                    self.Gatom.add_node(inode,atom="R",valence=10,nh=2)
                else:
                    self.Gatom.add_node(inode,atom="R",valence=10,nh=3)
                self.start=inode
            elif atomlink[0]=="=":
                nodename==re.sub("=","",atomlink)
                candids=get_nodes(self.Gatom,atom=nodename)
                if len(candids)>0:
                    self.start=candids[0]
                    self.Gatom.nodes[self.start]["atom"]=link
                else:
                    self.start=0                
            else:                
                candids=get_nodes(self.Gatom,atom=atomlink)
                if len(candids)>0:
                    self.start=candids[0]   
                else:
                    self.start=0
        else:
            self.start=0
    
            
    
    def nbh_graph(self):
        """
        compute the number of free Hydrogens in the atom graph
        """
        nbh=0
        for n in self.Gatom.nodes():            
            nbh+=self.Gatom.nodes[n]["nh"]
        return nbh
    
    def connectHeavy(self,nodestart):
        """        
        """
        heavy_nodes=[]
        
        for node in self.Gatom.nodes():
            if  self.Gatom.nodes[node]["valence"]>Composition.light and node!=self.start:                 
                    heavy_nodes.append(node)
                
        for n in heavy_nodes:
            self.bind_simple(nodestart,n)
            nodestart=n
    
    def bind_simple(self,node1,node2):   
        """
        """
        self.Gatom.nodes[node1]["nh"]-=1
        self.Gatom.nodes[node2]["nh"]-=1
        self.Gatom.add_edge(node1,node2)   
        #print((node1,node2))
    
    def connectLight(self,n=None):
        """
        """
        if n==None:    
            petit=sorted(self.Gatom.nodes(),key=lambda n:self.Gatom.nodes[n]["valence"])
            for node in petit:
                if  self.Gatom.nodes[node]["valence"]<=Composition.light and nx.degree(self.Gatom,node)==0:
                    self.connectLight(node)
        
        else:
            nodeH=None
            gros=sorted(self.Gatom.nodes(),key=lambda n:self.Gatom.nodes[n]["valence"],reverse=True)
            gros.remove(n)
            for node in gros:
                
                if  self.Gatom.nodes[node]["nh"]>=1:                    
                    nodeH=node
                    break
            if nodeH!=None:
                self.bind_simple(nodeH,n)
    
    def fill_nh(self,valence,gros=None):
        """
        tries to double bind atoms until the nb of H in formula is achieved
        """
        gnh=self.nbh_graph()
        if gnh>self.nH:
            if gros==None:
                if self.Gatom.nodes[self.start]["atom"]=="R":
                    self.fill_nh(valence,self.start)
                else:
                    for node in self.Gatom.nodes():
                        if self.Gatom.nodes[node]["valence"]==valence and self.Gatom.nodes[node]["nh"]>=1:
                            self.fill_nh(valence,node)
            else:
                neigh=list(self.Gatom[gros])                
                neigh=sorted(neigh,key=lambda n:self.Gatom.nodes[n]["valence"])
                if self.start in neigh and self.Gatom.nodes[self.start]["nh"]<2:
                    neigh.remove(self.start)
                ndb=(gnh-self.nH)//2
                if len(neigh)<=ndb:
                    for n in neigh:
                        if self.Gatom.nodes[n]["nh"]>=1:
                            self.Gatom.nodes[n]["atom"]="="+self.Gatom.nodes[n]["atom"]
                            self.Gatom.nodes[n]["nh"]-=1
                            self.Gatom.nodes[gros]["nh"]-=1
                    self.fill_nh(valence-1,None)
                else:
                    for i in range(ndb):
                        for n in neigh:                        
                            if self.Gatom.nodes[n]["nh"]>=1:
                                self.Gatom.nodes[n]["atom"]="="+self.Gatom.nodes[n]["atom"]
                                self.Gatom.nodes[n]["nh"]-=1
                                self.Gatom.nodes[gros]["nh"]-=1                        
                                break
   
       


class UncyclicSmiles:
    def __init__(self,smi):
        self.smi=smi
        self.G=smi2graph(smi)
        
    
    def jejardine(G,inode,lsnode):    
        """
        creates smiles branches (i.e. elements inside "(" and ")")
        """
        
        smi=G.nodes[inode]["atom"]
        lsnode.append(inode)
        neigh=list(set(G[inode])-set(lsnode))   
        
        
        if len(neigh)>0:        
            for n in neigh[:-1]:
                smi+="("+jejardine(G,n,lsnode)+")"
            smi+=jejardine(G,neigh[-1:][0],lsnode)
            
        return smi           
