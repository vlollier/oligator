#-*-coding:Utf-8-*-
__author__ ="Virginie Lollier"
__version__ = "1.0.1"
__license__ = "BSD"

import json
import re,os,string
from utils import Logger,Atom
from enum import Enum



class SingletonTopo:
   """
   """
   
   OSE_LIST = {}     #dictionnaire listant les oses et les substitutions avec comme cl� leur num�ro et l'objet OseModel comme valeur ou une chaine de caract�re pour les substitutions
   """
   "static" dictionary with ose model identifier as key and ose model object as value
   """
   OBOND_LIST = []       #dictionnaire listant les liaison avec comme cl� un tuple : les deux oses li�s et comme valeur un tuple avec les num�ro des deux carbones de la liaison
   """
   "static" list of OsidicBond objects
   """   
   C1FREE=[]
   """
   ose identifiers without osidic binding on anomeric carbon
   """   
   DEFAULT_LINK=[1,4]   
   DEFAULT_OSETYPE={"nct":6,"ncc":5,"startc":1,"name":""}
   
   
   
   def addOse(ose):
      """
      """
      # if ose comes from if design
      if ose.oseid==None:
         OseModel.NUM_OSE+=1
         ose.oseid=OseModel.NUM_OSE
      # if ose comes from graph (loaded format)
      else:
         if ose.oseid>OseModel.NUM_OSE:
            OseModel.NUM_OSE=ose.oseid
         
      SingletonTopo.OSE_LIST[OseModel.NUM_OSE]=ose      
      SingletonTopo.C1FREE.append(ose.oseid)
      
   
   def deledge(ob):      
      SingletonTopo.OBOND_LIST.remove(ob)                    
      if ob.child_carbon==SingletonTopo.OSE_LIST[ob.child_ose].startcycle:
         SingletonTopo.C1FREE.append(ob.child_ose)
      
   def addEdge(osefrom,carbfrom,oseto,carbto):
      """
      """      
      ob=OsidicBond(osefrom,oseto,carbfrom,carbto)
      SingletonTopo.OBOND_LIST.append(ob)
      SingletonTopo.OSE_LIST[osefrom].bind_carb(carbfrom)
      SingletonTopo.OSE_LIST[oseto].bind_carb(carbto)
      
      if osefrom in SingletonTopo.C1FREE:
         if carbfrom==SingletonTopo.OSE_LIST[osefrom].startcycle:
            SingletonTopo.C1FREE.remove(osefrom)      
         
      if oseto in SingletonTopo.C1FREE:
         if carbto==SingletonTopo.OSE_LIST[oseto].startcycle:
            SingletonTopo.C1FREE.remove(oseto) 
      
      return ob
   
   def get_parent(oseid):
      for ob in SingletonTopo.OBOND_LIST:
         if ob.child_ose==oseid:
            return ob.parent_ose
      return None
   
   def get_osidicbonds(oseid1,oseid2=None):
      """
      """      
      result=[]
      for ob in SingletonTopo.OBOND_LIST:
         if ob.contains(oseid1):
            if oseid2:
               if ob.contains(oseid2):
                  result.append(ob)
            else:
               result.append(ob)
      return result
      
      
   def remove_oid(oid):
      """
      """      
      bonds=SingletonTopo.get_osidicbonds(oid)
      om=SingletonTopo.OSE_LIST[oid]
      for bond in bonds:
         if oid==bond.parent_ose:
            SingletonTopo.OSE_LIST[bond.child_ose].unbind_carb(bond.child_carbon)
         else:
            SingletonTopo.OSE_LIST[bond.parent_ose].unbind_carb(bond.parent_carbon)
         SingletonTopo.OBOND_LIST.remove(bond)
      del SingletonTopo.OSE_LIST[oid]
      if oid in SingletonTopo.C1FREE:
         SingletonTopo.C1FREE.remove(oid)
      om=None
      
   def clear():
      """
      """      
      SingletonTopo.OSE_LIST={}
      SingletonTopo.OBOND_LIST.clear()
      SingletonTopo.C1FREE.clear()
      OseModel.NUM_OSE=0 
      OsidicBond.NUM=0
      
   def strbondlist():
      """
      """      
      lsb=[]
      for b in SingletonTopo.OBOND_LIST:
         lsb.append(b.getAttributString())
      return lsb
   
   def get_directedbond(parent):
      """
      """      
      bonds=[]
      for bond in SingletonTopo.OBOND_LIST:
         if bond.parent_ose==parent:
            bonds.append(bond)
      return bonds
   
   def topogrid():
      """
      assigns row and column numbers of oses into a grid according to carbon bindings
      """      
      grid=SingletonTopo.__basegrid__()
      conflict=[]      
      oids=list(SingletonTopo.OSE_LIST.keys())
      
      for i in range(1,OseModel.NUM_OSE+1):
         for j in range(i+1,OseModel.NUM_OSE+1):
            if i in oids and j in oids:
               if i in grid and j in grid:
                  coordi=grid[i]
                  coordj=grid[j]
                  if coordi==coordj:
                     conflict.append([i,j])
      if len(conflict)>0:         
         for collision in conflict:
            print(("collision imgs:",collision))
            group1=[]
            group2=[]
      
            oid1=collision[0]
            oid2=collision[1]
            
            ancetres1=SingletonTopo.__parseTopo__([oid1],oid1,"end")
            ancetres2=SingletonTopo.__parseTopo__([oid2],oid2,"end")
      
            shared=set(ancetres1).intersection(set(ancetres2))
      
            if len(shared)>0:
               for oid in shared:
                  ancetres1.remove(oid)
                  ancetres2.remove(oid)
               group1=SingletonTopo.__parseTopo__([],ancetres1[-1:][0],"start")
               group1.append(ancetres1[-1:][0])
               group2=SingletonTopo.__parseTopo__([],ancetres2[-1:][0],"start")      
               group2.append(ancetres2[-1:][0])
               #print(ancetres1[-1:][0])
               #print(ancetres2[-1:][0])
               coords_ancetre1=grid[ancetres1[-1:][0]]
               coords_ancetre2=grid[ancetres2[-1:][0]]
               if coords_ancetre1[0]>coords_ancetre2[0]:
                  for elt in group1:
                     grid[elt][0]+=1
                  for elt in group2:
                     grid[elt][0]-=1
                     
               else:
                  for elt in group1:
                     grid[elt][0]-=1
                  for elt in group2:
                     grid[elt][0]+=1
         
      #print(grid)
      return grid
         
   
   def __parseTopo__(branch,om,direction,limit=None):
      """
      """      
      if limit==None or len(branch)<limit:
         
         edge_next=[]
         
         bonds=SingletonTopo.get_osidicbonds(om)
         for bond in bonds:
            if direction=="start" and bond.parent_ose==om:
               edge_next.append(bond.child_ose)
            if direction=="end" and bond.child_ose==om:
               edge_next.append(bond.parent_ose)

         if len(edge_next)>0:              
            bb=[]                
            for e in edge_next:
               bb+=SingletonTopo.__parseTopo__(branch+[e],e,direction,limit)
                             
            return bb     
         else:            
            return branch
      else:
         return branch   
   
         
   def __basegrid__(bond=None,grid=None):
      """
      """      
      bonds=[]
      if bond:
         print(bond.getAttributString())
      if grid==None:
         grid={}         
         start=SingletonTopo.C1FREE[0]         
         grid[start]=[0,0]
         bonds=SingletonTopo.get_directedbond(start)
      else:
         row=grid[bond.parent_ose][0]
         col=grid[bond.parent_ose][1]-1
         if bond.parent_carbon>=5:
            row-=1
         elif bond.parent_carbon in [2,3]:
            row+=1
            if bond.parent_carbon == 2:
               col+=1
         grid[bond.child_ose]=[row,col]
         
         bonds=SingletonTopo.get_directedbond(bond.child_ose)
      
      if len(bonds)>0:         
         for b in bonds:
            SingletonTopo.__basegrid__(b,grid)   
         
      return grid
   
   
      
class OseModel:
   """
   """
   NUM_OSE = 0  
   """
   numbering of ose instances (unique identifier)
   """
   
   def __init__(self,ncc=None,nct=None,startc=None):   
      self.oseid=None
      self.name=""
      if startc:
         self.startcycle=startc
      else:
         self.startcycle=SingletonTopo.DEFAULT_OSETYPE["startc"]
      if ncc:
         self.ncc=ncc
      else:
         self.ncc=SingletonTopo.DEFAULT_OSETYPE["ncc"]
      if nct:         
         self.nct=nct
      else:
         self.nct=SingletonTopo.DEFAULT_OSETYPE["nct"]
      self.modifs=[]
      iso=[]
      mods=[]      
      #init carbon description table
      for icarb in range(self.nct):
         iso.append("")
         if icarb==self.ncc-1:
            mods.append(-1)
         else:   
            mods.append(SubstitutionLibrary.NOSUBID)
      self.modifs.append(iso)
      self.modifs.append(mods)
      self.anhydro=False
   
   def setname(self,name):
      self.name=name
      
   def in_cycle(self,cnum):
      """
      """      
      return cnum in range(self.startcycle,self.ncc+self.startcycle)
   
   def bind_carb(self,cnum):
      """
      """      
      self.modifs[1][cnum-1]=-1
   
   def unbind_carb(self,cnum):
      """
      """      
      self.modifs[1][cnum-1]=SubstitutionLibrary.NOSUBID
   
   def set_modcarb(self,cnum,subid):
      """
      """      
      self.modifs[1][cnum-1]=subid
      
   def get_modcarb(self,cnum):
      """
      """      
      return self.modifs[1][cnum-1]
      
   def get_carbmod(self,idsub):
      carbs=[]
      for icarb in range(len(self.modifs)):
         if self.modifs[1][icarb]==idsub:
            carbs.append(icarb+1)
      return carbs
   
   def set_isocarb(self,cnum,iso):
      """
      """      
      self.modifs[0][cnum-1]=iso
   
   def get_isocarb(self,cnum):
      return self.modifs[0][cnum-1]
      
   def rm_carbs(self,carbs):
      """
      """      
      for cnum in carbs:
         self.modifs[1].pop(cnum)
   
   def add_carbs(self,nbcarbs,mods=None):
      """
      """      
      for icarb in range(nbcarbs):
         if mods:
            self.modifs[0].append(mods[icarb][0])
            self.modifs[1].append(mods[icarb][1])
         else:
         
            self.modifs[0].append("")
            self.modifs[1].append(SubstitutionLibrary.NOSUBID)
      
   def set_anhydrobond(self,bind):
      """
      bind: boolean for now (text like 3,6 can be later)
      """
      self.anhydro=bind
   
   def get_endc(self):
      return self.ncc+self.startcycle-1
   
   def striso(self):
      iso=""
      for isocarb in self.modifs[0]:
         if isocarb in ["D","L"]:
            iso+=isocarb
         else:
            iso+="-"
      return iso
   
   def get_epimer(self):
      isotxt=""
      endc=self.startcycle+self.ncc-1
      
      # not chiral
      if self.nct==self.ncc:
         endc-=1
         
      # neuraminic acid as abnormal example ....
      if self.nct>6:
         endc=self.nct-1
         
      for cnum in range(self.startcycle+1,endc):
         iso=self.modifs[0][cnum-1]
         if iso in ["D","L"] :
            isotxt+=iso
         else:
            isotxt+="-"
      return isotxt

         

class Substitution:         
   """
   """
   UKN=0
   """
   naming of unreferenced substitution
   """
   
   PATFORM="[A-Z][a-z]?[0-9]{1,2}"
   CID=0
   """
   numbering of substitution instances (as unique identifier)
   """   
   
   def __init__(self,formula,link=None,name=None,smiles=None):    
      """
      """
      if name==None:
         self.name="ukn"+str(Substitution.UKN)
         Substitution.UKN+=1           
      else:
         self.name=name
      
      self.identifier=None      
      self.formula=formula
      self.smiles=smiles     
      self.link=""
      
      if link==None:
         if smiles!=None:
            self.link=re.match("^=?[A-Z][a-z]?",smiles).group()
         elif re.search("O",formula):
            self.link="O"
         else:
            self.link=re.match("^=?[A-Z][a-z]?",re.sub("H[0-9]+","",formula)).group()
        
      else:
         self.link=link
            
   # ose carbon + delta   
   def get_delta(self):
      """
      Mass of the formula minus the mass of the replaced OH (and 1H if linkage is on ose C)
      :return: the round difference of Dalton mass between OH and elements in the substitution
      :rtype: float
      """
      delta=self.massSubstitution()
          
      o=Atom.mass("O")
      h=Atom.mass("H")
         
      delta=delta-o-h
      if re.match("^=",self.smiles):
         delta-=h
      linkage=self.smiles
      while re.match("^\(",linkage):
         delta-=h
         if re.match("^\(=",self.smiles):
            delta-=h           
         linkage=linkage[linkage.find(")")+1:]
          
            
      #if re.match("^=",self.link):      
         #delta-=h
      
      ### linkage on carbon not at Oxygen place
      ## !!!! ne convient pas pour substit=desoxy, cas particulier ou` link=""
      ## ok si 2 atomes liés au C de l'ose      
      #if self.link=="" and self.formula!="H1":         
         ## remove an H to the carbon , add the mass of substituent formula minus OH
         #delta-=h
         
      return round(delta,3)
   
   def equals(self,compar):
      """
      """
      formula_ref=sorted(re.findall("[A-Z][a-z]?[0-9]+",self.formula))
      txtref=""
      for a in formula_ref:
         txtref+=a   
      return compar==txtref
   
 
   
      
   def massSubstitution(self): 
      """
      """
      return Substitution.__compute_mass__(self.formula)
   
   
   def __compute_mass__(formula):
      """
      """
      m=0
     
      f=re.findall(Substitution.PATFORM,formula)
      for ab in f:
         #m+=SubstitutionLibrary.massAtom(re.sub("[0-9]","",ab))*int(re.sub("[a-zA-Z]","",ab))
         m+=Atom.mass(re.sub("[0-9]","",ab))*int(re.sub("[a-zA-Z]","",ab))
      return m       
   

         
    
class SubstitutionLibrary:
   """
   """
   
       
   SUBSTITUTIONS=[]  
   """
   store the list of substitutions
   """
   NOSUBID=0
   """
   store the identifier of no substitution (OH on carbon)
   """
   
   def __init__(self):
      """
      """
      
         #name,formula,linkage,smiles=None
      SubstitutionLibrary.create_substit("oxydation","O2H1","","(=O)O")
      #SubstitutionLibrary.create_substit("keto","C1O2H4","","(CO)O")
      SubstitutionLibrary.create_substit("desoxy","H1","","")
      
      SubstitutionLibrary.create_substit("double_bond","O1",None,"=O")         
      SubstitutionLibrary.create_substit("ferulic_acid","C10H9O4","O","OC(=O)C=CC1=CC=C(O)C(OC)=C1") 
      
   def sort_by_name():     
      """
      Sort the list of substitutions according to their name
      """
      
      SubstitutionLibrary.SUBSTITUTIONS=sorted(SubstitutionLibrary.SUBSTITUTIONS,key=lambda s:s.name)  
      
                      
   @staticmethod
   def get_data():
      """
      """
      data=[]      
      for s in SubstitutionLibrary.SUBSTITUTIONS:      
            data.append({"identifier":s.identifier,"name":s.name,"formula":s.formula,"link":s.link,"smiles":s.smiles,"mass":s.massSubstitution()})
      
      return data
  
   def to_json():
      non=["oxydation","keto","desoxy","double_bond","ferulic_acid"]
      data={}  
      for s in SubstitutionLibrary.SUBSTITUTIONS:      
         if s.name not in non:
            data[s.name]={"formula":s.formula,"smiles": s.smiles}
            
      return json.dumps(data)
      
      
   def create_substit(name,formula,link,smiles=None):
      """
      :return: a new substitution if formula not found in internal ressource
      """
      
      query= SubstitutionLibrary.get_subformul(formula)
      
      if not query:         
         substit=Substitution(formula,link,name,smiles)         
         SubstitutionLibrary.add_substit(substit)
         return substit
      elif link!=query.link:
         substit=Substitution(formula,link,name,smiles)         
         SubstitutionLibrary.add_substit(substit)      
         return substit
      else:
         return query
      
   def add_substit(substit):
      """
      add the substitution object to the internal catalog
      """
      SubstitutionLibrary.SUBSTITUTIONS.append(substit)      
      Substitution.CID+=1
      substit.identifier=Substitution.CID 
      
     
   @staticmethod
   def get_subname(identifier):
      """
      :param identifier: indification number
      :type: int
      
      :return: the substituent name according to the identifier
      :rtype: string
      """
      s=SubstitutionLibrary.getSub(identifier)
      if s:
         return s.name
         
      Logger.debug("identifier not found: %i"%identifier,1)
      return ""
      
  
   
   def get_subfromname(name):      
      """
      :rtype: Substitution
      """
      for s in SubstitutionLibrary.SUBSTITUTIONS:
         if s.name.lower()==name.lower():
            return s   
      return None      
   
   @staticmethod
   def get_subformul(formula):
      """
      """
      for s in SubstitutionLibrary.SUBSTITUTIONS:
         if sorted(s.formula)==sorted(formula):
            return s
      return None

  
   

   def mod_id(**kwargs):
      """
      """
      for modsub in SubstitutionLibrary.SUBSTITUTIONS:
         if modsub.identifier==int(kwargs["identifier"]):
            modsub.name=kwargs["name"]
            modsub.formula=kwargs["formula"]
            modsub.smiles=kwargs["smiles"]
            modsub.link=kwargs["linkage"]
      
   
  
   
   def rm_id(identifier):   
      """
      """
      matchsub=None
      for substit in SubstitutionLibrary.SUBSTITUTIONS:
         if substit.identifier==identifier:
            matchsub=substit
      if matchsub:
         SubstitutionLibrary.SUBSTITUTIONS.remove(matchsub)
         
   def check_smiles(identifier):
      """
      """
      substit=SubstitutionLibrary.getSub(identifier)
      if substit!=None:
         smiles=substit.smiles
         compo=substit.formula
         if smiles=="" and compo!="H1":
            return False
         if smiles=="O" and compo not in ["O1H1","H1O1"]:
            return False
         if smiles=="=":
            return False
         return True
         
      else:
         return False
     
   
   def get_non_smiles():
      """
      """
      non_smi={}
      for s in SubstitutionLibrary.SUBSTITUTIONS:
         if s.smiles==None and s.formula!="H1":
            non_smi[s.name]=s.massSubstitution()
      if len(non_smi)>0:
         return non_smi
      else:
         return None
   
   
   def getSub(identifier):
      """
      """
      for s in SubstitutionLibrary.SUBSTITUTIONS:
         if s.identifier==identifier:
            return s
      return None
   
   # formula is dictionary atom:count
   @staticmethod
   def get_substitution(formula):
      """
      :param formula: dictionary atom:count
      """
      txtformula=""
      for atom,nb in sorted(formula.items()):
         if nb>0:
            txtformula+=atom+str(nb)
      
     
      for substit in SubstitutionLibrary.SUBSTITUTIONS:
         if substit.equals(txtformula):
            return substit
      return None
    
   
         

# small class to store binding direction from reducing end
class OsidicBond:
   NUM=0
   def __init__(self,parent_ose,child_ose,parent_carbon,child_carbon):   
      """
      directed binding where parent is on the side of the reducing end
      """
      self.parent_ose=parent_ose
      self.child_ose=child_ose
      self.parent_carbon=parent_carbon
      self.child_carbon=child_carbon
      self.onode=None # node number of O in gatom
      OsidicBond.NUM+=1
      self.identifier=OsidicBond.NUM

   def __inverse__(self):
      """      
      """
      o1=self.parent_ose
      c1=self.parent_carbon

      o2=self.child_ose
      c2=self.child_carbon

      self.parent_ose=o2
      self.parent_carbon=c2

      self.child_ose=o1
      self.child_carbon=c1
   
   def contains(self,oid,cnum=None):
      """
      """
      ose=False
      if oid==self.parent_ose or oid==self.child_ose:
         ose=True
         
      if ose:
         if cnum:
            if oid==self.parent_ose and cnum==self.parent_carbon:
               return True
            elif oid==self.child_ose and cnum==self.child_carbon:
               return True
            else:
               return False
         else:
            return True
      return False
               
            
   def setParent(self,osenum):
      """
      """
      if osenum==self.parent_ose:
         return self.child_ose
      elif osenum==self.child_ose:
         self.__inverse__()            
         return self.child_ose
      else:
         return None

   def set_opos(self,node_number):
      """
      """
      self.onode=node_number
      
   def getAttributString(self):
      """
      """
      bound="%i (%i+%i) %i "%(self.parent_ose,self.parent_carbon,self.child_carbon,self.child_ose)
      return bound
   
   # ose1,carb1,carb2,ose2
   @staticmethod   
   def getBondFromAttribute(attr):      
      """
      :param attr: ose(carb+carb)ose
      :type: string
      
      :return: list of int
      """
      result=[]
      b=attr.split('+')
      if len(b)==2:
         parent=b[0].split('(')
         child=b[1].split(')')
         b=parent+child
         if len(b)==4:
            for i in b:
               result.append(i.strip())                   
      return result
   
   
      
# at least one instance from default      
SubLib=SubstitutionLibrary()

with open ("Ressources"+os.path.sep+"substitution_list.json", 'r') as f :
   JSON = json.load(f)
   
   for k,v in JSON.items():
      if "smiles" in v:
         if v["smiles"][0]=="(":
            link=""
         else:
            link=None
         
        
         #name,formula,linkage,smiles=None        
         s=SubstitutionLibrary.create_substit(k,v["formula"],link,v["smiles"])
         if s.name=="hydroxy":
            SubstitutionLibrary.NOSUBID=s.identifier
            
"""
Create a pre-formed ose model from a short list of most frequent oses 
"""
class OseFactory:
   
   
   _code_={}
   _code_["D-Glcp"]=("-DLDD-",5,1)
   _code_["D-Manp"]=("-LLDD-",5,1)
   _code_["L-Galp"]=("-DLLL-" ,5,1)      
   _code_["D-Galp"]=("-DLLD-" ,5,1)      
   _code_["D-Xylp"]=("-DLD-",5,1)
   _code_["L-Araf"]=("-LDL-",4,1)
   
   _code_["D-Neu5Ac"]=("---LDLDD-",5,2,["1ox","3desox","5nac"])
   _code_["L-Fucp"]=("-DLLL-",5,1,["6desox"])
   _code_["L-Rhap"]=("-LLDL-",5,1,["6desox"])
   _code_["L-IdoAp"]=("-LDLL-",5,1,["6ox"])
   _code_["L-GulAp"]=("-DDLL-",5,1,["6ox"])
   _code_["D-GalAp"]=("-DLLD-" ,5,1,["6ox"])   
   _code_["D-GalNAcp"]=("-DLLD-" ,5,1,["2nac"])   
   _code_["D-GlcNAcp"]=("-DLDD-",5,1,["2nac"])

   _dicsubcode_={"nac":SubstitutionLibrary.get_subfromname("n-acetyl").identifier,"ox":SubstitutionLibrary.get_subfromname("oxydation").identifier,"desox":SubstitutionLibrary.get_subfromname("desoxy").identifier}
   #_dicsubcode_={"nac":SubstitutionLibrary.get_subfromname("n-acetyl").identifier,"ac":SubstitutionLibrary.get_subfromname("acetyl").identifier,"ox":SubstitutionLibrary.get_subfromname("oxydation").identifier,"desox":SubstitutionLibrary.get_subfromname("desoxy").identifier}
     
   def get_ordered_templates(order=None):
      if order=="alpha":
         return sorted(OseFactory._code_.keys())
      elif order=="iso":
         return sorted(OseFactory._code_.keys(),key=lambda k:OseFactory._code_[k][0])
      elif order=="stem":
         return sorted(OseFactory._code_.keys(),key=lambda k:re.sub("[DL]-","",k))
      else:         
         return OseFactory._code_.keys()
      
   def getcolor(name):      
         
      key=re.sub("^[DL]-","",name)
      key=key[0:3]
      
      #stereochemistry color code (not compliant with snfg system for fuc and neu5ac)
      _color = {}   
      _color["Glc"] = "blue"
      _color["Man"] = "green"   
      _color["Ido"] = "saddlebrown" # brown
      _color["Gal"] = "gold"   
      _color["Fuc"] = "gold"
      _color["Xyl"] = "darkorange"
      #_color["Fuc"] = "red"
      #_color["Neu"] = "purple"    
      _color["Neu"] = None          
      _color["Rha"] = "green"
      _color["Gul"] = "darkorange"
      _color["Ara"] = "green"
      
      return _color.get(key)
   
   def substituted(name):      
      return len(OseFactory._code_.get(name))==4
     
      
   
   
   def getose(name):
      om=None      
      
      code=OseFactory._code_[name]
      om=OseModel(code[1],len(code[0]),code[2])
      subs=[]
      for icarb in range(0,len(code[0])):
         if code[0][icarb] in ["D","L"]:
            om.set_isocarb(icarb+1,code[0][icarb])         
      if len(code)==4:
         subs=code[3]
      for subcode in subs:
         cnum=int(re.sub("[^0-9]","",subcode))
         sname=re.sub("[0-9]","",subcode).lower()
         om.set_modcarb(cnum,OseFactory._dicsubcode_.get(sname))
            
      
      om.setname(name)
         
      return om      
   
   
   def checkose(om):
      
      name=""
      #iso="-"+om.get_epimer()+"-"
      
      iso=om.striso()
      for ref in OseFactory._code_.keys():
         code=OseFactory._code_.get(ref)
         pattiso="^"+re.sub("-",".",code[0])+"$"
         if re.match(pattiso,iso):            
            if len(code)==4:
               alt=ref
               for subcode in code[3]:
                  cnum=int(re.sub("[^0-9]","",subcode))
                  sname=re.sub("[0-9]","",subcode).lower()                  
                  if om.get_modcarb(cnum)!=OseFactory._dicsubcode_[sname]:
                     alt=None                     
                     break
               if alt!=None:
                  name=alt
                  return name
            else:
               name=ref
     
      
      return name
   
               