#-*-coding:Utf-8-*-
__author__ ="Virginie Lollier"
__version__ = "1.0.1"
__license__ = "BSD"

from PyQt5.QtWidgets import *
from PyQt5.QtGui import QPen, QColor,QDoubleValidator, QIntValidator,QFont,QPixmap,QCursor,QPen
from PyQt5.QtCore import Qt, pyqtSignal,QAbstractTableModel, QVariant, pyqtSlot,  QPersistentModelIndex, QModelIndex,QRectF,QPointF,QSizeF,QObject
from PyQt5.uic import *



from math import sqrt
from GlycanTopology import SubstitutionLibrary,Substitution
from Conversion import IupacName,UncyclicSmiles
from operator import itemgetter
import re,os,numpy as np,time,random,string
from utils import Logger,Atom

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt


import networkx as nx

# number of carbons according to monosacch size name
NCT={"hex":6,"pen":5,"tetra":4,"tri":3,"dec":10,"non":9,"oct":8}
# number of carbons according to monosacch cycle form     
NCF={"P":5,"F":4}     



class Graphics( QGraphicsView ):
    
    def __init__(self,scene):
        QGraphicsView.__init__(self,scene)        

        self.setDragMode(QGraphicsView.RubberBandDrag) 

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)   
        #self.setAutoFillBackground(False)
        self.setStyleSheet("background:white")
        self.Trash=[]
        self.grid=False
        self.gridItem=None
        #self.identifier=0
        #self.signal=SignalObject()
        
            
   
    def applyTopoGrid(self,grid):
        """
        """
        imglist=[]
        imgfree=[]
        padding=10
        
        for elt in self.scene().items():
            if isinstance(elt,OseImage):                
                if len(elt.lines())>0:
                    imglist.append(elt)        
                else:
                    imgfree.append(elt)
        
        
        
        x=0
        y=0
    
        for img in imglist:   
            if img.identifier!=None: 
                coords=grid[img.identifier]                
                img.setPos(x+coords[1]*(OseImage.WIDTH*1.5),y+coords[0]*(OseImage.HEIGHT*1.5))
    
        
        if len(imgfree)>0:
            for foi in imgfree:
                foi.setPos(*self.position_oseimage())
        
        
    def position_oseimage(self):
        """
        """
        x=0
        y=0
        w=OseImage.WIDTH*1.5
        h=OseImage.HEIGHT
        
        
        for elt in self.scene().items():
            if isinstance(elt,OseImage):   
                if elt.pos().x()<x:
                    x=elt.pos().x()
                if elt.pos().y()>y:
                    y=elt.pos().y()
        
        x-=w
        
        return x,y
           
    def drawGrid(self,nrow,ncol):
        pen = QPen(Qt.green)
        coords=np.array(coords)
        side = OseView.WIDTH*20
        
        for i in range(nrow):
            for j in range(ncol):
                self.gridItem = QRectF(QPointF(i*side, j*side), QSizeF(side, side))
                self.scene().addRect(self.gridItem, pen)
        
        return gridItem
        
    #def drawGrid(self,coords,x,y):    
        #pen = QPen(Qt.green)
        #coords=np.array(coords)
        #side = OseView.WIDTH*20
        
        #for i in range(coords.shape[0]):
            #for j in range(coords.shape[1]):
                #self.gridItem = QRectF(QPointF(i*side, j*side), QSizeF(side, side))
                #self.scene().addRect(self.gridItem, pen)
        
        #return gridItem

    def hideGrid(self):
        self.scene().removeItem(self.gridItem)
        self.gridItem=None
        
    #def keyPressEvent(self, event):
        #"""
        #"""
        #if event.key() == Qt.Key_Delete :
            #self.signal.suppr_press()
            
    
    
    def emptyTrash(self):  
        """
        """
        #print(self.Trash)
        if len(self.Trash)>0:
            for elt in self.Trash:                
                self.scene().removeItem(elt)
        self.Trash.clear()        
        self.viewport().update()
        
    

class edgeType(QDialog):
    
    
    
    def __init__(self,fromto, parent=None):
        
        QDialog.__init__(self, parent)
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        gridlayout = QGridLayout()
        
        label = QLabel("Default binding between : ")
        gridlayout.addWidget(label, 0, 0, 1, 3)        
        
        
        self.carb1,self.carb2=fromto
        self.c1ComboBox = QComboBox()        
        self.c1ComboBox.addItem(str(self.carb1),self.carb1)
        self.c1ComboBox.addItem("2")
        self.c1ComboBox.setCurrentText(str(self.carb1)) 
        gridlayout.addWidget(self.c1ComboBox, 1, 0, 1, 1)        
        
        label2 = QLabel("->")
        gridlayout.addWidget(label2, 1, 1, 1, 1)
        
        self.num2 = QLineEdit(str(self.carb2))          
        gridlayout.addWidget(self.num2, 1, 2, 1, 1)
        
        buttonBox = QDialogButtonBox(self)
        buttonBox.setLayoutDirection(Qt.LeftToRight)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

        gridlayout.addWidget(buttonBox, 2, 1, 1, 2)
        self.setLayout(gridlayout)
    
        self.setWindowTitle("Default bond")
        self.setWindowModality(Qt.ApplicationModal)
        
    def accept(self):
        """
        """
        if str(self.c1ComboBox.currentText()).isdigit() and self.num2.text().isdigit():
            if  int(self.num2.text())>6 or int(self.num2.text())<1:        
                QMessageBox.warning(self, "Warning", "Please enter a correct carbon number")
            else :                
                self.carb1 = self.c1ComboBox.currentData()  
                if re.match("^[0-9]+$",self.num2.text()):
                    self.carb2 = int(self.num2.text())
                    QDialog.accept(self)
        else :
            QMessageBox.warning(self, "Warning", "Please enter a correct carbon number")
            
        
class oseType(QDialog):    
    
    ncc=5
    nct=6
    name=""
    
    def __init__(self,current, templates, parent=None):
        
        QDialog.__init__(self, parent)
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        gridlayout = QGridLayout()
        # --------
        groupType=QGroupBox("Ose size")
        osebox = QVBoxLayout()
        gridlayout.addWidget(groupType,1,0,1,2)
                
        
        self.hexoseButton = QRadioButton("Hexose")          
        osebox.addWidget(self.hexoseButton)
        
        self.pentoseButton = QRadioButton("Pentose")          
        osebox.addWidget(self.pentoseButton)
        groupType.setLayout(osebox)
        
        # --------
        groupForm=QGroupBox("Cycle form")
        formbox=QVBoxLayout()
        
        self.pyraneButton = QRadioButton("Pyrane")          
        formbox.addWidget(self.pyraneButton)
        self.furaneButton = QRadioButton("Furane")          
        formbox.addWidget(self.furaneButton)
        gridlayout.addWidget(groupForm,2,0,1,2)
        groupForm.setLayout(formbox)
        
        # --------         
        groupStd=QGroupBox("OR select a common ose")
        stdbox = QVBoxLayout()
        gridlayout.addWidget(groupStd,3,0,1,2)
        groupStd.setLayout(stdbox) 
        self.cmbStd=QComboBox()
        self.cmbStd.addItems([""]+templates)       
        self.cmbStd.currentTextChanged.connect(self.selectname)
        
        
        stdbox.addWidget(self.cmbStd)
        
        # --------
        buttonBox = QDialogButtonBox(self)
        buttonBox.setLayoutDirection(Qt.LeftToRight)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)


        gridlayout.addWidget(buttonBox, 5, 0, 1, 2)
        self.setLayout(gridlayout)
    
        self.setWindowTitle("Default ose")
        self.setWindowModality(Qt.ApplicationModal)
        
        self.ncc=current["ncc"]
        self.nct=current["nct"]
        self.cmbStd.setCurrentText(current["name"])
        self.ose()
        
    
    def ose(self):
        """
        """
        
        if self.nct==6:
            self.hexoseButton.setChecked(True)
        elif self.nct==5:
            self.pentoseButton.setChecked(True)
        
        if self.ncc==5:
            self.pyraneButton.setChecked(True)
        elif self.ncc==4:
            self.furaneButton.setChecked(True)
        
    
            
    def selectname(self):
        self.name=self.cmbStd.currentText()
    
    
    
    def accept(self):
        """
        """
        
        if self.furaneButton.isChecked():
            self.ncc=4
        else:
            self.ncc=5
        
        if self.pentoseButton.isChecked():
            self.nct=5
        else:
            self.nct=6
        QDialog.accept(self)


class SubGlycanDialog(QDialog):
    def __init__(self,add=True,parent=None):
        QDialog.__init__(self, parent)
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        self.resize(300, 150)
        gridlayout = QGridLayout(self)
        self.setLayout(gridlayout)
        if add:
            self.setWindowTitle("Replace all free glycosidic OH")
        else:
            self.setWindowTitle("Remove a substitution type on glycan")
        self.setWindowModality(Qt.ApplicationModal)     
        self.subid=None
        self.comboSub=QComboBox()
        self.comboSub.addItem("",None)
        for s in SubstitutionLibrary.SUBSTITUTIONS:
            if s.identifier!=SubstitutionLibrary.NOSUBID:
                self.comboSub.addItem(s.name,s.identifier)  
            #if s.name!="keto":
                #self.modifCombobox.addItem(s.name,s.identifier)  
        
        #self.comboSub.setCurrentText(SubstitutionLibrary.get_subname(SubstitutionLibrary.NOSUBID))         
        self.comboSub.currentTextChanged.connect(self.selectSub)
        gridlayout.addWidget(QLabel("Select a substituent from list:"),0,0,1,1)
        gridlayout.addWidget(self.comboSub,2,0,1,1)
        
        buttonBox = QDialogButtonBox()
        buttonBox.setLayoutDirection(Qt.LeftToRight)
        buttonBox.setOrientation(Qt.Horizontal)
    
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)        
        gridlayout.addWidget(buttonBox, 3, 0, 1, 1)        
                           
    def selectSub(self):
        self.subid=self.comboSub.currentData()
        
    def cancel(self):
        self.subid=None
        
    def accept(self):
        """
        """
        self.selectSub()
        QDialog.accept(self)        

class EdgeDialog(QDialog):
    def __init__(self,c1,c2,free,parent=None):
        self.carbfrom=c1
        self.carbto=c2
        QDialog.__init__(self, parent)
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        self.resize(300, 150)
        gridlayout = QGridLayout(self)
        label = QLabel("carbon selection:")
        gridlayout.addWidget(label, 0, 0, 1, 1)    
        
        self.ToComboBox = QComboBox()        
        self.ToComboBox.addItem("1 (aldose)",1)
        self.ToComboBox.addItem("2 (ketose)",2)
        #self.ToComboBox.setCurrentText(str(c2)) 
        self.ToComboBox.setCurrentIndex(self.ToComboBox.findData(c2))
        gridlayout.addWidget(self.ToComboBox, 0, 1, 1, 1)
        
        label2 = QLabel(" -> ")
        gridlayout.addWidget(label2, 0, 2, 1, 1)
        self.FromComboBox = QComboBox() 
        for cnum in free:
            self.FromComboBox.addItem(str(cnum),cnum)
        gridlayout.addWidget(self.FromComboBox, 0, 3, 1, 1)      
        self.FromComboBox.setCurrentText(str(c1)) 
        
        buttonBox = QDialogButtonBox()
        buttonBox.setLayoutDirection(Qt.LeftToRight)
        buttonBox.setOrientation(Qt.Horizontal)

        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)        
        gridlayout.addWidget(buttonBox, 2, 2, 1, 2)
        
        self.setLayout(gridlayout)
        self.setWindowTitle("Binding properties")
        self.setWindowModality(Qt.ApplicationModal)        
        
    def accept(self):
        """
        """
        self.carbFrom=self.FromComboBox.currentData()
        self.carbTo=self.ToComboBox.currentData()
        QDialog.accept(self)
    
   

# formulaire de creation ou de modification d'une substitution chimique du OH
class SubstituantForm(QDialog):
    """
    formulaire de creation ou de modification d'une substitution chimique du OH
    """
    
    VALIDSMI=["C","O","N","S","P","F","Na","Cl","Se","Br","K","I"]
    
             
    def __init__(self, window,**kwargs):
        
        QDialog.__init__(self, None)
        self.isdefault=len(kwargs)==0
        
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        self.window = window
        
        gridlayout = QGridLayout()
        
        self.txtname = QLineEdit()
        self.txtsmi = QLineEdit()
        self.btnShowG=QPushButton("Refresh graph view")        
        self.btnShowG.clicked.connect(self.plotG)
        
        self.chkDble=QCheckBox("Double bond")        
        self.chkDble.clicked.connect(self.doubleEvent)
        #self.chkDble.stateChanged.connect(self.doubleEvent)
        
        self.fig=Figure()        
        self.canvas=FigureCanvas(self.fig)
        
        self.lblCompo=QLabel()
        self.lblMass=QLabel()
        
        gridlayout.addWidget(QLabel("Name: "),0,0,1,1)
        lblCanonsmi=QLabel("Canonical\nSMILES: ")
        lblCanonsmi.setAlignment(Qt.AlignLeft)
        gridlayout.addWidget(lblCanonsmi,1,0,1,1)
        
        gridlayout.addWidget(self.txtname,0,1,1,1)      
        gridlayout.addWidget(self.txtsmi,1,1,1,2)                
        gridlayout.addWidget(self.chkDble,2,1,1,1)    
        gridlayout.addWidget(self.btnShowG,2,2,1,1)   
        gridlayout.addWidget(self.canvas,3,0,1,3)                
        gridlayout.addWidget(self.lblCompo,4,0,1,1)     
        gridlayout.addWidget(self.lblMass,4,2,1,1,alignment=Qt.AlignRight)     
        
        buttonBox = QDialogButtonBox(self)
        buttonBox.setLayoutDirection(Qt.LeftToRight)
        buttonBox.setOrientation(Qt.Horizontal)
        buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)        
        self.btnok=buttonBox.button(QDialogButtonBox.Ok)                
        gridlayout.addWidget(buttonBox,5,2,1,1) 
        
        self.setLayout(gridlayout)        
        
        self.valmin=1
        self.mass=0
        self.compo=""
        self.row = None
        self.identifier=None        
        
        if self.isdefault:            
            self.setWindowTitle("Add substituent")
            # temporary identifier
            self.identifier=-len(SubstitutionLibrary.SUBSTITUTIONS)  
            #random name picked in ascii letters
            randname=""
            for i in range(random.randrange(4,8)):
                randname+=string.ascii_letters[random.randrange(52)]
            self.txtname.setText(randname)
            
            self.smiles=""
            self.G=None
            self.txtsmi.textChanged.connect(self.smiChanged) 
            #self.connector.addItem(SubstituantForm.directlink)
        else :
            self.setWindowTitle("Modify substituent")
            self.row = kwargs["row"]
            if "identifier" in kwargs.keys():
                self.identifier = kwargs["identifier"]
            
            self.txtname.setText(kwargs["name"])
            self.smiles=kwargs["smiles"]            
            if self.smiles!="":                
                self.txtsmi.setText(self.smiles)     
                self.chkDble.setChecked(self.smiles[0]=="=" or self.smiles[0:2]=="(=") 
                if self.check_notation(self.smiles):
                    self.plotG()            
                    self.set_formula()                               
            self.txtsmi.textChanged.connect(self.smiChanged)   
                
        

        self.setWindowModality(Qt.ApplicationModal)
        self.exec_()   
          
        
    
            
    def check_notation(self,smi):   
        atoms=re.findall("[A-Z][a-z]?",smi)
        if re.match("^[A-Za-z\(\)=0-9]+$",smi) and len(atoms)>0:
            nbParO=len(re.findall(re.escape("("),smi))
            nbParF=len(re.findall(re.escape(")"),smi))
            if nbParF==nbParO:
                nbcycle_pair=len(re.findall("[0-9]",smi))%2
                if nbcycle_pair==0:                    
                    for atom in atoms:
                        if not atom in Atom.dicovalence.keys():
                            return False
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
            
    def smiChanged(self):
        self.changed=False
        txt=self.txtsmi.text()
        if self.txtsmi.hasFocus():                        
            if self.check_notation(txt):
                self.plotG()
                if  txt[0]=="=" or txt[0:2]=="(=":                    
                    if self.G.nodes()[min(self.G.nodes())]["valence"]>=self.valmin:                    
                        self.chkDble.setChecked(True)                                                          
                    else:
                        self.chkDble.setChecked(False)                
                    
                else:
                    self.smiles=txt    
                    self.chkDble.setChecked(False)                
                
                self.set_formula()                                                  
                self.changed=True                
                
            else:
                self.chkDble.setChecked(False)  
        else:
            self.changed=True           
            if txt!="":                
                self.chkDble.setChecked(txt[0]=="=" or txt[0:2]=="(=")    
                self.set_formula()
            
        
    def set_formula(self):
        dico_compte={}
        nH=0
        compo=""
        
        self.compo=""
        self.mass=0
        
        if self.smiles=="":
            self.lblCompo=""
            self.lblMass=""
        else:
            for symbol in Atom.dicovalence.keys():
                lsatom=re.findall("[A-Z][a-z]?",self.smiles)
                ns=0
                if symbol!="":
                    for atom in lsatom:
                        if atom==symbol:
                            ns+=1
                       
                    if ns>0:
                        lsatom.remove(symbol)
                    #dico_compte[symbol]=len(re.findall(symbol,self.smiles))
                    dico_compte[symbol]=ns 
                    nH+=dico_compte[symbol]*Atom.dicovalence.get(symbol)
            
            if self.G!=None :
                paths=list(nx.connected_components(self.G))
                # connexion(s) à l'ose                
                nH-=len(paths)
                if self.chkDble.isChecked():
                    nH-=1
                    
                for edge in self.G.edges():
                    nH-=self.G.edges[edge]["retireH"]*2
                
                
                
                
            dico_compte["H"]=nH
            
            for symbol in dico_compte:
                if dico_compte.get(symbol)>0:
                    self.compo+=symbol+str(dico_compte.get(symbol))
                    compo+=symbol+"<sub>"+str(dico_compte.get(symbol))+"</sub>"
                    self.mass+=dico_compte.get(symbol)*Atom.dicomass.get(symbol)
                    
            self.lblCompo.setText(compo)
            self.lblMass.setText("%0.3f Da"%self.mass)
        
            
        
    def get_atom(self,node):
        return self.G.nodes[node]["atom"]
    
    def smi2graph(self,smi,G=None,starti=0):
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
            
            self.node_cycles={}
            
            i=starti
            for elt in atoms:
                if elt!="Oz":
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
                            if cycnum in self.node_cycles.keys():                        
                                self.node_cycles[cycnum].append(i)
                            else:
                                self.node_cycles[cycnum]=[i]   
                        
                i+=1
            for cyc in self.node_cycles.keys():
                edge=self.node_cycles.get(cyc)
                #print("cycle "+str(edge[0])+" "+str(edge[1]))
                G.add_edge(*edge,close_cyc=True)                        

                
            self.smi2graph(smi,G)
            
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
                self.smi2graph(smi,G)
            
            elif len(re.findall("[0-9]+",smi))>0:
                nodes=re.findall("[0-9]+",smi)
                for n in range(0,len(nodes)-1):                    
                    G.add_edge(int(nodes[n]),int(nodes[n+1]))        
              
        return G    
    
    
    def egdes_width(self):
        weights=[]
        dble_width=5
        sple_width=1
        for e in self.G.edges():            
            e1=e[0]
            e2=e[1]
            atom1=self.G.nodes[e1]["atom"]
            atom2=self.G.nodes[e2]["atom"]
            if re.search("[0-9]",atom1) and re.search("[0-9]",atom2):
                cyc1=set(re.findall("[0-9]",atom1))
                cyc2=set(re.findall("[0-9]",atom2))
                if len(cyc1.intersection(cyc2))==0:
                    weights.append(dble_width)    
                    self.G.nodes[e1]["valence"]-=1
                    self.G.nodes[e2]["valence"]-=1
                    self.G[e1][e2]["retireH"]=2
                else:
                    self.G[e1][e2]["retireH"]=1
                    weights.append(sple_width)            
            else:
                if re.search("=",atom2):
                    weights.append(dble_width)
                    self.G.nodes[e1]["valence"]-=1
                    self.G.nodes[e2]["valence"]-=1     
                    self.G[e1][e2]["retireH"]=2                    
                else:
                    weights.append(sple_width)
                    self.G[e1][e2]["retireH"]=1
            
        return weights
    
    
    def doubleEvent(self):  
        checked=self.chkDble.isChecked()   
        txt=self.txtsmi.text()
        if checked:
            if txt[0]=="(":
                self.txtsmi.setText("(="+txt[1:])
            else:
                self.txtsmi.setText("="+txt)                
        else:
            if txt[0]=="=":
                self.txtsmi.setText(txt[1:])
            elif self.txtsmi.text()[0:2]=="(=":
                self.txtsmi.setText("("+txt[2:])   
                 
    
    def plotG(self):
        if self.check_notation(self.txtsmi.text()):
            self.smiles=self.txtsmi.text()
            if self.smiles[0]=="(":
                self.G=self.smi2graph("Oz"+self.smiles)
                self.G.remove_node(0)
            else:                
                self.G=self.smi2graph(self.smiles)
           
            self.weights=self.egdes_width()
            self.smiout=self.smiles         
            firstnode=min(self.G.nodes())
            valence=self.G.nodes()[firstnode]["valence"]-len(self.G[firstnode])
            atom=self.G.nodes()[firstnode]["atom"]
            if atom[0]=="=":
                self.chkDble.setEnabled(True) 
            else:                
                self.chkDble.setEnabled(valence>=self.valmin+1) 
            
           
            self.fig.clear()             
            self.ax=self.fig.add_subplot(111)  
            self.ax.set_title('click on atom node to select binder')
            
            w=self.weights
            pos=nx.spring_layout(self.G)
            nx.draw(self.G,pos,self.ax,node_color="white")
            colors= ["white"]*len(self.G.nodes())
            colors[0]="green"
            self.node_binder=min(self.G.nodes())
            nodes = nx.draw_networkx_nodes(self.G, pos,ax=self.ax, node_color=colors,node_size=300, alpha=0.2,
            edgecolors='white')             
            copie=self.G.copy()
           
            for n in copie.nodes():
                copie.nodes[n]['atom']=re.sub("[0-9=]","",self.G.nodes[n]['atom'])
            nx.draw_networkx_labels(copie,pos,nx.get_node_attributes(copie,'atom'), ax=self.ax)        
            nx.draw_networkx_edges(self.G,pos,width=w, ax=self.ax)        
            
            # afficher la figure (au lieu de plt.show())
            self.fig.canvas.draw()        
            
            nodes.set_picker(5)       
            self.fig.canvas.mpl_connect('pick_event', lambda event: self.onpick(event,self.G.nodes()))               
                      
                     
        else:
            print("invalid SMILES")
            self.txtsmi.setText(self.smiles)
        
    
      
    def g2smi(self,G,start,prev,done): 
        """
        creates recursively the smiles notation from a start node of gatom
        """
        neigh=list(G[start])          
        for node in done:
            if node in neigh:
                neigh.remove(node)
        done.append(start)   
        smi=re.sub("=","",G.nodes[start]['atom'])
        if prev!=None:               
            if G[prev][start]["retireH"]==2:
                smi="="+smi        
        
        if len(neigh)>0:
            neigh=sorted(neigh,key=lambda n:G.nodes[n]["valence"])             
                     
            for i in range(len(neigh)-1):  
                smi+="("+self.g2smi(G,neigh[i],start,done)+")"
            smi+=self.g2smi(G,neigh[-1],start,done)
        
        return smi
        
    
    def onpick(self,event,nodes):
        self.btnShowG.setFocus()
        if isinstance(event.artist, PathCollection):
            all_nodes=event.artist
            
            ind = event.ind[0]  # event.ind is a single element array.
            inode=list(nodes)[int(ind)]
            valence=nodes[inode]["valence"]-len(self.G[inode])
            # debug
            #print("Obj,Node: %i,%i" % (int(ind),list(nodes)[int(ind)]))
            #print("Atom: ", nodes[inode]["atom"])
            #print("valence", valence)
            
            self.chkDble.setEnabled(valence>=self.valmin+1)
            
            self.chkDble.setChecked(False) 
            if valence>=self.valmin:
                g=self.G.copy()
                for e in self.G.edges():
                    if "close_cyc" in self.G[e[0]][e[1]]:
                        g.remove_edge(*e)
                paths=list(nx.connected_components(g))
                if len(paths)>1:     
                    pp={}
                    for path in paths:
                        if inode in path:
                            pp=path
                            break                    
                    
                    self.smiout="("+self.g2smi(g.subgraph(pp),inode,None,[])+")"
                    paths.remove(pp)
                    if len(paths)>1:
                        for path in paths[:-1]:
                            start=min(path)
                            self.smiout+="("+self.g2smi(g.subgraph(path),start,None,[])+")"
                    pp=paths[-1:][0]
                    start=min(pp)
                    self.smiout+=self.g2smi(g.subgraph(pp),start,None,[])                      
                    
                    
                else:
                    self.smiout=self.g2smi(g,inode,None,[])
                #print(self.smiout)
                self.txtsmi.setText(self.smiout)
            
            # Set the colours for all the nodes, highlighting the picked node with
            # a different colour:
            colors = [(1, 1, 1)] * len(nodes)
            if valence>=self.valmin:
                colors[ind] = (0, 0.9, 0)
            else:
                colors[ind] = (1, 0, 0)
            all_nodes.set_facecolors(colors)
        
            ## Update the plot to show the change:
            self.fig.canvas.draw()  
            self.node_binder=inode
             
            
        return True  
        
    
    def accept(self):
        """
        """
        if not self.txtsmi.text():
            QMessageBox.about(self, "Warning", "Please fill SMILES notation")          
        else :
            if self.txtsmi.text()[0]=="(":
                link=""
            else:
                link=self.G.nodes[self.node_binder]["atom"]
            
            
                if self.chkDble.isChecked():                
                    link="="+link
            
            name=self.txtname.text()
            #formula=""            
            #for group in re.findall("[A-Z][a-z]?[0-9]*",self.formula.text()):
                #atom=re.sub("[0-9]","",group)
                #nb=re.sub("[^0-9]","",group)
                #if len(nb)==0:
                    #nb="1"
                #formula+=atom+nb
            
            smiles=self.txtsmi.text()
            
            
                
            # attendu  name, formula, linkage,smiles, mass,row=None
            self.window.addModif(self.identifier,name, self.compo, link,smiles,self.mass, self.row)
            QDialog.accept(self)  

           
class SubstituantView(QWidget):

    def __init__(self, *args):
        QWidget.__init__(self,parent=None)   
        self.setGeometry(100, 100, 500, 400)
        self.setWindowTitle("List of the substitutions")        
        self.addSubstitution = []
        self.modifSubstitution = []
        self.deleteSubstitution = []
        
        self.ctrl=args[0]
        data=args[1].copy()
        #for d in args[1]:
            #if d["name"]=="keto":
                #data.remove(d)
            
        self.tableView = TableView(data, len(data), len(data[0]))
        self.tableView.itemSelectionChanged.connect(self.selectedItem)
        self.gridlayout = QGridLayout(self)
        self.gridlayout.addWidget(self.tableView, 0, 0, 5, 3)      
        
        self.addButton = QPushButton("Add")
        self.gridlayout.addWidget(self.addButton, 0, 3, 1, 1)
        self.addButton.clicked.connect(lambda button="add" : self.dialog("add"))
        
        self.modifButton = QPushButton("Modify")
        self.gridlayout.addWidget(self.modifButton, 1, 3, 1, 1)
        self.modifButton.clicked.connect(lambda button="modify" : self.dialog("modify"))
        self.modifButton.setEnabled(False)       
        
        self.supprButton = QPushButton("Remove")
        self.gridlayout.addWidget(self.supprButton, 2, 3, 1, 1)
        self.supprButton.clicked.connect(self.deleteModif)
        self.supprButton.setEnabled(False) 
        
        self.okButton = QPushButton("Ok", self)
        self.okButton.clicked.connect(self.accept)
        self.gridlayout.addWidget(self.okButton, 7, 2, 1, 1)
        
        self.cancelButton = QPushButton("Cancel")
        self.cancelButton.clicked.connect(self.reject)
        self.gridlayout.addWidget(self.cancelButton, 7, 3, 1, 1)
        
        self.setLayout(self.gridlayout)
        
        
    
    def selectedItem(self):
        """
        """        
        indexesSelection = self.tableView.selectedIndexes()
        row = self.tableView.currentRow()
        name=self.tableView.item(row,0).text()
        if len(indexesSelection)==0 :
            self.modifButton.setEnabled(False)
            self.supprButton.setEnabled(False)
        else :
            if re.match("^desoxy$",name,re.IGNORECASE) or re.match("^hydroxy$",name,re.IGNORECASE):
                self.modifButton.setEnabled(False)
                self.supprButton.setEnabled(False)
                
            else:
                self.modifButton.setEnabled(True)
                self.supprButton.setEnabled(True)       

    def dialog (self, button):
        """
        """
        if button=="modify":
            row = self.tableView.currentRow()
            rowdata=self.tableView.getRow(row)
            rowdata["row"]=row
            
            if int(rowdata["identifier"])>=len(self.tableView.data):
                for iline in range(0,len(self.addSubstitution)):                    
                    line=self.addSubstitution[iline]
                    #print((line["identifier"],rowdata["identifier"]))
                    if line["identifier"]==rowdata["identifier"]:
                        self.addSubstitution.remove(line)
            
            SubstituantForm(self,**rowdata)
        else :
            SubstituantForm(self)        
    
    #@pyqtSlot(str, str, str, int)
    def addModif(self, identifier,name, formula, linkage,smiles, mass,row=None):
        """
        """
        # update table of substituents
        identifier=self.tableView.setRow(row,name=name,formula=formula,linkage=linkage,mass=mass,smiles=smiles,identifier=identifier)                    
        
        # store values if accept list edition
        if row == None :        
            self.addSubstitution.append({"identifier":identifier,"name":name, "formula":formula,"linkage": linkage,"smiles":smiles,"mass":mass})
        else :    
            modadd=None
            for add in self.addSubstitution:
                if add["identifier"]==identifier:
                    modadd=add
            if modadd:
                self.addSubstitution.remove(modadd)
                self.addSubstitution.append({"identifier":identifier,"name":name, "formula":formula,"linkage": linkage,"smiles":smiles,"mass":mass})                
            else:
                self.modifSubstitution.append({"identifier":identifier,"name":name, "formula":formula,"linkage": linkage,"smiles":smiles,"mass":mass})
            
       
    
    def deleteModif(self):
        """
        """
        row = self.tableView.currentRow()
        
        identifier=int(self.tableView.item(row,5).text())
        if identifier>0:
            self.deleteSubstitution.append(identifier)
        self.tableView.removeRow(row)
    
    def updateModification(self):
        """
        """
        #ajoute dans le dictionnaire SUBSTITUTIONS les substitutions contenues dans la liste addSubstitution
        if len(self.addSubstitution)!=0 :
            for i in range(len(self.addSubstitution)):
                
                self.ctrl.createSub(self.addSubstitution[i])
                
        else :
            pass
        #supprime dans le dictionnaire SUBSTITUTIONS les substitutions contenues dans la liste deleteSubstitution
        if len(self.deleteSubstitution)!=0 :
            for i in self.deleteSubstitution:
                self.ctrl.removeSub(int(i))
            
        else :
            pass
        #Modifie dans le dictionnaire SUBSTITUTIONS les substitutions contenues dans la liste modifSubstitution
        if len(self.modifSubstitution)!=0 :            
            for i in range(len(self.modifSubstitution)):                
                self.ctrl.modSub(self.modifSubstitution[i])
                
        else :
            pass
    
    def accept(self):
        """
        """
        self.updateModification()
        self.close()
        
    def reject(self):
        """
        """
        self.close()

# numeric sort of table view items
class QTableDoubleItem(QTableWidgetItem):
    """
    class used for numeric sort of table view items
    """
    def __init__ (self, value):
        super(QTableDoubleItem, self).__init__(str('%.3f' % value))    
        
    def __lt__(self, other):
        if (isinstance(other, QTableDoubleItem)):
            selfDataValue  = float(self.data(Qt.EditRole))
            otherDataValue = float(other.data(Qt.EditRole))
            return selfDataValue < otherDataValue
        else:
            return QtGui.QTableWidgetItem.__lt__(self, other)        
                    
class TableView(QTableWidget):
    def __init__(self, data, *args):
        QTableWidget.__init__(self, *args)
        self.data = data
        self.setData()       
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        self.setEditTriggers(QAbstractItemView.NoEditTriggers) #desactive le mode editable
        self.setSelectionBehavior(QAbstractItemView.SelectRows) #permet de selectionner toute la ligne
        self.setSortingEnabled(True)
    
    def setRow(self,row=None,**kwargs):
        """
        """
        self.setSortingEnabled(False)      
        num=None
        if len(kwargs)>0 and kwargs["formula"]!="":
            if row==None:
                rownum=self.rowCount()
                self.insertRow(rownum)
            else:
                rownum=row            
                
            for colname,value in kwargs.items():
                if colname=="name":
                    itm=QTableWidgetItem(value)
                    self.setItem(rownum,0,itm)
                elif colname=="formula":
                    itm=QTableWidgetItem(value)
                    self.setItem(rownum,1,itm)
                elif colname=="mass":
                    itm=QTableDoubleItem(float(value))                    
                    self.setItem(rownum,2,itm)
                elif colname=="smiles":
                    itm=QTableWidgetItem(value)
                    self.setItem(rownum,3,itm)   
               
                elif colname=="linkage":
                    itm=QTableWidgetItem(value)
                    self.setItem(rownum,4,itm)   
                elif colname=="identifier" :
               
                    num=str(value)
                    itm=QTableWidgetItem(num)                    
                    self.setItem(rownum,5,itm)                 
                    self.setColumnHidden(5,True)
                
        self.resizeRowsToContents()
        self.setSortingEnabled(True)
        return num
    
    def getRow(self,rownum):
        """
        """
        rowdata={}
        if rownum in range(self.rowCount()):  
            rowdata["name"]=self.item(rownum,0).text()
            rowdata["formula"]=self.item(rownum,1).text()
            rowdata["mass"]=float(self.item(rownum,2).text())
            rowdata["smiles"]=self.item(rownum,3).text()
            
            rowdata["linkage"]=self.item(rownum,4).text()
            if self.item(rownum,5)!=None:
                rowdata["identifier"]=self.item(rownum,5).text()
            else:
                rowdata["identifier"]=None
        return rowdata
        
    def setData(self):
        """
        """
        header = ["Name", "Composition", "Mass","SMILES","Linker","Id"]
        self.setHorizontalHeaderLabels(header)
        self.verticalHeader().hide()
        
        
        for substit in self.data:
            self.setRow(name=substit["name"],
                        formula=substit["formula"],
                        mass=substit["mass"],
                        smiles=substit["smiles"],
                        linkage=substit["link"],
                        identifier=substit["identifier"])
        
       

class OseDetail(QDialog):
    """
    Definition of ose (form,size,isomery and substitutions)
    """
    DCT_NBC={"5":"pen","6":"hex","7":"hep","8":"oct","9":"non","10":"dec"}    
    ISOCOD={2:"D",1:"L",0:""}
    CANHY=(3,6)
    
    
    def __init__(self,datain, iuipac,anhydro=False,parent=None):
        
        QDialog.__init__(self, parent)   
        self.setWindowFlags(Qt.WindowCloseButtonHint)
        
        row=0
        col=1
        
        
        self.datain=np.array(datain)
        self.dataout=[]
        self._iupac=iuipac
        self.obond=[]
       
        nbC=self.datain.shape[0]
        
        #print(self.datain)
        self.startc=int(np.where(self.datain[:,2]==True)[0].min())+1
        self.endc=int(np.where(self.datain[:,2]==True)[0].max())+1
        
        form=self.endc-self.startc+1
       
        self.gridlayout = QGridLayout(self)
        
        
        self.typeLabel = QLabel("Number of carbons : ")
        self.gridlayout.addWidget(self.typeLabel, row, col, 1, 1)
        
        self.typeComboBox = QComboBox()
        
        for itemnbc in range(5,11):
            self.typeComboBox.addItem(str(itemnbc),itemnbc)
        
        self.typeComboBox.setCurrentText(str(nbC))
        self.typeComboBox.currentTextChanged.connect(self.nbcarbChange)  #modifie le type d ose en fonction du nombre de carbone choisit
        self.gridlayout.addWidget(self.typeComboBox, row, col+1, 1, 1)
        
        self.linecyc=None
        self.startLabel = QLabel("Ring start : ")        
        self.gridlayout.addWidget(self.startLabel, row, col+2, 1, 1)
        
        
        self.startCmbBox = QComboBox()
        
        self.startCmbBox.addItem("1",1)
        self.startCmbBox.addItem("2",2) 
        self.startCmbBox.setCurrentText(str(self.startc))
        self.startCmbBox.currentTextChanged.connect(self.startChange) 
        self.gridlayout.addWidget(self.startCmbBox , row, col+3, 1, 1)             
        
        
        
        #Forme : furane / pyrane
        self.formLabel = QLabel("Form : ")
        self.gridlayout.addWidget(self.formLabel, row, col+4, 1, 1)
        
        self.formCombobox = QComboBox()
        
        self.formCombobox.addItem("pyrane",5)
        self.formCombobox.addItem("furane",4) 
        self.formCombobox.setCurrentIndex(self.formCombobox.findData(form))        
        
        self.formCombobox.currentTextChanged.connect(self.display_cycleline) 
        self.gridlayout.addWidget(self.formCombobox , row, col+5, 1, 1)  
        # chgmt de form impossible si C4 lié?
        if self.datain[3][1]<0:
            self.formCombobox.setDisabled(True)
        
        
        # NEW LINE
        row+=1
        
        # Intra-cyclic bond= anhydro function
        self.chkAnhydro=QCheckBox("3,6-anhydro bond")        
        self.chkAnhydro.stateChanged.connect(self.anhydro_ctrl)                
        self.gridlayout.addWidget(self.chkAnhydro,row,col+5,1,1)           
        
        
        #isomers
        #clockwise: S/L, anticlockwise:R
        self.isomerLbl = QLabel("Isomers:")
        self.gridlayout.addWidget(self.isomerLbl, row, col, 1, 1)        
        
        #Modification  
        self.modifLabel = QLabel("Substituents : ")
        self.gridlayout.addWidget(self.modifLabel, row, col+2, 1, 1) 

        
        
        row+=1
        self.nameLabel = QLabel()
        
            
        #creer les listes descriptives des carbones [[label,iso,mods]]        
        self.combos=[]
        self.addComboS(0,nbC,row,col)
        
        for cmbindex in range(0,len(self.combos)):
            #iso
            iso=self.datain[cmbindex,0]               
            self.combos[cmbindex][1].setCurrentText(OseDetail.ISOCOD[iso])
            
            
            #mods
            mod=self.datain[cmbindex,1]                                                                  
            if mod<0:
                
                self.combos[cmbindex][2].setCurrentText(SubstitutionLibrary.get_subname(SubstitutionLibrary.NOSUBID))     
                
                self.combos[cmbindex][2].setDisabled(True)   
                
                if cmbindex+1 in self.CANHY:
                    self.chkAnhydro.setEnabled(True)
                if cmbindex+1!=self.endc:
                    self.obond.append(cmbindex)
            else:
                self.combos[cmbindex][2].setCurrentText(SubstitutionLibrary.get_subname(mod))     
                if mod!=SubstitutionLibrary.NOSUBID and cmbindex+1 in self.CANHY:
                    self.chkAnhydro.setEnabled(False)
                    
        
        self.enable_formchange()                          
        self.display_cycleline()
        # has to be after comboboxes in order to disable substit for carb 3 and 6
        self.chkAnhydro.setChecked(anhydro)
        self.enable_anhydro()
        
        row+=nbC+1
        
        
        self.gridlayout.addWidget(self.nameLabel,row,col+1,1,1)     
        
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setLayoutDirection(Qt.LeftToRight)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.gridlayout.addWidget(self.buttonBox, row, col+3, 1, 2)
        
        self.setWindowTitle("Ose description")
        self.setWindowModality(Qt.ApplicationModal)
       
    
    def addComboS(self,start,stop,row,col):
        """
        """
        for i in range(start,stop):       
            row+=1
            self.carbLabel = QLabel(self)
            self.carbLabel.setText("C"+str(i+1))
            self.gridlayout.addWidget(self.carbLabel, row, col, 1, 1)    
            self.isoCombobox = QComboBox()
            
            for k,v in OseDetail.ISOCOD.items():
                self.isoCombobox.addItem(v,k)
            if i==stop-1: 
            
                self.isoCombobox.setCurrentText("")
                self.isoCombobox.setDisabled(True)
            self.gridlayout.addWidget(self.isoCombobox, row, col+1, 1, 1)  
            self.isoCombobox.setObjectName('Iso%d' % i)    #nomme chacun des combobox
            self.isoCombobox.currentTextChanged.connect(self.iupac_name)    #modifie le nom IUPAC si l utilisateur a fait un changement dans les combobox
            
                       
            
            self.modifCombobox = QComboBox()            
            
            for s in SubstitutionLibrary.SUBSTITUTIONS:
                self.modifCombobox.addItem(s.name,s.identifier)  
                #if s.name!="keto":
                    #self.modifCombobox.addItem(s.name,s.identifier)  
            
            self.modifCombobox.setCurrentText(SubstitutionLibrary.get_subname(SubstitutionLibrary.NOSUBID))      
            self.modifCombobox.setObjectName('Modif%d' % i)   #nomme chacun des combobox            
            if i+1==self.CANHY[1]:
                self.modifCombobox.currentTextChanged.connect(self.enable_anhydro)
            
            
            
            self.gridlayout.addWidget(self.modifCombobox, row, col+3, 1, 1)    #ajoute le combobox dans la grille  
                
            self.combos.append([self.carbLabel,self.isoCombobox,self.modifCombobox])    

    def display_cycleline(self): 
        """
        """
        if self.linecyc:
            self.gridlayout.removeWidget(self.linecyc)
        else:
            self.linecyc=QFrame()
            self.linecyc.setFrameShape(QFrame.VLine)                   
            self.linecyc.setFrameShadow(QFrame.Sunken)            
        column=0
        
        
        startcycle=self.combos[self.startCmbBox.currentData()-1][0]        
        indexstart=self.gridlayout.indexOf(startcycle)
        posistart=self.gridlayout.getItemPosition(indexstart)
        
        endcycle=self.combos[(self.startCmbBox.currentData()+self.formCombobox.currentData()-1)-1][0] 
        
        indexend=self.gridlayout.indexOf(endcycle)
        posiend=self.gridlayout.getItemPosition(indexend)        
        end=posiend[0]-posistart[0]+1
        self.gridlayout.addWidget(self.linecyc,posistart[0],column,end,1)    
       
        for icarb in range(len(self.combos)):
            self.combos[icarb][1].setDisabled(False)
            self.combos[icarb][2].setDisabled(False)
            # carbon ending cycle
            if icarb==(self.startCmbBox.currentData()+end-1)-1:
                self.combos[icarb][2].setDisabled(True)
                self.combos[icarb][2].setCurrentText(SubstitutionLibrary.get_subname(SubstitutionLibrary.NOSUBID))
          
          
            # not ending cycle               
            elif icarb< self.datain.shape[0] :
                if self.datain[icarb,1]<0 :
                    if icarb==self.endc-1:
                        self.combos[icarb][2].setCurrentText(SubstitutionLibrary.get_subname(SubstitutionLibrary.NOSUBID))
                    # keep osidic binding status
                    else:
                        self.combos[icarb][2].setDisabled(True)
                else:
                    self.combos[icarb][2].setCurrentText(SubstitutionLibrary.get_subname(self.datain[icarb,1]))
                
        
        self.combos[icarb][1].setDisabled(True)        
        self.iupac_name()
        self.chkAnhydro.setChecked(False)
        self.enable_anhydro()
        
    def startChange(self):
        """
        """        
        cycend=self.startCmbBox.currentData()+self.formCombobox.currentData()-1
        if cycend<=len(self.combos):
            self.display_cycleline()    
            print("start changed")
        else:
            self.startCmbBox.setCurrentText("1")

    def nbcarbChange(self, osetype): 
        """
        """          
        if self.typeComboBox.currentData()>=self.startCmbBox.currentData()+self.formCombobox.currentData()-1:                           
            #modifie les composants de la fenêtre de dialogue en fonction du type du nombre de carbone selectionne
            if self.typeComboBox.currentData()<len(self.combos):
                rm_tab=[]
                for icarb in range(self.typeComboBox.currentData(),len(self.combos)):
                    rm_tab.append(self.combos[icarb])
                    for elt in self.combos[icarb]:
                        self.gridlayout.removeWidget(elt)
                        elt.deleteLater()
                    
                    
                for tab in rm_tab:
                    self.combos.remove(tab)
                
                self.display_cycleline()
                self.enable_anhydro()
                
            elif self.typeComboBox.currentData() >len(self.combos):            
                lastlabel=self.combos[len(self.combos)-1][0]
                lastindex=self.gridlayout.indexOf(lastlabel)
                posistart=self.gridlayout.getItemPosition(lastindex)            
                
                
                self.gridlayout.removeWidget(self.nameLabel)
                self.gridlayout.removeWidget(self.buttonBox)
                self.addComboS(len(self.combos),self.typeComboBox.currentData(),posistart[0],posistart[1])
    
                row=self.gridlayout.rowCount()
                
                self.gridlayout.addWidget(self.nameLabel,row , posistart[1], 1, 2)
                self.gridlayout.addWidget(self.buttonBox,row , posistart[1]+2, 1, 2)
    
                if self.typeComboBox.currentData()==self.formCombobox.currentData():
                    self.startCmbBox.setCurrentText("1")
                self.display_cycleline()
                self.enable_anhydro()
            
        else:
            self.typeComboBox.setCurrentText(str(len(self.combos)))

        self.enable_formchange()
        self.enable_anhydro()
       
    
    
    def enable_formchange(self):
        curform=self.formCombobox.currentData()
        
        if curform==4:                      
            if self.typeComboBox.currentData()>=5:
                icarb=self.startCmbBox.currentData()+3
                try:
                    ok=(icarb in self.obond) or self.combos[icarb][2].currentData()!=SubstitutionLibrary.NOSUBID
                    self.formCombobox.setDisabled(ok)                
                except IndexError:
                    print((len(self.combos),icarb))
            else:                               
                self.formCombobox.setDisabled(False)
        elif curform==5:
            icarb=self.startCmbBox.currentData()+2
            
            ok=(icarb in self.obond) or self.combos[icarb][2].currentData()!=SubstitutionLibrary.NOSUBID
            self.formCombobox.setDisabled(ok)           
    
    
        
    def cancel(self):
        """
        """
        self.dataout=[]
    
    def accept(self):
        """
        """
        #recupere le nombre de carbone, le type et la forme de l'ose
        nbC = self.typeComboBox.currentData()
        form = self.formCombobox.currentData()
        istartcycle=self.startCmbBox.currentData()-1
        iendcycle=self.startCmbBox.currentData()-1+form
          
        
        #parcourt les listes deroulantes des substitutions et de l isomerie et modifie le tableau des modifications de l ose
        iso = []
        modif = []
        dataout=[]
        for icarb in range(nbC):
            incycle=False
            
            if icarb in range(istartcycle,iendcycle):
                incycle=True
            if icarb==iendcycle-1:
                dataout.append([self.combos[icarb][1].currentData(),-1,incycle])
            # keep osidic bond
            elif icarb in self.obond:
                dataout.append([self.combos[icarb][1].currentData(),-1,incycle])                
            else:
                dataout.append([self.combos[icarb][1].currentData(),self.combos[icarb][2].currentData(),incycle])
            
        
                
        if np.array_equal(self.datain,np.array(dataout)):
            self.dataout=[]            
        else:
            self.dataout=dataout.copy()
            
        QDialog.accept(self)
    
    def iupac_name(self):
        """
        """
        if len(self.combos)>0:
            name=""
            startc=self.startCmbBox.currentData()
            endc=startc+self.formCombobox.currentData()-1
            if endc==len(self.combos):
                endiso=endc-1
            else:
                endiso=endc
                
            if self.combos[startc-1][1].currentText()=="L":
                name+="b"
            elif self.combos[startc-1][1].currentText()=="D":
                name+="a"
            else:
                name+="x"
            name+='-'
            
            if self.combos[endiso-1][1].currentText()=="":
                name+="x"
            else:
                name+=self.combos[endiso-1][1].currentText().lower()        
            
            iso=""  
            for cnum in range(startc+1,endiso):               
                if self.combos[cnum-1][1].currentText()=="":
                    iso+="-"
                else:
                    iso+=self.combos[cnum-1][1].currentText()
            
            
            if iso in self._iupac:
                name+=self._iupac[iso]
            else:
                name=""
            
            
            
            self.nameLabel.setText(name)
    
  
    
    def enable_anhydro(self):                
        modfrom=self.combos[self.CANHY[0]-1][2].currentData()
        
        if len(self.combos)==self.CANHY[1]:
            modto=self.combos[self.CANHY[1]-1][2].currentData()
            if self.CANHY[1]-1 in self.obond or self.CANHY[0]-1 in self.obond:
                self.chkAnhydro.setEnabled(False)           
            elif modfrom==SubstitutionLibrary.NOSUBID and modto==SubstitutionLibrary.NOSUBID:
                self.chkAnhydro.setEnabled(True)
            else:
                self.chkAnhydro.setChecked(False)
                self.chkAnhydro.setEnabled(False)
        else:
            self.chkAnhydro.setChecked(False)
            self.chkAnhydro.setEnabled(False)
        
    def anhydro_ctrl(self):        
        if self.chkAnhydro.isChecked():            
            Logger.debug("connect C3 and C6")   
            self.combos[self.CANHY[0]-1][2].setEnabled(False)
            self.combos[self.CANHY[1]-1][2].setEnabled(False)                           
        else:
            Logger.debug("dis-connect C3 and C6")   
            self.combos[self.CANHY[0]-1][2].setEnabled(True)
            if len(self.combos)==self.CANHY[1]:
                self.combos[self.CANHY[1]-1][2].setEnabled(True)  

## next release preparation
#class DialogPeakList(QDialog):        
    #def __init__(self,name=None,number=None,precision=None):
        #super(DialogPeakList, self).__init__()
        #self.name=name
        #self.number=number
        #self.precision=precision
        #self.peaklist=()
        
class DialogMS(QDialog):  
    default_ion="+H"
    
    def __init__(self,random,ratio,dicoH=None,peaklist=None):
        """
        build Qt window from ressource file "SpectrometryIntensity_ion.ui"
        """
        super(DialogMS, self).__init__()
    
        loadUi("Ressources"+os.path.sep+'SpectrometryIntensity_ion.ui', self)
        
        self.cmbAdduct.addItem("+H",Atom.mass("H"))
        self.cmbAdduct.addItem("-H",-1*Atom.mass("H"))
        self.cmbAdduct.addItem("+Na",Atom.mass("Na"))
        self.cmbAdduct.addItem("+K",Atom.mass("K"))
        self.cmbAdduct.addItem("+Cs",Atom.mass("Cs"))
        self.cmbAdduct.addItem("+Li",Atom.mass("Li"))
        self.cmbAdduct.setCurrentText(self.default_ion)
        self.cmbAdduct.currentTextChanged.connect(self.change_adduct)
        
        
        self.Acount=0
        self.Xcount=0
        self.dataMZ=None
        self.ionlabels={}
        
        grid=self.gridSpectrum
       
        self.dico_Acontrols={0:(self.spinBox_a01,self.lblA_between,self.spinBox_a02,self.btnAplus,self.spinBox_Amin,self.spinBox_Amax,self.sliderA,self.lblRatioA)}
        self.dico_Xcontrols={0:(self.spinBox_x01,self.lblX_between,self.spinBox_x02,self.btnXplus,self.spinBox_Xmin,self.spinBox_Xmax,self.sliderX,self.lblRatioX)}
        
        self.chkRatio.toggled.connect(self.__toogleRatio__)
        self.chkRandom.toggled.connect(self.__toogleRandom__)
               
        self.btnAplus.clicked.connect(self.__add__A)   
        self.btnXplus.clicked.connect(self.__add__X) 
        
        self.lineRatioMax.textChanged.connect(self.__maxchanged__)
        
        self.lblRatioA.setText(str(self.sliderA.value()))
        self.lblRatioX.setText(str(self.sliderX.value()))
        self.lblRatioB.setText(str(self.sliderB.value()))
        self.lblRatioY.setText(str(self.sliderY.value()))
        self.lblRatioC.setText(str(self.sliderC.value()))
        self.lblRatioZ.setText(str(self.sliderZ.value()))
        
        self.connectLineControls(self.spinBox_Amin,self.spinBox_Amax,self.sliderA,self.lblRatioA,self.spH_A,self.spinBox_a01,self.spinBox_a02)
        self.connectLineControls(self.spinBox_Xmin,self.spinBox_Xmax,self.sliderX,self.lblRatioX,self.spH_X,self.spinBox_x01,self.spinBox_x02)
        self.connectLineControls(self.spinBox_Bmin,self.spinBox_Bmax,self.sliderB,self.lblRatioB,self.spH_B)
        self.connectLineControls(self.spinBox_Cmin,self.spinBox_Cmax,self.sliderC,self.lblRatioC,self.spH_C)
        self.connectLineControls(self.spinBox_Ymin,self.spinBox_Ymax,self.sliderY,self.lblRatioY,self.spH_Y)
        self.connectLineControls(self.spinBox_Zmin,self.spinBox_Zmax,self.sliderZ,self.lblRatioZ,self.spH_Z)        
        
        sptest=QSpinBox()  
        if dicoH==None:
            self.set_cid()
        else:
            
            if "A" in dicoH:
                self.spH_A.setValue(dicoH["A"])
            if "X" in dicoH:
                self.spH_X.setValue(dicoH["X"])        
            if "B" in dicoH:
                self.spH_B.setValue(dicoH["B"])                
            if "C" in dicoH:
                self.spH_C.setValue(dicoH["C"])                
            if "Y" in dicoH:
                self.spH_Y.setValue(dicoH["Y"])                                
            if "Z" in dicoH:
                self.spH_Z.setValue(dicoH["Z"])            
        
        
       
        
        self.random=random
        self.ratio=ratio
        
        
        
        fig = Figure(figsize=(9,2), dpi=96)
        
        self.canvas=FigureCanvas(fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        self.canvas.setParent(self.framePlot)
        self.ax = self.canvas.figure.add_subplot(111)      
        self.ax.set_title('spectrum preview')    
        on_move_id = self.canvas.mpl_connect('motion_notify_event', self.on_move)
                
    
    
    def change_adduct(self):
        self.adduct=Atom.mass(re.sub("[+-]","",self.cmbAdduct.currentText()))
        if re.sub("[^+-]","",self.cmbAdduct.currentText())=="-":
            self.adduct=self.adduct*-1
        self.set_cid()

        self.plotData()
    
    def set_cid(self):        
        self.apply_dicoH({"A":1,"X":1,"B":-1,"C":1,"Y":1,"Z":-1})
        
        
    def reset_nh(self):        
        self.apply_dicoH({"A":0,"X":0,"B":0,"C":0,"Y":0,"Z":0})
        
    def apply_dicoH(self,dicoH):                 
        self.spH_A.valueChanged.disconnect()
        self.spH_B.valueChanged.disconnect()
        self.spH_C.valueChanged.disconnect()
        self.spH_X.valueChanged.disconnect()
        self.spH_Y.valueChanged.disconnect()
        self.spH_Z.valueChanged.disconnect()        
        
        self.spH_A.setValue(dicoH["A"])        
        self.spH_X.setValue(dicoH["X"]) 
        self.spH_B.setValue(dicoH["B"]) 
        self.spH_C.setValue(dicoH["C"]) 
        self.spH_Y.setValue(dicoH["Y"])         
        self.spH_Z.setValue(dicoH["Z"])   
        self.__changeMass__()
        
        self.spH_A.valueChanged.connect(self.__changeMass__)
        self.spH_B.valueChanged.connect(self.__changeMass__)
        self.spH_C.valueChanged.connect(self.__changeMass__)
        self.spH_X.valueChanged.connect(self.__changeMass__)
        self.spH_Y.valueChanged.connect(self.__changeMass__)
        self.spH_Z.valueChanged.connect(self.__changeMass__)
        
    def connectLineControls(self,spinRangeMin,spinRangeMax,sliderRatio,lblRatio,spH=None,spinCfrom=None,spinCto=None):
        """
        """
        if spinCfrom!=None and spinCto!=None:
            spinCfrom.valueChanged.connect(lambda:self.__spinFromCarbons__(spinCfrom,spinCto))            
            spinCto.valueChanged.connect(lambda:self.__spinToCarbons__(spinCfrom,spinCto))   
        
        spinRangeMin.valueChanged.connect(lambda:self.__changeRangeMin__(spinRangeMin,spinRangeMax))
        spinRangeMax.valueChanged.connect(lambda:self.__changeRangeMax__(spinRangeMin,spinRangeMax))
        sliderRatio.valueChanged.connect(lambda:lblRatio.setText(str(sliderRatio.value())))
        sliderRatio.valueChanged.connect(self.__changeRatio__)
        
        if spH!=None:
            spH.valueChanged.connect(self.__changeMass__)
    
    def __update_filters__(self):
        """
        """
        aions=[]
        xions=[]
        self.random.clear()
        self.ratio.clear()
        if self.dataMZ!=None:
            # reset table
            for ion in self.dataMZ:
                ion["intensity"]=0        
        #spin,spin,btnMinus,spin,spin,slider,label
        for acount in range(0,self.Acount+1):
        
            sep=self.dico_Acontrols[acount][1].property("alternative")
            cfrom=self.dico_Acontrols[acount][0].value()
            cto=self.dico_Acontrols[acount][2].value()
            if sep=="-":                
                regexp="\\^\\{[%i%s5],[%i%s%i]\\}A"%(cfrom,sep,cfrom,sep,cto)
            else:
                regexp="\\^\\{%i%s%i\\}A"%(cfrom,sep,cto)            
            aions.append(regexp)
        for xcount in range(0,self.Xcount+1):
            sep=self.dico_Xcontrols[xcount][1].property("alternative")
            cfrom=self.dico_Xcontrols[xcount][0].value()
            cto=self.dico_Xcontrols[xcount][2].value()
            if sep=="-":                
                regexp="\\^\\{[%i%s5],[%i%s%i]\\}X"%(cfrom,sep,cfrom,sep,cto)
            else:
                regexp="\\^\\{%i%s%i\\}X"%(cfrom,sep,cto)
            xions.append(regexp)        
        
        if self.chkRandom.isChecked():            
            for acount in range(0,self.Acount+1):                
                self.random[aions[acount]]=[self.dico_Acontrols[acount][4].value(),self.dico_Acontrols[acount][5].value()]
            for xcount in range(0,self.Xcount+1):                
                self.random[xions[xcount]]=[self.dico_Xcontrols[xcount][4].value(),self.dico_Xcontrols[xcount][5].value()]                
        
            self.random["B"]=[self.spinBox_Bmin.value(),self.spinBox_Bmax.value()]
            self.random["Y"]=[self.spinBox_Ymin.value(),self.spinBox_Ymax.value()]
            self.random["C"]=[self.spinBox_Cmin.value(),self.spinBox_Cmax.value()]
            self.random["Z"]=[self.spinBox_Zmin.value(),self.spinBox_Zmax.value()]
        else:                   
            maxi=float(self.lineRatioMax.text())
            for acount in range(0,self.Acount+1):                
                self.ratio[aions[acount]]=self.dico_Acontrols[acount][6].value()*maxi/100
            for xcount in range(0,self.Xcount+1):                
                self.ratio[xions[xcount]]=self.dico_Xcontrols[xcount][6].value()*maxi/100            
            self.ratio["B"]=self.sliderB.value()*maxi/100
            self.ratio["Y"]=self.sliderY.value()*maxi/100
            self.ratio["C"]=self.sliderC.value()*maxi/100
            self.ratio["Z"]=self.sliderZ.value()*maxi/100    
            #print(self.ratio)        
        self.plotData()
        
    
    def accept(self): 
        """
        """
        QDialog.accept(self)
    
    def reject(self):
        """
        """
        QDialog.reject(self)
        
    def setMZ(self,data):
        """
        """
        self.dataMZ=data
        
        self.__update_filters__()

   
    def on_move(self,event): 
        if self.annotation!=None:
            self.annotation.remove()
            self.annotation=None
            
        for l in self.points:            
            if l.contains(event)[0]==True:     
                x=event.xdata
                y=event.ydata                
                point=None
                dist=100
                
                # get closest point
                for xyline in l._xy:                    
                    if xyline[1]>0 :
                        if abs(x-xyline[0])==dist:
                            if point==None:
                                dist=abs(x-xyline[0])
                                point=(xyline[0],xyline[1])
                            elif xyline[1]>point[1]:
                                dist=abs(x-xyline[0])
                                point=(xyline[0],xyline[1])                                
                                
                        elif abs(x-xyline[0])<dist:
                            dist=abs(x-xyline[0])
                            point=(xyline[0],xyline[1])
                            
                
                    
                if point !=None:  
                    lbl=""
                    alls=self.ionlabels.get(point[0])
                    maxlbl=3
                    for i in range(0,len(alls)):
                        if i==maxlbl:
                            lbl+="..."
                            break
                        else:
                            lbl+="$"+alls[i]+"$ "
                    
                        
                    self.annotation=self.ax.annotate("%0.2f\n%s" % (point[0],str(lbl)),xy=(point[0], point[1]), xycoords='data',xytext=(x + 10, y), 
                                               textcoords='data',
                                               horizontalalignment="left",
                                               arrowprops=dict(arrowstyle="->"),
                                               bbox=dict(boxstyle="round", facecolor="w",edgecolor="0.5", alpha=0.9)) 
                    self.annotation.set_visible(True)
                
        self.canvas.draw() 
         
        
            
    def plotData(self):    
        """
        """

        adduct=self.cmbAdduct.currentData()
        self.ionlabels.clear()
        self.annotation=None
        
        if self.dataMZ!=None:
            xydata=[]
            for ion in self.dataMZ:            
                if self.chkRatio.isChecked():                
                    for regexpion,ratio in self.ratio.items():                      
                        if re.search(regexpion,ion["name"])!=None:                          
                            ion["intensity"]=ratio
                        
                else:                               
                    for regexpion,rangeI in self.random.items():                        
                        if re.search(regexpion,ion["name"])!=None:                  
                            ion["intensity"]=random.uniform(rangeI[0],rangeI[1])                        
                          
            for ion in self.dataMZ:
                typion=re.sub("[^AXZYBC]","",ion["name"])
                if typion in self.dicoH:
                    mz=adduct+ion["mz"]+self.dicoH[typion]*Atom.mass("H")
                else:
                    mz=adduct+ion["mz"]
                xydata.append([mz-0.01,0])
                xydata.append([mz,ion["intensity"]])
                xydata.append([mz+0.01,0])  
                if ion["intensity"]>0:
                    # to display ion name(s) on mouse move over displayed chart
                    if mz in self.ionlabels:
                        self.ionlabels[mz].append(ion["name"])
                    else:
                        self.ionlabels[mz]=[ion["name"]]
                        
            xydata=np.array(xydata)
            xydata=xydata[xydata[:,0].argsort()]
            
            
            self.ax.clear()
            self.points=self.ax.plot(xydata[:,0],xydata[:,1],'g-')                         
            self.canvas.draw() 
        else:
            Logger.debug("MS/MS preview: nothing to plot")
        
    def __maxchanged__(self):
        """
        """
        if self.chkRatio.isChecked():
            self.__update_filters__()
    
    def __toogleRandom__(self):  
        """
        """
        self.chkRatio.setChecked(not self.chkRandom.isChecked())           
        self.__update_filters__()
            
    def __toogleRatio__(self):  
        """
        """
        self.chkRandom.setChecked(not self.chkRatio.isChecked())
        self.__update_filters__()
    
    def __spinFromCarbons__(self,spinFrom,spinTo):
        """
        """
        if spinFrom.value()>=spinTo.value():
            spinFrom.setValue(spinTo.value()-1)
        self.__update_filters__()
        
    def __spinToCarbons__(self,spinFrom,spinTo):   
        """
        """
        if spinTo.value()<=spinFrom.value():
            spinTo.setValue(spinFrom.value()+1)
        self.__update_filters__()
    
    def __changeRangeMin__(self,spinMin,spinMax):
        """
        """
        
        if self.chkRandom.isChecked():
            self.__update_filters__()
            
    def __changeRangeMax__(self,spinMin,spinMax):
        """
        """
       
        if self.chkRandom.isChecked():
            self.__update_filters__()    
    
    def __changeMass__(self):
        """
        """        
        self.dicoH={}
        self.dicoH["A"]=self.spH_A.value()
        self.dicoH["X"]=self.spH_X.value()
        self.dicoH["B"]=self.spH_B.value()
        self.dicoH["C"]=self.spH_C.value()
        self.dicoH["Y"]=self.spH_Y.value()
        self.dicoH["Z"]=self.spH_Z.value()
        self.plotData()
    
    def __changeRatio__(self):
        """
        """        
        if self.chkRatio.isChecked():
            self.__update_filters__()    
        
    def __add__A(self):        
        """
        """        
        self.Acount+=1        
        blurb1=QSpinBox()
        
        blurb1.setSingleStep(self.dico_Acontrols[0][0].singleStep())
        blurb1.setMinimum(self.dico_Acontrols[0][0].minimum())
        blurb1.setMaximum(self.dico_Acontrols[0][0].maximum())
        
        self.gridA.addWidget(blurb1, self.Acount, 2, 1, 1)
        lblbetween=QLabel("&")
        lblbetween.setProperty("alternative",",")
        self.gridA.addWidget(lblbetween, self.Acount, 3, 1, 1)
        
        blurb2=QSpinBox()
        blurb2.setSingleStep(self.dico_Acontrols[0][2].singleStep())
        blurb2.setMinimum(self.dico_Acontrols[0][2].minimum())
        blurb2.setMaximum(self.dico_Acontrols[0][2].maximum())  
        blurb2.setValue(blurb2.maximum())
        
        self.gridA.addWidget(blurb2, self.Acount, 4, 1, 1)
        
        
        
        blurb3=QSpinBox()
        blurb3.setSingleStep(self.dico_Acontrols[0][4].singleStep())
        blurb3.setMinimum(self.dico_Acontrols[0][4].minimum())
        blurb3.setMaximum(self.dico_Acontrols[0][4].maximum())          
        
        self.gridA_rand.addWidget(blurb3, self.Acount, 0, 1, 1)
        
        blurb4=QSpinBox()
        blurb4.setSingleStep(self.dico_Acontrols[0][5].singleStep())
        blurb4.setMinimum(self.dico_Acontrols[0][5].minimum())
        blurb4.setMaximum(self.dico_Acontrols[0][5].maximum())         
        blurb4.setValue(blurb4.maximum())
        
        self.gridA_rand.addWidget(blurb4, self.Acount, 1, 1, 1)        
        
        blurb5=QSlider()
        blurb5.setOrientation(self.dico_Acontrols[0][6].orientation())
        blurb5.setSingleStep(self.dico_Acontrols[0][6].singleStep())        
        blurb5.setMinimum(self.dico_Acontrols[0][6].minimum())
        blurb5.setMaximum(self.dico_Acontrols[0][6].maximum())
        blurb5.setSliderPosition(self.dico_Acontrols[0][6].sliderPosition())
        blurb5.setTickPosition(self.dico_Acontrols[0][6].tickPosition())
        blurb5.setTickInterval(self.dico_Acontrols[0][6].tickInterval())
        
        self.gridA_ratio.addWidget(blurb5, self.Acount, 0, 1, 1)    
        
        blurb6=QLabel(str(blurb5.value()))
        self.gridA_ratio.addWidget(blurb6, self.Acount, 1, 1, 1)
        
        
        btnMinus=QPushButton("-")                          
        btnMinus.setFixedSize(self.btnAplus.width(),self.btnAplus.height())
        self.gridA.addWidget(btnMinus,self.Acount,5,1,1)
        
        self.dico_Acontrols[self.Acount]=(blurb1,lblbetween,blurb2,btnMinus,blurb3,blurb4,blurb5,blurb6)
        self.connectLineControls(blurb3,blurb4,blurb5,blurb6,blurb1,blurb2)
        btnMinus.clicked.connect(lambda:self.__remove__(self.Acount,True))        
    
    def __add__X(self):        
        """
        """        
        self.Xcount+=1
        
        blurb1=QSpinBox()
        
        blurb1.setSingleStep(self.dico_Xcontrols[0][0].singleStep())
        blurb1.setMinimum(self.dico_Xcontrols[0][0].minimum())
        blurb1.setMaximum(self.dico_Xcontrols[0][0].maximum())
        
        self.gridX.addWidget(blurb1, self.Xcount, 2, 1, 1)
        
        lblbetween=QLabel("&")
        lblbetween.setProperty("alternative",",")
        
        self.gridX.addWidget(lblbetween, self.Xcount, 3, 1, 1)        
        
        
        blurb2=QSpinBox()
        blurb2.setSingleStep(self.dico_Xcontrols[0][2].singleStep())
        blurb2.setMinimum(self.dico_Xcontrols[0][2].minimum())
        blurb2.setMaximum(self.dico_Xcontrols[0][2].maximum())   
        blurb2.setValue(blurb2.maximum())
        
        self.gridX.addWidget(blurb2, self.Xcount, 4, 1, 1)
        
        blurb3=QSpinBox()
        blurb3.setSingleStep(self.dico_Xcontrols[0][4].singleStep())
        blurb3.setMinimum(self.dico_Xcontrols[0][4].minimum())
        blurb3.setMaximum(self.dico_Xcontrols[0][4].maximum())          
        
        self.gridX_rand.addWidget(blurb3, self.Xcount, 0, 1, 1)
        
        blurb4=QSpinBox()
        blurb4.setSingleStep(self.dico_Xcontrols[0][5].singleStep())
        blurb4.setMinimum(self.dico_Xcontrols[0][5].minimum())
        blurb4.setMaximum(self.dico_Xcontrols[0][5].maximum())         
        blurb4.setValue(blurb4.maximum())
        self.gridX_rand.addWidget(blurb4, self.Xcount, 1, 1, 1)        
        
        blurb5=QSlider()
        blurb5.setOrientation(self.dico_Xcontrols[0][6].orientation())
        blurb5.setSingleStep(self.dico_Xcontrols[0][6].singleStep())        
        blurb5.setMinimum(self.dico_Xcontrols[0][6].minimum())
        blurb5.setMaximum(self.dico_Xcontrols[0][6].maximum())
        blurb5.setSliderPosition(self.dico_Xcontrols[0][6].sliderPosition())
        blurb5.setTickPosition(self.dico_Xcontrols[0][6].tickPosition())
        blurb5.setTickInterval(self.dico_Xcontrols[0][6].tickInterval())
        
        self.gridX_ratio.addWidget(blurb5, self.Xcount, 0, 1, 1)    
        
        blurb6=QLabel(str(blurb5.value()))
        self.gridX_ratio.addWidget(blurb6, self.Xcount, 1, 1, 1)        
        
        
        btnMinus=QPushButton("-")
        btnMinus.setFixedSize(self.btnAplus.width(),self.btnAplus.height())
        self.gridX.addWidget(btnMinus,self.Xcount,5,1,1)
        
        self.dico_Xcontrols[self.Xcount]=(blurb1,lblbetween,blurb2,btnMinus,blurb3,blurb4,blurb5,blurb6)        
        self.connectLineControls(blurb3,blurb4,blurb5,blurb6,blurb1,blurb2)
        btnMinus.clicked.connect(lambda:self.__remove__(self.Xcount,False))        
    
    
    
    def __remove__(self,row,isA):
        """
        """        
        items=()
        if isA:
            items=self.dico_Acontrols[row]
            for i in range(0,4):
                self.gridA.removeWidget(items[i])            
                         
            self.gridA_rand.removeWidget(items[i])
            self.gridA_rand.removeWidget(items[i+1])                
            self.gridA_ratio.removeWidget(items[i+2]) 
            self.gridA_ratio.removeWidget(items[i+3]) 
            self.dico_Acontrols.pop(row)
            self.Acount-=1              
        else:
            items=self.dico_Xcontrols[row]
            for i in range(0,4):
                self.gridX.removeWidget(items[i])            
                         
            self.gridX_rand.removeWidget(items[i])
            self.gridX_rand.removeWidget(items[i+1])                
            self.gridX_ratio.removeWidget(items[i+2]) 
            self.gridX_ratio.removeWidget(items[i+3]) 
            self.dico_Xcontrols.pop(row)
            self.Xcount-=1                          
            
        for w in items:
            w.deleteLater()
        
        self.__update_filters__()
    

class EdgeLine(QGraphicsLineItem):
   
    
    def __init__(self,x1,y1,x2,y2,identifier,parent=None):
        QGraphicsLineItem.__init__(self, parent)
        self.setPen(QPen(QColor("black"), 2))  
        self.x1=x1
        self.x2=x2
        self.y1=y1
        self.y2=y2
        self.setLine(self.x1, self.y1, self.x2, self.y2) 
        self.signal=SignalObject()
        self.identifier=identifier
    
    def setLine(self,x1,y1,x2,y2):
        QGraphicsLineItem.setLine(self,x1,y1,x2,y2)
        self.x1=x1
        self.x2=x2
        self.y1=y1
        self.y2=y2        
        
    def mouseDoubleClickEvent(self,event): 
        """
        """
        self.signal.dble_click(self.identifier)
       
    def notify(self, dx, dy, position):   
        """
        """
        self.prepareGeometryChange()
        
        if position==1:
            self.x1+=dx
            self.y1+=dy
        elif position==2:
            self.x2+=dx
            self.y2+=dy
            
        self.setLine(self.x1, self.y1, self.x2, self.y2)    
        



class SignalObject(QObject):
    """    
    """
    itemDoubleClicked=pyqtSignal()    
    itemContextMenu=pyqtSignal()  
    itemSupprKey=pyqtSignal()
    
    value=0
    
    def dble_click(self,i):
        """
        """
        SignalObject.value=i
        #print("------------------------------")
        #print(i)
        
        self.itemDoubleClicked.emit()
        #print("------------------------------")
        
    def right_click(self,i):
        """
        """
        SignalObject.value=i
        #print("--right click-----------------")
        #print(i)
        
        self.itemContextMenu.emit()
        #print("------------------------------")
    
    
    def suppr_press(self):
        """
        """
        SignalObject.value=None
        self.itemSupprKey.emit()
        
      
class OseImage(QGraphicsPixmapItem):
    IMG_DIR="Image"+os.path.sep
    WIDTH=80
    HEIGHT=80
    DEFAULT_COLOR="gray"
    
    def __init__(self,identifier,nct,ncc ,startc=1,anhydro=False,parent=None):
                
        QGraphicsPixmapItem.__init__(self, parent)
        
        self.setFlags(QGraphicsItem.ItemIsSelectable |
                          QGraphicsItem.ItemIsMovable |
                          QGraphicsItem.ItemSendsGeometryChanges)
        
        
        
        self.setAcceptHoverEvents(True)  # hover events are used to change mouse cursor                              
        
        self.signal=SignalObject()
        self.shadowcolor=OseImage.DEFAULT_COLOR
        self.identifier=identifier
        self.nct=nct
        self.ncc=ncc
        self.startc=startc
        self.position = self.pos()
        
        self.deltas=[]
        self.connectors={}
        
        self.anhydro=anhydro           
        self.oseImage()   
    
    def setcolor(self,color=None):
        if color==None:
            self.shadowcolor=OseImage.DEFAULT_COLOR
        else:
            self.shadowcolor=color
        
        if self.isSelected():
            self.itemChange(QGraphicsItem.ItemSelectedChange,True)
        
    def oseImage(self):        
        """
        """
        path = OseImage.imagePath(self.nct,self.ncc,self.anhydro)
        self.image = QPixmap(path).scaled(OseImage.WIDTH, OseImage.HEIGHT, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.setPixmap(self.image)
        
        
        furanePosition={1:(OseImage.WIDTH,OseImage.HEIGHT*2/3),2:(OseImage.WIDTH*3/4,OseImage.HEIGHT),3:(OseImage.WIDTH*1/4,OseImage.HEIGHT),4:(0,OseImage.HEIGHT*2/3),5:(0,OseImage.HEIGHT*1/3)}        
        pyranePosition={1:(OseImage.WIDTH,OseImage.HEIGHT*2/3),2:(OseImage.WIDTH*3/4,OseImage.HEIGHT),3:(OseImage.WIDTH*1/4,OseImage.HEIGHT),4:(0,OseImage.HEIGHT*2/3),5:(OseImage.WIDTH*1/4,OseImage.HEIGHT*1/3),6:(OseImage.WIDTH*1/4,0)}                
        if self.ncc==4:
            self.carbonpos=furanePosition
        else:
            self.carbonpos=pyranePosition
        
        for i in range(self.nct):
            self.create_txtdelta(i+1)        
    
        
    
    # number of backbone carbons, number of cycle carbones
    def imagePath(nct,ncc,anh) :        
        """
        selects the image file to diplay according to number of carbons 
        
        :param nct: total number of carbons
        :param ncc: number of carbons in cycle
        :type nct: int
        :type:ncc: int
        """
        img=OseImage.IMG_DIR
        
        if nct>=6 and ncc==5 :
            if nct==6 and anh:
                img+="hexosep_anhydro"
            else:                
                img+="hexosep"
        elif  nct>=6 and ncc==4:
            img+="hexosef"
        elif nct==5 and ncc==5:
            img+="pentosep"
        elif nct==5 and ncc==4 :
            img+="pentosef"
        else :
            img+="hexosep"        
        img+=".png"
        if os.path.exists(img):
            return img
        else:
            return None
    
    def itemChange(self, change, value):
        """
        """
        if change == QGraphicsItem.ItemSelectedChange:
            
            if value :    
                
                #shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=4, yOffset=4)    
                #shadow.setColor(QColor(self.shadowcolor))
                effect=QGraphicsColorizeEffect()
                effect.setColor(QColor(self.shadowcolor))
                self.setGraphicsEffect(effect)                
                
                for d in self.deltas:  
                    f=float(d.toPlainText())
                    if f!=0.0:                           
                        #d.setDefaultTextColor(Qt.black)
                        d.setDefaultTextColor(QColor(self.shadowcolor))
                        #font = d.font()
                        #font.setBold(True)
                        #d.setFont(font)
                        
                        
            else : 
                self.setGraphicsEffect(None)   
                for d in self.deltas:  
                    f=float(d.toPlainText())
                    if f>0.0:                           
                        d.setDefaultTextColor(Qt.darkGreen)      
                    elif f<0: 
                        d.setDefaultTextColor(Qt.red) 
                    else:
                        d.setDefaultTextColor(Qt.transparent) 
                    #font = d.font()
                    #font.setBold(False)
                    #d.setFont(font)
                    
                    
        if change == QGraphicsItem.ItemPositionHasChanged:
            transformX = value.x() - self.position.x()
            transformY = value.y() - self.position.y()    
            self.position = value
            self.notify_lines(transformX,transformY)
           
            
        return QGraphicsItem.itemChange(self, change, value)       
      
    
         
    def deltavalue(self,cnum,delta):        
        """
        """
        if delta < 0 :
            color = Qt.red
        elif delta>0:
            color = Qt.darkGreen
        else:
            color=Qt.transparent
    
        qtxt=self.deltas[cnum-1]    
    
        qtxt.setDefaultTextColor(color) 
    
        x, y = self.get_carbPos(cnum,True)
        middle=self.pos().x()+OseImage.WIDTH/2
       
        if cnum>=self.ncc/2:
            qtxt.setHtml("<p align='%s'>%.2f</font>"%("right",delta))            
        else:            
            qtxt.setHtml("<p align='%s'>%.2f</font>"%("left",delta))            
        
    def mouseDoubleClickEvent(self, event):
        """
        """
        
        self.signal.dble_click(self.identifier)
        

    def mousePressEvent(self,event):
        """        
        """
        
        Logger.debug("image identifier: %i" %self.identifier)
    
    def contextMenuEvent(self,event):
        """
        """
        
        self.signal.right_click(self.identifier)
        
    def create_txtdelta(self, cnum):
        """
        """
        txtdelta=QGraphicsTextItem(self)
        
        txtdelta.setDefaultTextColor(Qt.transparent) 
                
        
        txtdelta.setHtml("<p align='%s'>%.2f</font>"%("left",0))
        
        self.position_delta(txtdelta,cnum)
        self.deltas.append(txtdelta)   
        
    def position_delta(self,txtitem,cnum):
        """
        """
        w=52
        h=22
        padding=0    
        middle=OseImage.WIDTH/2
        txtitem.setTextWidth(w)          
        
        # when cycle does not start at carbon1 (ketoses) => see get_carbPos function to shift carbon connector       
        
        x, y = self.get_carbPos(cnum,True)
        if x<middle:            
            x-=w
        
        
       
        if cnum>=(self.startc+self.ncc):
            y=y-((cnum-self.startc-self.ncc+1)*h)
        elif cnum<self.startc:            
            y=y+((cnum-self.startc)*h)
       
        
        else:
            y-=h/2 

        txtitem.setPos(x,y)        
    
    def get_carbPos(self,cnum,local=False):
        """
        """
       
        if local:
            x=0
            y=0            
        else:
            x=self.position.x()        
            y=self.position.y()                    
       
        first=list(self.carbonpos.keys())[0]
        last=list(self.carbonpos.keys())[-1]
        if cnum>=(self.startc+self.ncc):            
            x +=self.carbonpos.get(last)[0]
            y +=self.carbonpos.get(last)[1]                                                
        elif cnum<=self.startc:            
            x +=self.carbonpos.get(first)[0]
            y +=self.carbonpos.get(first)[1]          
        # position of carbons from 1 to 4 independant of ose type                                                                      
        elif cnum in self.carbonpos.keys():           
            x +=self.carbonpos.get(cnum+first-self.startc)[0]
            y +=self.carbonpos.get(cnum+first-self.startc)[1]           
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("There is no carbon number%s" %cnum)
            msg.setWindowTitle("Warning")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec()             
        
        return x, y    
    
    
    def notify_lines(self,dx,dy):
        """
        """
        for line,conn in self.connectors.items():
            line.notify(dx,dy,conn[1])
            
        
       
    def notify_deltas(self):
        for cnum in range(1,self.nct+1):
            txt=self.deltas[cnum-1]            
            self.position_delta(txt,cnum)    
     
    def connexion_reset(self):
        """
        """
        self.connectors={}
    
    def add_connexion(self,line,connector,orientation):
        """
        """
        Logger.debug("connexion ose %i,cnum %i,sens %i"%(self.identifier,connector,orientation))
        if connector in self.carbonpos.keys():
            self.connectors[line]=(connector,orientation)
    
    def remove_connexion(self,line):
        """
        """
        if line  in self.connectors.keys():
            del self.connectors[line]
    
    def lines(self):
        """
        """
        return list(self.connectors.keys())

    def getconnector(self,line):
        """
        """
        for k,l in self.connectors.items():
            if k==line:
                return line[0]
        return None
    
    
    
    def getconnector_position(self,line):
        """
        """
        connector=self.getconnector(line)
        if connector:
            return self.carbonpos[connector]
        
    def getline(self,cnum):
        """
        """
        for k,l in self.connectors.items():
            if l[0]==cnum:
                return k
        return None
    
    def change_image(self,ncc,nct,anhydro=False):
        """
        """
        self.ncc=ncc
        self.nct=nct
        self.anhydro=anhydro
        self.oseImage()
        self.notify_deltas()
       
        
