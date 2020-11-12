#-*-coding:Utf-8-*-

__author__ ="Virginie Lollier"
__version__ = "1.0.0"
__license__ = "BSD"

import sys,platform
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QMenu, QWidget, QFileDialog
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon


from GlycanQt5 import *
from MassSpectrometry import Spectrum
from Conversion import *
from GlycanTopology import SingletonTopo,OseModel

import utils

class Window( QMainWindow ):
     
    
    """
    main window
    """
    
    def __init__ (self) :  
        """
        create menus ...
        """
        QMainWindow.__init__(self,None)
        self.setWindowTitle('Oligator v.1.0.0')
        self.resize(900, 600)
        
        
        self.scene=QGraphicsScene()
        self.view = Graphics(self.scene)
        self.setCentralWidget(self.view) 
        
        self.maincontroller=TopoControler(self.scene)
        self.view.signal.itemSupprKey.connect(self.deleteSelection)
        
        # toolbar
        toolbar=self.addToolBar("Bar")
        addButton = QAction(QIcon("Image/plus.png"),"add ose",self)
        addButton.triggered.connect(lambda: self.inserOse(SingletonTopo.DEFAULT_OSETYPE))
        toolbar.addAction(addButton)   
        
        delAction = QAction(QIcon("Image/del.png"),"delete selection",self)
        delAction.triggered.connect(self.deleteSelection)
        toolbar.addAction(delAction)        
        
        gridAction= QAction(QIcon("Image/grid.png"),"organize as grid",self)
        gridAction.triggered.connect(self.grid)
        toolbar.addAction(gridAction)     
        
        
        importSmilesBtn = QAction("Import SMILES notation", self)
        importSmilesBtn.triggered.connect( self.importSmiles )  
        
        importWurcsBtn = QAction("Import WURCS notation", self)
        importWurcsBtn.triggered.connect( self.importWurcs )          
        
        
        smilesButton=QAction("Export SMILES notation", self)
        smilesButton.triggered.connect( self.exportSMILES )  
                
        wurcsButton=QAction("Export WURCS notation", self)
        wurcsButton.triggered.connect( self.exportWURCS )              
        
        
        msmsButton=QAction("Export theoretical MS/MS", self)
        msmsButton.triggered.connect( self.exportSpectrum )          
        
        msintensityButton=QAction("Preview spectrum", self)
        msintensityButton.triggered.connect(self.load_SpectroWin)
        
        exitButton = QAction("Exit", self)
        exitButton.triggered.connect( self.close )    
               
        modifButton = QAction("Edit list", self)
        modifButton.triggered.connect( self.updateModif )              
        
        clearButton = QAction("Clear all", self)
        clearButton.triggered.connect( self.clearScene )
        
        edgeButton = QAction("Default binding", self)
        edgeButton.triggered.connect( self.edgeValues )        
        
        oseButton = QAction("Default ose", self)
        oseButton.triggered.connect( self.oseValues )   

        menuBar = self.menuBar()
        
        notation=menuBar.addMenu("Notation")
        design=menuBar.addMenu("Design")
        substituents=menuBar.addMenu("Substituents")
        
        msms=menuBar.addMenu("MS/MS")
        exit=menuBar.addMenu("Exit")
        
        notation.addAction(importSmilesBtn)
        notation.addAction(importWurcsBtn)
        notation.addSeparator()

        notation.addAction(smilesButton)
        notation.addAction(wurcsButton)
       
        design.addAction(oseButton)
        design.addAction(edgeButton)
        design.addAction(clearButton)
        
        substituents.addAction(modifButton)
        
        msms.addAction(msintensityButton)
        msms.addAction(msmsButton)
        msms.addSeparator()        
        
        
        exit.addAction(exitButton)
        self.reset_msms_options()
        
    def reset_msms_options(self):
        self.RatioIon={}
        self.RandomIon={}
        self.dicoH=Spectrum.default
        self.adduct=""        
        
    def grid(self):    
        """
        """
        grrr=SingletonTopo.topogrid()
        if len(SingletonTopo.C1FREE)==1:
            self.view.applyTopoGrid(grrr)
            
        
    #fournit la liste des oses et celle des liaisons à partir de la description d'une structure de glycane au format GlycoCT
    def importGlycoCT(self):        
        """        
        """
        fileName = QFileDialog.getOpenFileName(self, 'Import GlycoCT', "", "All Files (*);;Text Files (*.txt)") #récupère le fichier
        
        res=[]
        lin=[]
        
        if fileName[0]:
            if len(SingletonTopo.OSE_LIST)>0:        
                self.clearScene()
            transfer=0
            #lecture du fichier
            with open(fileName[0], 'r') as f:
                for line in f:
                    line=line.strip()
                    if line=="RES":
                        transfer = 1
                    elif line=="LIN":
                        transfer = 2
                    elif line=="":
                        pass
                    elif transfer == 1:
                        res.append(line)
                    elif transfer == 2:  
                        lin.append(line[line.find(":")+1:])
                
               
                f.close()
        
            if len(res)>0 and len(lin)>0:
                check=GlycoCTDecoder.checkformat(res,lin)
                if check==None:
                    decoder=GlycoCTDecoder(res,lin)
                    
                    if decoder!=None:                    
                        if decoder.check==None:
                            Logger.debug("nb oses from GlycoCT:"+str(len(decoder.gose.nodes())),0) 
                            self.maincontroller.loadmodel(decoder.gatom,decoder.gose)                   
                            self.maincontroller.display_topo()
                else:
                    errdial=QMessageBox()
                    errdial.setIcon(QMessageBox.Critical)
                    errdial.setText("malformed GlycoCT")
                    errdial.setWindowTitle("format error")
                    errdial.setDetailedText(check)                    
                    errdial.exec_()
                
          
            else:
                errdial=QMessageBox()
                errdial.setIcon(QMessageBox.Critical)
                errdial.setText("not a GlycoCT file ")
                errdial.setWindowTitle("format error")
                errdial.exec_()                
        
        
    
    #exporte la structure dessinée au format GlycoCT
    def exportGlycoCT(self):
        """
        """
        
        res=[]
        lin=[]
        fileName = QFileDialog.getSaveFileName(self, 'Export', "", "All Files (*);;Text Files (*.txt)")
        if len(SingletonTopo.OSE_LIST)>0: 
            conversion=None
            
            TopoControler.updateGraph()
            conversion=GlycoCTEncoder(TopoControler.TopoG.Gatom,TopoControler.TopoG.Gose)
            if conversion:            
                if fileName[0]:
                    with open(fileName[0], 'w') as f:
                        f.write("RES"+os.linesep)
                        for r in conversion.reslist:
                            f.write(r+os.linesep)
                        f.write("LIN"+os.linesep)
                        for l in conversion.linlist:
                            f.write(l+os.linesep)
                    f.close()
                    
    def close( self ):
        """
        """
        sys.exit()
    
    #ajoute un ose dans le modèle de donnée et sur la surface de dessin
    def inserOse(self, ose_type):          
        """
        """
        self.maincontroller.createose(self.view.position_oseimage())
        
        
    
    #instancie l'objet SubstituantView      
    def updateModif(self):        
        """
        """        
        self.dialog = SubstituantView(self,SubstitutionLibrary.get_data())        
        self.dialog.setWindowModality(Qt.ApplicationModal)        
        self.dialog.show()
        
    
    def createSub(self,data):
        """
        """
        SubstitutionLibrary.create_substit(data["name"],data["formula"],data["linkage"],data["smiles"])    
    
    def modSub(self,data):
        """
        """
        identifier=data["identifier"]     
        sub=SubstitutionLibrary.getSub(int(data["identifier"]))
        mass1=sub.massSubstitution()
        if int(identifier)>0:            
            SubstitutionLibrary.mod_id(**data)     
            mass2=sub.massSubstitution()
            if mass1!=mass2: 
                Logger.debug("propaget mass change to deltas on ose images")
                self.maincontroller.subchanged(sub)
                
        else:            
            self.createSub(data)
            
    def removeSub(self,identifier):
        """
        """
        SubstitutionLibrary.rm_id(identifier)
    
    def deleteSelection(self):        
        """
        """
       
        
        for elt in self.maincontroller.clearSelection() :
            if isinstance(elt,OseImage):                    
                elt.setSelected(False)
                elt.setVisible(False)
                elt.identifier=None
                elt.setPos(0,0)
            else:
                self.scene.removeItem(elt)
                
             
        
    def clearScene(self):
        """
        """
        SingletonTopo.clear()
        self.maincontroller.lsoseimg={}
        self.scene.clear()        
        self.maincontroller.TopoG=None
        self.reset_msms_options()
        
        
   
    
    #fenêtre de dialogue des paramètres par défaut des liaisons
    def edgeValues(self):
        """
        """
        dialog = edgeType(SingletonTopo.DEFAULT_LINK)
        if dialog.exec_(): 
            SingletonTopo.DEFAULT_LINK[0]=dialog.carb1
            SingletonTopo.DEFAULT_LINK[1]=dialog.carb2
    
    #fenêtre de dialogue des paramètres par défaut des oses    
    def oseValues(self):
        """
        """
        dialog = oseType(SingletonTopo.DEFAULT_OSETYPE)
        if dialog.exec_():
            SingletonTopo.DEFAULT_OSETYPE["nct"]=dialog.nct
            SingletonTopo.DEFAULT_OSETYPE["ncc"]=dialog.ncc
    
       
    
    def importWurcs(self):
        """
        """
        fileName = QFileDialog.getOpenFileName(self, 'Import WURCS', "", "All Files (*);;Text Files (*.txt)") #récupère le fichier
        notation=""
        if fileName[0]:
            if len(SingletonTopo.OSE_LIST)>0:
            #if len(Topology.OSE_LIST)>0:
                self.clearScene()
            #lecture du fichier
            with open(fileName[0], 'r') as f:
                for line in f:
                    notation+=line.strip()
                f.close()
        if notation=="":            
            errdial=QMessageBox()
            errdial.setIcon(QMessageBox.Critical)
            errdial.setText("empty notation")
            errdial.setWindowTitle("format error")
            errdial.setDetailedText(check)                    
            errdial.exec_()                
        else:            
            if re.match("^WURCS=2.0.*$",notation):
                try:
                    decoder=WurcsDecoder(notation)  
                    if decoder==None:
                        Logger.debug("empty graph error ",2)                          
                    else:
                        self.maincontroller.loadmodel(decoder.Gatom,decoder.Gose)
                        self.maincontroller.display_topo()                                          
                except :
                    Logger.debug("empty graph error ",2)                    
            else:
                errdial=QMessageBox()
                errdial.setIcon(QMessageBox.Critical)
                errdial.setText("not a WURCS notation")
                errdial.setWindowTitle("format error")                
                errdial.exec_()                 
            
    def importSmiles(self):
        """
        """
        fileName = QFileDialog.getOpenFileName(self, 'Import SMILES', "", "All Files (*);;Text Files (*.txt)") #récupère le fichier
        notation=""
        if fileName[0]:
            if len(SingletonTopo.OSE_LIST)>0:    
                self.clearScene()
            #lecture du fichier
            with open(fileName[0], 'r') as f:
                for line in f:
                    notation+=line.strip()
                f.close()
        
            if SmilesDecoder.check_smiles(notation):
                try:
                    decoder=SmilesDecoder(notation)                     
                    if decoder==None:                    
                        Logger.debug("empty graph error ",2)    
                    else:
                        Logger.debug("nb oses from smiles:"+str(len(decoder.Gose.nodes())),0) 
                        
                        self.maincontroller.loadmodel(decoder.Gatom,decoder.Gose)
                        self.maincontroller.display_topo()
                except:
                    Logger.debug("empty graph error ",2)  
                    errdial=QMessageBox()
                    errdial.setIcon(QMessageBox.Critical)
                    errdial.setText("error on notation")
                    errdial.setWindowTitle("SMILES format error")                
                    errdial.exec_()                     
            else:
                errdial=QMessageBox()
                errdial.setIcon(QMessageBox.Critical)
                errdial.setText("not a SMILES notation")
                errdial.setWindowTitle("format error")                
                errdial.exec_()               
            
    def exportSMILES(self):
        """
        """
        fileName = QFileDialog.getSaveFileName(self, 'Export SMILES', "", "All Files (*);;Text Files (*.txt)")
        if len(SingletonTopo.OSE_LIST)>0:            
        
            TopoControler.updateGraph()
            conversion=SmilesEncoder(TopoControler.TopoG.Gatom)
                
            text=conversion.graph2smi()
            
            
            if fileName[0]:
                self.__write_lines__([text],fileName[0])
     
    
                    
    def exportWURCS(self):
        """
        """
        
        if len(SingletonTopo.OSE_LIST)>0:   
            TopoControler.updateGraph()
            if TopoControler.is_wurcs_ok():                    
                fileName = QFileDialog.getSaveFileName(self, 'Export WURCS', "", "All Files (*);;Text Files (*.txt)")            
            
                if fileName[0]:                
                
                    conversion=WurcsEncoder(TopoControler.TopoG.Gatom,TopoControler.TopoG.Gose)    
                    text=conversion.notation                
                    self.__write_lines__([text],fileName[0])                    
            else:
                
                errdial=QMessageBox()
                errdial.setIcon(QMessageBox.Critical)
                errdial.setText("WURCS notation with cyclic substituants not available")
                errdial.setWindowTitle("export error")
                errdial.setDetailedText("The structure may contain one or more anhydro bonds or aromatic substituants. \nThis should be exported as SMILES notation instead")                    
                errdial.exec_()                    
                    
    def exportSpectrum(self):   
        """
        """
        
        if len(SingletonTopo.OSE_LIST)>0:  
            fileName = QFileDialog.getSaveFileName(self, 'Export theoretical MS/MS', "", "All Files (*);;Text Files (*.txt)")
        
            print(SubstitutionLibrary.get_non_smiles())
            TopoControler.updateGraph()
            conversion=SmilesEncoder(TopoControler.TopoG.Gatom)                            
            smi=conversion.graph2smi()            
            s=Spectrum(TopoControler.TopoG.Gatom,TopoControler.TopoG.Gose,TopoControler.TopoG.get_ostart())        
            
                
            ions=[]
            if len(self.RandomIon)>0:                
                ions=s.get_peaks(self.RandomIon,True)
            elif len(self.RatioIon)>0:
                ions=s.get_peaks(self.RatioIon,False)
            else:
                ions=s.get_ions()
            
            if len(self.dicoH)>0    :
                ions=Spectrum.__dicomode__(self.dicoH,ions,"H")        
            
            if self.adduct!="":
                ions=Spectrum.__adduct__(self.adduct,ions)
                
            if fileName[0] and len(ions)>0:
                lines=[]
                lines.append("#SMILES:%s"%smi)
                lines.append("#ION:%s"%self.adduct)
                
                for ion in ions:
                    lines.append(ion.replace(s.getseparator(),"\t"))
                self.__write_lines__(lines,fileName[0])
                
       
    def __write_lines__(self,lines,fname):    
        """
        """
        if len(lines)>0:
            if os.path.exists(os.path.dirname(fname)):
                with open(fname,'w',encoding="utf-8") as fd:
                    if len(lines)>1:
                        for i in lines:   
                            if platform.system()=="Windows":
                                fd.write(i+"\n")
                            else:
                                fd.write(i+os.linesep)
                    else:
                        fd.write(lines[0])
                    fd.close()  
                return True
            else:            
                Logger.debug("error path not valid:%s"%os.path.basename(fname),2)
                return False      
        else:
            Logger.debug("nothing to write here",1)
            
        
            
            
    def load_SpectroWin(self):
        """
        open the window of spectrum preview
        """
        random=self.RandomIon.copy()
        ratio=self.RatioIon.copy()    
        default_H=Spectrum.default.copy()
        dialMS=DialogMS(random,ratio)
               
        if len(SingletonTopo.OSE_LIST)>0:          
    
            TopoControler.updateGraph()
            s=Spectrum(TopoControler.TopoG.Gatom,TopoControler.TopoG.Gose,TopoControler.TopoG.get_ostart(),SubstitutionLibrary.get_non_smiles())
            
            dialMS.setMZ(s.ions)
    
                    
        isOk=dialMS.exec_()
        if isOk:
            self.RandomIon=random.copy()
            self.RatioIon=ratio.copy()
            self.dicoH=dialMS.dicoH   
            self.adduct=dialMS.cmbAdduct.currentText()
        
class TopoControler:
    """
    Class for synchronisation model<->view
    """    
    TopoG=Topograph()
    
    def __init__(self,scene):
        self.lsoseimg={}
        self.scene=scene
        SubstitutionLibrary.sort_by_name()
        
    # from display to graph
    def updateGraph():
        """
        """
        TopoControler.TopoG.topo2graph()
    
    def subchanged(self,sub):        
        # get ose en cnum having this substition
        for oid,om in SingletonTopo.OSE_LIST.items():
            carbs=om.get_carbmod(sub.identifier)
            if len(carbs)>0:
                for cnum in carbs:                    
                    delta=sub.get_delta()                    
                    oimg=self.lsoseimg[oid]
                    oimg.deltavalue(cnum,delta)
        
    # from graph to display
    def loadmodel(self,gatom,gose):
        """
        from graph to display
        """
        TopoControler.TopoG.Gatom=gatom.copy()
        TopoControler.TopoG.Gose=gose.copy()
        SingletonTopo.clear()
        
        oses,edges=TopoControler.TopoG.graph2topo() 
        for o in oses:
            SingletonTopo.addOse(o)
        for b in edges:
            SingletonTopo.addEdge(b[0],b[1],b[2],b[3])
        
        
    def dialose(self,oid):
        """
        """
       
        data=[]
        om=SingletonTopo.OSE_LIST[oid]
        codiso={"D":2,"L":1,"":0}
        for i in range(om.nct):
            incycle=False
            if i+1 in range(om.startcycle,om.startcycle+om.ncc):
                incycle=True
            if om.modifs[0][i] in codiso.keys():
                data.append([codiso[om.modifs[0][i]],om.modifs[1][i],incycle])
            else:
                data.append([0,om.modifs[1][i],incycle])
        out=[]
        anhydro=om.anhydro
        dial=OseDetail(data,IupacName._iupac,om.anhydro)
        
        if dial.exec_():            
            out=dial.dataout
            anhydro=dial.chkAnhydro.isChecked()
        
        if len(out)>0 :            
            print("ose has changed")            
            print(out)
            nbcarbs=len(out)
            startcycle=0
            endcycle=nbcarbs
            incycle=False
            
            nctchanged=nbcarbs!=om.nct
            
            #apply to model
            if nctchanged:
                if nbcarbs<om.nct:
                    om.rm_carbs(range(nbcarbs,om.nct))
                    om.nct=nbcarbs
                elif nbcarbs>om.nct:               
                    om.add_carbs(nbcarbs-om.nct)
                    om.nct=nbcarbs            
            
            for icarb in range(nbcarbs):
                
                om.set_isocarb(icarb+1,OseDetail.ISOCOD[out[icarb][0]])
               
                om.set_modcarb(icarb+1,out[icarb][1])
                if out[icarb][2]:
                    if incycle==False:
                        startcycle=icarb+1
                        incycle=True
                else:
                    if incycle:
                        endcycle=icarb
                        incycle=False
                    
                
            formchanged=om.ncc!=(endcycle-startcycle+1)
            starchanged=om.startcycle!=startcycle
            if formchanged:
                om.ncc=endcycle-startcycle+1
                
            if starchanged:
                om.startcycle=startcycle
            
                      
            oimg=self.lsoseimg[oid]
            if starchanged:
                oimg.startc=startcycle
                
            if nctchanged or formchanged:
                oimg.change_image(om.ncc,nbcarbs,anhydro)
            
            
            if om.anhydro!=anhydro:
                oimg.change_image(om.ncc,nbcarbs,anhydro)
                om.set_anhydrobond(anhydro)
            
                
            for icarb in range(nbcarbs):                
                if icarb <oimg.nct:
                    if out[icarb][1]!=SubstitutionLibrary.getSub(SubstitutionLibrary.NOSUBID) and out[icarb][1]>0:
                        delta=SubstitutionLibrary.getSub(out[icarb][1]).get_delta()
                        oimg.deltavalue(icarb+1,delta)
                        if starchanged:
                            oimg.position_delta(oimg.deltas[icarb],icarb+1)
                else:
                    oimg.create_txtdelta(icarb+1)
                    if out[icarb][1]!=SubstitutionLibrary.getSub(SubstitutionLibrary.NOSUBID)and out[icarb][1]>0:
                        delta=SubstitutionLibrary.getSub(out[icarb][1]).get_delta()
                        oimg.deltavalue(icarb+1,delta)                       
    
        
    
        elif om.anhydro!=anhydro:
            oimg=self.lsoseimg[oid]
            oimg.change_image(om.ncc,om.nct,anhydro)
            om.set_anhydrobond(anhydro)            
            
    def dialedge(self,identifier):
        """
        """
        bond=None
        
        for ob in SingletonTopo.OBOND_LIST:
            if ob.identifier==identifier:
                bond=ob
                break
        if bond:
            img1=self.lsoseimg[bond.parent_ose]
            img2=self.lsoseimg[bond.child_ose]
            om1=SingletonTopo.OSE_LIST[bond.parent_ose]
            om2=SingletonTopo.OSE_LIST[bond.child_ose]
            line=img1.getline(bond.parent_carbon)
            print("selection %i-%i binding, %s"%(bond.parent_ose,bond.child_ose,bond.getAttributString()))
            if line!=None:
                free=[bond.parent_carbon]
                for icarb in range(0,om1.nct):            
                    if om1.modifs[1][icarb]==SubstitutionLibrary.NOSUBID and icarb+1!=om1.ncc : 
                        if om1.anhydro and icarb+1 in OseDetail.CANHY:
                            pass
                        else:
                            free.append(icarb+1)
                        
                            
                   
                dial=EdgeDialog(bond.parent_carbon,bond.child_carbon,free)
                changed=False
                if dial.exec_():                   
                    if dial.carbFrom!=bond.parent_carbon:                    
                        img1.remove_connexion(line)
                        img1.add_connexion(line,dial.carbFrom,1)
                        om1.unbind_carb(bond.parent_carbon)
                        
                        bond.parent_carbon=dial.carbFrom                                           
                        om1.bind_carb(dial.carbFrom)
                        
                        changed=True
                        
                    if dial.carbTo!=bond.child_carbon:
                        img2.remove_connexion(line)
                        img2.add_connexion(line,dial.carbTo,2)                         
                        om2.unbind_carb(bond.child_carbon)
                        bond.child_carbon=dial.carbTo
                        om2.bind_carb(dial.carbTo)   
                        
                        changed=True
                    
                    if changed:
                        
                        x1,y1=img1.get_carbPos(dial.carbFrom)
                        x2,y2=img2.get_carbPos(dial.carbTo)
                        line.setLine(x1,y1,x2,y2)
                    
                    
        
    # removes on model and returns corresponding and de-referenced graphic elements                
    def clearSelection(self):
        """
        """
        selection= self.scene.selectedItems()
        size=len(selection)
        lines=[]
        
        for i in range(size):           
            oid=selection[i].identifier
            SingletonTopo.remove_oid(oid)
        
            del self.lsoseimg[oid]            
            for l in selection[i].lines():
                lines.append(l)
                
        #unreference lines in remaining images
        for i,img in self.lsoseimg.items():                
            for l in lines:                
                img.remove_connexion(l)
                
        
        return selection+lines
        
    def binding(self,oid):    
        """
        decides what to do on image right click
        """
        contextualMenu = QMenu()
        selection= self.scene.selectedItems()
        alone=(len(selection)==1 and selection[0]==self) or len(selection)==0
        nothingtodo=False
        
        if len(selection)==1:        
            oimg=selection[0]
           
            if (len(oimg.lines())>0):
                bindButton = contextualMenu.addAction("unbind")
                bindButton.triggered.connect(lambda:self.unbind_ose(oimg))                    
            else:
                nothingtodo=True
           
        elif len(selection)>1:
            x=[]
            for img in selection:
           
                x.append(img.x())
                
            
            # sort from right to left           
            trix=sorted(x)
            trioimg=selection.copy()
            for i in trix[::-1]:
                trioimg[trix.index(i)]=selection[x.index(i)]
            
            
            bindButton = contextualMenu.addAction("bind")
            bindButton.triggered.connect(lambda:self.bind_selection(trioimg))               
                
            bindButton = contextualMenu.addAction("unbind")
            bindButton.triggered.connect(lambda:self.unbind_selection(selection))              
        else:
            nothingtodo=True
            
        if nothingtodo:            
            self.lsoseimg[oid].setSelected(False)
        else:
            contextualMenu.exec_(QCursor.pos())
            Logger.debug("reducing end = %s"%str(SingletonTopo.C1FREE))        
        
    def unbind_selection(self, imgs):
        """
        """
        for i in range(len(imgs)):
            for j in range(i+1,len(imgs)):
                obs=SingletonTopo.get_osidicbonds(imgs[i].identifier,imgs[j].identifier)
                if len(obs)>0:
                    for ob in obs:
                        # unbind in model
                        om_right=SingletonTopo.OSE_LIST[ob.parent_ose]
                        om_left=SingletonTopo.OSE_LIST[ob.child_ose]
                        om_right.unbind_carb(ob.parent_carbon)
                        om_left.unbind_carb(ob.child_carbon)            
                        # unbind in view                        
                        line=list(set(self.lsoseimg[ob.parent_ose].lines()) & set(self.lsoseimg[ob.child_ose].lines()))
                        if len(line)==1:
                            self.scene.removeItem(line[0])
                            self.lsoseimg[ob.parent_ose].remove_connexion(line[0])
                            self.lsoseimg[ob.child_ose].remove_connexion(line[0])
                            
                        else:
                            print("pas normal ça")
                        SingletonTopo.OBOND_LIST.remove(ob)                    
        
    
    def get_nearest_carbon(self,img1,img2):
        """
        """
        listdist = []
        carb1=SingletonTopo.OSE_LIST[img1.identifier].modifs[1]
        carb2=SingletonTopo.OSE_LIST[img2.identifier].modifs[1]
        
        for icarb1 in range(len(carb1)):
            if carb1[icarb1]==SubstitutionLibrary.NOSUBID:
                for icarb2 in range(len(carb2)):
                    if carb2[icarb2]==SubstitutionLibrary.NOSUBID:
                        x1, y1 =  img1.get_carbPos(icarb1+1)
                        x2, y2 =  img2.get_carbPos(icarb2+1)                        
                        distance = sqrt(((x2-x1)**2)+((y2-y1)**2))
                        listdist.append([icarb1+1, icarb2+1, distance])                        
        
        
        if listdist!=[]:
            distance_sorted = sorted(listdist, key=lambda colonnes: colonnes[2])
            cnum1 = distance_sorted[0][0]
            cnum2 = distance_sorted[0][1]
            return cnum1, cnum2  
        else :
            return 0, 0        
    
    def bind_selection(self,selection):        
        """
        """
        for i in range(len(selection)-1):
            j=i+1
            oidi=selection[i].identifier
            oidj=selection[j].identifier
            if len(SingletonTopo.get_osidicbonds(oidi,oidj))==0:
                carbright=SingletonTopo.DEFAULT_LINK[1]
                carbleft=SingletonTopo.DEFAULT_LINK[0]
                if selection[i].x()>selection[j].x():
                    right=selection[i]
                    left=selection[j]
                else:
                    right=selection[j]
                    left=selection[i] 
                    
                # furane
                if right.ncc==carbright:
                    carbright-=1                
                    
                if SingletonTopo.OSE_LIST[right.identifier].modifs[1][carbright-1]==SubstitutionLibrary.NOSUBID:
                    if SingletonTopo.OSE_LIST[left.identifier].modifs[1][carbleft-1]!=SubstitutionLibrary.NOSUBID:                        
                        carb1,carbleft=self.get_nearest_carbon(right,left)
                else:
                    if SingletonTopo.OSE_LIST[left.identifier].modifs[1][carbleft-1]==SubstitutionLibrary.NOSUBID:                        
                        carbright,carb2=self.get_nearest_carbon(right,left)
                    else:
                        carbright,carbleft=self.get_nearest_carbon(right,left)
                        
                if carbleft!=0 and carbright!=0:    
                    ob=SingletonTopo.addEdge(right.identifier,carbright,left.identifier,carbleft)
                    posright=right.get_carbPos(carbright)
                    posleft=left.get_carbPos(carbleft)
                    line=EdgeLine(posright[0],posright[1],posleft[0],posleft[1],ob.identifier)
                    self.scene.addItem(line)
                    right.add_connexion(line,carbright,1)
                    left.add_connexion(line,carbleft,2)
                    line.signal.itemDoubleClicked.connect(lambda:self.dialedge(line.signal.value)) 
                    
                else:
                    print("binding impossible")

    def unbind_ose(self,oimage): 
        """
        """
        lines=[]
        # edge lines to remove 
        for lo in oimage.lines():
            lines.append(lo)
        
        # remove references to image lines in bound images
        for i,img in self.lsoseimg.items():
            for lb in img.lines():
                if lb in lines:
                    img.remove_connexion(lb)
        
        # remove graphic lines from scene
        for lo in lines:
            self.scene.removeItem(lo)
        
        # init table of lines in the selected image 
        oimage.connexion_reset()   
       
        oid=oimage.identifier
        
        rm_ob=SingletonTopo.get_osidicbonds(oid)
        #update model
        for ob in rm_ob:
            om_right=SingletonTopo.OSE_LIST[ob.parent_ose]
            om_left=SingletonTopo.OSE_LIST[ob.child_ose]
            om_right.unbind_carb(ob.parent_carbon)
            om_left.unbind_carb(ob.child_carbon)            
            SingletonTopo.OBOND_LIST.remove(ob)
                
    def createose(self,position):
        """
        """
        # generaliser ce bloc
        omdl=OseModel()
        SingletonTopo.addOse(omdl)
        oimg=OseImage(omdl.oseid,omdl.nct,omdl.ncc ,omdl.startcycle)        
        oimg.signal.itemDoubleClicked.connect(lambda:self.dialose(omdl.oseid)) 
        oimg.signal.itemContextMenu.connect(lambda:self.binding(omdl.oseid))
        self.scene.addItem(oimg) 
        self.lsoseimg[omdl.oseid]=oimg
        # fin bloc        
        
        if len(self.lsoseimg)>0:           
            oimg.setPos(*position)
        
    
    
    def display_topo(self): 
        """
        """
        grid={}
        if len(SingletonTopo.C1FREE)==1:                        
            grid=SingletonTopo.topogrid()
        x=0 
        y=0
        padding=20
        for oid,coords in grid.items():                      
            x=coords[1]*(OseImage.WIDTH*1.5)
            y=coords[0]*(OseImage.HEIGHT*1.5)
            
            omdl=SingletonTopo.OSE_LIST[oid]
            oimg=OseImage(oid,omdl.nct,omdl.ncc ,omdl.startcycle,omdl.anhydro)
            oimg.signal.itemDoubleClicked.connect(lambda:self.dialose(oimg.signal.value)) 
            oimg.signal.itemContextMenu.connect(lambda:self.binding(oimg.signal.value))
            self.scene.addItem(oimg) 
            oimg.setPos(x,y)
            self.lsoseimg[oid]=oimg
            
            for icarb in range(omdl.nct):
                if omdl.modifs[1][icarb]!=SubstitutionLibrary.NOSUBID and omdl.modifs[1][icarb]>0:
                    oimg.deltavalue(icarb+1,SubstitutionLibrary.getSub(omdl.modifs[1][icarb]).get_delta())
         
                    
        for ob in SingletonTopo.OBOND_LIST:
            oimg_right_x,oimg_right_y=self.lsoseimg[ob.parent_ose].get_carbPos(ob.parent_carbon)
            oimg_left_x,oimg_left_y=self.lsoseimg[ob.child_ose].get_carbPos(ob.child_carbon)
            
            
            
            line=EdgeLine(oimg_right_x,oimg_right_y,oimg_left_x,oimg_left_y,ob.identifier)
            self.scene.addItem(line)  
            self.lsoseimg[ob.parent_ose].add_connexion(line,ob.parent_carbon,1)
            self.lsoseimg[ob.child_ose].add_connexion(line,ob.child_carbon,2)
            line.signal.itemDoubleClicked.connect(lambda:self.dialedge(line.signal.value)) 
            
            
        
    def is_wurcs_ok():        
        if TopoControler.TopoG.get_nbcycles()>TopoControler.TopoG.get_nboses():
            return False
        else:
            return True
        
            
            
if __name__ == "__main__" :
    app = QApplication(sys.argv)
    utils.Logger.LEVEL=2
    window = Window()
    window.show()
    
    sys.exit( app.exec_() )