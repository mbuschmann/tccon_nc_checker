from __future__ import print_function, division
import os, sys, time
from PyQt5.QtGui import QIcon
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.ticker as mtick
import numpy as np
import datetime as dt
import xarray as xr

class ftsreader():
    '''Python class to interact with FTS files.\n\n
    Version 2019-08-15\n\n
    Usage:
        ftsreaderobject = ftsreader(path, verbose=False, getspc=False, getifg=False)

    returns an instance of this class, with access to the header and data blocks of the file.

    Example:

        import matplotlib.pyplot as plt

        ftsreaderobject = ftsreader(path, verbose=False, getspc=True)

        ftsreaderobject.print_header()

        if ftsreaderobject.has_block('Data Block SpSm'):
            plt.figure()
            plt.title('Spectrum of '+ftsreaderobject.filename)
            plt.plot(ftsreaderobject.spcwvn, ftsreaderobject.spc, 'k-')
            plt.show()
    '''

    def search_header_par(self, par):
        '''search the header for parameter <par> and return datablock designation '''
        pars = []
        for i in list(self.header.keys()):
            for j in list(self.header[i].keys()):
                if par == j:
                    pars.append(i)
        if len(pars)==1:
            return pars[0]
        elif len(pars)>1:
            if self.verbose: print('Found parameter in multiple datablocks')
            return pars
        else:
            if self.verbose: print('Parameter', par, 'not found in header!')
            return None

    def get_header_par(self, par):
        try:
            return self.header[self.search_header_par(par)][par]
        except:
            print('Parameter not found in header ...')
            return None

    def read_structure(self):
        '''Read the structure of the file and write to ftsreader.fs'''
        # known blocks so far, there is always a block zero, that is still unidentified
        self.__blocknames =    {'160': 'Sample Parameters',
                        '23': 'Data Parameters',
                        '96': 'Optic Parameters',
                        '64': 'FT Parameters',
                        '48': 'Acquisition Parameters',
                        '32': 'Instrument Parameters',
                        '7':  'Data Block',
                        '0':  'something'}
        self.__blocknames2 = {'132': ' ScSm', # another declaration to differentiate blocks between ifg, spc, etc.
                        '4': ' SpSm',
                        '8': ' IgSm',
                        '20': ' TrSm',
                        '12': ' PhSm',
                        b'\x84': ' SpSm/2.Chn.', # some weird stuff going on with python3 decoding here, use binary representation
                        b'\x88': ' IgSm/2.Chn.'}
        self.fs = {}
        with open(self.path, 'rb') as f:
            f.seek(0)
            self.log.append('Reading structure of file')
            # read beginning of file to assert magic number, total number of blocks and first offset
            # some unidentified numbers in between, do not seem to be necessary for header, spc or ifg blocks
            (magic, something, something, offset1, something, numberofblocks) = struct.unpack('6i', f.read(struct.calcsize('6i')))
            f.seek(offset1) # start at first offset
            for i in range(numberofblocks): # go through all blocks and save all found blocks in self.fs
                s = f.read(struct.calcsize('2BH2i'))
                #read beginning of block, with infos on block types, something yet unidentified/unimportant of size 'H' for now, length and gobal offset of the block
                (blocktype, blocktype2, something, length, offset2) = struct.unpack('2BH2i',s)
                blocktype = str(blocktype)
                blocktype2 = str(blocktype2)
                if blocktype in self.__blocknames.keys():
                    hdrblockname = self.__blocknames[blocktype]
                else:
                    hdrblockname = '[unknown block '+blocktype+']'
                if blocktype2 in self.__blocknames2.keys():
                    hdrblockname += self.__blocknames2[blocktype2]
                else: pass
                self.log.append('Found block '+str(blocktype)+', '+str(blocktype2)+' and identified as '+hdrblockname)
                if blocktype == '0' or blocktype not in self.__blocknames.keys():
                    hdrblockname += ' len %3i' % (length)
                else:
                    pass
                #print(hdrblockname, type(hdrblockname))
                self.fs[hdrblockname] = {'blocktype': blocktype, 'blocktype2': blocktype2, 'length': length, 'offset': offset2}

    def getparamsfromblock(self, offset, length, full=False):
        '''Read all parameters in a block at binary <length> and <offset> and return as dictionary. On request also include binary length and offset of that parameter.'''
        params = {}
        i=0
        test = True
        fullblock = []
        with open(self.path, 'rb') as f:
            while test:
                f.seek(offset+i) # goto block offset
                s = f.read(8) # read 8 bytes
                para, thistype, length = struct.unpack('4s2H', s) # unpack to get info on how to unpack block
                if full:
                    fullblocktmp = [para, thistype, length, offset+i]
                i+=8
                if struct.unpack('4c', para)[-1]==b'\x00': #get null terminating string
                    para=para[:-1]
                else: pass
                if para[:3] != b'END' and length>0: # if not empty block
                    f.seek(offset+i)
                    data = f.read(2*length)
                    i+=2*length
                    try:
                        if thistype == 0:
                            val = struct.unpack('%1ii'%(len(data)/4), data)[0]
                        elif thistype == 1:
                            val = struct.unpack('%1id'%(len(data)/8), data)[0]
                        elif thistype >= 2 and thistype <=4:
                            t = struct.unpack('%1is'%(2*length), data)[0].decode('ISO-8859-1')
                            t2 = ''
                            for ji in t: # deal with zeros in byte array
                                if ji!='\x00' and type(ji)==str: # in python2 you might want to add ... or type(ji)=='unicode'):
                                    t2 += ji
                                else:
                                    break
                            val=t2
                        else:
                            val= '[read error]'
                        params[para.decode()] = val
                        if full:
                            fullblocktmp.append(val)
                            fullblock.append(fullblocktmp)
                    except Exception as e:
                        print('Exception in getparamsfromblock')
                        self.log.append(e)
                        print (e)
                else:
                    test = False
        if full:
            return fullblock
        else:
            return params

    def read_header(self):
        '''Read the header and return as a dictionary.'''
        self.log.append('Reading Header ...')
        self.read_structure()
        self.header = {}
        for block in self.fs.keys():
            if block[:10]!='Data Block' and self.fs[block]['length']>0: # if not data block and not empty, try reading header info
                if 'unknown' in block or 'something' in block:
                    pass
                else:
                    try:
                        self.log.append('Reading Header Block: '+block)
                        self.header[block] = self.getparamsfromblock(self.fs[block]['offset'], self.fs[block]['length'], full=False)
                    except Exception as e:
                        print(e)
                        self.log.append(e)
            else: pass
        return 0

    def fwdifg(self):
        if self.header['Instrument Parameters']['GFW']==1:
            return self.ifg[:len(self.ifg)/2]
        else:
            return None

    def bwdifg(self):
        if self.header['Instrument Parameters']['GBW']==1:
            return self.ifg[len(self.ifg)/2:][::-1]
        else:
            return None

    def print_header(self, getlist=False):
        '''Print a nice representation of the header including the names assigned to the header parameters (not complete). Return list of this if requested via <getlist=True>.'''
        headernames = {'Data Parameters': {
            'DPF': 'Data Point Format',
            'FXV': 'Frequency of First Point',
            'LXV': 'Frequency of Last Point',
            'DAT': 'Date of Measurement',
            'TIM': 'Time of Measurement'},
        'Acquisition Parameters': {
            'AQM': 'Acquisition Mode',
            'HFW': 'Wanted High Frequency Limit',
            'LFW': 'Wanted Low Frequency Limit',
            'NSS': 'Sample Scans',
            'RES': 'Resolution'},
        'FT Parameters': {
            'APF': 'Apodization Function',
            'PHR': 'Phase Resolution',
            'ZFF': 'Zero Filling Factor'},
        'Optic Parameters': {
            'APT': 'Aperture Setting',
            'BMS': 'Beamsplitter Setting',
            'CHN': 'Measurement Channel',
            'DTC': 'Detector Setting',
            'HPF': 'High Pass Filter',
            'LPF': 'Low Pass Filter',
            'OPF': 'Optical Filter Setting',
            'PGN': 'Preamplifier Gain',
            'SRC': 'Source Setting',
            'VEL': 'Scanner Velocity'},
        'Sample Parameters': {},
        'Instrument Parameters': {
            'HFL': 'High Folding Limit',
            'LFL': 'Low Folding Limit',
            'LWN': 'Laser Wavenumber',
            'GFW': 'Number of Good FW Scans',
            'GBW': 'Number of Good BW Scans',
            'BFW': 'Number of Bad FW Scans',
            'BBW': 'Number of Bad BW Scans',
            'PKA': 'Peak Amplitude',
            'PKL': 'Peak Location'},
        }
        headerlist = []
        for i in self.header.keys():
            print(i)
            for j in self.header[i].keys():
                if i in headernames.keys() and j in headernames[i].keys():
                    print('  %3s %030s %030s'%(j, headernames[i][j], self.header[i][j]))
                    headerlist.append((i, j, headernames[i][j], self.header[i][j]))
                else:
                    print('  %3s '%(j)+' '*30+'%030s'%(self.header[i][j]))
                    headerlist.append((i, j, ' ', self.header[i][j]))
        if getlist:
            return headerlist
        else: pass

    def get_block(self, pointer, length):
        '''Get data block from ftsreader.path at <pointer> with length <length>.'''
        self.log.append('Getting data block at '+str(pointer)+' with length '+str(length))
        with open(self.path, 'rb') as f:
            f.seek(pointer)
            dat = np.array(struct.unpack('%1if'%(length), f.read(length*4)))
        return dat

    def get_datablocks(self, block):
        '''Read a datablock named <block> and retrieve x- and y-axis np.arrays from it.'''
        self.log.append('Getting data blocks')
        yax = np.array(self.get_block(self.search_block(block)['offset'], self.search_block(block)['length']))
        #print(block)
        if block == 'Data Block IgSm':
            self.log.append('Getting ifg data block')
            # crude estimate of opd axis, only for illustratiion purposes, zpd's not included in calculation, and triangular apod. assumption -> 0.9
            xax = np.linspace(0,2*0.9/float(self.header['Acquisition Parameters']['RES']), len(yax))
        if block == 'Data Block SpSm':
            self.log.append('Getting spc data block')
            # calculate wavenumber axis for spectrum from frequencies of first and last point stored in header
            xax = np.linspace(self.header['Data Parameters SpSm']['FXV'], self.header['Data Parameters SpSm']['LXV'], len(yax))
        if block == 'Data Block ScSm':
            self.log.append('Getting spc data block')
            xax = np.linspace(self.header['Data Parameters ScSm']['FXV'], self.header['Data Parameters ScSm']['LXV'], len(yax))
        if block == 'Data Block TrSm':
            self.log.append('Getting trm data block')
            xax = np.linspace(self.header['Data Parameters TrSm']['FXV'], self.header['Data Parameters TrSm']['LXV'], len(yax))
        if block == 'Data Block PhSm':
            self.log.append('Getting pha data block')
            xax = np.linspace(self.header['Data Parameters PhSm']['FXV'], self.header['Data Parameters PhSm']['LXV'], len(yax))
        return xax, yax

    def get_slices(self, path):
        '''First attempt to implement concatinated slices from automated measurement routines. Probably only works for Uni Bremen setup currently.'''
        self.slices = {}
        self.slices_headers = {}
        slice_list = os.listdir(os.path.join(self.path, 'scan'))
        slice_list.sort()
        good_slice_list = []
        for i in slice_list:
            if i[-5:]!='.info':
                try:
                    self.filename = i
                    self.folder = os.path.join(path, 'scan')
                    self.path = os.path.join(path, 'scan', i)
                    #print('testing file', i)
                    self.test_if_ftsfile()
                    if self.status:
                        #print('read header', self.path)
                        self.read_header()
                        if self.has_block('Data Block IgSm'):
                            opd, ifg = self.get_datablocks('Data Block IgSm')
                            self.read_header()
                            self.slices_headers[i[1:9]+'_header'] = self.header
                            self.slices[i[1:9]] = ifg
                            good_slice_list.append(i)
                        else: pass
                    else: pass
                except: pass
            else: pass
        if len(good_slice_list)>0:
            #print(self.slices.keys())
            self.filename = good_slice_list[0]
            self.folder = os.path.join(self.path, 'scan')
            self.path = os.path.join(path, 'scan', good_slice_list[0])
            self.read_header()
            ifg = np.array([])
            for i in good_slice_list:
                if i[-5:]!='.info' and i[1:9] in self.slices.keys():
                    ifg = np.concatenate([ifg, self.slices[i[1:9]]])
                else: pass
            self.ifg = ifg
            self.opd = np.linspace(0,2*0.9/float(self.header['Acquisition Parameters']['RES']), len(self.ifg))
        else:
            print('Error loading slices from ', path)
            self.status = False
        return 0

    def test_if_ftsfile(self):
        '''Check the initialized filename for FTS magic number.'''
        self.log.append('testing if FTS file')
        # same 4-byte binary representation found on all valid FTS files ... must be magic
        ftsmagicval = b'\n\n\xfe\xfe'
        try:
            with open(self.path, 'rb') as f:
                f.seek(0)
                magic = f.read(4)
            if magic==ftsmagicval:
                if self.verbose:
                    self.log.append('Identified '+self.path+' as FTS file ...')
                self.status=True
                self.isftsfile = True
            else:
                self.log.append('Bad Magic found in '+self.path)
                print('Bad Magic in ', self.path)
                self.status=False
                self.isftsfile = False
        except Exception as e:
            self.log.append(e)
            self.status=False
            self.isftsfile = False

    def search_block(self, blockname):
        '''Searches a <blockname> within the identifies FTS file structure. Returns dictionary entry of the block <blockname>.'''
        #ipdb.set_trace()
        if blockname in list(self.fs.keys()):
            #print(blockname)
            return self.fs[blockname]
        else:
            self.log.append('Could not find '+str(blockname)+' in self.fs.keys()')

    def print_fs(self):
        '''Printing the structure of the FTS file. This includes found data blocks, their binary lengths and offsets.'''
        for i in self.fs.keys():
            print(i, '\n\toffset =', self.fs[i]['offset'], '\n\tlength =', self.fs[i]['length'])

    def print_log(self):
        '''Printing the log of everything that has happened to the class object to std out'''
        for i in self.log:
            print(i)

    def compare_fts_header(self, header2, verbose=True):
        '''Compare this instances header with another <header2>. If <verbose=False> only differences are shown.'''
        S = ' this header           the other header \n'
        for i in self.header.keys():
            if i in header2.keys():
                if verbose:
                    S += '\n'+str(i)+'\n'
                else: pass
                for j in self.header[i].keys():
                    try:
                        a, b = self.header[i][j], header2[i][j]
                        if a==b and verbose:
                            s = j+' '*67+'\n'
                            s = s[:21]+'identical'+s[30:]
                        elif a!=b:
                            s = j+' '*67+'\n'
                            s = s[:8]+str(a)+s[8+len(str(a)):]
                            s = s[:32]+str(b)+s[32+len(str(b)):]
                        else:
                            s = ''
                    except:
                        s = j+' '*67+'\n'
                        s = s[:18]+'problem with key'+s[34:]
                    S += s
            else:
                S += '\n'+str(i)+' missing in other header \n'
        return S

    def has_block(self, blockname):
        '''Check if <blockname> is present in ftsreader.fs'''
        if blockname in self.fs.keys():
            return True
        else:
            return False

    def __init__(self, path, verbose=False, getspc=False, getifg=False, gettrm=False, getpha=False, getslices=False):
        self.log = []
        self.status = True
        self.verbose = verbose
        self.path = path
        if self.verbose:
            print('Initializing ...')
        self.log.append('Initializing')
        try:
            if path.rfind('/')>0:
                self.folder = path[:path.rfind('/')]
                self.filename = path[path.rfind('/')+1:]
            else:
                self.folder = './'
                self.filename = path
            if not getslices:
                self.test_if_ftsfile()
            if self.status:
                if not getslices:
                    self.read_header()
                else: pass
                # get spc if requested
                if getspc and self.has_block('Data Block SpSm'):
                    self.spcwvn, self.spc = self.get_datablocks('Data Block SpSm')
                elif getspc and self.has_block('Data Block ScSm'):
                    self.log.append('Setting self.spc tp ScSm instead of SpSm')
                    self.spcwvn, self.spc = self.get_datablocks('Data Block ScSm')
                else:
                    self.log.append('No Spectrum requested or not found ... skipping.')
                # get transmission spc if requested
                if gettrm and self.has_block('Data Block TrSm'):
                    self.trmwvn, self.trm = self.get_datablocks('Data Block TrSm')
                else:
                    self.log.append('No Transmissionspectrum requested or not found ... skipping.')
                # get ifg if requested
                if getpha and self.has_block('Data Block PhSm'):
                    self.phawvn, self.pha = self.get_datablocks('Data Block PhSm')
                else:
                    self.log.append('No Phasespectrum requested or not found ... skipping.')
                # get ifg if requested
                if getifg and self.has_block('Data Block IgSm'):
                    self.ifgopd, self.ifg = self.get_datablocks('Data Block IgSm')
                else:
                    self.log.append('No Interferogram requested or not found ... skipping.')
                # try getting slices if requested
                if getslices:
                    self.get_slices(path)
                else: self.log.append('No slices requested or not found ... skipping.')
                if self.verbose and self.status:
                    self.log.append('Finished initializing FTS object.\n\n')
                    print('\n\tFinished initializing ftsreader object.')
            else: raise(ValueError('Does not seem to be an FTS file ... skipping'))
            if self.verbose and not self.status:
                self.log.append('An error occured.')
                print('An error occured.')
        except Exception as e:
            self.log.append('Problem with '+str(e))
            print('Error while processing '+path+' ... check self.log or do self.print_log()')

class TcconCheck(QtWidgets.QMainWindow):
    def __init__(self, ncfile):
        super().__init__()
        self.ncfile = ncfile
        self.current_var = 'xluft'
        self.current_var2 = ''
        self.fixedxlims = False
        self.xlims = None
        self.flagged = False
        self.printerrorbars = False
        self.savefname = 'spc_to_remove_'+dt.datetime.now().strftime('%Y%m%d%H%M%S')+'.dat'
        self.spcname = ''
        #
        #self.spcdir = '/procdata/125HR_Bremen/Bremen_Solar/'
        self.spcdir = '/procdata/ggg2020_spectra/125HR_Bremen/Bremen_Solar_lse_corrected/'
        #
        self.vars2 = ['flag', 'year', 'day', 'hour', 'run', 'lat', 'long', 'zobs', 'zmin', 'solzen',
                      'azim', 'osds', 'opd', 'fovi', 'amal', 'graw', 'tins', 'pins', 'tout', 'pout',
                      'hout', 'sia', 'fvsi', 'wspd', 'wdir', 'tmod', 'pmod', 'h2o_dmf_out', 'h2o_dmf_mod']
        #
        self.load_ncfile()
        self.spcdict = {}
        self.find_all_spc()
        #
        #
        self.title = 'TCCON private.nc Checker'
        self.setWindowTitle(self.title)
        #
        self.initUI()
        #
        self._update_canvas()
        #

    def initUI(self):
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        self.gridlayout = QtWidgets.QGridLayout()
        self.gridlayout.setSpacing(10)
        self.setWindowTitle(self.title)
        left = 10
        top = 10
        width = 1500
        height = 900
        self.setGeometry(left,top,width,height)
        self.statusBar().showMessage('Loading')
        #
        mainMenu=self.menuBar()
        fileMenu=mainMenu.addMenu('File')
        #
        #openButton = QtWidgets.QAction(QIcon('open24.png'), 'Open', self)
        #openButton.setShortcut('Ctrl+O')
        #openButton.setStatusTip('Open directory')
        #openButton.triggered.connect(self.getfolder)
        #fileMenu.addAction(openButton)
        #
        exitButton=QtWidgets.QAction(QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)
        #
        #savebutton=QtWidgets.QPushButton('Save list',self)
        #savebutton.setToolTip("Saving a list of all spectra marked with the 'mark spectrum' button as a textfile")
        #self.gridlayout.addWidget(savebutton, 0, 4, 1 ,1, QtCore.Qt.AlignRight)
        #savebutton.clicked.connect(self.savelist)
        #
        #markbutton=QtWidgets.QPushButton('Mark spectrum',self)
        #markbutton.setToolTip("Mark a spectrum to be put in the list")
        #self.gridlayout.addWidget(markbutton, 0, 1, 1 ,1, QtCore.Qt.AlignRight)
        #markbutton.clicked.connect(self.appendlist)
        #
        #nextbutton=QtWidgets.QPushButton('Next',self)
        #nextbutton.setToolTip("Go to the next spectrum in the list")
        #self.gridlayout.addWidget(nextbutton, 0, 0, 1 ,1, QtCore.Qt.AlignRight)
        #nextbutton.clicked.connect(self.nextspc)
        #
        #
        self.cb = QtWidgets.QComboBox()
        for k in self.vars:
            self.cb.addItem(k)
        self.cb.currentIndexChanged.connect(self.selectionchange)
        self.gridlayout.addWidget(self.cb, 0,0)
        #self.cb2 = QtWidgets.QComboBox()
        #for k in self.vars2:
        #     self.cb2.addItem(k)
        #self.cb2.currentIndexChanged.connect(self.selectionchange2)
        #self.gridlayout.addWidget(self.cb2, 0,1)
        #
        #
        #self.textboxw = QtWidgets.QLabel(self) #QLineEdit(self)
        #self.textboxw.setText("")
        #self.gridlayout.addWidget(self.textboxw, 0,2)
        #
        #self.textbox = QtWidgets.QLabel(self) #QLineEdit(self)
        #self.textbox.setText("")
        #self.gridlayout.addWidget(self.textbox, 3,0)
        #
        self.dirlistwidget = QtWidgets.QListWidget()
        #self.dirlistwidget.itemClicked.connect(self.listclick)
        self.gridlayout.addWidget(self.dirlistwidget, 4,0,2,3, QtCore.Qt.AlignRight)
        #
        markbutton=QtWidgets.QPushButton('Delete spectrum',self)
        markbutton.setToolTip("Put a spectrumname on the list to filter out")
        self.gridlayout.addWidget(markbutton, 0, 5, 1 ,1, QtCore.Qt.AlignRight)
        markbutton.clicked.connect(self.remspc)
        #
        #checkBox = QtWidgets.QCheckBox("persistent x-axis?")
        ##        if self.quickplot: checkBox.toggle()
        #checkBox.stateChanged.connect(self.togglexlims)
        #self.gridlayout.addWidget(checkBox, 0, 2, QtCore.Qt.AlignRight)
        showallbutton=QtWidgets.QPushButton('Show all',self)
        self.gridlayout.addWidget(showallbutton, 1, 0, 1 ,1, QtCore.Qt.AlignRight)
        showallbutton.clicked.connect(self.showall)
        #
        previousbutton=QtWidgets.QPushButton('Previous day',self)
        self.gridlayout.addWidget(previousbutton, 1, 1, 1 ,1, QtCore.Qt.AlignRight)
        previousbutton.clicked.connect(self.previousday)
        #
        nextbutton=QtWidgets.QPushButton('Next day',self)
        self.gridlayout.addWidget(nextbutton, 1, 2, 1 ,1, QtCore.Qt.AlignRight)
        nextbutton.clicked.connect(self.nextday)
        #
        self.textboxday = QtWidgets.QLabel(self)
        try:
            self.textboxday.setText(self.currentday.strftime('%Y-%m-%d'))
        except: pass
        self.gridlayout.addWidget(self.textboxday, 1,3)
        #
        #
        self.checkBox3 = QtWidgets.QCheckBox("day by day?")
        self.checkBox3.stateChanged.connect(self.nextday)
        self.gridlayout.addWidget(self.checkBox3, 0, 2, QtCore.Qt.AlignRight)
        #
        self.checkBox = QtWidgets.QCheckBox("exclude flagged?")
        self.checkBox.stateChanged.connect(self.excludeflagged)
        self.gridlayout.addWidget(self.checkBox, 0, 3, QtCore.Qt.AlignRight)
        #
        self.checkBox2 = QtWidgets.QCheckBox("show errorbars?")
        self.checkBox2.stateChanged.connect(self.toggleerrorbars)
        self.gridlayout.addWidget(self.checkBox2, 0, 4, QtCore.Qt.AlignRight)
        #
        ##matplotlib integration from:
        ##https://matplotlib.org/gallery/user_interfaces/embedding_in_qt_sgskip.html#sphx-glr-gallery-user-interfaces-embedding-in-qt-sgskip-py
        self.dynamic_canvas = FigureCanvas(Figure(figsize=(6, 5)))
        self.gridlayout.addWidget(self.dynamic_canvas, 2,0,2,6)
        self.addToolBar(QtCore.Qt.BottomToolBarArea, NavigationToolbar(self.dynamic_canvas, self))
        self._dynamic_ax1 = self.dynamic_canvas.figure.subplots(1)
        self._dynamic_ax1.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x,p:x.strftime('%Y-%m-%d %H:%M:%S')))
        #self._dynamic_ax2 = self._dynamic_ax1.twinx()
        self.dynamic_canvas.figure.autofmt_xdate()
        self.dynamic_canvas.figure.canvas.mpl_connect('pick_event', self.mplonpick)
        #
        #
        self.dynamic_canvas2 = FigureCanvas(Figure(figsize=(6, 5)))
        self.gridlayout.addWidget(self.dynamic_canvas2, 4,4,2,2)
        self.addToolBar(QtCore.Qt.BottomToolBarArea, NavigationToolbar(self.dynamic_canvas2, self))
        #self._dynamic2_ax1, self._dynamic2_ax2 = self.dynamic_canvas2.figure.subplots(2)
        self._dynamic2_ax1 = self.dynamic_canvas2.figure.subplots(1)
        #
        #
        #self.dirlistwidget = QtWidgets.QListWidget()
        #self.make_listwidget()
        #import ipdb; ipdb.set_trace()
        #self.dirlistwidget.itemClicked.connect(self.listclick)
        #self.gridlayout.addWidget(self.dirlistwidget, 1,3,3,2, QtCore.Qt.AlignRight)
        #
        self.gridlayout.setColumnStretch(4, 4)
        self._main.setLayout(self.gridlayout)
        self.show()

    #def listclick(self, item):
    #    self.filename = item.text().replace('  checked', '').replace('  marked', '')
    #    self.i = self.dirlist.index(self.filename)
    #    self.dirlist2[self.i] = self.dirlist[self.i]+'  checked'
    #    item.setText(item.text()+'  checked')

    def remspc(self):
        if self.spcname!='':
            with open(self.savefname, 'a', encoding='utf8') as f:
                f.write(self.spcname+'\n')
            print(self.spcname, 'marked for deletion')

    def find_all_spc(self):
        for root, dirs, files in os.walk(self.spcdir, topdown=False):
            for name in files:
                self.spcdict[name] = os.path.join(root, name)

    def showall(self):
        self.xcon = self.dates>np.min(self.dates)
        self._dynamic_ax1.set_xlim(self.currentday-dt.timedelta(hours=2), self.currentday+dt.timedelta(hours=26))
        self._update_canvas()

    def setday(self):
        self.xcon = (self.dates > self.currentday) & (self.dates<=self.currentday+dt.timedelta(days=1))

    def previousday(self):
        self.currentday = self.currentday-dt.timedelta(days=1)
        self.setday()
        while not np.any(self.xcon):
            self.currentday = self.currentday-dt.timedelta(days=1)
            self.setday()
        self.textboxday.setText(self.currentday.strftime('%Y-%m-%d'))
        self._update_canvas()

    def nextday(self):
        self.currentday = self.currentday+dt.timedelta(days=1)
        self.setday()
        while not np.any(self.xcon):
            self.currentday = self.currentday+dt.timedelta(days=1)
            self.setday()
        self.textboxday.setText(self.currentday.strftime('%Y-%m-%d'))
        self._update_canvas()

    def find_spc_path(self,spcname):
        try:
            return self.spcdict[spcname]
        except:
            return ''

    def mplonpick(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        points = tuple(zip(xdata[ind], ydata[ind]))
        ll = ''
        lll = []
        for l in points[:10]:
            self.spcname = self.data['spectrum'][self.dates == l[0]].values[0]
            ll+=l[0].strftime('%Y-%m-%d %H:%M:%S \t')+self.spcname+'\t'+str(self.data['pout'][self.dates == l[0]].values[0])+' hPa\t'+'\n'
            lll.append(l[0].strftime('%Y-%m-%d %H:%M:%S \t')+self.spcname+'\t'+str(self.data['pout'][self.dates == l[0]].values[0])+' hPa\t')
        print(ll)
        #item = QtWidgets.QListWidgetItem(item)
        self.dirlistwidget.clear()
        for l in lll:
            self.dirlistwidget.addItem(l)
        #self.textbox.setText(ll)
        self._dynamic2_ax1.clear()
        #self._dynamic2_ax2.clear()
        print('Looking for', self.spcname)
        spcpath = self.find_spc_path(self.spcname)
        print('Found', spcpath)
        if spcpath!='':
            self._dynamic2_ax1.set_title(self.spcname)
            o = ftsreader(spcpath, getspc=True, getifg=False)
            try:
                self._dynamic2_ax1.plot(o.spcwvn, o.spc, linewidth=1)
                self._dynamic2_ax1.figure.canvas.draw()
            except AttributeError:
                print('No SPC found')
            #try:
            #    self._dynamic2_ax2.plot(o.ifg, linewidth=1)
            #    self._dynamic2_ax2.figure.canvas.draw()
            #except AttributeError:
            #    print('No IFG found')

    def toggleerrorbars(self):
        if self.checkBox2.isChecked():
            self.printerrorbars = True
        else:
            self.printerrorbars = False
        self._update_canvas()

    def excludeflagged(self):
        if self.checkBox.isChecked():
            self.flagged = True
        else:
            self.flagged = False
        self._update_canvas()

    #def togglexlims(self):
    #    if self.fixedxlims:
    #        self.fixedxlims = False
    #    else:
    #        self.fixedxlims = True
    #        self.xlims = self._dynamic_ax1.get_xlim()
    #    print('NOT IMPLEMENTED YET ... fixedxlims =', self.fixedxlims)

    def selectionchange2(self,i):
        self.current_var2 = self.vars2[i]
        self._update_canvas()

    def selectionchange(self,i):
        self.current_var = self.vars[i]
        self._update_canvas()
        #print("Items in the list are :")
        #for count in range(self.cb.count()):
        #    print(self.cb.itemText(count))
        #    print("Current index",i,"selection changed ",self.cb.currentText())

    def load_ncfile(self):
        print('Reading', self.ncfile)
        self.data = xr.open_dataset(self.ncfile, decode_times=False)
        #c = np.arange(self.data['time'].shape[0]) < 100
        #self.data = self.data[c]
        datetime_from_tccon = lambda d: dt.datetime(1970,1,1)+dt.timedelta(seconds=int(d))
        #ddd =dt.datetime.now()
        #self.dates = np.array([datetime_from_tccon(self.data['time'][i]) for i in range(self.data.dims['time'])])
        #print(dt.datetime.now()-ddd)
        #ddd =dt.datetime.now()
        self.dates = np.array(list(map(lambda d: dt.datetime(1970,1,1)+dt.timedelta(seconds=int(d)), self.data['time'])))
        #print(dt.datetime.now()-ddd)
        print(len(self.dates[self.dates<dt.datetime(2009,1,1)]))
        #
        #
        #np.datetime64(dt.date.today()) - df.index.values.astype('datetime64[D]')
        #
        #
        self.vars = []
        for k in self.data.variables.keys():
            if k.startswith('x') and not '_' in k: #.endswith('_error'):
                self.vars.append(k)
        #self.vars.sort()
        self.vars = self.vars+self.vars2
        self.xcon = self.dates>np.min(self.dates)
        d = np.min(self.dates)
        self.currentday = dt.datetime(d.year, d.month, d.day)
        print(self.data['xluft'].shape[0], 'data points found')


    def _update_canvas(self):
        #print('Plotting ', self.filename)
        self._dynamic_ax1.clear()
        self._dynamic_ax1.set_title(self.ncfile.split('/')[-1])
        self._dynamic_ax1.set_ylabel(self.data[self.current_var].attrs['long_name']+' ['+self.data[self.current_var].attrs['units']+']')
        #self._dynamic_ax1.set_ylim(np.min(self.data[self.current_var]), np.max(self.data[self.current_var]))
        d1 = np.min(self.dates[self.xcon])
        d2 = np.max(self.dates[self.xcon])
        self._dynamic_ax1.set_xlim(dt.datetime(d1.year, d1.month, d1.day)-dt.timedelta(hours=2), dt.datetime(d2.year, d2.month, d2.day)+dt.timedelta(hours=26))
        #if self.fixedxlims:
        #    self._dynamic_ax1.set_xlim(self.xlims)
        if self.flagged:
            cond1 = self.xcon.copy()
            cond1 = cond1 & (self.data['flag']==0)
            if self.printerrorbars and self.current_var+'_error' in self.data.variables.keys():
                self._dynamic_ax1.errorbar(self.dates[cond1], self.data[self.current_var][cond1], yerr=self.data[self.current_var+'_error'][cond1], c='k', marker='None', linestyle='None', ecolor='gray')
            self._dynamic_ax1.plot(self.dates[cond1], self.data[self.current_var][cond1], c='k', marker='.', linestyle='None', picker=True, pickradius=5)
            #if self.current_var2!='':
            #    self._dynamic_ax2.plot(self.dates[cond1], self.data[self.current_var2][cond1], c='b', marker='.', linestyle='None')
            #    self._dynamic_ax2.set_ylim(np.min(self.data[self.current_var2]), np.max(self.data[self.current_var2]))
        else:
            cond1 = self.xcon.copy()
            cond1 = cond1 & (self.data['flag']==0)
            cond2 = self.xcon.copy()
            cond2 = cond2 & (self.data['flag']>0)
            if self.printerrorbars and self.current_var+'_error' in self.data.variables.keys():
                self._dynamic_ax1.errorbar(self.dates[cond1], self.data[self.current_var][cond1], yerr=self.data[self.current_var+'_error'][cond1], c='k', marker='None', linestyle='None', ecolor='gray')
                self._dynamic_ax1.errorbar(self.dates[cond2], self.data[self.current_var][cond2], yerr=self.data[self.current_var+'_error'][cond2], c='r', marker='None', linestyle='None', ecolor='orange')
            self._dynamic_ax1.plot(self.dates[cond1], self.data[self.current_var][cond1], c='k', marker='.', linestyle='None', picker=True, pickradius=5)
            self._dynamic_ax1.plot(self.dates[cond2], self.data[self.current_var][cond2], c='r', marker='.', linestyle='None', picker=True, pickradius=5)
            #if self.current_var2!='':
            #    self._dynamic_ax2.plot(self.dates[cond1], self.data[self.current_var2][cond1], c='b', marker='.', linestyle='None')
            #    self._dynamic_ax2.plot(self.dates[cond2], self.data[self.current_var2][cond2], c='g', marker='.', linestyle='None')
            #    self._dynamic_ax2.set_ylim(np.min(self.data[self.current_var2]), np.max(self.data[self.current_var2]))
        self._dynamic_ax1.figure.autofmt_xdate()
        self._dynamic_ax1.figure.canvas.draw()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    if len(sys.argv)==2:
        ex = TcconCheck(sys.argv[1])
    else:
        sys.exit('Run like this:\n\tpython3 tccon_nc_checker.py path/to/tccon_nc_file.nc')
    sys.exit(app.exec_())
