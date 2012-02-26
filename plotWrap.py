# DO-NOT-DELETE revisionify.begin() 
#
#   Copyright (c) 2007-2009 Lawrence Livermore National Security,
#   LLC. Produced at the Lawrence Livermore National Laboratory (Nathan
#   Barton <barton22@llnl.gov>) CODE-OCEC-08-104.
#   
#   Please also read the file NOTICES.
#   
#   This file is part of the mdef package (version 0.2) and is
#   free software: you can redistribute it and/or modify it under the
#   terms of the GNU Lesser General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#   
#   A copy of the GNU Lesser General Public License may be found in the
#   file NOTICES. If this file is missing, see
#   <http://www.gnu.org/licenses/>.
#
# DO-NOT-DELETE revisionify.end()
'''
***
For now, plotWrap is hardwired for TkAgg
this is not great, but also not a high priority to fix for now;
PlotWinP might be the start of a decent fix

but plotWrap was built with the idea that you could have a plot without a window, and pyplot relies on the new_figure_manager functions in the backends, which always make a window; 
what we need is something to make the figure and the canvas without the figure manager
... but FigureCanvasMac

... how does one get a drawable figure that is not necessarily drawn?!



'''

import matplotlib
matplotlib.use('TkAgg')

# from matplotlib.backend_bases import Event
# class MyResizeEvent(Event):
#     def __init__(self, name, canvas, width, height):
#         Event.__init__(self, name, canvas)
#         self.width  = width
#         self.height = height
#         return

def autoTicks(x, n, format='%0.2e'):
    import numpy as num
    xmin = num.array(x).min()
    xmax = num.array(x).max()
    delta = (xmax-xmin)/n
    xticks = num.arange(xmin, xmax+delta*0.5, delta)
    if format is not None:
        xticksFormatted = num.array(map(lambda x: float(format % (x)), xticks))
        return xticksFormatted
    return xticks

def argToPW(arg):
    local = True
    if hasattr(arg, 'callXY'): # isinstance(arg, plotWrap.PlotWrap):
        retval = arg
        local = False
    elif hasattr(arg, 'getNextAxes'): # isinstance(arg, plotWrap.PlotWin):
        retval = PlotWrap(window=arg)
    elif arg is None:
        retval = PlotWrap()
    else:
        raise RuntimeError, 'do not know what to do with arg of type : '+str(type(arg))
    return retval, local

class PlotWin:
    __provideToolbarDflt = True
    __debug = False
    __maxAutoFigSize = (7,7)
    __softDestroy = False
    def __init__(self,
                 numRows=1, numCols=-1,
                 title='PlotWin window', 
                 figure=None,
                 relfigsize=(3,3),
                 axesList=None,
                 noAutoAxesList=False,
                 dpi=100
                 ):
        '''
        If pass negative numCols, then numRows is the number of plots
        and the layout is done automatically
        '''
        
        self.dpi = dpi
        self.f = None
        self.iaCur = 0
        if axesList is not None:
            self.axesList = axesList
        else:
            self.axesList = []
        if self.__debug:
            print 'len(axesList): %g' % (len(self.axesList))
        #self.autoAxes = False # not necessary
        self.pwList = [ None for iaCur in range(len(self.axesList)) ]
        self.iaCur = len(self.axesList)

        if numCols <= 0:
            if numRows <= 0:
                #self.autoAxes = True
                self.nc = 0
                self.nr = 0
            else:
                self.__autoNRC(numRows)
            self.figsize = (
                min(max(1,self.nc)*relfigsize[0],self.__maxAutoFigSize[0]),
                min(max(1,self.nr)*relfigsize[1],self.__maxAutoFigSize[1])
                )
        else:
            self.nr = numRows
            self.nc = numCols
            self.figsize = (
                max(1,self.nc)*relfigsize[0],
                max(1,self.nr)*relfigsize[1]
                )
        #
        self.title = title
        self.root = None # dead = True
        
        self.provideToolbar = self.__provideToolbarDflt

        self.f = figure
        self.__checkWin()
        
        if self.__debug:
            print 'nr, nc, len(axesList): %g %g %g' % (self.nr, self.nc, len(self.axesList))
        if len(self.axesList) == 0 and not noAutoAxesList:
            ia = 0
            for ir in range(1,self.nr+1):
                for ic in range(1,self.nc+1):
                    ia += 1
                    a = self.f.add_subplot(self.nr, self.nc, ia)
                    a.clear()
                    a.axis('off')
                    self.axesList.append(a)
                    self.pwList.append(None)
        self.__checkAutoNRAxes()

        return
    def __autoNRC(self, n, resizeFig=False):
        import math
        self.nc = int(math.ceil(math.sqrt(float(n))))
        self.nr = int(math.ceil(float(n)/self.nc))
        if resizeFig and hasattr(self.root, 'wm_geometry'):
            import Tkinter as Tk
            self.__checkWin()
            wInch = max(1,self.nc)*3
            hInch = max(1,self.nr)*3
            w = wInch*self.dpi
            h = hInch*self.dpi
            self.f.set_size_inches((wInch,hInch))
            curGeom = self.root.wm_geometry()
            hxw, xp, yp = curGeom.split('+')
            newGeom = str(w)+'x'+str(h)+'+'+xp+'+'+yp
            self.root.wm_geometry(newGeom)
            self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
            self.canvas.show()
            
            # wInch = max(1,self.nc)*3
            # hInch = max(1,self.nr)*3
            # # self.f.set_size_inches((wInch,hInch))
            # #self.canvas.resize(wInch*self.dpi, hInch*self.dpi)
            # #self.canvas.show() 
            # myResizeEvent = MyResizeEvent('autoNRC_resize', self.canvas, 
            #                               wInch*self.dpi, hInch*self.dpi)
            # self.canvas.resize(myResizeEvent)
        return
    def __checkAutoNRAxes(self):
        if self.__debug:
            print 'nr, nc, len(axesList): %g %g %g' % (self.nr, self.nc, len(self.axesList))
        if len(self.axesList) == 0: return
        assert len(self.axesList) <= self.nr*self.nc, \
            'axesList is wrong length'
        for ia in range(len(self.axesList)):
            axes = self.axesList[ia]
            if hasattr(axes,'change_geometry'):
                axes.change_geometry(self.nr, self.nc, ia+1)
        return
    def __checkWin(self):
        '''
        This is pretty ugly in how it is hardwired for tk;
        consider going to something like PlotWinP
        '''
        import sys
        
        if self.root is not None: return  # if not self.dead: return
        
        useTkAgg = True
        if 'matplotlib.backends' in sys.modules:
            if matplotlib.get_backend() != 'TkAgg':
                useTkAgg = False
        if not useTkAgg:
            #self.root = matplotlib.get_backend()
            # assert len(self.axesList) == 1, 'plotWrap case note coded, axesList len : %d' % (len(self.axesList))
            # self.root = PlotWinP(axes=self.axesList, figsize=self.figsize, dpi=self.dpi, 
            #                     title=self.title)
            assert self.f is None, 'self.f is not None'
            self.root = PlotWinP(axes=self.axesList, figsize=self.figsize, dpi=self.dpi)
            self.canvas = self.root.getCanvas()
            self.f = self.root.f
        else:
            # matplotlib.use('TkAgg') # moved this back to above
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
            
            import Tkinter as Tk
    
            self.root = root = Tk.Tk() # pops up a little window
            root.wm_title(self.title)
            root.wm_resizable(width=True, height=True)
            #self.dead = False
            
            if self.__softDestroy:
                def destroy(e): 
                    self.root = None # self.dead = True
                root.bind("<Destroy>", destroy)
            
            # a tk.DrawingArea
            if self.f is None: # if figure is None:
                from matplotlib.figure import Figure
                self.f = Figure(figsize=self.figsize,
                                dpi=self.dpi)

            self.canvas = canvas = FigureCanvasTkAgg(self.f, master=root)
            #if self.showByDefault: 
            canvas.show()
            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1) # actually displays something now!
            
            #if self.showByDefault: canvas.show() # updates display
            canvas.show() # updates display
            
            if self.provideToolbar:
                toolbar = NavigationToolbar2TkAgg( self.canvas, root )
                toolbar.update()
                canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        return
    def fillPW(self):
        for iPW in range(len(self.pwList)):
            if self.pwList[iPW] is None:
                self.pwList[iPW] = PlotWrap(window=self)
        return
    def haveXLabels(self):
        self.f.subplots_adjust(bottom=0.2,hspace=0.3) 
        return
    def haveYLabels(self):
        self.f.subplots_adjust(left=0.2,wspace=0.3) 
        return
    def destroy(self):
        if self.root is not None: # not self.dead:
            if not self.__softDestroy:
                self.axesList = []
                # for pw in self.pwList:
                #     if pw is not None:
                #         pw.destroy()
                self.pwList = []
            if hasattr(self.root,'destroy'):
                self.root.destroy()
    def save(self, **keyArgs):
        if keyArgs.has_key('filename'):
            filename = keyArgs.pop('filename')
            #self.canvas.print_figure(filename, **keyArgs)
            self.f.savefig(filename, **keyArgs)
        else:
            raise RuntimeError, 'need filename entry in keyArgs'
        return
    def isDead(self):
        return self.root is None # self.dead
    def getCanvas(self):
        return self.canvas
    def show(self):
        self.__checkWin()
        if hasattr(self.f, 'show'):
            self.f.show()
        elif hasattr(self.canvas, 'show'):
            self.canvas.show()
        else:
            raise RuntimeError, 'do not know how to do show'
    def getFigure(self):
        return self.f
    def getNextAxes(self, rect=None, attachPW=None, **axprops):
        aCur = self.getAxes(self.iaCur, rect=rect, **axprops)
        self.pwList[self.iaCur] = attachPW
        self.iaCur += 1
        return aCur
    def getAxes(self, plotNum, rect=None, withPW=False, **axprops):
        '''careful: plotNum is from 0, not 1 as is the case for subplot
        axprops is useful for this like setting sharex and sharey'''
        #if plotNum+1 > self.nr*self.nc:
        if plotNum+1 > len(self.axesList) or len(axprops) > 0:
            #if not self.autoAxes:
            #    raise RuntimeError, 'not more axes available'
            if rect is not None:
                aCur = self.f.add_axes(rect, **axprops)
                self.axesList.append(aCur)
                self.pwList.append(None)
            else:
                if plotNum < len(self.axesList):
                    'replacing existing' 
                    aCur = self.f.add_subplot(self.nr, self.nc, plotNum+1, **axprops)
                else:
                    'just add more axes, instead of fussing!' 
                    self.__autoNRC(plotNum+1, resizeFig=True)
                    aCur = self.f.add_subplot(self.nr, self.nc, self.iaCur+1, **axprops)
                    self.axesList.append(aCur)
                    self.pwList.append(None)
                self.__checkAutoNRAxes()
            if self.__debug:
                print 'length of axesList is now %d' % len(self.axesList)
            retval = aCur
        else:
            aCur = self.axesList[plotNum]
            pw = self.pwList[plotNum]
            if withPW:
                retval = aCur, pw
            else:
                retval = aCur
        return retval

class PlotWinLite:
    '''
    Lightweight PlotWin substitute for when windows are being controlled by code outside of plotWrap
    '''
    def __init__(self, canvas, figure, axes):
        self.c = canvas
        self.f = figure
        self.a = axes
        return
    def getAxes(self, plotNum):
        assert plotNum == 0, \
               'plotNum must be 0 for PlotWinLite instances'
        return self.axes
    def haveXLabels(self):
        'may want to turn off any functionality in this method'
        self.f.subplots_adjust(bottom=0.2,hspace=0.3) 
        return
    def haveYLabels(self):
        'may want to turn off any functionality in this method'
        self.f.subplots_adjust(left=0.2,wspace=0.3) 
        return
    def destroy(self):
        raise RuntimeError, 'should not call destroy on PlotWinLite instances'
        return
    
class PlotWinP:
    '''
    Just wrap pyplot
    '''
    def __init__(self, axes=None, as3D=False, 
                 title=None, 
                 **kwargs):
        import matplotlib.pyplot as plt
        
        self.iaCur = 0
        
        figsize = None
        dpi = None
        if kwargs.has_key('figsize'):
            figsize = kwargs.pop('figsize')
        if kwargs.has_key('dpi'):
            dpi = kwargs.pop('dpi')
        self.f = plt.figure(figsize=figsize, dpi=dpi)
        self.c = self.f.canvas
        
        if axes is None:
            if as3D:
                import mpl_toolkits.mplot3d.axes3d as p3
                self.a = p3.Axes3D(self.f, **kwargs)
            else:
                self.a = [self.f.add_subplot(1,1,1, **kwargs)]
                self.iaCur += 1
        else:
            if hasattr(axes,'__len__'):
                self.a = axes
            else:
                self.a = [axes]
        self.pwList = [None for a in self.a]
        
        'for now, punt on setting a title'
        # if title is not None:
        #     'not sure how to set title on window, so settle for axes title'
        #     self.a.set_title(title)
        
        return
    def getNextAxes(self, rect=None, attachPW=None, **axprops):
        if self.iaCur < len(self.a):
            aCur = self.a[self.iaCur]
            self.pwList[self.iaCur] = attachPW
        else:
            aCur = self.f.add_axes(rect, **axprops)
            self.pwList.append(attachPW)
        self.iaCur += 1
        return aCur
    def getCanvas(self):
        return self.c
    def getFigure(self):
        return self.f
    def getAxes(self, plotNum):
        assert plotNum < len(self.a), \
               'plotNum must be 0 for PlotWinLite instances'
        return self.a[plotNum]
    def haveXLabels(self):
        'may want to turn off any functionality in this method'
        self.f.subplots_adjust(bottom=0.2,hspace=0.3) 
        return
    def haveYLabels(self):
        'may want to turn off any functionality in this method'
        self.f.subplots_adjust(left=0.2,wspace=0.3) 
        return
    def destroy(self):
        import matplotlib.pyplot as plt
        plt.close(self.f)
        self.a = None
        self.pwList = None
        self.c = None
        self.f = None
        return
    
    
        
class PlotWrap(object):
    __debug = False
    __keyArgDict =  {
        'title' : None,
        'winTitle' : None,
        'xlabel' : None,
        'xticks' : None,
        'ylabel' : None,
        'yticks' : None,
        'ylog' : False, 
        'xlog' : False, 
        'ybound' : (None,None),
        'xbound' : (None,None),
        'legend' : None,
        'style' : None,
        'aging' : None,
        'alphaMin' : 0.02,
        'agingCull' : True,
        'bball' : None,
        'accum' : False,
        'showByDefault' : True,
        'makeWinByDefault': True,
        'axes' :  None,
        'axesRect' :  None,
        'axprops' :  {},
        'figsize' : (5,4),
        'dpi' : 100,
        'window' :  None,
        'as3D' : False,
        }
    def __init__(self,
                 **keyArgs
                 ):
        
        # defaults
        for parm, val in self.__keyArgDict.iteritems():
            self.__setattr__(parm, val)
        #
        # parse keyword args
        for keyArg in self.getKeyArgList():
            if keyArgs.has_key(keyArg):
                self.__setattr__(keyArg, keyArgs.pop(keyArg))
        if len(keyArgs) > 0:
            raise RuntimeError, 'unparsed keyword args : '+str(keyArgs)
        #
        window = self.window
        del self.window
        #
        self.x_accum = self.accum
        del self.accum
        #
        axes = self.axes
        del self.axes
        
        self.x_prev = None
        self.x_store = [[],[]]

        self.win = None
        self.a   = None
        self.canvas = None
        
        self.ownCanvas = True
        if window is None:
            self.ownCanvas = True
            'checking self.showByDefault here causes trouble because it is hard to attach a figure to a window later for a general backend ***'
            if self.showByDefault or self.makeWinByDefault:
                # 'go ahead and make a window, using PlotWinP to increase backend flexibility'
                # self.win    = ...PlotWinP(as3D=self.as3D, figsize=self.figsize, dpi=self.dpi)
                # self.a      = self.win.getAxes(0)
                # self.f      = self.win.getFigure()
                # self.canvas = self.win.getCanvas()
                # self.a.set_autoscale_on(True)
                self.__checkWin()
                self.f = self.win.f
            else:
                from matplotlib.figure import Figure
                self.f = Figure(figsize=self.figsize, dpi=self.dpi)
            if self.as3D:
                import mpl_toolkits.mplot3d.axes3d as p3
                self.a = p3.Axes3D(self.f)
            else:
                if self.axesRect is None:
                    self.a = self.f.add_subplot(1, 1, 1, **self.axprops)
                else:
                    self.a = self.f.add_axes(self.axesRect, **self.axprops)
                
            if axes is not None:
                raise RuntimeError, 'do not specify axes when have not passed a window'
        else:
            self.ownCanvas = False
            self.win = window
            self.canvas = self.win.getCanvas()
            self.f = self.win.getFigure()
            if axes is None:
                if self.__debug: print 'using axprops: '+str(self.axprops)
                self.a = self.win.getNextAxes(rect=self.axesRect, attachPW=self, **self.axprops)
            elif isinstance(axes, int):
                self.a = self.win.getAxes(axes, **self.axprops)
            else:
                self.a = axes
        #
        self.clear() # sets labels and so forth

        # self.style = style # can be None
        # self.bball = bball # can be None
        #
        #self.aging = aging
        #self.agingCull = agingCull
        #self.alphaMin = alphaMin
        self.agingNumAge = -1
        if self.aging is not None:
            from math import log
            self.agingNumAge = max(int(log(self.alphaMin)/log(self.aging)),1)
        self.plotProps = {}
        # self.lineListList = [] # not needed, can, and should, use a.lines
        
        self.asImIJ = False
        # self.transform = None
        self.axIm = None
        self.__colorbar = None

        return
    def add_axes(self, arect, kwAddAxes={}, **kwarg):
        aNew = self.f.add_axes(arect, **kwAddAxes)
        # self.a.set_position((0.05, 0.05, 0.45, 0.9))
        new = self.__class__(window=self.win, axes=aNew, **kwarg)
        return new
    @classmethod
    def popKeyArgs(cls, kwargs):
        retval = {}
        for key in cls.getKeyArgList():
            if kwargs.has_key(key):
                retval[key] = kwargs.pop(key)
        return retval
    @classmethod
    def getKeyArgList(cls):
        return cls.__keyArgDict.keys()
    def __checkCanvas(self):
        'leave dead True, not checking for window, just for canvas'
        #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        if self.win is not None: return
        if self.canvas is not None: return

        from matplotlib.backends.backend_tkagg import FigureCanvasAgg
        
        self.canvas = FigureCanvasAgg(self.f)
        if hasattr(self.a,'mouse_init'):
            self.a.mouse_init()
        return

    def __checkWin(self):
        if self.win is not None: return
        winTitle = self.winTitle
        if winTitle is None : winTitle = self.title
        if winTitle is None : winTitle = 'plotWin'
        if hasattr(self,'f'):
            self.win = PlotWin(1,1,title=winTitle, relfigsize=self.figsize, figure=self.f, axesList=[self.a])
        else:
            self.win = PlotWin(1,1,title=winTitle, relfigsize=self.figsize, axesList=[self.a])
        #self.win = PlotWinP(title=winTitle, figure=self.f, axes=self.a)
        self.canvas = self.win.getCanvas()
        return
    def removeLines(self):
        lines = self.a.get_lines()
        for line in lines:
            line.remove()
        if self.showByDefault:
            self.show() # self.canvas.show()
        return
    def show(self):
        self.__checkWin()
        if self.legend is not None:
            pass # do not yet have this right
            # if isinstance(self.legend, str) or isinstance(self.legend, int): #  or hasattr(self.legend,'__len__')
            #     self.a.legend(loc=self.legend)
            # elif isinstance(self.legend, bool):
            #     self.a.legend()
            # elif isinstance(self.legend, list):
            #     self.a.legend(*self.legend)
            # elif isinstance(self.legend, dict):
            #     self.a.legend(**self.legend)
            # else:
            #     raise RuntimeError, 'do not know what to do with legend specification: '+str(self.legend)
        if hasattr(self.a,'mouse_init'):
            self.a.mouse_init()
        if hasattr(self.f, 'show'):
            self.f.show()
        elif hasattr(self.canvas, 'show'):
            self.canvas.show()
        else:
            raise RuntimeError, 'do not know how to do show'
        return
    def clear(self):
        self.asImIJ = False
        self.axIm = None
        # self.transform = None
        self.x_prev = None
        self.__colorbar = None
        self.x_store = [[],[]]
        if self.showByDefault:
            self.__checkWin()
        if self.a is not None:
            if hasattr(self.a,'mouse_init'):
                'doing a.clear() seems like it causes trouble with plotting collections?' 
                self.a.clear() # self.a.cla()
                self.a.set_autoscale_on(True)
            else:
                self.a.clear() # self.a.cla()
                self.a.set_autoscale_on(True)
        if self.xticks is not None:
            self.a.set_xticks(self.xticks)
        if self.title is not None:
            self.a.set_title(self.title)
        if self.xlabel is not None:
            self.a.set_xlabel(self.xlabel)  # r'X axis label $\mu\mbox{lentz}$'
            if self.win is not None:
                'make some room for the label'
                self.win.haveXLabels()
        if self.yticks is not None:
            self.a.set_yticks(self.yticks)
        if self.xlog: self.a.set_xscale('log')
        if self.ylog: self.a.set_yscale('log')
        #self.a.set_ybound(self.ybound)
        self.a.set_ylim(ymin=self.ybound[0], ymax=self.ybound[1])
        self.a.set_xlim(xmin=self.xbound[0], xmax=self.xbound[1])
        if self.ylabel is not None:
            self.a.set_ylabel(self.ylabel)
            if self.win is not None:
                'make some room for the label'
                self.win.haveYLabels()
        if self.legend is not None:
            pass # do not yet have this right
            # if isinstance(self.legend, str) or isinstance(self.legend, int): #  or hasattr(self.legend,'__len__')
            #     self.a.legend(loc=self.legend)
            # elif isinstance(self.legend, bool):
            #     self.a.legend()
            # elif isinstance(self.legend, list):
            #     self.a.legend(*self.legend)
            # elif isinstance(self.legend, dict):
            #     self.a.legend(**self.legend)
            # else:
            #     raise RuntimeError, 'do not know what to do with legend specification: '+str(self.legend)
        if self.a is not None:
            self.a.set_autoscale_on(True)
        if self.showByDefault:
            self.show() # self.canvas.show()
        return
    def setLegend(self, val=True):
        # if self.legend is None:
        #     self.legend = True
        # if val is not None:
        #     self.legend = val
        self.legend = val
        return
    def save(self, **keyArgs):
        '''make hardcopy of the current figure;
        specify filename or prefix keyword argument
        '''
        import numpy as num
        
        filename = None
        prefix = None
        if keyArgs.has_key('filename'):
            filename = keyArgs.pop('filename')
        if keyArgs.has_key('prefix'):
            prefix = keyArgs.pop('prefix')
            filename = prefix+'.pdf' # .eps
        
        if prefix is not None:
            'export data'
            if len(self.x_store[0]) > 0:
                dataFilename = prefix+'.data'
                import arrayUtil
                try:
                    arrayUtil.writeArray(dataFilename, num.array(self.x_store))
                except:
                    import sys
                    print 'problem writing to '+dataFilename+' : '+str(self.x_store)
                    sys.exit(1)
        
        if not self.ownCanvas:
            if self.__debug:
                print 'skipping print_figure because this plot does not own its canvas'
        else:
            self.__checkCanvas()
            if filename is None:
                raise RuntimeError, 'need filename or prefix entry in keyArgs'
            #self.canvas.print_figure(filename, **keyArgs)
            self.f.savefig(filename, **keyArgs)
        
        return
    def destroy(self):
        'does not clean up self.a, just kills the window if this plot owns the window'
        if self.ownCanvas and self.win is not None:
            self.win.destroy()
        return
    def drawBBox(self, bbox, **kwargs):
        import numpy as num
        bbox_x = bbox[0]
        bbox_y = bbox[1]
        xBBox = num.array([ bbox_x[0], bbox_x[1], bbox_x[1], bbox_x[0], bbox_x[0] ])
        yBBox = num.array([ bbox_y[0], bbox_y[0], bbox_y[1], bbox_y[1], bbox_y[0] ])
        'go through call, instead of callXY, in case of asImIJ'
        self.__call__(xBBox, yBBox, **kwargs)
        return
    def __call__(self, *args, **keyArgs):
        import numpy as num
        noShow = False
        retVal = None
        # if keyArgs.has_key('transform'):
        #     'pop transform so that can keep it for other (overlaid) plots'
        #     self.transform = keyArgs.pop('transform')
        if keyArgs.has_key('noShow'):
            noShow = keyArgs.pop('noShow')
        if len(args) == 2:
            alphas = None
            if self.asImIJ:
                x = args[1]
                y = args[0]
            else:
                x = args[0]
                y = args[1]
            if keyArgs.has_key('alphas'):
                alphas = keyArgs.pop('alphas')
                assert len(args) == 2, 'len(args) != 2'
                assert len(alphas) == len(y), 'len(alphas) != len(y)'
                self.x_prev = None
                for iX in range(len(x)):
                    self.callXY([x[iX]],[y[iX]],alpha=alphas[iX],**keyArgs)
            else:
                self.callXY(x, y, **keyArgs)
        elif len(args) == 3:
            X = args[0]; Y = args[1]; data = args[2];
            cont = self.callContour(X, Y, data, **keyArgs)
            retVal = cont
        elif len(args) == 1 and isinstance(args[0],str):
            'interpret string as name of a file to plot in axes'
            filename = args[0]
            im = self.callImage(filename, **keyArgs)
            retVal = im
        elif len(args) == 1 and isinstance(args[0],num.ndarray):
            im = self.callIm(args[0], **keyArgs)
            retVal = im
        else:
            raise RuntimeError, 'do not know what to do with args'
        self.a.set_ylim(ymin=self.ybound[0], ymax=self.ybound[1])
        self.a.set_xlim(xmin=self.xbound[0], xmax=self.xbound[1])
        if not noShow:
            if self.showByDefault:
                self.show()
        return retVal
    def callImage(self, filename, 
                  **keyArgs):
        import Image
        im = Image.open(filename)
        s = im.tostring() # convert PIL image -> string
        import numpy as num
        rgb = num.fromstring(s, dtype=num.uint8).astype(num.float)/255.0 # convert string -> array of floats
        rgb = num.resize(rgb, (im.size[1], im.size[0], 3)) # resize to RGB array
        retval = self.callIm(rgb, **keyArgs)
        return retval
    def callIm(self, im, 
               interpolation='nearest', aspect='equal', 
               ijAsXY=False,
               clear=True, **keyArgs):
        if clear: 
            self.clear()
        self.a.axis('off')
        self.a.set_autoscale_on(True)
        if ijAsXY:
            self.asImIJ = False
            'imshow does not yet really support transform'
            if len(im.shape) == 3:
                imT = im.transpose(1,0,2)
            else:
                imT = im.T
            axIm = self.a.imshow(imT, 
                                 interpolation=interpolation, aspect=aspect, origin='lower', 
                                 # transform=self.transform,
                                 **keyArgs)
            self.a.format_coord = lambda x,y: 'i=%d; j=%d; val=%s' % \
                (round(x), round(y), str(im[round(x),round(y)]))
        else:
            self.asImIJ = True
            'imshow does not yet really support transform'
            axIm = self.a.imshow(im, 
                                 interpolation=interpolation, aspect=aspect, 
                                 # transform=self.transform,
                                 **keyArgs)
            self.a.format_coord = lambda x,y: 'i=%d; j=%d; val=%s' % \
                (round(y), round(x), str(im[round(y),round(x)]))
        'turn off autoscale so that axis limits do not get reset on replots'
        self.axIm = axIm
        self.mappable = axIm
        self.a.set_autoscale_on(False)
        return axIm
    def callContour(self, X, Y, data, 
                    interpolation=None, aspect=None,
                    **keyArgs):
        pp = {}
        pp.update(self.plotProps)
        pp.update(keyArgs) # plotProps
        #cont = self.a.contourf(X, Y, data, 200, **keyArgs)
        'imshow does not yet really support transform'
        self.a.set_autoscale_on(True)
        cont = self.a.imshow(data, origin='lower',
                             extent=(X[0,0],X[0,-1],Y[0,0],Y[-1,0]), 
                             interpolation=interpolation,
                             aspect=aspect,
                             # transform=self.transform,
                             **keyArgs)
        self.a.set_autoscale_on(False)
        self.mappable = cont
        return cont
    def discontXY(self):
        self.x_prev = None
        return
    def callXY(self, x, y, style=None, **keyArgs):
        assert len(x) == len(y), \
               'x and y must be same length'

        xUse = x
        yUse = y
        if len(x) == 1:
            if self.x_prev is not None:
                xUse = []; xUse.append(self.x_prev[0]); xUse.append(x)
                yUse = []; yUse.append(self.x_prev[1]); yUse.append(y)

        pp = {}
        pp.update(self.plotProps)
        pp.update(keyArgs) # plotProps
        if self.bball is not None:
            'first, get rid of last bball'
            lenCur = len(self.a.lines)
            if lenCur>0:
                self.a.lines = self.a.lines[0:lenCur-1]
        if self.aging is not None:
            #for lineList in self.lineListList[:-1]:
            #    for line in lineList:
            lenCur = len(self.a.lines)
            for line in self.a.lines[max(0,lenCur-self.agingNumAge):lenCur]:
                alphaCur = line.get_alpha()
                line.set_alpha(alphaCur*self.aging)
            if self.agingCull:
                if lenCur > self.agingNumAge:
                    self.a.lines = self.a.lines[lenCur-self.agingNumAge:lenCur]
        if style is not None:
            'passing transform of None can cause trouble'
            lines = self.a.plot(xUse, yUse, style,
                                # transform=self.transform,
                                **pp)
        elif self.style is not None:
            #self.lineListList.append(self.a.plot(x,y,self.style,**pp))
            lines = self.a.plot(xUse, yUse, self.style,
                                # transform=self.transform,
                                **pp)
        else:
            #self.lineListList.append(self.a.plot(x,y,**pp))
            lines = self.a.plot(xUse, yUse, 
                                # transform=self.transform,
                                **pp)
        if self.bball is not None:
            'add bball at end'
            self.a.plot([x[-1]],[y[-1]], self.bball)
            # transform=self.transform
        if len(x) == 1:
            if self.x_accum:
                self.x_prev = [x, y]
            self.x_store[0] = self.x_store[0] + list(x)
            self.x_store[1] = self.x_store[1] + list(y) 
        else:
            self.x_prev = None
            if len(x) == len(self.x_store[0]):
                'assume called with the same x values'
                self.x_store.append(list(y))
            else:
                self.x_store[0] = list(x)
                self.x_store[1] = list(y) 
        return
    def setVMM(self, vMM):
        if type(vMM) == int:
            iCol = vMM
            vMM = self.a.collections[iCol].get_clim()
        for col in self.a.collections:
            col.set_clim(vmin=vMM[0],vmax=vMM[1])
        return
    def colorbar(self, rect=(0.80,0.1,0.05,0.8), adjustPos=True, thing=None, **kwargs):
        '''
        if set rect to None, then colorbar steals self.a
        '''
        self.__checkWin()
        w = self.win
        f = self.f
        if len(self.a.collections) > 0:
            self.setVMM(0)
        if thing is None:
            if len(self.a.collections) > 0:
                thing = self.a.collections[0]
            elif hasattr(self, 'mappable'):
                thing = self.mappable
            else:
                raise RuntimeError, 'do not know what to colobar, pass in a thing'
        if rect is None:
            self.__colorbar = f.colorbar(thing, ax=self.a)
        else:
            cax = w.getNextAxes(rect=rect)
            self.__colorbar = f.colorbar(thing, cax=cax)
            if adjustPos:
                'adjust existing position to make room for colorbar'
                bbox = self.a.get_position().get_points()
                arect = (bbox[0,0], # left
                         bbox[0,1], # bottom
                         rect[0]-bbox[0,0]-0.02, # width
                         bbox[1,1]-bbox[0,1], # height
                         ) 
                self.a.set_position(arect) 
        self.show()
        return

def hist2D(xVals, yVals, bins, hRange = None, weights = None,
           **kwArgs
           ):
    '''
    Plot 2D histogram of data yVals versus xVals, with number of bins given by bins (int or 2-tuple)
    '''
    xedges, yedges, H = makeHist2D(xVals, yVals, bins, hRange=hRange, weights=weights)
    window, a, cax = plotHist2D(xedges, yedges, H, hRange=hRange, **kwArgs)
    return window, a, cax

def plotHist2D(
    xedges, yedges, H,
    hRange=None, 
    logScale = False,
    minCount=1,
    win=None,
    xlabel=None,
    ylabel=None,
    xformat=None,
    yformat=None,
    nXTics=0,
    nYTics=0,
    winArgs = {}
    ):
    '''
    winArgs can include things like title, relfigsize, dpi
    '''
    if win is None:
        window = PlotWin(-1, **winArgs) # dummy=True
        a      = window.getNextAxes(rect=(0.1,0.15,0.6,0.8))
        cax    = window.getNextAxes(rect=(0.80,0.15,0.05,0.8))
    else:
        window = win
        a      = window.getNextAxes()
        cax    = window.getNextAxes()
    if xlabel is not None:
        a.set_xlabel(xlabel)
    if ylabel is not None:
        a.set_ylabel(ylabel)
    if nXTics:
        assert hRange is not None, 'need hRange to do nXTics'
        a.set_xticks(autoTicks(hRange[0],nXTics,format=xformat))
    if nYTics:
        assert hRange is not None, 'need hRange to do nYTics'
        a.set_yticks(autoTicks(hRange[1],nYTics,format=yformat))
    if logScale:
        #pColor = a.pcolor(xedges, yedges, num.log(H+1))
        pColor = a.pcolor(xedges, yedges, H.T, norm=matplotlib.colors.LogNorm(vmin=minCount))
    else:
        pColor = a.pcolor(xedges, yedges, H.T)
    f = window.getFigure()
    f.colorbar(pColor, cax=cax)
    #win.haveXLabels()
    #win.haveYLabels()

    window.show()
    #win.save(filename=fname+'.eps')
    #win.save(filename=fname+'.png')
    return window, a, cax

def makeHist2D(xVals, yVals, bins, hRange = None, weights = None):
    import numpy as num

    if hRange is None:
        hR = ((xVals.min(), xVals.max()), (yVals.min(), yVals.max()))
    else:
        hR = hRange
    H, xedges, yedges = num.histogram2d(
        xVals,
        yVals,
        weights=weights,
        bins=bins, range=hR)
    return xedges, yedges, H
    
def main():
    import plotWrap
    import numpy as num
    import math
    inBits = True
    withPoints = False
    interactive = False

    t = num.arange(0.0,2.0+0.01,0.01)
    s = num.sin(2*math.pi*t)
    c = num.cos(2*math.pi*t)

    if withPoints:
        pw = plotWrap.PlotWrap(style='ro',aging=0.75)
        if inBits:
            for iSub in range(len(t)):
                pw([t[iSub]],[s[iSub]])
        else:
            pw(t, s)
    else:
        pw = plotWrap.PlotWrap(style='r-',aging=0.75,alphaMin=0.2,agingCull=False,bball='ro')
        if inBits:
            #for iSub in range(len(t)-1):
                #pw(t[iSub:iSub+2],s[iSub:iSub+2])
            for iSub in range(len(t)):
                pw([t[iSub]],[s[iSub]])
        else:
            pw(t, s)
    
    pWin = plotWrap.PlotWin(2,1,title='test PlotWin')
    p1 = plotWrap.PlotWrap(window=pWin,ylabel='s')
    p2 = plotWrap.PlotWrap(window=pWin,ylabel='c')
    p1(t,s,label='sin')
    p2(t,c,label='cos')
    
    if interactive: 
        import pylab as p
        p.ion()
        p.show()
    else:
        pw.save(filename='pw_test_1.eps')
        pw.destroy()
        pWin.save(filename='pw_test_2.eps')
        pWin.destroy()

if __name__ == '__main__':
    main()


