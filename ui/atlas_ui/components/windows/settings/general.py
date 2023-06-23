import wx

from . import SettingsWindow
from atlas_ui.components.fields import *

class GeneralSettingsWindow (SettingsWindow):
    def __init__ (self, *args, **kwargs):
        super (GeneralSettingsWindow, self).__init__ ( *args, **kwargs )
        self.Build()
        
    def Build (self):
        self.sizer = wx.FlexGridSizer ( cols = 2, hgap = 20, vgap = 20 )
        self.sizer.SetFlexibleDirection(wx.VERTICAL)
        self.sizer.AddGrowableCol (idx = 0, proportion = 1)
        self.sizer.AddGrowableCol (idx = 1, proportion = 1)
        
        self.sizer.AddSpacer(15)
        self.sizer.AddSpacer(15)
        
        self.isday = wx.CheckBox ( parent = self, label = "Day measurements?" )
        self.isday.SetValue ( True )
        self.sizer.Add ( self.isday )

        self.newdata = wx.CheckBox ( parent = self, label = "New data?" )
        self.newdata.SetValue ( True )
        self.sizer.Add ( self.newdata )
        
        self.visualize = wx.CheckBox ( parent = self, label = "Create plots?" )
        self.visualize.SetValue ( True )
        self.sizer.Add ( self.visualize )
        
        self.sizer.AddStretchSpacer()
        
        process_qa_options = {
            "ray": "Rayleigh Fit",
            "tlc": "Telecover Test",
            "pcb": "Polarization Calibration"
        }
        
        self.process_qa = MultipleChoiceField ( self, label = "Process QA Tests", choices = process_qa_options )
        self.sizer.Add ( self.process_qa, 1, wx.EXPAND )
        
        process_qcx_options = {
            "ray": "Rayleigh Fit",
            "tlc": "Telecover Test",
            "pcb": "Polarization Calibration"
        }
        self.process_qcx = MultipleChoiceField ( self, label = "Process Quicklooks", choices = process_qcx_options )
        self.sizer.Add ( self.process_qcx, 1, wx.EXPAND )
        
        self.SetSizer ( self.sizer )