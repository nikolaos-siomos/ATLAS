import wx

from . import SettingsWindow
from atlas_ui.components.fields import *

class PreprocessorSettingsWindow (SettingsWindow):
    def __init__ (self, *args, **kwargs):
        super (PreprocessorSettingsWindow, self).__init__ ( *args, **kwargs )
        self.Build()
        
    def Build (self):
        self.sizer = wx.GridBagSizer (vgap = 20, hgap = 20)
        
        self.vertical_trimming = AtlasUIInputField.FromSetting ( self._settings.vertical_trimming, parent = self, label = "Vertical trimming?" )
        self.sizer.Add ( self.vertical_trimming, pos = wx.GBPosition (0, 0), flag = wx.EXPAND )
        
        self.vertical_limit = AtlasUIInputField.FromSetting ( self._settings.vertical_limit, parent = self, label = "Vertical limit" )
        self.sizer.Add ( self.vertical_limit, pos = wx.GBPosition (0, 1), flag = wx.EXPAND )
        
        self.channels = AtlasUIInputField.FromSetting ( self._settings.channels, parent = self, label = "Only process channels", choices = self._channels )
        self.sizer.Add ( self.channels, pos = wx.GBPosition (1, 0), span = wx.GBSpan (1, 2), flag = wx.EXPAND )
        
        self.exclude_telescope_type = AtlasUIInputField.FromSetting ( self._settings.exclude_telescope_type, parent = self, label = "Exclude telescope type" )
        self.sizer.Add ( self.exclude_telescope_type, pos = wx.GBPosition (2, 0), flag = wx.EXPAND )
        
        exclude_detection_mode_choices = {
            "a": "Analog",
            "p": "Photon counting"
        }
        self.exclude_acquisition_mode = AtlasUIInputField.FromSetting ( self._settings.exclude_acquisition_mode, parent = self, label = "Exclude acquisition mode" )
        self.sizer.Add (self.exclude_acquisition_mode, pos = wx.GBPosition (2, 1), flag = wx.EXPAND)

        self.exclude_scattering_type = AtlasUIInputField.FromSetting ( self._settings.exclude_scattering_type, parent = self, label = "Exclude acquisition mode" )
        self.sizer.Add ( self.exclude_scattering_type, pos = wx.GBPosition (3, 0), flag = wx.EXPAND )
        
        self.exclude_channel_subtype = AtlasUIInputField.FromSetting ( self._settings.exclude_channel_subtype, parent = self, label = "Exclude acquisition mode" )
        self.sizer.Add (self.exclude_channel_subtype, pos = wx.GBPosition (3, 1), flag = wx.EXPAND)
        
        self.sizer.SetFlexibleDirection(wx.VERTICAL)
        self.sizer.AddGrowableCol (idx = 0, proportion = 1)
        self.sizer.AddGrowableCol (idx = 1, proportion = 1)
        
        self.SetSizer ( self.sizer )