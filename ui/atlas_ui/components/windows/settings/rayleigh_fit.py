import wx

from . import SettingsWindow
from atlas_ui.components.fields import *

class RayleighFitSettingsWindow (SettingsWindow):
    def __init__ (self, *args, **kwargs):
        super (RayleighFitSettingsWindow, self).__init__ ( *args, **kwargs )
        self.Build()
        
    def Build (self):
        self.sizer = wx.GridBagSizer (vgap = 20, hgap = 20)
        
        self.channels = AtlasUIInputField.FromSetting ( self._settings.channels, parent = self, label = "Only process channels", choices = self._channels )
        self.sizer.Add ( self.channels, pos = wx.GBPosition (1, 0), span = wx.GBSpan (1, 2), flag = wx.EXPAND )
        
        self.exclude_telescope_type = AtlasUIInputField.FromSetting ( self._settings.exclude_telescope_type, parent = self, label = "Exclude telescope type")
        self.sizer.Add ( self.exclude_telescope_type, pos = wx.GBPosition (2, 0), flag = wx.EXPAND )
        
        self.exclude_acquisition_mode = AtlasUIInputField.FromSetting ( self._settings.exclude_acquisition_mode, parent = self, label = "Exclude detection mode")
        self.sizer.Add (self.exclude_acquisition_mode, pos = wx.GBPosition (2, 1), flag = wx.EXPAND)

        self.exclude_scattering_type = AtlasUIInputField.FromSetting ( self._settings.exclude_scattering_type, parent = self, label = "Exclude scattering type")
        self.sizer.Add ( self.exclude_scattering_type, pos = wx.GBPosition (3, 0), flag = wx.EXPAND )

        self.exclude_channel_subtype = AtlasUIInputField.FromSetting ( self._settings.exclude_channel_subtype, parent = self, label = "Exclude channel subtype")
        self.sizer.Add (self.exclude_channel_subtype, pos = wx.GBPosition (3, 1), flag = wx.EXPAND)
        
        self.limits_sizer = wx.StaticBoxSizer ( parent = self, orient = wx.VERTICAL, label = "Plotting" )
        limits_row_1 = wx.BoxSizer ( orient = wx.HORIZONTAL )
        limits_row_2 = wx.BoxSizer ( orient = wx.HORIZONTAL )
        limits_row_3 = wx.BoxSizer ( orient = wx.HORIZONTAL )
        
        x_inputs = wx.BoxSizer ( orient = wx.VERTICAL )
        
        self.x_lims = AtlasUIInputField.FromSetting ( self._settings.x_lims, parent = self.limits_sizer.GetStaticBox(), label = "X Limits [bin]")
        self.x_tick = AtlasUIInputField.FromSetting ( self._settings.x_tick, parent = self.limits_sizer.GetStaticBox(), label = "X Tick [bin]")
        
        x_inputs.Add ( self.x_lims, 0, wx.EXPAND )
        x_inputs.Add ( self.x_tick, 0, wx.EXPAND )
        
        y_inputs = wx.BoxSizer ( orient = wx.VERTICAL )
        
        self.y_lims = AtlasUIInputField.FromSetting ( self._settings.y_lims, parent = self.limits_sizer.GetStaticBox(), label = "Y Limits [km]")
        
        y_inputs.Add ( self.y_lims, 0, wx.EXPAND )
        
        self.dpi = AtlasUIInputField.FromSetting ( self._settings.dpi, parent = self.limits_sizer.GetStaticBox(), label = "Image DPI")
        
        self.color_reduction = AtlasUIInputField.FromSetting ( self._settings.color_reduction, parent = self.limits_sizer.GetStaticBox(), label = "Color reduction?")
        self.use_range = AtlasUIInputField.FromSetting ( self._settings.use_range, parent = self.limits_sizer.GetStaticBox(), label = "Use distance?")
        self.normalization_region = AtlasUIInputField.FromSetting ( self._settings.normalization_region, parent = self.limits_sizer.GetStaticBox(), label = "Normalization region [km]")
        self.use_lin_scale = AtlasUIInputField.FromSetting ( self._settings.use_lin_scale, parent = self.limits_sizer.GetStaticBox(), label = "Use linear scale?")
        
        limits_row_1.Add ( self.use_range )
        limits_row_1.Add ( self.color_reduction )
        limits_row_1.Add ( self.use_lin_scale )
        
        limits_row_2.Add ( x_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        limits_row_2.Add ( y_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        
        limits_row_3.Add ( self.normalization_region, 1, wx.ALIGN_CENTER_VERTICAL )
        limits_row_3.Add ( self.dpi, 1, wx.ALIGN_CENTER_VERTICAL )
        
        self.limits_sizer.Add ( limits_row_1, 0, wx.EXPAND )
        self.limits_sizer.Add ( limits_row_2, 0, wx.EXPAND )
        self.limits_sizer.Add ( limits_row_3, 0, wx.EXPAND )
        
        self.sizer.Add ( self.limits_sizer, pos = wx.GBPosition (4, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        
        self.smoothing_sizer = wx.StaticBoxSizer ( parent = self, orient = wx.HORIZONTAL, label = "Smoothing" )
        
        self.smooth = AtlasUIInputField.FromSetting ( self._settings.smooth, parent = self.smoothing_sizer.GetStaticBox(), label = "Smoothing enabled?")
        self.smoothing_exponential = AtlasUIInputField.FromSetting ( self._settings.smoothing_exponential, parent = self.smoothing_sizer.GetStaticBox(), label = "Exponentian smoothing?")
        self.smoothing_range = AtlasUIInputField.FromSetting ( self._settings.smoothing_range, parent = self.smoothing_sizer.GetStaticBox(), label = "Smoothing range [km]")
        self.smoothing_window = AtlasUIInputField.FromSetting ( self._settings.smoothing_range, parent = self.smoothing_sizer.GetStaticBox(), label = "Smoothing window [m]")
        
        v_box_sizer = wx.BoxSizer (wx.VERTICAL)
        v_box_sizer.Add ( self.smooth )
        v_box_sizer.Add ( self.smoothing_exponential )
        
        self.smoothing_sizer.Add ( v_box_sizer, 0, wx.ALIGN_CENTER_VERTICAL )
        self.smoothing_sizer.Add ( self.smoothing_range, 1 )
        self.smoothing_sizer.Add ( self.smoothing_window, 1 )
        
        self.sizer.Add ( self.smoothing_sizer, pos = wx.GBPosition (5, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        
        self.sizer.SetFlexibleDirection(wx.VERTICAL)
        self.sizer.AddGrowableCol (idx = 0, proportion = 1)
        self.sizer.AddGrowableCol (idx = 1, proportion = 1)
        
        self.SetSizer ( self.sizer )