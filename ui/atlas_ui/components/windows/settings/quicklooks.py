import wx

from . import SettingsWindow
from atlas_ui.components.fields import *

class QuicklooksSettingsWindow (SettingsWindow):
    def __init__ (self, *args, **kwargs):
        super (QuicklooksSettingsWindow, self).__init__ ( *args, **kwargs )
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
        
        t_inputs = wx.BoxSizer ( orient = wx.VERTICAL )
        
        self.upper_t_lim = AtlasUIInputField.FromSetting ( self._settings.upper_t_lim, parent = self.limits_sizer.GetStaticBox(), label = "Time axis start [hh:mm]")
        self.lower_t_lim = AtlasUIInputField.FromSetting ( self._settings.lower_t_lim, parent = self.limits_sizer.GetStaticBox(), label = "Time axis end [hh:mm]")
        self.t_tick = AtlasUIInputField.FromSetting ( self._settings.t_tick, parent = self.limits_sizer.GetStaticBox(), label = "Time axis tick [hh:mm]")
        
        t_inputs.Add ( self.lower_t_lim, 0, wx.EXPAND )
        t_inputs.Add ( self.upper_t_lim, 0, wx.EXPAND )
        t_inputs.Add ( self.t_tick, 0, wx.EXPAND )
        
        y_inputs = wx.BoxSizer ( orient = wx.VERTICAL )
        
        self.y_lims = AtlasUIInputField.FromSetting ( self._settings.y_lims, parent = self.limits_sizer.GetStaticBox(), label = "Y Limits [km]")
        self.y_tick = AtlasUIInputField.FromSetting ( self._settings.y_tick, parent = self.limits_sizer.GetStaticBox(), label = "Y Tick [km]")
        
        y_inputs.Add ( self.y_lims, 0, wx.EXPAND )
        y_inputs.Add ( self.y_tick, 0, wx.EXPAND )
        
        z_inputs = wx.BoxSizer ( orient = wx.VERTICAL )
        
        self.z_lims = AtlasUIInputField.FromSetting ( self._settings.z_lims, parent = self.limits_sizer.GetStaticBox(), label = "Z Limits")
        
        z_inputs.Add ( self.z_lims, 0, wx.EXPAND )
        
        z_max_inputs = wx.BoxSizer ( orient = wx.VERTICAL )
        
        self.z_max_zone = AtlasUIInputField.FromSetting ( self._settings.z_max_zone, parent = self.limits_sizer.GetStaticBox(), label = "Z max Zone")
        
        z_max_inputs.Add ( self.z_max_zone, 0, wx.EXPAND )
        
        self.dpi = AtlasUIInputField.FromSetting ( self._settings.dpi, parent = self.limits_sizer.GetStaticBox(), label = "Image DPI")
        
        self.color_reduction = AtlasUIInputField.FromSetting ( self._settings.color_reduction, parent = self.limits_sizer.GetStaticBox(), label = "Color reduction?")
        self.use_range = AtlasUIInputField.FromSetting ( self._settings.use_range, parent = self.limits_sizer.GetStaticBox(), label = "Use distance?")
        self.use_log_scale = AtlasUIInputField.FromSetting ( self._settings.use_log_scale, parent = self.limits_sizer.GetStaticBox(), label = "Use log scale?")
        
        limits_row_1.Add ( self.use_log_scale )
        # limits_row_1.AddSpacer (10)
        limits_row_1.Add ( self.use_range )
        # limits_row_1.AddSpacer (10)
        limits_row_1.Add ( self.color_reduction )
        
        limits_row_2.Add ( x_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        # limits_row_2.AddSpacer (10)
        limits_row_2.Add ( t_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        # limits_row_2.AddSpacer (10)
        limits_row_2.Add ( y_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        
        limits_row_3.Add ( z_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        # limits_row_3.AddSpacer (10)
        limits_row_3.Add ( z_max_inputs, 1, wx.ALIGN_CENTER_VERTICAL )
        # limits_row_3.AddSpacer (10)
        limits_row_3.Add ( self.dpi, 1, wx.ALIGN_CENTER_VERTICAL )
        
        self.limits_sizer.Add ( limits_row_1, 0, wx.EXPAND )
        # self.limits_sizer.AddSpacer (10)
        self.limits_sizer.Add ( limits_row_2, 0, wx.EXPAND )
        # self.limits_sizer.AddSpacer (10)
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
        
        self.use_log_scale.OnUpdate ( self.OnLogScaleChanged )
        
        
    def OnLogScaleChanged ( self, value ):
        print (f"Here: {value} {value is True} {value==True}")
        if value:
            try:
                self.z_lims.SetIncrement ( 0.00001 )
                self.z_lims.SetValue ( (0.00005, 1.) )
                self._settings.z_lims.SetDefault ( (0.00005, 1.) )
            except Exception as e:
                # print (e)
                pass
        else:
            self.z_lims.SetValue ( (0., 1.) )
            self._settings.z_lims.SetDefault ( (0., 1.) )
            self.z_lims.SetIncrement ( 0.01 )