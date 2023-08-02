import wx
from atlas_ui.components.fields import *

from atlas_ui.settings import ATLASChannelSettings

typeEVT_CHANNEL_DELETE = wx.NewEventType ()
EVT_CHANNEL_DELETE = wx.PyEventBinder(typeEVT_CHANNEL_DELETE, 1)

typeEVT_CHANNEL_CHANGE = wx.NewEventType ()
EVT_CHANNEL_CHANGE = wx.PyEventBinder(typeEVT_CHANNEL_CHANGE, 1)

class ChannelDeleteEvent ( wx.PyCommandEvent ):
    def __init__(self, id):
        wx.PyCommandEvent.__init__(self, typeEVT_CHANNEL_DELETE, id)
        window = None
        
    def SetWindow ( self, window ):
        self.window = window
        
    def GetWindow ( self ):
        return self.window
        
class ChannelChangeEvent ( wx.PyCommandEvent ):
    def __init__(self, id):
        wx.PyCommandEvent.__init__(self, typeEVT_CHANNEL_CHANGE, id)
        window = None
        
    def SetWindow ( self, window ):
        self.window = window
        
    def GetWindow ( self ):
        return self.window

class ChannelData:
    def __init__ (
        self,
        telescope_type,
        channel_type,
        channel_subtype,
        dead_time,
        daq_trigger_offset,
        recorder_channel_id,
        laser,
        scc_channel_id,
        emitted_wavelength,
        detected_wavelength,
        channel_bandwidth,
        background_low_bin,
        background_high_bin,
        acquisition_mode,
        bins,
        data_acquisition_range,
        adc_resolution,
        range_resolution,
        laser_repetition_rate
    ):
        pass

class ChannelWindow (wx.ScrolledWindow):
    def __init__ (self, settings, *args, **kwargs):
        wx.ScrolledWindow.__init__(self, *args, **kwargs)
        self.SetScrollbars(20, 20, 55, 40)
        
        self.SetSettings ( settings )
        
    def SetSettings ( self, settings ):
        self._settings = settings
        self.Build()
        
    def Build ( self ):
        self.panel = wx.Panel (self)
        self.sizer = wx.GridBagSizer ( hgap = 10, vgap = 10 )
        
        self.smallhbox = wx.BoxSizer (wx.HORIZONTAL)
        
        self.telescope_type = AtlasUIInputField.FromSetting ( self._settings.telescope_type, parent = self.panel, label = "Telescope type" )
        self.smallhbox.Add ( self.telescope_type, 0, flag = wx.EXPAND )
        self.smallhbox.AddSpacer (5)

        self.channel_type = AtlasUIInputField.FromSetting ( self._settings.channel_type, parent = self.panel, label = "Channel type" )
        # self.sizer.Add ( self.channel_type, pos = wx.GBPosition(0, 1), flag = wx.EXPAND )
        self.smallhbox.Add ( self.channel_type, 0, flag = wx.EXPAND )
        self.smallhbox.AddSpacer (5)
        
        self.acquisition_mode = AtlasUIInputField.FromSetting ( self._settings.acquisition_mode, parent = self.panel, label = "Acquisition mode" )
        # self.optional_box_sizer.Add ( self.acquisition_mode, flag = wx.EXPAND )
        self.smallhbox.Add ( self.acquisition_mode, 0, flag = wx.EXPAND )
        self.smallhbox.AddSpacer (5)
        
        self.channel_subtype = AtlasUIInputField.FromSetting ( self._settings.channel_subtype, parent = self.panel, label = "Channel subtype" )
        # self.sizer.Add ( self.channel_subtype, pos = wx.GBPosition(0, 2), flag = wx.EXPAND )
        self.smallhbox.Add ( self.channel_subtype, 1, flag = wx.EXPAND )
        
        self.sizer.Add (self.smallhbox, pos = wx.GBPosition(0, 0), span = wx.GBSpan(1, 3), flag = wx.EXPAND)
        
        self.dead_time = AtlasUIInputField.FromSetting ( self._settings.dead_time, parent = self.panel, label = "Dead time [ns]" )
        self.sizer.Add ( self.dead_time, pos = wx.GBPosition(1, 0), flag = wx.EXPAND )
        
        self.daq_trigger_offset = AtlasUIInputField.FromSetting ( self._settings.daq_trigger_offset, parent = self.panel, label = "Trigger delay [bins]" )
        self.sizer.Add ( self.daq_trigger_offset, pos = wx.GBPosition(1, 1), flag = wx.EXPAND )
        
        self.recorder_channel_id = AtlasUIInputField.FromSetting ( self._settings.recorder_channel_id, parent = self.panel, label = "Channel ID" )
        self.sizer.Add ( self.recorder_channel_id, pos = wx.GBPosition(2, 0), flag = wx.EXPAND )
        
        self.laser = AtlasUIInputField.FromSetting ( self._settings.laser, parent = self.panel, label = "Laser" )
        self.sizer.Add ( self.laser, pos = wx.GBPosition(2, 1), flag = wx.EXPAND )
        
        self.scc_channel_id = AtlasUIInputField.FromSetting ( self._settings.scc_channel_id, parent = self.panel, label = "SCC ID" )
        self.sizer.Add ( self.scc_channel_id, pos = wx.GBPosition(2, 2), flag = wx.EXPAND )
        
        self.emitted_wavelength = AtlasUIInputField.FromSetting ( self._settings.emitted_wavelength, parent = self.panel, label = "Emitted wavelength [nm]" )
        self.sizer.Add ( self.emitted_wavelength, pos = wx.GBPosition(3, 1), flag = wx.EXPAND )
        
        self.detected_wavelength = AtlasUIInputField.FromSetting ( self._settings.detected_wavelength, parent = self.panel, label = "Detected wavelength [nm]" )
        self.sizer.Add ( self.detected_wavelength, pos = wx.GBPosition(3, 2), flag = wx.EXPAND )
        
        self.channel_bandwidth = AtlasUIInputField.FromSetting ( self._settings.channel_bandwidth, parent = self.panel, label = "Channel bandwidth [nm]" )
        self.sizer.Add ( self.channel_bandwidth, pos = wx.GBPosition(4, 0), flag = wx.EXPAND )
        
        self.background_range = AtlasUIInputField.FromSetting ( self._settings.background_range, parent = self.panel, label = "Background range [bins]" )
        self.sizer.Add ( self.background_range, pos = wx.GBPosition(4, 2), flag = wx.EXPAND )
        
        self.optional_box = wx.StaticBox (self.panel, label = "Optional Settings")
        self.optional_box_sizer = wx.GridBagSizer( hgap = 10, vgap = 10 )
        
        # self.optional_box_sizer.AddSpacer (10)
        
        self.total_bins = AtlasUIInputField.FromSetting ( self._settings.bins, parent = self.optional_box, label = "Total bins" )
        self.optional_box_sizer.Add ( self.total_bins, flag = wx.EXPAND | wx.TOP, border = 25, pos = wx.GBPosition(0, 0) )
        
        self.daq_range = AtlasUIInputField.FromSetting ( self._settings.data_acquisition_range, parent = self.optional_box, label = "DAQ range [mV]" )
        self.optional_box_sizer.Add ( self.daq_range, flag = wx.EXPAND, pos = wx.GBPosition(1, 1) )
        
        self.adc_resolution = AtlasUIInputField.FromSetting ( self._settings.analog_to_digital_resolution, parent = self.optional_box, label = "ADC resolution [bits]" )
        self.optional_box_sizer.Add ( self.adc_resolution, flag = wx.EXPAND, pos = wx.GBPosition(2, 0) )
        
        self.range_resolution = AtlasUIInputField.FromSetting ( self._settings.range_resolution, parent = self.optional_box, label = "Range resolution [m]" )
        self.optional_box_sizer.Add ( self.range_resolution, flag = wx.EXPAND, pos = wx.GBPosition(2, 1) )
        
        bottom_border = 30 if wx.GetOsVersion()[0] & wx.OS_UNIX else 10
        
        self.laser_repetition_rate = AtlasUIInputField.FromSetting ( self._settings.laser_repetition_rate, parent = self.optional_box, label = "Laser repetition rate [Hz]" )
        self.optional_box_sizer.Add ( self.laser_repetition_rate, flag = wx.EXPAND | wx.BOTTOM, border = bottom_border, pos = wx.GBPosition(3, 1) )
        
        # self.optional_box_sizer.AddSpacer(50)
        self.optional_box_sizer.AddGrowableCol(idx = 0, proportion = 1)
        self.optional_box_sizer.AddGrowableCol(idx = 1, proportion = 1)
        
        self.optional_box.SetSizer ( self.optional_box_sizer )
        self.sizer.Add ( self.optional_box, pos = wx.GBPosition (5, 0), span = wx.GBSpan ( 4, 3 ), flag = wx.EXPAND )
        
        self.sizer.AddGrowableCol ( 0, 1 )
        self.sizer.AddGrowableCol ( 1, 1 )
        self.sizer.AddGrowableCol ( 2, 1 )
        
        self.panel.SetSizer (self.sizer)
        self.window_sizer = wx.BoxSizer (wx.VERTICAL)
        
        self.delete_button = wx.Button (self, label = "Remove channel")
        self.delete_button.Bind ( wx.EVT_BUTTON, self.OnDeleteBtnClicked )
        
        self.window_sizer.Add ( self.delete_button, border = 10, flag = wx.ALL | wx.CENTER )
        self.window_sizer.Add (self.panel, border = 10, flag = wx.ALL | wx.EXPAND)
        
        self.SetSizer (self.window_sizer)
        
        self.Bind (wx.EVT_SPINCTRL, self.OnDataChanged)
        self.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnDataChanged)
        self.Bind (wx.EVT_TEXT, self.OnDataChanged)
        self.Bind (wx.EVT_CHECKBOX, self.OnDataChanged)
        self.Bind (wx.EVT_CHOICE, self.OnDataChanged)
        
        self.acquisition_mode.OnUpdate ( self.OnAcquisitionModeChanged )
        self.channel_type.OnUpdate ( self.OnChannelTypeChanged )
        
        self.OnChannelTypeChanged ( self.channel_type.GetValue() )
        self.OnAcquisitionModeChanged ( self.acquisition_mode.GetValue() )
        
    def OnDeleteBtnClicked ( self, e ):
        event = ChannelDeleteEvent ( self.GetId() )
        event.SetWindow ( self )
        self.GetEventHandler().ProcessEvent (event)
        
    def OnDataChanged ( self, e ):
        event = ChannelChangeEvent ( self.GetId() )
        event.SetWindow ( self )
        self.GetEventHandler().ProcessEvent (event)
        
    def GetChannelName ( self ):
        am = 'a' if self._settings.acquisition_mode() == '0' else 'p'
        return f"{self._settings.detected_wavelength():04.0f}{self._settings.telescope_type()}{self._settings.channel_type()}{am}{self._settings.channel_subtype()}"
        
    def GetChannelCode ( self ):
        return ""
        
    def GetSettings ( self ):
        return self._settings
        
    def Settings ( self ):
        return self.GetSettings()
        
    def OnAcquisitionModeChanged ( self, value ):
        if value == "0": # Analog
            self.daq_range.Enable()
            self.adc_resolution.Enable()
            self.dead_time.Disable()
        else: # Photon
            self.daq_range.Disable()
            self.adc_resolution.Disable()
            self.dead_time.Enable()
        
    def OnChannelTypeChanged ( self, value ):
        if value in [ "p", "c", "o", "x" ]:
            choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS"
            }
        elif value in [ "t" ]:
            choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS",
                "x": "Unspecified"
            }
        elif value in [ "v" ]:
            choices = {
                "n": "N2 Raman line",
                "o": "O2 Raman line",
                "w": "H2O Raman line",
                "c": "CH4 Raman line"
            }
        elif value in [ "r" ]:
            choices = {
                "h": "High Rotational Raman",
                "l": "Low Rotational Raman",
                "x": "Unspecified"
            }
        elif value in [ "a" ]:
            choices = {
                "a": "Mie (aerosol) HSRL signal",
                "m": "Molecular HSRL signal"
            }
        elif value in [ "f" ]:
            choices = {
                "b": "Broadband Fluorescence",
                "s": "Spectral Fluorescence"
            }
        else:
            choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS",
                "n": "N2 Raman line",
                "o": "O2 Raman line",
                "w": "H2O Raman line",
                "c": "CH4 Raman line",
                "h": "High J Rotational Raman",
                "l": "Low J Rotational Raman",
                "a": "Mie (aerosol) HSRL signal",
                "m": "Molecular HSRL signal",
                "b": "Broadband Fluorescence",
                "s": "Spectral Fluorescence",
                "x": "Unspecified"
            }
            
        self.channel_subtype.SetChoices ( choices )
