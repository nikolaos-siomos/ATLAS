import wx

from . import SettingsWindow
from atlas_ui.components.fields import *
from atlas_ui.settings import ChannelPair

import copy

class ChannelPairDialog(wx.Dialog, wx.ScrolledWindow):
    def __init__ (self, parent, available_channels, title="Edit channel pair", channel_pair = None):
        super(ChannelPairDialog, self).__init__(None, title = title)
        
        if len(available_channels) < 2:
            raise ValueError("Cannot have a channel pair when there are less than 2 channels.")
            
        if channel_pair is None:
            channel_pair = ChannelPair (available_channels[0], available_channels[1])
            
        self.initial_channel_pair = copy.deepcopy(channel_pair)
        self.channel_pair = channel_pair
        
        panel = wx.Panel(self)
        sizer = wx.GridBagSizer (vgap = 5, hgap = 5)
        
        self.reflected_channel = SingleChoiceField( parent = panel, choices = available_channels, value = self.channel_pair.ReflectedChannel(), label = "Reflected channel" )
        self.reflected_channel.choice.Bind (wx.EVT_CHOICE, self.OnReflectedChannel)
        sizer.Add ( self.reflected_channel, pos = wx.GBPosition (0, 0), flag = wx.EXPAND )
        
        self.transmitted_channel = SingleChoiceField( parent = panel, choices = available_channels, value = self.channel_pair.TransmittedChannel(), label = "Transmitted channel" )
        self.transmitted_channel.choice.Bind (wx.EVT_CHOICE, self.OnTransmittedChannel)
        sizer.Add ( self.transmitted_channel, pos = wx.GBPosition (0, 1), flag = wx.EXPAND )
        
        self.k = NumberField( parent = panel, increment = .01, value = self.channel_pair.K(), label = "K value" )
        self.k.spinner.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnK)
        sizer.Add ( self.k, pos = wx.GBPosition (1, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        
        self.g_r = NumberField( parent = panel, increment = .01, value = self.channel_pair.G_R(), label = "G value (Reflected)" )
        self.g_r.spinner.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnG_R)
        sizer.Add ( self.g_r, pos = wx.GBPosition (2, 0), flag = wx.EXPAND )
        
        self.g_t = NumberField( parent = panel, increment = .01, value = self.channel_pair.G_T(), label = "G value (Transmitted)" )
        self.g_t.spinner.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnG_T)
        sizer.Add ( self.g_t, pos = wx.GBPosition (2, 1), flag = wx.EXPAND )
        
        self.h_r = NumberField( parent = panel, increment = .01, value = self.channel_pair.H_R(), label = "H value (Reflected)" )
        self.h_r.spinner.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnH_R)
        sizer.Add ( self.h_r, pos = wx.GBPosition (3, 0), flag = wx.EXPAND )
        
        self.h_t = NumberField( parent = panel, increment = .01, value = self.channel_pair.H_T(), label = "H value (Transmitted)" )
        self.h_t.spinner.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnH_T)
        sizer.Add ( self.h_t, pos = wx.GBPosition (3, 1), flag = wx.EXPAND )
        
        self.rt_t_ratio = NumberField( parent = panel, increment = .01, value = self.channel_pair.RT_T_Ratio() )
        self.rt_t_ratio.spinner.Bind (wx.EVT_SPINCTRLDOUBLE, self.OnRT_T_Ratio)
        sizer.Add ( self.rt_t_ratio, pos = wx.GBPosition (4, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        
        panel.SetSizer(sizer)
        
        BottomHBox=wx.BoxSizer(wx.HORIZONTAL)
        
        self.applyBtn = wx.Button(self, label="OK")
        self.applyBtn.Bind(wx.EVT_BUTTON, self.OnApplyBtnClicked)
        self.cancelBtn = wx.Button(self, label="Cancel")
        self.cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancelBtnClicked)
        
        BottomHBox.Add(self.applyBtn, proportion=0, flag=wx.RIGHT, border=5)
        BottomHBox.Add(self.cancelBtn, proportion=0, flag=wx.LEFT, border=5)
        
        self.sizer = wx.BoxSizer (wx.VERTICAL)
        self.sizer.Add ( panel, proportion = 1, flag=wx.EXPAND|wx.BOTTOM, border = 10 )
        self.sizer.Add ( BottomHBox, flag=wx.CENTER|wx.ALL, border = 10 )
        
        self.SetSizer(self.sizer)
        
    def OnReflectedChannel (self, e):
        self.channel_pair.ReflectedChannel ( self.reflected_channel.choice.GetString( self.reflected_channel.choice.GetSelection() ) )
        
    def OnTransmittedChannel (self, e):
        self.channel_pair.TransmittedChannel ( self.transmitted_channel.choice.GetString( self.transmitted_channel.choice.GetSelection() ) )
        
    def OnK (self, e):
        self.channel_pair.K ( self.k.spinner.GetValue() )
        
    def OnG_R (self, e):
        self.channel_pair.G_R ( self.g_r.spinner.GetValue() )
        
    def OnG_T (self, e):
        self.channel_pair.G_T ( self.g_t.spinner.GetValue() )
        
    def OnH_R (self, e):
        self.channel_pair.H_R ( self.h_r.spinner.GetValue() )
        
    def OnH_T (self, e):
        self.channel_pair.H_T ( self.h_t.spinner.GetValue() )
        
    def OnRT_T_Ratio (self, e):
        self.channel_pair.RT_T_Ratio ( self.rt_t_ratio.spinner.GetValue() )
        
    def OnApplyBtnClicked (self, e):
        self.EndModal(wx.ID_OK)
        
    def OnCancelBtnClicked (self, e):
        self.channel_pair = self.initial_channel_pair
        self.Destroy()

class ChannelPairsDialog(wx.Dialog):
    def __init__(self, parent, available_channels, title="Select channel pairs", channel_pairs = []):
        super(ChannelPairsDialog, self).__init__(None, title = title)
        
        self.available_channels = available_channels
        self.initial_channel_pairs = channel_pairs
        self.channel_pairs = channel_pairs
        
        panel = wx.Panel(self)
        panelHBox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.list = wx.ListBox ( panel, style=wx.LB_SINGLE )
        self.list.Bind(wx.EVT_LISTBOX, self.OnPairSelected)
        panelHBox.Add (self.list, 1, wx.EXPAND)
        
        # ADD LIST OF CHANNEL PAIRS
        
        buttonVBox = wx.BoxSizer(wx.VERTICAL)
        
        self.addChannelPairBtn = wx.Button ( panel, label = "Add new pairing" )
        self.addChannelPairBtn.Bind(wx.EVT_BUTTON, self.OnAddPair)
        buttonVBox.Add ( self.addChannelPairBtn, 1, wx.EXPAND )
        
        self.editChannelPairBtn = wx.Button ( panel, label = "Edit pairing" )
        self.editChannelPairBtn.Bind(wx.EVT_BUTTON, self.OnEditPair)
        buttonVBox.Add ( self.editChannelPairBtn, 1, wx.EXPAND )
        
        self.removeChannelPairBtn = wx.Button ( panel, label = "Remove pairing" )
        self.removeChannelPairBtn.Bind(wx.EVT_BUTTON, self.OnRemovePair)
        buttonVBox.Add ( self.removeChannelPairBtn, 1, wx.EXPAND )
        
        panelHBox.Add (buttonVBox)
        
        panel.SetSizer (panelHBox)
        
        BottomHBox=wx.BoxSizer(wx.HORIZONTAL)
        
        self.applyBtn = wx.Button(self, label="OK")
        self.applyBtn.Bind(wx.EVT_BUTTON, self.OnApplyBtnClicked)
        self.cancelBtn = wx.Button(self, label="Cancel")
        self.cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancelBtnClicked)
        
        BottomHBox.Add(self.applyBtn, proportion=0, flag=wx.RIGHT, border=5)
        BottomHBox.Add(self.cancelBtn, proportion=0, flag=wx.LEFT, border=5)
        
        self.sizer = wx.BoxSizer (wx.VERTICAL)
        self.sizer.Add ( panel, proportion = 1, flag=wx.EXPAND|wx.BOTTOM, border = 10 )
        self.sizer.Add ( BottomHBox, flag=wx.CENTER|wx.ALL, border = 10 )
        
        self.SetSizer (self.sizer)
        
        self.PopulateListBox()
        
    def RefreshActionButtons ( self ):
        item = self.list.GetSelection()
        
        if item == wx.NOT_FOUND:
            self.editChannelPairBtn.Disable()
            self.removeChannelPairBtn.Disable()
        else:
            self.editChannelPairBtn.Enable()
            self.removeChannelPairBtn.Enable()
        
    def PopulateListBox (self):
        items = [ f"{cp.ReflectedChannel()} - {cp.TransmittedChannel()}" for cp in self.channel_pairs ]
        self.list.Set (items)
        
        for index in self.list.GetSelections():
            self.list.Deselect(index)
            
        self.RefreshActionButtons()
            
    def OnPairSelected (self, e):
        self.RefreshActionButtons()
        
    def OnAddPair (self, e):
        dlg = ChannelPairDialog(self, self.available_channels, channel_pair = None)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.channel_pairs.append(dlg.channel_pair)
        
        self.PopulateListBox()
        
        dlg.Destroy()
        
    def OnEditPair (self, e):
        item = self.list.GetSelection()
        
        if item == wx.NOT_FOUND:
            raise ValueError ("Trying to edit non-existing channel pair!")
            
        dlg = ChannelPairDialog(self, self.available_channels, channel_pair = self.channel_pairs[item])
        
        if dlg.ShowModal() == wx.ID_OK:
            self.channel_pairs[item] = dlg.channel_pair
        
        self.PopulateListBox()
        dlg.Destroy()
        
    def OnRemovePair (self, e):
        item = self.list.GetSelection()
        
        if item == wx.NOT_FOUND:
            raise ValueError ("Trying to delete non-existing channel pair!")
            
        self.channel_pairs.pop (item)
        self.PopulateListBox()
        
    def OnApplyBtnClicked (self, e):
        self.EndModal(wx.ID_OK)
        
    def OnCancelBtnClicked (self, e):
        self.channel_pairs = self.initial_channel_pairs
        self.Destroy()

class ChannelPairsField (AtlasUIInputField):
    def __init__ ( self, available_channels, *args, **kwargs ):
        AtlasUIInputField.__init__ ( self, *args, **kwargs )
        
        self.available_channels = available_channels
        
        self.input_sizer = wx.BoxSizer (orient = wx.HORIZONTAL)
        
        self.text = wx.TextCtrl ( self.input_window, style = wx.TE_READONLY )
        
        self.btn = wx.Button ( self.input_window, label = "Select pairings..." )
        self.btn.Bind ( wx.EVT_BUTTON, self.OnBtnClicked )
        
        self.input_sizer.Add (self.text, 1, wx.EXPAND)
        self.input_sizer.AddSpacer (5)
        self.input_sizer.Add (self.btn, 0, wx.EXPAND)
        
        self.input_window.SetSizer (self.input_sizer)
        
        self.SetValue ( self.value )
        self.Select ( self.IsSelected() )
        
    def OnBtnClicked (self, event):
        dlg = ChannelPairsDialog ( self, self.available_channels, channel_pairs = self.value)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.SetValue ( dlg.channel_pairs )
        else:
            self.SetValue ( [] )
            
        dlg.Destroy()
        
    def GetValue ( self ):
        return self.value
        
    def SetValue (self, value = None):
        final_value = []
        
        if value is None or not isinstance(value, list):
            pass
        else:
            for channel_pair in value:
                try:
                    if channel_pair.ReflectedChannel() in self.available_channels and channel_pair.TransmittedChannel() in self.available_channels:
                        final_value.append ( channel_pair )
                except Exception:
                    continue
            
        self.value = final_value
        
        self.text.SetValue(", ".join([ f"{cp.ReflectedChannel()} - {cp.TransmittedChannel()}" for cp in self.value ]))
            
        if self.setting:
            self.setting ( final_value )

class PolarizationCalibrationSettingsWindow (SettingsWindow):
    def __init__ (self, *args, **kwargs):
        super (PolarizationCalibrationSettingsWindow, self).__init__ ( *args, **kwargs )
        self.Build()
        
    def Build (self):
        self.sizer = wx.GridBagSizer (vgap = 20, hgap = 20)
        
        self.channel_pair_field = ChannelPairsField (parent = self, available_channels = list(self._channels.keys()), setting = self._settings.channel_pairs, label = "Channel pairs")
        self.sizer.Add ( self.channel_pair_field, pos = wx.GBPosition(0, 0), span = wx.GBSpan (1, 2), flag = wx.EXPAND )
        
        # self.calibration_height = AtlasUIInputField.FromSetting ( self._settings.calibration_height, parent = self, label = "Calibration height [km]" )
        # self.sizer.Add ( self.calibration_height, pos = wx.GBPosition (1, 0), flag = wx.EXPAND )
        
        # self.rayleigh_height = AtlasUIInputField.FromSetting ( self._settings.rayleigh_height, parent = self, label = "Rayleigh height [km]" )
        # self.sizer.Add ( self.rayleigh_height, pos = wx.GBPosition (1, 1), flag = wx.EXPAND )
        
        # self.half_calibration_window = AtlasUIInputField.FromSetting ( self._settings.half_calibration_window, parent = self, label = "Calibration half window [m]" )
        # self.sizer.Add ( self.half_calibration_window, pos = wx.GBPosition (2, 0), flag = wx.EXPAND )
        
        # self.half_rayleigh_window = AtlasUIInputField.FromSetting ( self._settings.half_rayleigh_window, parent = self, label = "Rayleigh half window [m]" )
        # self.sizer.Add ( self.half_rayleigh_window, pos = wx.GBPosition (2, 1), flag = wx.EXPAND )
        
        self.calibration_region = AtlasUIInputField.FromSetting ( self._settings.calibration_region, parent = self, label = "Calibration region [km]" )
        self.sizer.Add ( self.calibration_region, pos = wx.GBPosition (1, 0), flag = wx.EXPAND )
        
        self.rayleigh_region = AtlasUIInputField.FromSetting ( self._settings.rayleigh_region, parent = self, label = "Rayleigh region [km]" )
        self.sizer.Add ( self.rayleigh_region, pos = wx.GBPosition (1, 1), flag = wx.EXPAND )
        
        self.x_lims_calibration = AtlasUIInputField.FromSetting ( self._settings.x_lims_calibration, parent = self, label = "Calibration X Lims [m]" )
        self.sizer.Add ( self.x_lims_calibration, pos = wx.GBPosition (2, 0), flag = wx.EXPAND )
        
        self.x_lims_rayleigh = AtlasUIInputField.FromSetting ( self._settings.x_lims_rayleigh, parent = self, label = "Rayleigh X Lims [m]" )
        self.sizer.Add ( self.x_lims_rayleigh, pos = wx.GBPosition (2, 1), flag = wx.EXPAND )
        
        self.x_tick_calibration = AtlasUIInputField.FromSetting ( self._settings.x_tick_calibration, parent = self, label = "Calibration X Tick [km]" )
        self.sizer.Add ( self.x_tick_calibration, pos = wx.GBPosition (3, 0), flag = wx.EXPAND )
        
        self.x_tick_rayleigh = AtlasUIInputField.FromSetting ( self._settings.x_tick_rayleigh, parent = self, label = "Rayleigh X Tick [km]" )
        self.sizer.Add ( self.x_tick_rayleigh, pos = wx.GBPosition (3, 1), flag = wx.EXPAND )
        
        self.y_lims_calibration = AtlasUIInputField.FromSetting ( self._settings.y_lims_calibration, parent = self, label = "Calibration Y Lims [km]" )
        self.sizer.Add ( self.y_lims_calibration, pos = wx.GBPosition (4, 0), flag = wx.EXPAND )
        
        self.y_lims_rayleigh = AtlasUIInputField.FromSetting ( self._settings.y_lims_rayleigh, parent = self, label = "Rayleigh Y Lims [km]" )
        self.sizer.Add ( self.y_lims_rayleigh, pos = wx.GBPosition (4, 1), flag = wx.EXPAND )
        
        self.color_reduction = AtlasUIInputField.FromSetting ( self._settings.color_reduction, parent = self, label = "Color reduction?")
        self.sizer.Add ( self.color_reduction, pos = wx.GBPosition (5, 0), flag = wx.EXPAND )
        self.dpi = AtlasUIInputField.FromSetting ( self._settings.dpi, parent = self, label = "Image DPI")
        self.sizer.Add ( self.dpi, pos = wx.GBPosition (5, 1), flag = wx.EXPAND )
        
        # self.smoothing_sizer_calibration = wx.StaticBoxSizer ( parent = self, orient = wx.HORIZONTAL, label = "Smoothing (Calibration)" )
        
        # self.smooth_calibration = AtlasUIInputField.FromSetting ( self._settings.smooth_calibration, parent = self.smoothing_sizer_calibration.GetStaticBox(), label = "Smoothing enabled?")
        # self.smoothing_exponential_calibration = AtlasUIInputField.FromSetting ( self._settings.smoothing_exponential_calibration, parent = self.smoothing_sizer_calibration.GetStaticBox(), label = "Exponential smoothing?")
        # self.smoothing_range_calibration = AtlasUIInputField.FromSetting ( self._settings.smoothing_range_calibration, parent = self.smoothing_sizer_calibration.GetStaticBox(), label = "Smoothing range [km]")
        # self.smoothing_window_calibration = AtlasUIInputField.FromSetting ( self._settings.smoothing_window_calibration, parent = self.smoothing_sizer_calibration.GetStaticBox(), label = "Smoothing window [m]")
        
        # v_box_sizer_calibration = wx.BoxSizer (wx.VERTICAL)
        # v_box_sizer_calibration.Add ( self.smooth_calibration )
        # v_box_sizer_calibration.Add ( self.smoothing_exponential_calibration )
        
        # self.smoothing_sizer_calibration.Add ( v_box_sizer_calibration, 0, wx.ALIGN_CENTER_VERTICAL )
        # self.smoothing_sizer_calibration.Add ( self.smoothing_range_calibration, 1 )
        # self.smoothing_sizer_calibration.Add ( self.smoothing_window_calibration, 1 )
        
        # self.smoothing_sizer_rayleigh = wx.StaticBoxSizer ( parent = self, orient = wx.HORIZONTAL, label = "Smoothing (Rayleigh)" )
        
        # self.smooth_rayleigh = AtlasUIInputField.FromSetting ( self._settings.smooth_rayleigh, parent = self.smoothing_sizer_rayleigh.GetStaticBox(), label = "Smoothing enabled?")
        # self.smoothing_exponential_rayleigh = AtlasUIInputField.FromSetting ( self._settings.smoothing_exponential_rayleigh, parent = self.smoothing_sizer_rayleigh.GetStaticBox(), label = "Exponential smoothing?")
        # self.smoothing_range_rayleigh = AtlasUIInputField.FromSetting ( self._settings.smoothing_range_rayleigh, parent = self.smoothing_sizer_rayleigh.GetStaticBox(), label = "Smoothing range [km]")
        # self.smoothing_window_rayleigh = AtlasUIInputField.FromSetting ( self._settings.smoothing_window_rayleigh, parent = self.smoothing_sizer_rayleigh.GetStaticBox(), label = "Smoothing window [m]")
        
        # v_box_sizer_rayleigh = wx.BoxSizer (wx.VERTICAL)
        # v_box_sizer_rayleigh.Add ( self.smooth_rayleigh )
        # v_box_sizer_rayleigh.Add ( self.smoothing_exponential_rayleigh )
        
        # self.smoothing_sizer_rayleigh.Add ( v_box_sizer_rayleigh, 0, wx.ALIGN_CENTER_VERTICAL )
        # self.smoothing_sizer_rayleigh.Add ( self.smoothing_range_rayleigh, 1 )
        # self.smoothing_sizer_rayleigh.Add ( self.smoothing_window_rayleigh, 1 )
        
        self.smoothing_sizer = wx.StaticBoxSizer ( parent = self, orient = wx.HORIZONTAL, label = "Smoothing" )
        
        self.smooth = AtlasUIInputField.FromSetting ( self._settings.smooth, parent = self.smoothing_sizer.GetStaticBox(), label = "Smoothing enabled?")
        self.smoothing_exponential = AtlasUIInputField.FromSetting ( self._settings.smoothing_exponential, parent = self.smoothing_sizer.GetStaticBox(), label = "Exponential smoothing?")
        self.smoothing_range = AtlasUIInputField.FromSetting ( self._settings.smoothing_range, parent = self.smoothing_sizer.GetStaticBox(), label = "Smoothing range [km]")
        self.smoothing_window = AtlasUIInputField.FromSetting ( self._settings.smoothing_window, parent = self.smoothing_sizer.GetStaticBox(), label = "Smoothing window [m]")
        
        v_box_sizer = wx.BoxSizer (wx.VERTICAL)
        v_box_sizer.Add ( self.smooth )
        v_box_sizer.Add ( self.smoothing_exponential )
        
        self.smoothing_sizer.Add ( v_box_sizer, 0, wx.ALIGN_CENTER_VERTICAL )
        self.smoothing_sizer.Add ( self.smoothing_range, 1 )
        self.smoothing_sizer.Add ( self.smoothing_window, 1 )
        
        # self.sizer.Add ( self.smoothing_sizer_calibration, pos = wx.GBPosition (6, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        # self.sizer.Add ( self.smoothing_sizer_rayleigh, pos = wx.GBPosition (7, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        self.sizer.Add ( self.smoothing_sizer, pos = wx.GBPosition (6, 0), span = wx.GBSpan(1, 2), flag = wx.EXPAND )
        
        self.sizer.SetFlexibleDirection(wx.VERTICAL)
        self.sizer.AddGrowableCol (idx = 0, proportion = 1)
        self.sizer.AddGrowableCol (idx = 1, proportion = 1)
        
        self.SetSizer ( self.sizer )