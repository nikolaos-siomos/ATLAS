APP_TITLE = "ATLAS GUI"

import copy
import pickle
import wx

import atlas_ui.components.windows.settings

from atlas_ui.components.windows.settings.general import *
from atlas_ui.components.windows.settings.converter import *
from atlas_ui.components.windows.settings.preprocessor import *
from atlas_ui.components.windows.settings.quicklooks import *
from atlas_ui.components.windows.settings.rayleigh_fit import *
from atlas_ui.components.windows.settings.telecover import *
from atlas_ui.components.windows.settings.polarization_calibration import *
from atlas_ui.components.fields import *

from atlas_ui.components.windows.channel import *

from atlas_ui.settings import ATLASSettings, ATLASConfiguration
        
class SettingsDialog(wx.Dialog):
    def __init__(self, parent, settings, channels = {}, title='Settings'):
        super(SettingsDialog, self).__init__(None, title=title)
        
        self.SetSize((800, 700))
        
        self.panel = wx.Panel(self)
        
        self.tabs = wx.Notebook(self.panel, style = wx.NB_TOP)
        
        self.settings = settings
        self.initial_settings = settings
        
        self.converter_settings_window = ConverterSettingsWindow ( self.tabs, settings = self.settings.converter, channels = channels )
        self.preprocessor_settings_window = PreprocessorSettingsWindow ( self.tabs, settings = self.settings.preprocessor, channels = channels )
        self.quicklooks_settings_window = QuicklooksSettingsWindow(self.tabs, settings = self.settings.quicklooks, channels = channels)
        self.telecover_settings_window = TelecoverSettingsWindow(self.tabs, settings = self.settings.telecover, channels = channels)
        self.rayleighfit_settings_window = RayleighFitSettingsWindow(self.tabs, settings = self.settings.rayleigh_fit, channels = channels )
        self.polarizationcalibration_settings_window = PolarizationCalibrationSettingsWindow (self.tabs, settings = self.settings.polarization_calibration, channels = channels)
        
        page_index = 0
        
        # self.tabs.InsertPage ( index = 0, page = GeneralSettingsWindow(self.tabs), text = "General", select = True )
        self.tabs.InsertPage ( index = page_index, page = self.converter_settings_window, text = "Converter", select = True )
        page_index += 1
        self.tabs.InsertPage ( index = page_index, page = self.preprocessor_settings_window, text = "Preprocessor", select = True )
        page_index += 1
        self.tabs.InsertPage ( index = page_index, page = self.quicklooks_settings_window, text = "Quicklooks", select = True )
        page_index += 1
        self.tabs.InsertPage ( index = page_index, page = self.rayleighfit_settings_window, text = "Rayleigh Fit", select = True )
        page_index += 1
        self.tabs.InsertPage ( index = page_index, page = self.telecover_settings_window, text = "Telecover", select = True )
        page_index += 1
        self.tabs.InsertPage ( index = page_index, page = self.polarizationcalibration_settings_window, text = "Polarization calibration", select = True )
        page_index += 1
        
        self.tabs.SetSelection (0)
        
        sizer = wx.BoxSizer (wx.HORIZONTAL)
        sizer.Add (self.tabs, 1, wx.EXPAND)
        self.panel.SetSizer (sizer)
        
        BottomHBox=wx.BoxSizer(wx.HORIZONTAL)
        
        self.applyBtn = wx.Button(self, label="OK")
        self.applyBtn.Bind(wx.EVT_BUTTON, self.OnApplyBtnClicked)
        self.cancelBtn = wx.Button(self, label="Cancel")
        self.cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancelBtnClicked)
        
        BottomHBox.Add(self.applyBtn, proportion=0, flag=wx.RIGHT, border=5)
        BottomHBox.Add(self.cancelBtn, proportion=0, flag=wx.LEFT, border=5)
        
        VBox=wx.BoxSizer(wx.VERTICAL)
        VBox.Add(self.panel, proportion=1, flag=wx.EXPAND|wx.ALL, border=10)
        VBox.Add(BottomHBox, flag=wx.CENTER|wx.ALL, border=10)
        
        self.SetSizer(VBox)
        
    def OnApplyBtnClicked(self, e):
        '''
        "Apply" button event handler.
        
        Keyword arguments:
        e -- the event.
        '''        
        self.EndModal(wx.ID_OK)
        
    def OnCancelBtnClicked(self, e):
        '''
        "Apply" button event handler.
        
        Keyword arguments:
        e -- the event.
        '''        
        self.settings = self.initial_settings
        self.EndModal(wx.ID_CANCEL)
        
    def GetSettings (self):
        return self.settings
        
class ConfigurationWindow(wx.Window):
    def __init__ (self, *args, **kwargs):
        wx.Window.__init__(self, *args, **kwargs)
        
        self.panel = None
        self.SetSettings (kwargs.pop("settings", ATLASConfiguration()))
        
    def SetSettings ( self, settings ):
        self.settings = settings
        self.Build()
        
    def Build (self):
        self.DestroyChildren()
        
        self.panel = wx.Panel(self)
        self.sizer = wx.GridBagSizer (vgap = 10, hgap = 10)

        self.lidar_name = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.lidar_name, label = "Lidar name" )
        self.station_id = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.station_id, label = "EARLINET DB Lidar ID" )
        
        self.sizer.Add ( self.lidar_name, pos = wx.GBPosition (0, 0), span = wx.GBSpan (1, 1), flag = wx.EXPAND )
        self.sizer.Add ( self.station_id, pos = wx.GBPosition (0, 1), span = wx.GBSpan (1, 1), flag = wx.EXPAND )
                
        self.sizer.Add ( self.location, pos = wx.GBPosition(0, 2), span = wx.GBSpan(1, 1), flag = wx.EXPAND )
        
        self.latitude = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.latitude, label = "Latitude [°]" )
        self.longitude = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.longitude, label = "Longitude [°]" )
        self.altitude = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.altitude, label = "Altitude [m]" )
        
        self.sizer.Add ( self.latitude, pos = wx.GBPosition(1, 0), flag = wx.EXPAND )
        self.sizer.Add ( self.longitude, pos = wx.GBPosition(1, 1), flag = wx.EXPAND )
        self.sizer.Add ( self.altitude, pos = wx.GBPosition(1, 2), flag = wx.EXPAND )
        
        self.azimuth_angle = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.azimuth_angle, label = "Azimuth Angle [°]" )
        self.zenith_angle = AtlasUIInputField.FromSetting ( parent = self.panel, setting = self.settings.lidar.zenith_angle, label = "Zenith Angle [°]" )
        
        self.sizer.Add ( self.azimuth_angle, pos = wx.GBPosition(2, 0), flag = wx.EXPAND )
        self.sizer.Add ( self.zenith_angle, pos = wx.GBPosition(2, 1), flag = wx.EXPAND )
        
        self.sizer.AddGrowableCol (idx = 0, proportion = 1)
        self.sizer.AddGrowableCol (idx = 1, proportion = 1)
        self.sizer.AddGrowableCol (idx = 2, proportion = 1)
        
        self.panel.SetSizer (self.sizer)
        
        self.window_sizer = wx.FlexGridSizer (cols = 1)
        self.window_sizer.Add (self.panel, 1, flag = wx.ALL | wx.EXPAND, border = 10)
        
        self.add_channel_btn = wx.Button ( self, label = "Insert new channel" )
        self.add_channel_btn.Bind(wx.EVT_BUTTON, self.OnInsertNewChannel)
        
        self.window_sizer.Add ( self.add_channel_btn )
        
        self.tabs = wx.Notebook(self, style = wx.CHB_TOP)
        # self.tabs.InsertPage ( index = 0, page = ChannelWindow(self.tabs), text = "1", select = True )
        # self.tabs.InsertPage ( index = 1, page = ConfigurationWindow(self.tabs), text = "Configuration", select = True )
        
        self.window_sizer.Add (self.tabs, 1, wx.EXPAND)
        
        self.window_sizer.AddGrowableCol (0, 1)
        self.window_sizer.AddGrowableRow (2, 1)

        self.SetSizerAndFit (self.window_sizer)
        
        self.Refresh()
        
    def RefreshChannelList( self ):
        self.settings.channels = [ self.tabs.GetPage(index).Settings() for index in range(self.tabs.GetPageCount()) ]
        
    def InsertChannel ( self, channel_settings = None ):
        if channel_settings is None:
            channel_settings = ATLASChannelSettings()
            
        window = ChannelWindow ( parent = self.tabs, settings = channel_settings )
        window.Bind ( EVT_CHANNEL_DELETE, self.OnDeleteChannel )
        window.Bind ( EVT_CHANNEL_CHANGE, self.OnChangeChannel )
        
        i = self.tabs.GetPageCount()
        self.tabs.InsertPage ( index = self.tabs.PageCount, page = window, text = f"{i+1}: {window.GetChannelName()}", select = True )
        
    def InitChannelTabs (self):
        for settings in self.settings.channels:
            self.InsertChannel ( channel_settings = settings )
            
        self.RefreshChannelList()
        
    def OnInsertNewChannel (self, e):
        self.InsertChannel ()
        self.RefreshChannelList()
        
    def OnDeleteChannel ( self, e ):
        index = self.tabs.FindPage ( e.GetWindow() )
        self.settings.channels.remove ( e.GetWindow().Settings() )
        self.tabs.DeletePage ( index )
        
        for i in range (self.tabs.GetPageCount()):
            self.tabs.SetPageText ( i, f'{i+1}: {self.tabs.GetPage(i).GetChannelName()}' )
            
        # self.RefreshChannelList()
            
    def OnChangeChannel ( self, e ):
        i = self.tabs.FindPage ( e.GetWindow() )
        self.tabs.SetPageText ( i, f'{i+1}: {self.tabs.GetPage(i).GetChannelName()}' )
        
    def GetChannels ( self ):
        return self.settings.channels

        
class DataWindow(wx.Window):
    def OnParentFolderChanged (self, event):
        print (event.GetString())
        
    def __init__ (self, *args, **kwargs):
        wx.Window.__init__(self, *args, **kwargs)
        
        sizer = wx.BoxSizer (wx.VERTICAL)
        
        sizer.AddSpacer (15)
        
        self.parent_folder = FolderField(parent = self, label = "Parent folder", optional = True, selected = False)
        sizer.Add ( self.parent_folder, 0, wx.EXPAND)
        self.parent_folder.Bind (wx.EVT_TEXT, self.OnParentFolderChanged)
        
        sizer.AddSpacer (15)
        
        self.dark_folder = FolderField(parent = self, label = "Dark folder")
        sizer.Add ( self.dark_folder, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.rayleigh_folder = FolderField(parent = self, label = "Rayleigh folder")
        sizer.Add ( self.rayleigh_folder, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.telecover_sectors_folder = FolderField(parent = self, label = "Telecover sectors folder")
        sizer.Add ( self.telecover_sectors_folder, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.telecover_rings_folder = FolderField(parent = self, label = "Telecover rings folder")
        sizer.Add ( self.telecover_rings_folder, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.pcb_cal_p45 = FolderField(parent = self, label = "Polarization calibration +45 folder")
        sizer.Add ( self.pcb_cal_p45, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.pcb_cal_m45 = FolderField(parent = self, label = "Polarization calibration -45 folder")
        sizer.Add ( self.pcb_cal_m45, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.pcb_cal_stc = FolderField(parent = self, label = "Polarization calibration single calibrator folder")
        sizer.Add ( self.pcb_cal_stc, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.radiosonde_folder = FolderField(parent = self, label = "Radiosonde folder")
        sizer.Add ( self.radiosonde_folder, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.converter_out_folder = FolderField(parent = self, label = "Converter output folder")
        sizer.Add ( self.converter_out_folder, 0, wx.EXPAND)
        
        sizer.AddSpacer (15)
        
        self.preprocessor_out_folder = FolderField(parent = self, label = "Preprocessor output folder")
        sizer.Add ( self.preprocessor_out_folder, 0, wx.EXPAND)
        
        self.SetSizer ( sizer )

class ATLASCommandFrame (wx.Window):
    def __init__ (self, *args, **kwargs):
        wx.Window.__init__(self, *args, **kwargs)
        
        self.tabs = wx.Notebook(self, style = wx.NB_TOP)
        # self.tabs.InsertPage ( index = 0, page = DataWindow(self.tabs), text = "Data", select = True )
        
        self.configuration_window = ConfigurationWindow(self.tabs)
        self.tabs.InsertPage ( index = 0, page = self.configuration_window, text = "Configuration", select = True )
        
        self.tabs.SetSelection (0)
        
        sizer = wx.BoxSizer (wx.HORIZONTAL)
        sizer.Add (self.tabs, 1, wx.EXPAND)
        self.SetSizer (sizer)
    
class ATLASTerminalFrame (wx.Window):
    def __init__ (self, *args, **kwargs):
        wx.Window.__init__(self, *args, **kwargs)
        self.SetBackgroundColour (wx.Colour (0, 0, 0))

class MainFrame (wx.Frame):
    def __init__ (self):
        wx.Frame.__init__(self, None, title=APP_TITLE)
        
        self.SetSize (800, 700)
        
        self.settings = ATLASSettings ()
        self.command_frame = ATLASCommandFrame ( self, size = wx.Size(100, 400) )
        
        bs = wx.BoxSizer (wx.HORIZONTAL)
        bs.Add ( self.command_frame, 1, wx.EXPAND )
        # bs.Add ( ATLASTerminalFrame ( self, size = wx.Size(300, 400) ), 0, wx.EXPAND )
        
        self.toolbar = self.CreateToolBar (name = "Toolbar", style = wx.TB_HORZ_TEXT)
        
        self.loadTool = self.toolbar.AddTool ( wx.ID_ANY, 'Save configuration', wx.ArtProvider.GetBitmap (wx.ART_FILE_SAVE) )
        self.Bind ( wx.EVT_TOOL, self.OnSaveConfigurationBtnClicked, id = self.loadTool.Id )
        
        self.loadTool = self.toolbar.AddTool ( wx.ID_ANY, 'Load configuration', wx.ArtProvider.GetBitmap (wx.ART_FILE_OPEN) )
        self.Bind ( wx.EVT_TOOL, self.OnLoadConfigurationBtnClicked, id = self.loadTool.Id )
        
        self.settingsTool = self.toolbar.AddTool ( wx.ID_ANY, 'Settings', wx.ArtProvider.GetBitmap( wx.ART_FIND_AND_REPLACE ) )
        self.Bind (wx.EVT_TOOL, self.OnSettingsButtonClicked, id=self.settingsTool.Id)
        
        self.toolbar.AddStretchableSpace()
        
        self.writeConfigTool = self.toolbar.AddTool ( wx.ID_ANY, 'Write config_file.ini', wx.ArtProvider.GetBitmap (wx.ART_FILE_SAVE_AS) )
        self.Bind (wx.EVT_TOOL, self.OnConfigurationFileBtnClicked, id = self.writeConfigTool.Id)
        
        self.writeSettingsTool = self.toolbar.AddTool ( wx.ID_ANY, 'Write settings_file.ini', wx.ArtProvider.GetBitmap (wx.ART_FILE_SAVE_AS) )
        self.Bind (wx.EVT_TOOL, self.OnSettingsFileBtnClicked, id = self.writeSettingsTool.Id)
        
        self.toolbar.Realize()
        
        self.SetSizer (bs)
        
        self.Bind (wx.EVT_CLOSE, self.OnClose)
        
    def OnLoadConfigurationBtnClicked ( self, e ):
        # Code to load configuration goes here!
        with wx.FileDialog(
            self, "Load progress",
            wildcard = "Pickle files (*.pkl)|*.pkl",
            style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        ) as file_dialog:
        
            if file_dialog.ShowModal() == wx.ID_CANCEL:
                return # User canceled
                
            path = file_dialog.GetPath ()
            
            try:
                with open( path, 'rb' ) as file:
                    data = pickle.load (file)
                    
                self.command_frame.configuration_window.settings = data["configuration"]
                self.settings = data["settings"]
                
                # self.command_frame.configuration_window.DestroyChildren()
                self.command_frame.configuration_window.Build()
                self.command_frame.configuration_window.InitChannelTabs()
                self.command_frame.Layout()
                
            except Exception as e:
                print (e)
                wx.MessageDialog(self, "Could not load pickle file", "Error", style = wx.ICON_ERROR | wx.CENTRE).ShowModal()
                
    def OnSaveConfigurationBtnClicked (self, e):
        if self.SavePickle():
            e.Skip()
        
    def OnClose ( self, e ):
        if self.SavePickle():
            e.Skip()
        
    def SavePickle(self):
        with wx.FileDialog(
            self, "Save progress",
            wildcard="Pickle files (*.pkl)|*.pkl",
            defaultFile = "atlas.pkl",
            style = wx.FD_SAVE
        ) as file_dialog:
        
            if file_dialog.ShowModal() == wx.ID_CANCEL:
                return True # User canceled
                
            try:
                with open ( file_dialog.GetPath (), 'wb' ) as file:
                    data = {
                        "configuration": self.command_frame.configuration_window.settings,
                        "settings": self.settings
                    }
                    
                    pickle.dump( data, file )
                return True
            except Exception as e:
                # Could not save file!
                print (e)
                wx.MessageDialog(self, "Could not save progress", "Error", style = wx.ICON_ERROR | wx.CENTRE).ShowModal()
                
        return False
                
    def OnConfigurationFileBtnClicked (self, e):
        with wx.FileDialog(
            self, "Save configuration file",
            wildcard="INI files (*.ini)|*.ini",
            defaultFile = "config_file.ini",
            style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as file_dialog:
            if file_dialog.ShowModal() == wx.ID_CANCEL:
                return # User canceled
                
            try:
                success = self.command_frame.configuration_window.settings.save_to_file ( path = file_dialog.GetPath() )
            except Exception as e:
                print (e)
                success = False
            
        if success:
            wx.MessageDialog(self, "Configuration file saved successfully!", "File saved!", style = wx.ICON_INFORMATION | wx.CENTRE).ShowModal()
        else:
            wx.MessageDialog(self, "Configuration file could not be saved!", "File could not be saved!", style = wx.ICON_ERROR | wx.CENTRE).ShowModal()
        
    def OnSettingsFileBtnClicked (self, e):
        with wx.FileDialog(
            self, "Save settings file",
            wildcard="INI files (*.ini)|*.ini",
            defaultFile = "settings_file.ini",
            style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as file_dialog:
            if file_dialog.ShowModal() == wx.ID_CANCEL:
                return # User canceled
                
            try:
                success = self.settings.save_to_file ( path = file_dialog.GetPath() )
            except Exception:
                success = False
            
        if success:
            wx.MessageDialog(self, "Configuration file saved successfully!", "File saved!", style = wx.ICON_INFORMATION | wx.CENTRE).ShowModal()
        else:
            wx.MessageDialog(self, "Configuration file could not be saved!", "File could not be saved!", style = wx.ICON_ERROR | wx.CENTRE).ShowModal()
        
    def OnSettingsButtonClicked(self, e):
        temp_settings = copy.deepcopy ( self.settings )
        
        channels = { name : name for name in [ channel.GetChannelName() for channel in self.command_frame.configuration_window.GetChannels() ] }
        
        dialog = SettingsDialog(self, title="Settings", settings=temp_settings, channels = channels)
        
        if dialog.ShowModal() == wx.ID_OK:
            self.settings = dialog.GetSettings()
            
        dialog.Destroy()
        
class ATLASGui(wx.App):
    def OnInit(self):
        self.frame = MainFrame()
        self.SetTopWindow (self.frame)
        self.frame.Show()
        
        return True

ATLASGui().MainLoop()