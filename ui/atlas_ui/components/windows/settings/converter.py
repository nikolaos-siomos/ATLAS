import wx

from . import SettingsWindow
from atlas_ui.components.fields import *

class ConverterSettingsWindow (SettingsWindow):
    def __init__ (self, *args, **kwargs):
        super (ConverterSettingsWindow, self).__init__ ( *args, **kwargs )
        self.Build()
        
    def Build (self):
        self.sizer = wx.FlexGridSizer ( cols = 2, hgap = 10, vgap = 10 )
        self.sizer.SetFlexibleDirection(wx.VERTICAL)
        self.sizer.AddGrowableCol (idx = 0, proportion = 1)
        self.sizer.AddGrowableCol (idx = 1, proportion = 1)
        
        self.sizer.AddStretchSpacer()
        self.sizer.AddStretchSpacer()
        
        self.debug = AtlasUIInputField.FromSetting ( parent = self, label = "Debug?", setting = self._settings.debug )
        # self.debug.SetValue ( self._settings.debug() )
        self.sizer.Add ( self.debug, 1, wx.EXPAND )

        self.trim_overflows = AtlasUIInputField.FromSetting ( parent = self, label = "Trim overflows", setting = self._settings.trim_overflows )
        self.sizer.Add ( self.trim_overflows, 1, wx.EXPAND )
        
        self.files_per_sector = AtlasUIInputField.FromSetting ( parent = self, label = "Files per sector", setting = self._settings.files_per_sector )
        self.sizer.Add ( self.files_per_sector, 1, wx.EXPAND )
        
        self.files_per_ring = AtlasUIInputField.FromSetting ( parent = self, label = "Files per ring", setting = self._settings.files_per_ring )
        self.sizer.Add ( self.files_per_ring, 1, wx.EXPAND )
        
        self.rsonde_skip_header = AtlasUIInputField.FromSetting ( parent = self, label = "Radiosonde header lines to ignore", setting = self._settings.rsonde_skip_header )
        self.sizer.Add ( self.rsonde_skip_header, 1, wx.EXPAND )
        
        self.rsonde_skip_footer = AtlasUIInputField.FromSetting ( parent = self, label = "Radiosonde footer lines to ignore", setting = self._settings.rsonde_skip_footer )
        self.sizer.Add ( self.rsonde_skip_footer, 1, wx.EXPAND )
        
        self.rsonde_delimiter = AtlasUIInputField.FromSetting ( parent = self, label = "Radiosonde files column delimiter", setting = self._settings.rsonde_delimiter )
        self.sizer.Add ( self.rsonde_delimiter, 1, wx.ALIGN_CENTER_VERTICAL | wx.EXPAND )
        
        self.sizer.AddStretchSpacer()
        
        self.rsonde_format_box = wx.StaticBox (self, label = "Radiosonde file format")
        self.rsonde_format_sizer = wx.FlexGridSizer ( cols = 3, hgap = 0, vgap = 0 )
        self.rsonde_format_sizer.AddGrowableCol ( idx = 2, proportion = 1 )
        
        self.rsonde_format_sizer.Add (wx.StaticText ( self.rsonde_format_box, label = "Height data" ), flag = wx.ALIGN_CENTER_VERTICAL | wx.TOP, border = 25)
        self.rsonde_format_h_column = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Column", setting = self._settings.rsonde_h_column )
        self.rsonde_format_sizer.Add ( self.rsonde_format_h_column, flag = wx.TOP, border = 25 )
        self.rsonde_format_h_unit = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Unit", setting = self._settings.rsonde_h_unit )
        self.rsonde_format_sizer.Add ( self.rsonde_format_h_unit, 1, wx.EXPAND | wx.TOP, border = 25 )
        
        self.rsonde_format_sizer.Add (wx.StaticText ( self.rsonde_format_box, label = "Pressure data" ), flag = wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.TOP, border = 10)
        self.rsonde_format_p_column = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Column", setting = self._settings.rsonde_p_column )
        self.rsonde_format_sizer.Add ( self.rsonde_format_p_column, flag = wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.TOP, border = 10 )
        self.rsonde_format_p_unit = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Unit", setting = self._settings.rsonde_p_unit )
        self.rsonde_format_sizer.Add ( self.rsonde_format_p_unit, 1, flag = wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM | wx.TOP, border = 10 )
        
        self.rsonde_format_sizer.Add (wx.StaticText ( self.rsonde_format_box, label = "Temperature data" ), flag = wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM, border = 10)
        self.rsonde_format_t_column = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Column", setting = self._settings.rsonde_t_column )
        self.rsonde_format_sizer.Add ( self.rsonde_format_t_column, flag = wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM, border = 10 )
        self.rsonde_format_t_unit = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Unit", setting = self._settings.rsonde_t_unit )
        self.rsonde_format_sizer.Add ( self.rsonde_format_t_unit, 1, flag = wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM, border = 10 )
        
        self.rsonde_format_rh = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Relative Humidity Data", setting = self._settings.rsonde_rh )
        self.rsonde_format_sizer.Add (self.rsonde_format_rh, flag = wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM, border = 10)
        self.rsonde_format_rh_column = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Column", setting = self._settings.rsonde_rh_column )
        self.rsonde_format_sizer.Add ( self.rsonde_format_rh_column, flag = wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM, border = 10 )
        self.rsonde_format_rh_unit = AtlasUIInputField.FromSetting ( parent = self.rsonde_format_box, label = "Unit", setting = self._settings.rsonde_rh_unit )
        self.rsonde_format_sizer.Add ( self.rsonde_format_rh_unit, 1, flag = wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | wx.BOTTOM, border = 10 )

        self.rsonde_format_box.SetSizer ( self.rsonde_format_sizer )
        self.sizer.Add ( self.rsonde_format_box, 1, wx.EXPAND )
        
        self.rsonde_station_box = wx.StaticBox ( self, label = "Radiosonde station information" )
        self.rsonde_station_sizer = wx.BoxSizer (wx.VERTICAL)
        
        self.rsonde_geodata_sizer = wx.BoxSizer ( wx.HORIZONTAL )
        
        self.rsonde_geodata_latitude = AtlasUIInputField.FromSetting ( parent = self.rsonde_station_box, label = "Latitude [°]", setting = self._settings.rsonde_latitude )
        self.rsonde_geodata_sizer.Add ( self.rsonde_geodata_latitude, 1, flag = wx.ALIGN_CENTER_VERTICAL )
        self.rsonde_geodata_longitude = AtlasUIInputField.FromSetting ( parent = self.rsonde_station_box, label = "Longitude [°]", setting = self._settings.rsonde_longitude )
        self.rsonde_geodata_sizer.Add ( self.rsonde_geodata_longitude, 1, flag = wx.ALIGN_CENTER_VERTICAL )
        self.rsonde_geodata_altitude = AtlasUIInputField.FromSetting ( parent = self.rsonde_station_box, label = "Altitude [m]", setting = self._settings.rsonde_altitude )
        self.rsonde_geodata_sizer.Add ( self.rsonde_geodata_altitude, 1, flag = wx.ALIGN_CENTER_VERTICAL )
        
        self.rsonde_station_sizer.Add ( self.rsonde_geodata_sizer, 1, flag = wx.EXPAND | wx.TOP, border = 25 )
        
        bottom_border = 30 if wx.GetOsVersion()[0] & wx.OS_UNIX else 10
        
        self.rsonde_info_station_name = AtlasUIInputField.FromSetting ( parent = self.rsonde_station_box, label = "Radiosonde station name", setting = self._settings.rsonde_station_name )
        self.rsonde_station_sizer.Add ( self.rsonde_info_station_name, flag = wx.EXPAND | wx.BOTTOM, border = bottom_border )
        
        self.rsonde_info_wmo_number = AtlasUIInputField.FromSetting ( parent = self.rsonde_station_box, label = "Radiosonde WMO number", setting = self._settings.rsonde_wmo_number )
        self.rsonde_station_sizer.Add ( self.rsonde_info_wmo_number, flag = wx.EXPAND | wx.BOTTOM, border = bottom_border )
        
        self.rsonde_info_wban_number = AtlasUIInputField.FromSetting ( parent = self.rsonde_station_box, label = "Radiosonde WBAN number", setting = self._settings.rsonde_wban_number )
        self.rsonde_station_sizer.Add ( self.rsonde_info_wban_number, flag = wx.EXPAND | wx.BOTTOM, border = bottom_border )
        
        self.rsonde_station_box.SetSizer ( self.rsonde_station_sizer )
        self.sizer.Add ( self.rsonde_station_box, 1, wx.EXPAND )
        
        self.SetSizer ( self.sizer )
