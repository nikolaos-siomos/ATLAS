import wx
import copy

class SettingsWindow (wx.ScrolledWindow):
    def __init__ (self, *args, **kwargs):
        self._settings = kwargs.pop("settings", None)
        self._channels = kwargs.pop("channels", {})
        
        wx.ScrolledWindow.__init__(self, *args, **kwargs)
        self.SetScrollbars(20, 20, 55, 40)
        
    def Settings (self):
        return self._settings
        
    def Channels (self):
        return self._channels
        
    def SetSettings ( self, settings = None ):
        if settings:
            self._settings = settings
            self.Build()
        
    def Build ( self ):
        # To be implemented by subclasses
        pass