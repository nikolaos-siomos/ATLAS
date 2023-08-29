from math import ceil

class ATLASOption:
    def __init__ (
        self,
        value = None,
        selected = False,
        default = None,
        mandatory = False
    ):
        self.SetValue ( value )
        self.SetMandatory ( mandatory )
        
        self.default = default
        
        if selected or mandatory:
            self.Select()
            
            if not value:
                self.SetValue ( self.default )
        else:
            self.Deselect()
            
        if not self.value:
            self.SetValue ( self.default )
            
    def __call__ ( self, *args ):
        if len(args) > 0:
            self.SetValue ( args[0] )
        else:
            return self.GetValue ()
            
    def __getattribute__ ( self, name ):
        if name == 'optional' or name == 'pretty':
            return self.PrettyValue()
        elif name == 'required':
            return self.MandatoryValue()
            
        return object.__getattribute__ ( self, name )
        
    def SetValue ( self, value ):
        self.value = value
        
    def GetValue ( self ):
        if self.Selected():
            return self.value
            
        return None
        
    def SetDefault ( self, value ):
        self.default = value
        
    def PrettyValue ( self ):
        value = self.GetValue()
        
        if value is None:
            return ""
            
        return str(value)
        
    def MandatoryValue ( self ):
        value = self.GetValue()
        
        if value is None:
            return "_"
        
        return self.PrettyValue()
        
    def Selected ( self ):
        return self.selected
        
    def Select ( self, select = True ):
        self.selected = select
        
    def Deselect ( self ):
        self.selected = False
        
    def Mandatory ( self ):
        return self.mandatory
        
    def SetMandatory ( self, mandatory ):
        self.mandatory = True if mandatory is True else False
        
        
class ATLASSelectOption (ATLASOption):
    pass
        
class ATLASTextOption (ATLASOption):
    pass
    
class ATLASFolderOption (ATLASOption):
    pass
    
class ATLASChannelPairsSetting(ATLASOption):
    pass
        
class ATLASIntegerOption (ATLASOption):
    def __init__ (self, *args, **kwargs):
        self.min = kwargs.pop ("min", None)
        self.max = kwargs.pop ("max", None)
        self.inc = kwargs.pop ("increment", 1)
        ATLASOption.__init__ ( self, *args, **kwargs )
        
    def Min ( self ):
        return self.min
        
    def SetMin ( self, min ):
        self.min = min
        
    def Max ( self ):
        return self.max
        
    def SetMax ( self, max ):
        self.max = max
        
    def Inc ( self ):
        return self.inc
        
    def SetInc ( self, min ):
        self.inc = inc
        
    def PrettyValue ( self ):
        if self.GetValue() is None:
            return ""
            
        return str(self.GetValue())
        
class ATLASIntegerRangeOption (ATLASIntegerOption):
    def PrettyValue ( self ):
        try:
            value = self.GetValue()
            return f"{value[0]}, {value[1]}"
        except Exception as e:
            return ""
        
class ATLASDoubleOption (ATLASOption):
    def __init__ (self, *args, **kwargs):
        self.min = kwargs.pop ("min", None)
        self.max = kwargs.pop ("max", None)
        self.inc = kwargs.pop ("increment", .1)
        ATLASOption.__init__ ( self, *args, **kwargs )
        
    def Min ( self ):
        return self.min
        
    def SetMin ( self, min ):
        self.min = min
        
    def Max ( self ):
        return self.max
        
    def SetMax ( self, max ):
        self.max = max
        
    def Inc ( self ):
        return self.inc
        
    def SetInc ( self, min ):
        self.inc = inc
        
    def PrettyValue ( self ):
        value = self.GetValue()
        
        if value is None:
            return ""

        return f"{value:.2f}"
        
class ATLASDoubleRangeOption (ATLASDoubleOption):
    def PrettyValue ( self ):
        try:
            value = self.GetValue()
            return f"{value[0]:f}, {value[1]:f}"
        except Exception as e:
            return ""
        
class ATLASChoiceOption (ATLASOption):
    def __init__ (self, *args, **kwargs):
        self.choices = kwargs.pop ("choices", None)
        ATLASOption.__init__ ( self, *args, **kwargs )
    
    def Choices ( self ):
        return self.choices
        
    def SetChoices ( self, choices ):
        if choices is None or not isinstance(choices, dict):
            self.choices = {}
        else:
            self.choices = choices
            
        self.SetValue ( None )
        
    def SetValue ( self, value ):
        if value in self.choices.keys():
            self.value = value
        else:
            self.value = None
            
class ATLASMultipleChoiceOption (ATLASOption):
    def __init__ (self, *args, **kwargs):
        self.choices = kwargs.pop ("choices", None)
        ATLASOption.__init__ ( self, *args, **kwargs )
    
    def Choices ( self ):
        return self.choices
        
    def SetChoices ( self, choices ):
        if choices is None or not isinstance(choices, dict):
            self.choices = {}
        else:
            self.choices = choices
            
        self.SetValue ( [] )
        
    def SetValue ( self, value ):
        if not self.choices:
            # Accepted values list is dynamically generated
            # so we can't compare to a known list beforehand.
            self.value = value
            return
            
        self.value = []
        if value is not None and isinstance(value, list):
            for selected_value in value:
                if selected_value in self.choices.keys():
                    self.value.append( selected_value )
                    
    def PrettyValue ( self ):
        try:
            return ', '.join(self.GetValue())
        except Exception as e:
            return ""
                    
class BaseSettings:
    def __getitem__ ( self, key ):
        try:
            return getattr( self, key )()
        except AttributeError:
            raise KeyError ( "No such field!" )
            
    def __setitem__ ( self, key, item ):
        try:
            getattr( self, key )
            setattr ( self, key, item )
        except AttributeError:
            raise KeyError ( "No such field!" )
            
    def __contains__ ( self, key ):
        try:
            getattr ( self, key )
            return True
        except AttributeError:
            return False

class ATLASConverterSettings(BaseSettings):
    def __init__ (self):
        self.debug = ATLASSelectOption ( default = False )
        
        self.trim_overflows = ATLASChoiceOption (
            default = '0',
            mandatory = True,
            choices = {
                "0": "No action",
                "1": "Screen out overflows",
                "2": "Interpolate overflows",
                "3": "Include overflows (debug)"
            } 
        )
        
        self.files_per_sector = ATLASIntegerOption ( default = None )
        self.files_per_ring = ATLASIntegerOption ( default = None )
        
        self.rsonde_skip_header = ATLASIntegerOption ( default = 1, min = 0, max = 100 )
        self.rsonde_skip_footer = ATLASIntegerOption ( default = 0, min = 0, max = 100 )
        self.rsonde_delimiter = ATLASChoiceOption (
            default = 'S',
            choices = {
                "S": "Space",
                "C": "Comma"
            }
        )
        
        self.rsonde_h_column = ATLASIntegerOption ( default = 2, min = 1, max = 40, mandatory = True )
        self.rsonde_p_column = ATLASIntegerOption ( default = 1, min = 1, max = 40, mandatory = True )
        self.rsonde_t_column = ATLASIntegerOption ( default = 3, min = 1, max = 40, mandatory = True )
        self.rsonde_rh_column = ATLASIntegerOption ( default = None, min = 0, max = 40, mandatory = True )
        
        self.rsonde_h_unit = ATLASChoiceOption (
            default = 'm_asl',
            mandatory = True,
            choices = {
                "m_asl": "m_asl",
                "km_asl": "km_asl",
                "m_agl": "m_agl",
                "km_agl": "km_agl"
            }
        )
        
        self.rsonde_p_unit = ATLASChoiceOption (
            default = 'hPa',
            mandatory = True,
            choices = {
                "hPa": "hPa",
                "Pa": "Pa",
                "atm": "atm"
            }
        )
        
        self.rsonde_t_unit = ATLASChoiceOption (
            default = 'K',
            mandatory = True,
            choices = {
                "K": "K",
                "C": "°C",
                "Cx10": "°C x10"
            }
        )
        
        self.rsonde_rh_unit = ATLASChoiceOption (
            default = 'percent',
            mandatory = True,
            choices = {
                "percent": "Percent",
                "fraction": "Fraction"
            }
        )
        
        self.rsonde_rh = ATLASSelectOption ( default = False )
        self.rsonde_latitude = ATLASDoubleOption ( default = None, increment = .01, min = -90, max = 90 )
        self.rsonde_longitude = ATLASDoubleOption ( default = None, increment = .01, min = -90, max = 90 )
        self.rsonde_altitude = ATLASDoubleOption ( default = None, increment = .1, min = 0, max = 8000 )
        
        self.rsonde_station_name = ATLASTextOption ( default = None )
        self.rsonde_wmo_number = ATLASTextOption ( default = None )
        self.rsonde_wban_number = ATLASTextOption ( default = None )
        
class ATLASPreprocessorSettings(BaseSettings):
    def __init__ (self):
        self.vertical_trimming =  ATLASSelectOption ( default = False )
        self.vertical_limit =  ATLASDoubleOption ( default = 20, min = 0, max = 50 )
        self.channels =  ATLASMultipleChoiceOption ( default = [] )
        
        self.exclude_telescope_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "n": "Near range",
                "f": "Far range",
                "x": "Unspecified"
            }
        )
        self.exclude_scattering_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "p": "Co-polar linear analyzer",
                "c": "Cross-polar linear analyzer",
                "t": "Total (no depolarization)",
                "o": "Co-polar circular analyzer",
                "x": "Cross-polar circular analyzer",
                "v": "Vibrational Raman",
                "r": "Rotational Raman",
                "a": "Cabannes",
                "f": "Fluorescence"
            }
        )
        self.exclude_acquisition_mode =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "a": "Analog",
                "p": "Photon counting"
            }
        )
        self.exclude_channel_subtype =  ATLASMultipleChoiceOption (
            default = [],
                choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS",
                "n": "N2 Ramal line",
                "o": "O2 Ramal line",
                "w": "H2O Ramal line",
                "c": "CH4 Ramal line",
                "h": "High Rotational Raman",
                "l": "Low Rotational Raman",
                "a": "Mie (aerosol) HSRL signal",
                "m": "Molecular HSRL signal",
                "b": "Broadband Fluorescence",
                "s": "Spectral Fluorescence",
                "x": "Unspecified"
            }
        )
        
class ATLASQuicklooksSettings(BaseSettings):
    def __init__ (self):
        self.use_log_scale = ATLASSelectOption ( default = False )
        self.use_range =  ATLASSelectOption ( default = True )
        
        self.channels = ATLASMultipleChoiceOption ( default = [] )
        self.exclude_telescope_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "n": "Near range",
                "f": "Far range",
                "x": "Unspecified"
            }
        )
        
        self.exclude_scattering_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "p": "Co-polar linear analyzer",
                "c": "Cross-polar linear analyzer",
                "t": "Total (no depolarization)",
                "o": "Co-polar circular analyzer",
                "x": "Cross-polar circular analyzer",
                "v": "Vibrational Raman",
                "r": "Rotational Raman",
                "a": "Cabannes",
                "f": "Fluorescence"
            }
        )
        
        self.exclude_acquisition_mode =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "a": "Analog",
                "p": "Photon counting"
            }
        )
        
        self.exclude_channel_subtype =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS",
                "n": "N2 Ramal line",
                "o": "O2 Ramal line",
                "w": "H2O Ramal line",
                "c": "CH4 Ramal line",
                "h": "High Rotational Raman",
                "l": "Low Rotational Raman",
                "a": "Mie (aerosol) HSRL signal",
                "m": "Molecular HSRL signal",
                "b": "Broadband Fluorescence",
                "s": "Spectral Fluorescence",
                "x": "Unspecified"
            }
        )
        
        self.x_lims =  ATLASIntegerRangeOption ( default = None, min = 1 )
        self.x_tick =  ATLASIntegerOption ( default = None, min = 0 )
        
        self.upper_t_lim = ATLASTextOption ( default = '' )
        self.lower_t_lim = ATLASTextOption ( default = '' )
        
        # self.t_lims = ATLASOption ( default = None )
        self.t_tick =  ATLASTextOption ( default = None )
        self.y_lims = ATLASDoubleRangeOption ( default = (0, 14), min = 1 )
        self.y_tick =  ATLASDoubleOption ( default =  1 )
        self.z_lims =  ATLASDoubleRangeOption ( default = (0, 1) )
        # self.z_min_zone =  ATLASDoubleRangeOption ( default = (2., 10.) )
        self.z_max_zone =  ATLASDoubleRangeOption ( default = (.1, 2) )
        self.smooth =  ATLASSelectOption ( default = False )
        self.smoothing_range =  ATLASDoubleRangeOption ( default = (.05, 14) )
        self.smoothing_window =  ATLASDoubleRangeOption ( default = (50, 500) )
        self.smoothing_exponential =  ATLASSelectOption ( default = False )
        self.dpi =  ATLASIntegerOption ( default = 300, min = 100, max = 500 )
        self.color_reduction =  ATLASSelectOption ( default = False )
        
class ATLASRayleighFitSettings(BaseSettings):
    def __init__ (self):
        self.channels = ATLASMultipleChoiceOption ( default = [] )
        self.exclude_telescope_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "n": "Near range",
                "f": "Far range",
                "x": "Unspecified"
            }
        )
        
        self.exclude_scattering_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "p": "Co-polar linear analyzer",
                "c": "Cross-polar linear analyzer",
                "t": "Total (no depolarization)",
                "o": "Co-polar circular analyzer",
                "x": "Cross-polar circular analyzer",
                "v": "Vibrational Raman",
                "r": "Rotational Raman",
                "a": "Cabannes",
                "f": "Fluorescence"
            }
        )
        
        self.exclude_acquisition_mode =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "a": "Analog",
                "p": "Photon counting"
            }
        )
        
        self.exclude_channel_subtype =  ATLASMultipleChoiceOption (
            default = [],
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
        )
        
        self.x_lims =  ATLASDoubleRangeOption ( default = (0, 20), min = 1, increment = .1 )
        self.x_tick =  ATLASDoubleOption ( default = 1, min = 0, increment = .1 )
        
        self.y_lims = ATLASDoubleRangeOption ( default = None, min = 1, increment = .1 )

        self.smooth =  ATLASSelectOption ( default = False )
        self.smoothing_range =  ATLASDoubleRangeOption ( default = (.05, 15) )
        self.smoothing_window =  ATLASDoubleRangeOption ( default = (100, 100) )
        self.smoothing_exponential =  ATLASSelectOption ( default = False )
        self.dpi =  ATLASIntegerOption ( default = 300, min = 100, max = 500 )
        self.color_reduction =  ATLASSelectOption ( default = False )
        
        self.use_range =  ATLASSelectOption ( default = True )
        self.use_lin_scale =  ATLASSelectOption ( default = False )
        
        self.normalization_region =  ATLASDoubleRangeOption ( default = (8.5, 9.5) )
        
class ATLASTelecoverSettings(BaseSettings):
    def __init__ (self):
        self.channels = ATLASMultipleChoiceOption ( default = [] )
        self.exclude_telescope_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "n": "Near range",
                "f": "Far range",
                "x": "Unspecified"
            }
        )
        
        self.exclude_scattering_type =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "p": "Co-polar linear analyzer",
                "c": "Cross-polar linear analyzer",
                "t": "Total (no depolarization)",
                "o": "Co-polar circular analyzer",
                "x": "Cross-polar circular analyzer",
                "v": "Vibrational Raman",
                "r": "Rotational Raman",
                "a": "Cabannes",
                "f": "Fluorescence"
            }
        )
        
        self.exclude_acquisition_mode =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "a": "Analog",
                "p": "Photon counting"
            }
        )
        
        self.exclude_channel_subtype =  ATLASMultipleChoiceOption (
            default = [],
            choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS",
                "n": "N2 Ramal line",
                "o": "O2 Ramal line",
                "w": "H2O Ramal line",
                "c": "CH4 Ramal line",
                "h": "High Rotational Raman",
                "l": "Low Rotational Raman",
                "a": "Mie (aerosol) HSRL signal",
                "m": "Molecular HSRL signal",
                "b": "Broadband Fluorescence",
                "s": "Spectral Fluorescence",
                "x": "Unspecified"
            }
        )
        
        self.x_lims =  ATLASDoubleRangeOption ( default = (0, 2.4), min = 1 )
        self.x_tick =  ATLASDoubleOption ( default = 0.2, min = 0, increment = .1 )
        
        self.y_lims = ATLASDoubleRangeOption ( default = None, min = 1 )

        self.smooth =  ATLASSelectOption ( default = True )
        self.smoothing_range =  ATLASDoubleRangeOption ( default = (.5, 2.5) )
        self.smoothing_window =  ATLASDoubleRangeOption ( default = (100., 100.) )
        self.smoothing_exponential =  ATLASSelectOption ( default = False )
        self.dpi =  ATLASIntegerOption ( default = 300, min = 100, max = 500 )
        self.color_reduction =  ATLASSelectOption ( default = False )
        
        self.use_range =  ATLASSelectOption ( default = True )
        self.use_lin_scale =  ATLASSelectOption ( default = True )
        self.use_non_range_cor =  ATLASSelectOption ( default = False )
        
        self.use_last = ATLASSelectOption ( default = True )
        
        self.normalization_region =  ATLASDoubleRangeOption ( default = (1.8, 2.2) )
        
class ChannelPair:
    def __init__ (
        self, reflected_channel = None, transmitted_channel = None,
        k = 1., g_r = 1., g_t = 1., h_r = None, h_t = None, rt_t_ratio = 1.
    ):
        if h_r is None:
            if reflected_channel[5] in [ 'p', 'o' ]:
                h_r = -1
            else:
                h_r = 1
                
        if h_t is None:
            if transmitted_channel[5] in [ 'p', 'o' ]:
                h_t = -1
            else:
                h_t = 1
        
        self.ReflectedChannel(reflected_channel)
        self.TransmittedChannel(transmitted_channel)
        
        self.K(k)
        self.G_R(g_r)
        self.G_T(g_t)
        self.H_R(h_r)
        self.H_T(h_t)
        self.RT_T_Ratio(rt_t_ratio)
        
    # def ReflectedChannel(self):
        # return self.reflected_channel
        
    def ReflectedChannel(self, reflected_channel = None):
        if reflected_channel is None:
            return self.reflected_channel
            
        self.reflected_channel = reflected_channel
        
    # def TransmittedChannel(self):
        # return self.transmitted_channel
        
    def TransmittedChannel(self, transmitted_channel = None):
        if transmitted_channel is None:
            return self.transmitted_channel
            
        self.transmitted_channel = transmitted_channel
        
    # def K(self):
        # return self.k
        
    def K(self, k = None):
        if k is None:
            return self.k
            
        self.k = k
        
    # def G_R(self):
        # return self.g_r
        
    def G_R(self, g_r = None):
        if g_r is None:
            return self.g_r
            
        self.g_r = g_r
        
    # def G_T(self):
        # return self.g_t
        
    def G_T(self, g_t = None):
        if g_t is None:
            return self.g_t
            
        self.g_t = g_t
        
    # def H_R(self):
        # return self.h_r
        
    def H_R(self, h_r = None):
        if h_r is None:
            return self.h_r
            
        self.h_r = h_r
        
    # def H_T(self):
        # return self.h_t
        
    def H_T(self, h_t = None):
        if h_t is None:
            return self.h_t
            
        self.h_t = h_t
        
    # def RT_T_Ratio (self):
        # return self.rt_t_ratio
        
    def RT_T_Ratio (self, rt_t_ratio = None):
        if rt_t_ratio is None:
            return self.rt_t_ratio
            
        self.rt_t_ratio = rt_t_ratio
        
class ATLASPolarizationCalibrationSettings(BaseSettings):
    def __init__ (self):
        self.channel_pairs = ATLASChannelPairsSetting( default = [], mandatory = True )

        # self.calibration_height = ATLASDoubleOption ( default = 3., min = 0, max = 20 )
        # self.half_calibration_window = ATLASIntegerOption ( default = 500, min = 0, max = 5000 )
        
        self.calibration_region = ATLASDoubleRangeOption ( default = (2., 4.) )
        self.rayleigh_region = ATLASDoubleRangeOption ( default = (8.5, 9.5) )
        
        # self.rayleigh_height = ATLASDoubleOption ( default = 9., min = 0, max = 20 )
        # self.half_rayleigh_window = ATLASIntegerOption ( default = 500, min = 0, max = 5000 )
        
        self.x_lims_calibration = ATLASIntegerRangeOption ( default = (0, 5), min = 0 )
        self.x_lims_rayleigh = ATLASIntegerRangeOption ( default = (0, 10), min = 0 )
        
        self.x_tick_calibration = ATLASDoubleOption ( default = .5, min = 0 )
        self.x_tick_rayleigh = ATLASDoubleOption ( default = 1., min = 0 )
        
        self.y_lims_calibration = ATLASDoubleRangeOption ( default = None )
        self.y_lims_rayleigh = ATLASDoubleRangeOption ( default = None )
        
        self.use_range = ATLASSelectOption ( default = True )
        
        # self.smooth_calibration =  ATLASSelectOption ( default = True )
        # self.smoothing_range_calibration =  ATLASDoubleRangeOption ( default = (1., 15) )
        # self.smoothing_window_calibration =  ATLASDoubleRangeOption ( default = (100., 100.) )
        # self.smoothing_exponential_calibration =  ATLASSelectOption ( default = False )
        
        # self.smooth_rayleigh =  ATLASSelectOption ( default = True )
        # self.smoothing_range_rayleigh =  ATLASDoubleRangeOption ( default = (1, 14) )
        # self.smoothing_window_rayleigh =  ATLASDoubleRangeOption ( default = (100., 100.) )
        # self.smoothing_exponential_rayleigh =  ATLASSelectOption ( default = False )
        
        self.smooth =  ATLASSelectOption ( default = True )
        self.smoothing_range =  ATLASDoubleRangeOption ( default = (1, 14) )
        self.smoothing_window =  ATLASDoubleRangeOption ( default = (100., 100.) )
        self.smoothing_exponential =  ATLASSelectOption ( default = False )
        
        self.dpi =  ATLASIntegerOption ( default = 300, min = 100, max = 500 )
        self.color_reduction =  ATLASSelectOption ( default = True )
        
    def SetChannelPairs (self, channel_pairs):
        self.channel_pairs = channel_pairs
        
    def AddChannelPair (self, channel_pair):
        already_exists = False
        
        try:
            for existing_pair in self.channel_pairs:
                if existing_pair.ReflectedChannel() == channel_pair.ReflectedChannel() and \
                   existing_pair.TransmittedChannel() == channel_pair.TransmittedChannel():
                   
                    already_exists = True
                    break
        except Exception:
            already_exists = True

        if not already_exists:
            self.channel_pairs.append ( channel_pair )
            
    def RemoveChannelPair(self, channel_pair):
        remove_index = -1
        
        try:
            for index, existing_pair in enumerate(self.channel_pairs):
                if existing_pair.ReflectedChannel() == channel_pair.ReflectedChannel() and \
                   existing_pair.TransmittedChannel() == channel_pair.TransmittedChannel():
                   
                    remove_index = index
                    break
        except Exception:
            remove_index = -1
            
        if remove_index != -1:
            self.channel_pairs.pop (remove_index)

class ATLASSettings:
    def __init__ (self):
        self.converter = ATLASConverterSettings()
        self.preprocessor = ATLASPreprocessorSettings()
        self.quicklooks = ATLASQuicklooksSettings()
        self.rayleigh_fit = ATLASRayleighFitSettings()
        self.telecover = ATLASTelecoverSettings()
        self.polarization_calibration = ATLASPolarizationCalibrationSettings()
        
    def save_to_file ( self, path = "settings.ini" ):
        try:
            from mako.template import Template
            import os
            
            TEMPLATE_FILE_NAME = 'settings.ini.tpl'
            
            cwd = os.path.dirname(os.path.realpath(__file__) )
            
            with open ( os.path.join ( cwd, TEMPLATE_FILE_NAME ), 'r' ) as file:
                result = Template(file.read()).render (
                    converter = self.converter,
                    preprocessor = self.preprocessor,
                    quicklooks = self.quicklooks,
                    rayleigh_fit = self.rayleigh_fit,
                    telecover = self.telecover,
                    polarization_calibration = self.polarization_calibration
                )
                
                with open ( os.path.realpath ( path ), 'w' ) as output_file:
                    output_file.write ( result )
        except Exception as e:
            print (e)
            return False
            
        return True
        
        
class ATLASLidarSettings(BaseSettings):
    def __init__ (self):
        # Mandatory variables
        self.lidar_name = ATLASTextOption ( default = '', mandatory = True )
        self.station_id = ATLASTextOption ( default = '', mandatory = True )
        
        # Optional variables
        self.file_format = ATLASChoiceOption ( default = 'licel', choices = { 'licel': 'Licel', 'polly_xt': 'PollyXT' } )
        self.altitude = ATLASDoubleOption ( default = None, min = 0, max = 8000 )
        self.latitude = ATLASDoubleOption ( default = None, min = -90, max = 90, increment = .001 )
        self.longitude = ATLASDoubleOption ( default = None, min = -90, max = 90, increment = .001 )
        self.zenith_angle = ATLASDoubleOption ( default = None, min = 0, max = 90, increment = .001 )
        self.azimuth_angle = ATLASDoubleOption ( default = None, min = 0, max = 90, increment = .001 )
        
class ATLASChannelSettings(BaseSettings):
    def __init__ (self):
        # Required values
        self.telescope_type = ATLASChoiceOption (
            default = 'x',
            mandatory = True,
            choices = {
                'n': 'Near range',
                'f': 'Far range',
                'x': 'Single telescope'
            }
        )
        
        self.channel_type = ATLASChoiceOption (
            default = 'p',
            mandatory = True,
            choices = {
                "p": "Co-polar linear analyzer",
                "c": "Cross-polar linear analyzer",
                "t": "Total (no depolarization)",
                "o": "Co-polar circular analyzer",
                "x": "Cross-polar circular analyzer",
                "v": "Vibrational Raman",
                "r": "Rotational Raman",
                "a": "Cabannes",
                "f": "Fluorescence"
            }
        )
        
        self.acquisition_mode = ATLASChoiceOption (
            default = '0',
            choices = {
                '0': 'Analog',
                '1': 'Photon counting'
            },
            mandatory = True
        )
        
        self.channel_subtype = ATLASChoiceOption (
            default = 'r',
            mandatory = True,
            choices = {
                "r": "Signal Reflected from a PBS",
                "t": "Signal Transmitted through a PBS",
                "n": "N2 Ramal line",
                "o": "O2 Ramal line",
                "w": "H2O Ramal line",
                "c": "CH4 Ramal line",
                "h": "High Rotational Raman",
                "l": "Low Rotational Raman",
                "a": "Mie (aerosol) HSRL signal",
                "m": "Molecular HSRL signal",
                "b": "Broadband Fluorescence",
                "s": "Spectral Fluorescence",
                "x": "Unspecified"
            }
        )
        
        self.dead_time = ATLASDoubleOption ( default = 0, increment = .001, mandatory = True )
        self.daq_trigger_offset = ATLASIntegerOption ( default = 0, min = -4000, max = 4000, mandatory = True )
        self.recorder_channel_id = ATLASTextOption ( default = '' )
        self.laser = ATLASIntegerOption ( default = 1, min = 1, increment = 1 )
        self.scc_channel_id = ATLASIntegerOption ( default = 0 )
        
        self.emitted_wavelength = ATLASDoubleOption ( default = 355, mandatory = True, min = 1, max = 2000, increment = .01 )
        self.detected_wavelength = ATLASDoubleOption ( default = 355, mandatory = True, min = 1, max = 2000, increment = .01 )
        self.channel_bandwidth = ATLASDoubleOption ( default = .5, increment = .01 )
        self.background_low_bin = ATLASIntegerOption ( default = 0 )
        self.background_high_bin = ATLASIntegerOption ( default = 100 )
        self.background_range = ATLASIntegerRangeOption ( default = (0, 100), max = 1000 )
        
        # Optional variables
        self.bins = ATLASIntegerOption ( default = None )
        
        self.data_acquisition_range = ATLASIntegerOption ( default = None, min = 0 )
        self.analog_to_digital_resolution = ATLASIntegerOption ( default = None, min = 0 )
        self.range_resolution = ATLASDoubleOption ( default = None, min = 0 )
        self.laser_repetition_rate = ATLASIntegerOption ( default = None, min = 0 )
        
    def GetChannelName ( self ):
        am = 'a' if self.acquisition_mode() == '0' else 'p'
        return f"{self.detected_wavelength():04}{self.telescope_type()}{self.channel_type()}{am}{self.channel_subtype()}"
        
class ATLASConfiguration:
    def __init__ (self):
        self.lidar = ATLASLidarSettings()
        self.channels = []
        
    def save_to_file ( self, path = "config.ini" ):
        try:
            from mako.template import Template
            import os
            
            def print_channels_field ( channels, field ):
                try:
                    if len([ getattr(channel, field)() for channel in channels if getattr(channel, field)() is not None ]) > 0:
                        return ', '.join([ getattr(channel, field).required for channel in channels])
                    
                    return ''
                except Exception:
                    return '\n# FIELD ABOVE MIGHT HAVE ENCOUNTERED AN ERROR!'
            
            TEMPLATE_FILE_NAME = 'config.ini.tpl'
            
            cwd = os.path.dirname(os.path.realpath(__file__) )
            
            with open ( os.path.join ( cwd, TEMPLATE_FILE_NAME ), 'r' ) as file:
                result = Template(file.read()).render (
                    print_channels_field = print_channels_field,
                    lidar = self.lidar,
                    channels = self.channels
                )
                
                with open ( os.path.realpath ( path ), 'w' ) as output_file:
                    output_file.write ( result )
        except Exception as e:
            print (e)
            return False
            
        return True