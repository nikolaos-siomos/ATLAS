[converter]
#debug: If set to True, debugging files will be generated in the  <path_to_the_converter_folder>/debug folder. Note that the path to the converter folder Defaults to: <parent_folder>/netcdf but it can be different if manually provided by the user. These included the metadata gathered from the configuration file, the licel header, and the combination of the two. Defaults to: False
debug = ${converter.debug .pretty}

#trim_overflows: This options determines how overflow values in the raw input files will be treated if detected. Choose among:
#    • 0: the algorithm will stop and provide a diagnostic error (default)
#    • 1: the files containing at least one overflow value will be screened out
#    • 2: overflows will be interpolated (use with care and only only a few bins per profile have overflow values)
#    • 3: overflows will not be excluded, use this only for debugging purposes
trim_overflows = ${converter.trim_overflows.pretty}

#files_per_sector: The number of telecover files per sector (integer). If provided, the telecover files must be placed in a single folder (this should be <path_to_parent_folder>/tlc according to the default folder structure). An automated assignment of the telecover files in different sectors will be attempted serially assuming the following temporal sequence of sectors: north – east – south – west. Note that the telecover test can have more than 1 rounds.
files_per_sector = ${converter.files_per_sector.pretty}

#files_per_ring: The number of telecover files per ring (integer). If provided, the telecover files must be placed in a single folder (this should be <path_to_parent_folder>/tlc_rin according to the default folder structure). An automated assignment of the telecover files in different sectors will be attempted serially assuming the following temporal sequence of sectors: inner - outer. Note that the telecover test can have more than 1 rounds.
files_per_ring = ${converter.files_per_ring.pretty}

#rsonde_skip_header: Radiosonde parser option. Number of lines to skip at the beginning of the radiosonde ascii file. Defaults to: 1 (1 line reserved for header info)
rsonde_skip_header = ${converter.rsonde_skip_header.pretty}

#rsonde_skip_footer: Radiosonde parser option. Number of lines to skip at the end of the radiosonde file. Defaults to: 0 (no footer assumed)
rsonde_skip_footer = ${converter.rsonde_skip_footer.pretty}

#rsonde_delimiter: Radiosonde parser option. The delimiter that separates columns in the radiosonde file choose one of:
#    • S: space
#    • C: comma
#Defaults to: S.
rsonde_delimiter = ${converter.rsonde_delimiter.pretty}

#rsonde_column_index: Radiosonde parser option. The index (integer) of the column that contains the Height, Pressure, Temperature, and Relative Humidity (optional) information. For example, rsonde_columns = 1, 3, 2, 6 means: 
#    • Height: 1st column in the radiosonde file 
#    • Pressure: 3rd column 
#    • Temperature: 2nd column
#    • Relative Humidity: 6th column
#The relative humidity column is OPTIONAL and can be omitted. Defaults to: 2 1 3
rsonde_column_index = ${converter.rsonde_h_column.pretty}, ${converter.rsonde_p_column.pretty}, ${converter.rsonde_t_column.pretty}${f", {converter.rsonde_rh_column.pretty}" if converter.rsonde_rh() else ""}

#rsonde_column_units: Radiosonde parser option. The units of Height, Pressure, Temperature, and Relative Humidity (optional) columns in the radiosonde file. The number of values must be the same as in rsonde_column_index. Supported units for:
#    • height: m_asl (default), km_asl, m_agl, km_agl
#    • pressure: hPa (default), Pa, atm
#    • temperature: K (default), C, Cx10
#    • relative humidity: percent (default), fraction
#If the Height units are in agl (above ground level) then the station altitude must be provided in the rsonde_geodata
rsonde_column_units = ${converter.rsonde_h_unit.pretty}, ${converter.rsonde_p_unit.pretty}, ${converter.rsonde_t_unit.pretty}${f", {converter.rsonde_rh_unit.pretty}" if converter.rsonde_rh() else ""}

#rsonde_latitude: The radiosonde station latitude in degrees. Defaults to None.
rsonde_latitude = ${converter.rsonde_latitude.pretty}

#rsonde_longitude: The radiosonde station longitude in degrees. Defaults to None.
rsonde_longitude = ${converter.rsonde_longitude.pretty}

#rsonde_altitude: The radiosonde station altitude in meters. Defaults to None. Mandatory if the radiosonde_column_units of the altitude is either m_agl or km_agl.  Defaults to None. In that case the lowermost altitude reported in the radiosonde file is used instead.
rsonde_altitude = ${converter.rsonde_altitude.pretty}

#rsonde_station_name: The name of the station where the radiosonde measurement was performed. Defaults to: None.
rsonde_station_name = ${converter.rsonde_station_name.pretty}

#rsonde_wmo_number: The WMO number of the radiosonde station. Defaults to: None
rsonde_wmo_number = ${converter.rsonde_wmo_number.pretty}

#rsonde_wban_number: The WBAN number of the radiosonde station. Defaults to: None
rsonde_wban_number = ${converter.rsonde_wban_number.pretty}

[preprocessor]
#vertical_trimming: If set to True then bins above a certain altitude (20km by default) will be removed. Can speed up computations and save RAM. Defaults to: True
vertical_trimming = ${preprocessor.vertical_trimming.pretty}

#vertical_limit: The maximum altitude in km above which no calculations will be performed. Solar background calculations are performed prior to vertical signal trimming to enable background calculations up to the maximum signal range. Defaults to: 20km 
vertical_limit = ${preprocessor.vertical_limit.pretty}

#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels = ${preprocessor.channels.pretty}

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type = ${preprocessor.exclude_telescope_type.pretty}

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_scattering_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_scattering_type = ${preprocessor.exclude_scattering_type.pretty}

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode = ${preprocessor.exclude_acquisition_mode.pretty}

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. 
exclude_channel_subtype = ${preprocessor.exclude_channel_subtype.pretty}

#-------------------------------------------------------------------------------------------------
# Visualizer options following below. Please note:
#    • If new data is set to False, the converter and preprocessor sections will not be taken into
#      account as long as there are preprocessed output files. 
#    • A combination of the visualize, process, and process_qck options will determine which of
#      the following sections will be taken into account
#-------------------------------------------------------------------------------------------------

[quicklooks]
#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels = ${quicklooks.channels.pretty}

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type =${quicklooks.exclude_telescope_type.pretty}

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_scattering_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_scattering_type = ${quicklooks.exclude_scattering_type.pretty}

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode = ${quicklooks.exclude_acquisition_mode.pretty}

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. 
exclude_channel_subtype = ${quicklooks.exclude_channel_subtype.pretty}

#use_range: If set to False, the y axis units of the quicklook will correspond to the altitude in meters above sea level. If set to True, the y axis units corresponds to the distance from the system. This variable determines the units of the following variables: 
#    • y_lims 
#    • z_max_zone
#    • smoothing_range
#Defaults to: True
use_range = ${quicklooks.use_range.pretty}

#x_lims: The x axis limits (lower and upper). Use two integers corresponding to the first and last temporal samples (not date!) that will be plotted in the quicklooks. Use 1 to start from the first sample. If values below 1 or above the total number of samples are used, they will be ignored
x_lims = ${quicklooks.x_lims.pretty}

#x_tick: The x axis finest major tick in number of samples. Defaults to: automatic selection.
x_tick = ${quicklooks.x_tick.pretty}

#t_lims: The x axis limits in time units. Use the following format hh:mm for both limits (not meant to be used for day-long quicklooks). Defaults to: automatic selection.
t_lims = ${quicklooks.lower_t_lim.pretty} ${quicklooks.upper_t_lim.pretty}

#t_tick: The major tick for the time axis (same as x axis but with different units). Defaults to: automatic selection
t_tick = ${quicklooks.t_tick.pretty}

#y_lims: The y axis limits (altitude/distance) in km (lower and upper). If values below 0 or above the maximum altitude/distance are used, they will be ignored. Defaults to: 0., 14.
y_lims = ${quicklooks.y_lims.pretty}

#y_tick: The y axis finest major tick in km. Defaults to: 1km 
y_tick = ${quicklooks.y_tick.pretty}

#z_lims: The colorscale limits of the normalized range-corrected signal normalized with the mean value in the z_max_zone region. This setting is useful if for example the normalization takes place in a very strong aerosol layer making all other layers difficult to discern.  Defaults to: 
#    • z_lims = 0., 1. when use_log_scale = False and 
#    • z_lims = 1E-5, 1. when use_log_scale = True 
z_lims = ${quicklooks.z_lims.pretty}

#z_max_zone: Provide a zone (min and max altitude/distance) in km. The maximum range-corrected signal value encountered in the zone will be used for the normalization of the range-corrected signal for the quicklook. Particularly useful in order to avoid scaling the colors with a cloud. Defaults to: 0.1, 2.
z_max_zone = ${quicklooks.z_max_zone.pretty}

#smooth: If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. The smoothing_exponential, smoothing_exponential, and smoothing_window will be ignored if smooth is set to False. Defaults to: False
smooth = ${quicklooks.smooth.pretty}

#smoothing_range: Set the first and last altitude/distance boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 0.05, 14.
smoothing_range = ${quicklooks.smoothing_range.pretty}

#smoothing_window: The smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the same value twice to apply a constant window. Defaults to: smoothing_window = 50., 500.
smoothing_window = ${quicklooks.smoothing_window.pretty}

#smoothing_exponential: This variable is ignored if the upper and lower values of smoothing_window are the same. Choose one of:
#    • True: a smoothing window that exponentially increases with altitude/distance will be applied. 
#    • False:  a smoothing window that exponentially increases with altitude/distance will be applied
#Defaults to: False.
smoothing_exponential = ${quicklooks.smoothing_exponential.pretty}

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi = ${quicklooks.dpi.pretty}

#color_reduction: If set to True, and provided that Imagemagick is installed, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to False. A warning will be raised if Imagemagick is not installed.
color_reduction = ${quicklooks.color_reduction.pretty}

[rayleigh_fit]
#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels = ${rayleigh_fit.channels.pretty}

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type = ${rayleigh_fit.exclude_telescope_type.pretty}

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_scattering_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_scattering_type =${rayleigh_fit.exclude_scattering_type.pretty}

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode = ${rayleigh_fit.exclude_acquisition_mode.pretty}

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subty Defaults to: b, s, a, w, c
exclude_channel_subtype = ${rayleigh_fit.exclude_channel_subtype.pretty}

#use_range: If set to False, the x axis units of the Rayleigh fit will correspond to the altitude in meters above sea level. If set to True, the x axis units corresponds to the distance from the system. This variable determines the units of the following variables: 
#    • reference_height
#    • x_lims 
#    • smoothing_range
#Defaults to: True
use_range = ${rayleigh_fit.use_range.pretty}

#use_lin_scale: If set to True, a linear scale will be used for the y axis (signal axis). If set to False a logarithmic scale will be applied instead. Defaults to False:
use_lin_scale = ${rayleigh_fit.use_lin_scale.pretty}

#normalization_region: The lower and upper limits of the region used for normalizing the signal in the Rayleigh fit. If use_range is called, the limits correspond to distance. Defaults to: 8.5, 9.5
normalization_region = ${rayleigh_fit.normalization_region.pretty}

#x_lims: Set the x axis (altitude/distance) limits in km (lower and upper). If values below 0 or above the maximum altitude/distance are used, they will be ignored. Defaults to: 0., 14.
x_lims = ${rayleigh_fit.x_lims.pretty}

#x_lims: The x axis finest major tick in km. Defaults to: 1km 
x_tick = ${rayleigh_fit.x_tick.pretty}

#y_lims: The y axis (signal) limits (lower and upper) of the normalized rangecorrected signal in m-1 sr-1 (pseudo-units). It is recommended to skip this variable and use the automatic selection because channels in deifferent wavelengths have different optimal limits. Defaults to: automatic selection
y_lims = ${rayleigh_fit.y_lims.pretty}

#smooth: Refer to the smooth option in the quicklook section. Defaults to: False
smooth = ${rayleigh_fit.smooth.pretty}

#smoothing_exponential: Refer to the smooth option in the quicklook section. Defaults to: False.
smoothing_exponential = ${rayleigh_fit.smoothing_exponential.pretty}

#smoothing_range: Refer to the smooth option in the quicklook section. Defaults to: 1., 14.
smoothing_range = ${rayleigh_fit.smoothing_range.pretty}

#smoothing_window: Refer to the smooth option in the quicklook section. Defaults to: smoothing_window = 100., 100.
smoothing_window = ${rayleigh_fit.smoothing_window.pretty}

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi = ${rayleigh_fit.dpi.pretty}

#color_reduction: If set to True, and provided that Imagemagick is installed, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to False. A warning will be raised if Imagemagick is not installed.
color_reduction = ${rayleigh_fit.color_reduction.pretty}

[telecover]
#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels = ${telecover.channels.pretty}

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type = ${telecover.exclude_telescope_type.pretty}

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_scattering_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_scattering_type = ${telecover.exclude_scattering_type.pretty}

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode = ${telecover.exclude_acquisition_mode.pretty}

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. 
exclude_channel_subtype = ${telecover.exclude_channel_subtype.pretty}

#use_range: If set to False, the x axis units of the Telecover test will correspond to the altitude in meters above sea level. If set to True, the x axis units corresponds to the distance from the system. This variable determines the units of the following variables: 
#    • normalization_height
#    • x_lims 
#    • smoothing_range
#Defaults to: True
use_range = ${telecover.use_range.pretty}

#x_lims: Set the x axis (altitude/distance) limits in km (lower and upper). If values below 0 or above the maximum altitude/distance are used, they will be ignored. Defaults to: 0., 2.5
x_lims = ${telecover.x_lims.pretty}

#x_lims: The x axis finest major tick in km. Defaults to: 0.2km 
x_tick = ${telecover.x_tick.pretty}

#y_lims: The y axis (signal) limits (lower and upper) of the normalized rangecorrected signal in m-1 sr-1 (pseudo-units). It is recommended to skip this variable and use the automatic selection because channels in deifferent wavelengths have different optimal limits. Defaults to: automatic selection
y_lims = ${telecover.y_lims.pretty}

#use_non_rangecor: If set to True, the non range corrected signals will be used for the left subplot of the telecover test. If set to False the range corrected signals will be used instead. Defaults to: False
use_non_rangecor = ${telecover.use_non_range_cor.pretty}

#use_last: If set to True an additional black line with the difference between the first and the last sector is added in the 4th subplot (deviations) of an interleaved telecover test e.g. N2 – N1. If the telecover is not interleaved the difference line will be visible even if use_last is set to False.
use_last = ${telecover.use_last.pretty}

#normalization_region: The reference altitude/distance in km where the signals will be normalized for the Telecover test (relevant to middle and right subplots). It should be well above the suspected full overlap height. Defaults to: normalization_region = 1.8, 2.2
normalization_region = ${telecover.normalization_region.pretty}

#smooth: Refer to the smooth option in the quicklook section. Defaults to: True
smooth = ${telecover.smooth.pretty}

#smoothing_exponential: Refer to the smooth option in the quicklook section. Defaults to: False.
smoothing_exponential = ${telecover.smoothing_exponential.pretty}

#smoothing_range: Refer to the smooth option in the quicklook section. Defaults to: 0., 2.5
smoothing_range = ${telecover.smoothing_range.pretty}

#smoothing_window: Refer to the smooth option in the quicklook section. Defaults to: smoothing_window = 100., 100.
smoothing_window = ${telecover.smoothing_window.pretty}

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi = ${telecover.dpi.pretty}

#color_reduction: If set to True, and provided that Imagemagick is installed, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to False. A warning will be raised if Imagemagick is not installed.
color_reduction = ${telecover.color_reduction.pretty}

[polarization_calibration]
#ch_r: Provide here channel names from the available ones with the reflected channel_subtype (e.g. 0355xpar). The number of reflected channels must be the same as the number of the respective transmitted channels provided by ch_t and by the GHK parameters. By default, all available reflected channels will be matched to all available transmitted channels that share the same telescope type, detection mode, and wavelength. WARNING! The field is mandatory if non-default GHK values are applied. 
ch_r = ${', '.join([cp.ReflectedChannel() for cp in polarization_calibration.channel_pairs()])}
 
#ch_t: Provide here channel names from the available ones with the transmitted channel_subtype (e.g. 0355xcat). The number of reflected channels must be the same as the number of the respective transmitted channels provided by ch_r and by the GHK parameters. By default, all available reflected channels will be matched to all available transmitted channels that share the same telescope type, detection mode, and wavelength. WARNING! The field is mandatory if non-default GHK values are applied. 
ch_t = ${', '.join([cp.TransmittedChannel() for cp in polarization_calibration.channel_pairs()])}

#calibration_region: The lower and upper limits of the region used for Δ90 calibration. If use_range is called, the limits correspond to distance. Defaults to: 2., 4. km 
calibration_region = ${polarization_calibration.calibration_region.pretty}

#rayleigh_region: The lower and upper limits of the region used for the comparison with the Rayleigh atmosphere. Defaults to: 8.5, 9.5
rayleigh_region = ${polarization_calibration.rayleigh_region.pretty}

#x_lims_calibration: The x axis (altitude/distance) limits in km (lower and upper) for the pol. calibration subplot. If values below 0 or above the maximum altitude/distance are used, they will be ignored. Defaults to: 0, 5
x_lims_calibration = ${polarization_calibration.x_lims_calibration.pretty}

#x_lims_rayleigh: The x axis (altitude/distance) limits in km (lower and upper) for the Rayleigh comparison subplot. If values below 0 or above the maximum altitude/distance are used, they will be ignored. Defaults to: 0, 10
x_lims_rayleigh = ${polarization_calibration.x_lims_rayleigh.pretty}

#x_tick_calibration: The x axis finest major tick in km for the pol. calibration sunplot. Defaults to: 0.5
x_tick_calibration = ${polarization_calibration.x_tick_calibration.pretty}

#x_tick_rayleigh: The x axis finest major tick in Km for the Rayleigh plot. Defaults to: 1.
x_tick_rayleigh = ${polarization_calibration.x_tick_rayleigh.pretty}

#y_lims_calibration: The y axis limits (lower and upper) of the gain ratios at +-45. Used for the pol. calibration subplot. Defaults to: automatic selection
y_lims_calibration = ${polarization_calibration.y_lims_calibration.pretty}

#y_lims_rayleigh: The y axis limits (lower and upper) of the volume depolarization ratio. Used for the Rayleigh subplot. Defaults to: automatic selection
y_lims_rayleigh = ${polarization_calibration.y_lims_rayleigh.pretty}

#K: The K value for each channel pair. Defaults to: 1 for all channels 
K = ${', '.join([f"{cp.K():.2f}" for cp in polarization_calibration.channel_pairs()])}

#G_R: The G value for the reflected channel of the pair. Defaults to: 1 for all channels (assuming no correction for the emission and the receiver)
G_R = ${', '.join([f"{cp.G_R():.2f}" for cp in polarization_calibration.channel_pairs()])}

#G_T: The G value for the transmitted channel of the pair. Defaults to: 1 for all channels (assuming no correction for the emission and the receiver)
G_T = ${', '.join([f"{cp.G_T():.2f}" for cp in polarization_calibration.channel_pairs()])}

#H_R: The H value for the reflected channel of the pair. Defaults to: 
#    • 1 for all co-polar (p) reflected channels
#    • -1 for all  cross-polar (c) reflected channels
#(assuming no correction for the emission and the receiver)
H_R = ${', '.join([f"{cp.H_R():.2f}" for cp in polarization_calibration.channel_pairs()])}

#H_T: The H value for the transmitted channel of the pair. Defaults to: 1 or -1 for all co-polar (p) and cross-polar (c) transmitted channels, respectively (no receiver optics + emitted pcb. state correction)
#    • 1 for all co-polar (p) transmitted channels
#    • -1 for all  cross-polar (c) transmitted channels
#(assuming no correction for the emission and the receiver)
H_T = ${', '.join([f"{cp.H_T():.2f}" for cp in polarization_calibration.channel_pairs()])}

#R_to_T_transmission_ratio: The transmission ratio of the R to T channels per pair (TR/TT). TR = 1 and/or TT = 1 if no filter was applied during calibration. Defaults to: 1 for all pairs
R_to_T_transmission_ratio = ${', '.join([f"{cp.RT_T_Ratio():.2f}" for cp in polarization_calibration.channel_pairs()])}

#smooth: If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. Defaults to: True
smooth = ${polarization_calibration.smooth.pretty}

#smoothing_exponential: This variable is ignored if smooth = False and also if the upper and lower values of the smoothing_window variable are the same. Choose one of:
#    • True: a smoothing window that exponentially increases with altitude/distance will be applied. 
#    • False:  a smoothing window that exponentially increases with altitude/distance will be applied
#Defaults to: False.
smoothing_exponential = ${polarization_calibration.smoothing_exponential.pretty}

#smoothing_range: Set the first and last altitude/distance boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 1., 14.
smoothing_range = ${polarization_calibration.smoothing_range.pretty}

#smoothing_window: Half smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the same value twice to apply a constant window. Defaults to: smoothing_window = 100., 100.
smoothing_window = ${polarization_calibration.smoothing_window.pretty}

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi = ${polarization_calibration.dpi.pretty}

#color_reduction: If set to True, and provided that Imagemagick is installed, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True. A warning will be raised if Imagemagick is not installed.
color_reduction = ${polarization_calibration.color_reduction.pretty}

