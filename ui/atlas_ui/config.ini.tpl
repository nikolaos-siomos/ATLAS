[Lidar]
#-------------------------------------------------------------------------------------------------------------
#Mandatory Variables
#-------------------------------------------------------------------------------------------------------------
#lidar_name: The lidar name. This field will be reported in the plots. 
lidar_name = ${lidar.lidar_name.pretty}

#station_id: The lidar id according to the EARLINET DB. 
station_id = ${lidar.station_id.pretty}

#-------------------------------------------------------------------------------------------------------------
#Optional Variables
#-------------------------------------------------------------------------------------------------------------

#altitude: The altitude of the station above sea level (in meters). 
altitude = ${lidar.altitude.pretty}

#latitude: The station latitude (in degrees). 
latitude = ${lidar.latitude.pretty}

#longitude: The station latitude (in degrees). 
longitude = ${lidar.longitude.pretty}

#zenith_angle: Zenith angle of the lidar (in degrees, 0 at zenith, 90 at horizon). 
zenith_angle = ${lidar.zenith_angle.pretty}

#azimuth_angle: Azimuth angle of the lidar (in degrees, North = 0, E = 90).
azimuth_angle = ${lidar.azimuth_angle.pretty}

[Channels]
# All variables provided for this section should include exactly as many channels as the recorder_channel_id variable. 
# If a variable is not relevant for a specific channel (e.g. dead_time for analog channels) fill with "_". 
#-------------------------------------------------------------------------------------------------------------
#Mandatory Variables
#-------------------------------------------------------------------------------------------------------------

#telescope_type: The telescope configuration. Choose among:
#    •n: near range
#    • f: far range
#    • x: single telescope 
telescope_type = ${', '.join([channel.telescope_type.required for channel in channels])}

#channel_type: The channel type. Choose one among:
#    • p: co-polar linear analyzer
#    • c: cross-polar linear analyzer
#    • t: total (no depolarization)
#    • o: co-polar circular analyzer
#    • x: cross-polar circular analyzer
#    • v: vibrational Raman
#    • r: rotational Raman
#    • a: Cabannes
#    • f: fluorescence 
channel_type = ${', '.join([channel.channel_type.required for channel in channels])}

#channel_subtype: The channel subtype. Choose among:
#    • r: Signal Reflected from a PBS
#    • t: Signal Transmitted through a PBS
#    • n: N2 Ramal line
#    • o: O2 Ramal line
#    • w: H2O Ramal line
#    • c: CH4 Ramal line
#    • h: High Rotational Raman
#    • l: Low Rotational Raman
#    • a: Mie (aerosol) HSRL signal
#    • m: Molecular HSRL signal
#    • b: Broadband Fluorescence
#    • s: Spectral Fluorescence
#    • x: No specific subtype 
channel_subtype = ${', '.join([channel.channel_subtype.required for channel in channels])}

#dead_time: The dead time of the photon counting channels (in ns). (for analog channels set _).
dead_time = ${print_channels_field(channels, 'dead_time')}

#daq_trigger_offset: The trigger delay per channel in bins. Provide negative values if the recording starts before the Q-switch (pre-triggering) and positive values if the recording starts after the Q-switch. 
daq_trigger_offset = ${print_channels_field(channels, 'daq_trigger_offset')}

#-------------------------------------------------------------------------------------------------------------
#Mandatory Variables for Licel
#-------------------------------------------------------------------------------------------------------------
#recorder_channel_id: IDs of each channel according to the raw file header. For licel provide the licel id per channel that is going to be include. Currently only BT and BC channels are supported. For PollyXT systems all channels are included by default. Provide only if you want to process specific channels
recorder_channel_id = ${print_channels_field(channels, 'recorder_channel_id')}

#laser: The ascending laser number according to the licel manual. Use to link with the licel file. The recorder_channel_id is not a unique identifier for licel channels. A single channel can synchronize with more than one lasers. This variable is ignored when processing PollYXT files
laser = ${print_channels_field(channels, 'laser')}

#-------------------------------------------------------------------------------------------------------------
#Mandatory Variables in the Future
#-------------------------------------------------------------------------------------------------------------
#scc_channel_id: IDs of each channel according to the SCC configuration. In the future, linking with the HOI will be done through the scc_channel_id (currently optional). 
scc_channel_id = ${print_channels_field(channels, 'scc_channel_id')}

#-------------------------------------------------------------------------------------------------------------
#Partly Optional Variables
#-------------------------------------------------------------------------------------------------------------
#These variables take default values. It is highly recommended though to fill them in because  might not be valid for the system.
#dead_time_correction_type: The type of dead time correction. Choose among:
#    • 0: Non Paralyzable
#    • 1: Paralyzable
# Defaults to: 0 for photon channels.
dead_time_correction_type = ${print_channels_field(channels, 'dead_time_correction_type')}

#emitted_wavelength: The wavelength of the originally emitted laser beam per channel. A very rough guess process is applied by default using that detected channel information that it is valid ONLY for Nd:Yag lasers with elastic and vibrational Raman channels. Please provide the values explicitly if your system does not fall in this category. Defaults to:
#    • 355 for channels with detected wavelength <520nm
#    • 532 for channels with detected wavelength <1000nm 
#    • 1064 for channels with detected wavelength >=1000nm
#If accurate molecular depolarization ratio calculations are required this variable should be filled in with the actual emitted wavelength with sub-nanometer accuracy.
emitted_wavelength = ${print_channels_field(channels, 'emitted_wavelength')}

#detected_wavelength: The wavelength of the detected signal according to the applied Interference Filter of each channel (in nm). This variable is provided by default from the raw file metadata for both Licel and PollyXT. However, if accuracy is required for molecular calculations the exact central wavelength should be provided.
detected_wavelength = ${print_channels_field(channels, 'detected_wavelength')}

#channel_bandwidth: Interference filter Bandwidth (in nm, Defaults to: 1nm). The default value is a very rough guess and it can easily lead to high uncertainty in the molecular depolarization value. Exact values provided from the manufacturers should be used instead.
channel_bandwidth = ${print_channels_field(channels, 'channel_bandwidth')}

#background_low_bin: Starting bin of the background correction averaging range. Defaults to:
#    • The 600th bin before the end of the profile if the daq_trigger_offset < 400
#    • The 100th bin if the daq_trigger_offset >= 400 
#Using default values is currently risky. Please make sure that no spikes are appearing in the pretrigger range that would affect the background correction or manually provide this variable.
% if len([channel.background_range() for channel in channels if channel.background_range() is not None]) > 0:
background_low_bin = ${', '.join([f"{channel.background_range()[0]:.0f}" if channel.background_range() is not None else "_" for channel in channels])}
% else:
background_low_bin = 
% endif

#background_high_bin: ending bin of the background correction averaging range. Defaults to:
#    • The 100th bin before the end of the profile if the daq_trigger_offset < 400
#    • The 300th bin if the daq_trigger_offset >= 400 
#Using default values is currently risky. Please make sure that no spikes are appearing in the pretrigger range that would affect the background correction or manually provide this variable.
% if len([channel.background_range() for channel in channels if channel.background_range() is not None]) > 0:
background_high_bin = ${', '.join([f"{channel.background_range()[1]:.0f}" if channel.background_range() is not None else "_" for channel in channels])}
% else:
background_high_bin = 
% endif

#-------------------------------------------------------------------------------------------------------------
#Optional Variables (include only to override the raw file metadata)
#-------------------------------------------------------------------------------------------------------------
#acquisition_mode: The mode of the recorded signals per channel. Choose among:
#    • 0: analog
#    • 1: photon counting
acquisition_mode = ${print_channels_field(channels, 'acquisition_mode')}

#bins:  The total bins of the recorded signals per channel. 
bins = ${print_channels_field(channels, 'bins')}

#laser_shots: The number of acquired laser shots per channel. Not recommended to provide it manually. If fetched from the metadata, this variable is different per file which is more accurate than providing a constant value here.
laser_shots = ${print_channels_field(channels, 'laser_shots')}

#data_acquisition_range: The Data Acquisition Range of each analog channel [mV]. Use _ for photon channels
data_acquisition_range = ${print_channels_field(channels, 'data_acquisition_range')}

#analog_to_digital_resolution: The analog to digital resolution (in bits) of each analog channel. Use _ for photon channels.
analog_to_digital_resolution = ${print_channels_field(channels, 'analog_to_digital_resolution')}

#range_resolution: The range resolution of each channel [m]. It will be used to calculate the Sampling Rate.
range_resolution = ${print_channels_field(channels, 'range_resolution')}

#pmt_high_voltage: The high voltage [V] provided to the detection unit per channel.
pmt_high_voltage = ${print_channels_field(channels, 'pmt_high_voltage')}

#laser_repetition_rate: The Laser Repetition Rate [Hz]
laser_repetition_rate = ${print_channels_field(channels, 'laser_repetition_rate')}

