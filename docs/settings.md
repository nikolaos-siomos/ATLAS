# Settings INI reference (WIP)

*Generated from* `settings_file.ini` on 2026-02-11.

> This page is generated from the INI file comments so it stays in sync. Add/edit comments in the INI file and re-run the generator.

## File purpose

_TODO: Add a short narrative description here._

## Sections and options

### [converter]

| Key | Default / Example | Help |
|---|---|---|
| `debug` | `` | debug: If set to True, debugging files will be generated in the  <path_to_the_converter_folder>/debug folder. Note that the path to the converter folder Defaults to: <parent_folder>/netcdf but it can be different if manually provided by the user. These included the metadata gathered from the configuration file, the licel header, and the combination of the two. Defaults to: False |
| `trim_overflows` | `` | trim_overflows: This options determines how overflow values in the raw input files will be treated if detected. Choose among:<br>• 0: the algorithm will stop and provide a diagnostic error (default)<br>• 1: the files containing at least one overflow value will be screened out<br>• 2: overflows will be interpolated (use with care and only only a few bins per profile have overflow values)<br>• 3: overflows will not be excluded, use this only for debugging purposes |
| `files_per_sector` | `` | files_per_sector: The number of telecover files per sector (integer). If provided, the telecover files must be placed in a single folder (this should be <path_to_parent_folder>/tlc according to the default folder structure). An automated assignment of the telecover files in different sectors will be attempted serially assuming the following temporal sequence of sectors: north – east – south – west. Note that the telecover test can have more than 1 rounds. |
| `files_per_ring` | `` | files_per_ring: The number of telecover files per ring (integer). If provided, the telecover files must be placed in a single folder (this should be <path_to_parent_folder>/tlc_rin according to the default folder structure). An automated assignment of the telecover files in different sectors will be attempted serially assuming the following temporal sequence of sectors: inner - outer. Note that the telecover test can have more than 1 rounds. |
| `rsonde_skip_header` | `` | rsonde_skip_header: Radiosonde parser option. Number of lines to skip at the beginning of the radiosonde ascii file. Defaults to: 1 (1 line reserved for header info) |
| `rsonde_skip_footer` | `` | rsonde_skip_footer: Radiosonde parser option. Number of lines to skip at the end of the radiosonde file. Defaults to: 0 (no footer assumed) |
| `rsonde_delimiter` | `` | rsonde_delimiter: Radiosonde parser option. The delimiter that separates columns in the radiosonde file choose one of:<br>• S: space<br>• C: comma<br>Defaults to: S. |
| `rsonde_column_index` | `` | rsonde_column_index: Radiosonde parser option. The index (integer) of the column that contains the Height, Pressure, Temperature, and Relative Humidity (optional) information. For example, rsonde_columns = 1, 3, 2, 6 means:<br>• Height: 1st column in the radiosonde file<br>• Pressure: 3rd column<br>• Temperature: 2nd column<br>• Relative Humidity: 6th column<br>The relative humidity column is OPTIONAL and can be omitted. Defaults to: 2 1 3 |
| `rsonde_column_units` | `` | rsonde_column_units: Radiosonde parser option. The units of Height, Pressure, Temperature, and Relative Humidity (optional) columns in the radiosonde file. The number of values must be the same as in rsonde_column_index. Supported units for:<br>• height: m_asl (default), km_asl, m_agl, km_agl<br>• pressure: hPa (default), Pa, atm<br>• temperature: K (default), C, Cx10<br>• relative humidity: percent (default), fraction<br>If the Height units are in agl (above ground level) then the station altitude must be provided in the rsonde_geodata |
| `rsonde_latitude` | `` | rsonde_latitude: The radiosonde station latitude in degrees. Defaults to None. |
| `rsonde_longitude` | `` | rsonde_longitude: The radiosonde station longitude in degrees. Defaults to None. |
| `rsonde_altitude` | `` | rsonde_altitude: The radiosonde station altitude in meters. Defaults to None. Mandatory if the radiosonde_column_units of the altitude is either m_agl or km_agl.  Defaults to None. In that case the lowermost altitude reported in the radiosonde file is used instead. |
| `rsonde_station_name` | `` | rsonde_station_name: The name of the station where the radiosonde measurement was performed. Defaults to: None. |
| `rsonde_wmo_number` | `` | rsonde_wmo_number: The WMO number of the radiosonde station. Defaults to: None |
| `rsonde_wban_number` | `` | rsonde_wban_number: The WBAN number of the radiosonde station. Defaults to: None |

### [preprocessor]

| Key | Default / Example | Help |
|---|---|---|
| `vertical_trimming` | `` | vertical_trimming: If set to True then bins above a certain height above the lidar (20km by default) will be removed. Can speed up computations and save RAM. Defaults to: True |
| `vertical_limit` | `` | vertical_limit: The maximum height above the lidar in km above which no calculations will be performed. Solar background calculations are performed prior to vertical signal trimming to enable background calculations up to the maximum signal range. Defaults to: 20km |
| `channels` | `` | channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3. |
| `exclude_telescope_type` | `` | exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_type` | `` | exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:<br>• if isday = False --> no channel is excluded if this variable is not provided<br>• if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded) |
| `exclude_acquisition_mode` | `` | exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_subtype` | `` | exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. |

### [quicklooks]

-------------------------------------------------------------------------------------------------
Visualizer options following below. Please note:
• If new data is set to False, the converter and preprocessor sections will not be taken into
account as long as there are preprocessed output files.
• A combination of the visualize, process, and process_qck options will determine which of
the following sections will be taken into account
-------------------------------------------------------------------------------------------------

| Key | Default / Example | Help |
|---|---|---|
| `channels` | `` | channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3. |
| `exclude_telescope_type` | `` | exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_type` | `` | exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:<br>• if isday = False --> no channel is excluded if this variable is not provided<br>• if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded) |
| `exclude_acquisition_mode` | `` | exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_subtype` | `` | exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. |
| `use_range` | `` | use_range: If set to False, the y axis units of the quicklook will correspond to the height in km above the lidar. If set to True, the y axis units corresponds to the range above the lidar. This variable determines the units of the following variables:<br>• y_lims<br>• z_max_zone<br>• smoothing_range<br>Defaults to: True |
| `t_lims` | `` | t_lims: The x axis limits in time units. Use the following format: hhmm, hhmm for both limits (not meant to be used for daylong quicklooks). Defaults to: automatic selection. |
| `t_tick` | `` | t_tick: The major tick for the time axis (same as x axis but with different units). Defaults to: automatic selection |
| `y_lims` | `` | y_lims: The y axis limits (height/range) in km (lower and upper). Defaults to: 0., 14. |
| `y_tick` | `` | y_tick: The y axis finest major tick in km. Defaults to: 1km |
| `z_lims` | `` | z_lims: The colorscale limits of the normalized range-corrected signal normalized with the mean value in the z_max_zone region. This setting is useful if for example the normalization takes place in a very strong aerosol layer making all other layers difficult to discern.  Defaults to:<br>• z_lims = 0., 1. when use_log_scale = False and<br>• z_lims = 1E-5, 1. when use_log_scale = True |
| `z_max_zone` | `` | z_max_zone: Provide a zone (min and max height/range) in km. The maximum range-corrected signal value encountered in the zone will be used for the normalization of the range-corrected signal for the quicklook. Particularly useful in order to avoid scaling the colors with a cloud. Defaults to: 0.1, 2. |
| `smooth` | `` | smooth: If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. The smoothing_exponential, smoothing_exponential, and smoothing_window will be ignored if smooth is set to False. Defaults to: False |
| `smoothing_range` | `` | smoothing_range: Set the first and last height/range boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 0.05, 14. |
| `smoothing_window` | `` | smoothing_window: The smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the same value twice to apply a constant window. Defaults to: smoothing_window = 50., 500. |
| `smoothing_exponential` | `` | smoothing_exponential: This variable is ignored if the upper and lower values of smoothing_window are the same. Choose one of:<br>• True: a smoothing window that exponentially increases with height/range will be applied.<br>• False:  a smoothing window that exponentially increases with height/range will be applied<br>Defaults to: False. |
| `dpi` | `` | dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi |
| `color_reduction` | `` | color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True |

### [rayleigh_fit]

| Key | Default / Example | Help |
|---|---|---|
| `channels` | `` | channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3. |
| `exclude_telescope_type` | `` | exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_type` | `` | exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:<br>• if isday = False --> no channel is excluded if this variable is not provided<br>• if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded) |
| `exclude_acquisition_mode` | `` | exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_subtype` | `` | exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subty Defaults to: b, s, a, w, c |
| `use_range` | `` | use_range: If set to False, the x axis units of the Rayleigh fit will correspond to the height in km above the lidar. If set to True, the x axis units corresponds to the range above the lidar. This variable determines the units of the following variables:<br>• reference_height<br>• x_lims<br>• smoothing_range<br>Defaults to: True |
| `use_lin_scale` | `` | use_lin_scale: If set to True, a linear scale will be used for the y axis (signal axis). If set to False a logarithmic scale will be applied instead. Defaults to False: |
| `cross_check_lower_limit` | `` | cross_check_lower_limit: Lower limit applied for the cross-check criterion for the calculation of the Rayleigh fit mask. Use in case a high distance of full overlap is causing the cross-check test to fail. Defaults to 2km |
| `normalization_region` | `` | normalization_region: The lower and upper limits of the region used for normalizing the signal in the Rayleigh fit in km. If use_range is set to True, the limits correspond to distance. If not provided, the normalization region will be automatically identified (default by CARS) Example: normalization_region = 4., 6. |
| `fit_mask_region` | `` | fit_mask_region:  Lower and upper thresholds for the region where the Rayleigh fit mask will be applied in km. The cross criterion might might fail if the lower limit is below the distance of full overlap. Defaults to: 3 to 34 km Example: fit_mask_region = 4., 20. |
| `fit_mask_window` | `` | fit_mask_window: The size limits (min and max) of the window used for Rayleigh fit mask. The mask will be applied for window sizes ranging between the two provided values. Defaults to: 2., 6. km Example: fit_mask_window = 0.5, 4. |
| `fit_mask_window_step` | `` | fit_mask_window_step: The center of the window used for the Rayleigh fit mask will range between the fit_mask_window with a step equal to the value provided here in km. Keep in mind that smaller values will increase the resolution of the mask but will make processing slower. Defaults to: 0.2 km. Example: fit_mask_window_step = 0.1 |
| `rsem_threshold` | `` | rsem_threshold: The relative standard error of the mean threshold that is applied by the Rayleigh fit mask. Regions with relative SEM above this threshold will not be considered molecular. Defaults to 0.2 (20%). Example: rsem_threshold = 0.02 |
| `first_derivative_threshold` | `` | first_derivative_threshold: The standard deviation threshold applied by the Rayleigh fit mask. Regions where the first derivative is less than its uncertainty multiplied by the derivative_threshold will not be considered molecular. Defaults to 2. (2 * sigma). Example: first_derivative_threshold = 1. |
| `second_derivative_threshold` | `` | second_derivative_threshold: The standard deviation threshold applied by the Rayleigh fit mask. Regions where the second derivative is less than its uncertainty multiplied by the derivative_threshold will not be considered molecular. Defaults to 2. (2 * sigma). Example: second_derivative_threshold = 1. |
| `shapiro_wilk_threshold` | `` | shapiro_wilk_threshold: The p value threshold for the Shapiro-Wilk test applied by the Rayleigh fit mask. Regions where the Shapiro-Wilk test returns p values smaller than this threshold will not be considered molecular because the noise does not follow a Gaussian distribution. Defaults to 0.05. Example: shapiro_wilk_threshold = 0.10 |
| `cross_criterion_threshold` | `` | cross_criterion_threshold: The cross criterion threshold applied by the Rayleigh fit mask when checking if negative differences of the normalized rangecorrected signal and the molecular attenuated backscatter are significant. If even one region below the normalization region is found where the difference between the normalized rangecorected signal and the molecular profile is smaller than the negative uncertainty of the normalized rangecorrected signal multiplied by the cross_criterion_threshold means that the corresponding normalization region is not molecular. Defaults to 1. (1 * sigma). Example: cross_criterion_threshold = 1. |
| `durbin_watson_threshold` | `` | durbin_watson_threshold: Lower and upper thresholds for the Durbin-Watson test. Values close to 2. indicate no auto-correlation. Defaults to [1., 3.]. Example: cross_criterion_threshold = 1.5, 2.5 |
| `x_lims` | `` | x_lims: Set the x axis (height/range) limits in km (lower and upper). Defaults to: 0., 20. |
| `x_tick` | `` | x_tick: The x axis finest major tick in km. Defaults to: 2km |
| `y_lims` | `` | y_lims: The y axis (signal) limits (lower and upper) of the normalized rangecorrected signal in m-1 sr-1 (pseudo-units). It is recommended to skip this variable and use the automatic selection because channels in deifferent wavelengths have different optimal limits. Defaults to: automatic selection |
| `smooth` | `` | smooth: Refer to the smooth option in the quicklook section. Defaults to: False |
| `smoothing_exponential` | `` | smoothing_exponential: Refer to the smooth option in the quicklook section. Defaults to: False. |
| `smoothing_range` | `` | smoothing_range: Refer to the smooth option in the quicklook section. Defaults to: 0.05, 20. |
| `smoothing_window` | `` | smoothing_window: Refer to the smooth option in the quicklook section. Defaults to: smoothing_window = 100., 100. |
| `dpi` | `` | dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi |
| `color_reduction` | `` | color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True |

### [telecover]

| Key | Default / Example | Help |
|---|---|---|
| `channels` | `` | channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3. |
| `exclude_telescope_type` | `` | exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_type` | `` | exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:<br>• if isday = False --> no channel is excluded if this variable is not provided<br>• if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded) |
| `exclude_acquisition_mode` | `` | exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. |
| `exclude_channel_subtype` | `` | exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. |
| `use_range` | `` | use_range: If set to False, the x axis units of the Telecover test will correspond to the height in km above the lidar. If set to True, the x axis units corresponds to the range above the lidar. This variable determines the units of the following variables:<br>• normalization_height<br>• x_lims<br>• smoothing_range<br>Defaults to: True |
| `x_lims` | `` | x_lims: Set the x axis (height/range) limits for the near range subplots of the telecover test in km (lower and upper). Defaults to: 0., 2.5 for l, m, n, and x telescope types and: 0., 5. for the rest of the types |
| `x_tick` | `` | x_tick: The x axis finest tick for the near range subplots of the telecover test in km. Defaults to: 0.5 for l, m, n, and x telescope types and: 1. for the rest of the types |
| `y_lims` | `` | y_lims: The y axis (signal) limits (lower and upper) of the normalized rangecorrected signal in m-1 sr-1 (pseudo-units). It is recommended to skip this variable and use the automatic selection because channels in deifferent wavelengths have different optimal limits. Defaults to: automatic selection |
| `use_non_rangecor` | `` | use_non_rangecor: If set to True, the non range corrected signals will be used for the left subplot of the telecover test. If set to False the range corrected signals will be used instead. Defaults to: False |
| `use_last` | `` | use_last: While set to True an additional purple line will be added in all subplots of the telecover test. In the first 3 subplots the lines corresponds to the last sector (e.g. N2), while in the 4th subplot (deviations) it is the difference between the normalized signal of the last sector and the normalized mean signal of the correspondin sector e.g. N2 – N. Set to False to not visualize this line. Defaults to True. |
| `auto_fit` | `` | auto_fit: If set to True an automatic identification of the normalization region will be attempted. If the automatic procedure is successful, the normalization_region variable will be ignored. If the procedure is not successful or auto_fit is set to False, the manually-provided/default values will be used. Defaults to True |
| `normalization_region` | `` | normalization_region: The reference height/range in km where the signals will be normalized for the Telecover test (relevant to middle and right subplots). It should be well above the suspected full overlap height. Defaults to: normalization_region = 1.8, 2.2 |
| `smooth` | `` | smooth: Refer to the smooth option in the quicklook section. Defaults to: True |
| `smoothing_exponential` | `` | smoothing_exponential: Refer to the smooth option in the quicklook section. Defaults to: False. |
| `smoothing_window` | `` | smoothing_window: Refer to the smooth option in the quicklook section. Defaults to: smoothing_window = 100., 100. |
| `dpi` | `` | dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi |
| `color_reduction` | `` | color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True |

### [polarization_calibration]

| Key | Default / Example | Help |
|---|---|---|
| `ch_r` | `` | ch_r: Provide here channel names from the available ones with the reflected channel_subtype (e.g. 0355xpar). The number of reflected channels must be the same as the number of the respective transmitted channels provided by ch_t and by the GHK parameters. By default, all available reflected channels will be matched to all available transmitted channels that share the same telescope type, detection mode, and wavelength. WARNING! The field is mandatory if non-default GHK values are applied. |
| `ch_t` | `` | ch_t: Provide here channel names from the available ones with the transmitted channel_subtype (e.g. 0355xcat). The number of reflected channels must be the same as the number of the respective transmitted channels provided by ch_r and by the GHK parameters. By default, all available reflected channels will be matched to all available transmitted channels that share the same telescope type, detection mode, and wavelength. WARNING! The field is mandatory if non-default GHK values are applied. |
| `auto_fit` | `` | auto_fit: If set to True an automatic identification of the calibration and Rayleigh regions will be attempted. If the automatic procedure is successful, the calibration_region and rayleigh_region variables will be ignored. If the procedure is not successful or auto_fit is set to False, the manually-provided/default values will be used. Defaults to True |
| `calibration_region` | `` | calibration_region: The lower and upper limits of the region used for Δ90 calibration. If use_range is called, the limits correspond to distance. If auto_fit is set to True and the automatic identification is successful for a specific channel combination, the calibration_region variable will be ignored even if provided. Defaults to: 2., 4. km |
| `rayleigh_region` | `` | rayleigh_region: The lower and upper limits of the region used for the comparison with the Rayleigh atmosphere. If auto_fit is set to True and the automatic identification is successful for a specific channel combination, the rayleigh_region variable will be ignored even if provided. Defaults to: 8.5, 9.5 |
| `x_lims_calibration` | `` | x_lims_calibration: The x axis (height/range) limits in km (lower and upper) for the pol. calibration subplot. Defaults to: 0, 15 |
| `x_lims_rayleigh` | `` | x_lims_rayleigh: The x axis (height/range) limits in km (lower and upper) for the Rayleigh comparison subplot. Defaults to: 0, 15 |
| `x_tick_calibration` | `` | x_tick_calibration: The x axis finest major tick in km for the pol. calibration sunplot. Defaults to: 0.5 |
| `x_tick_rayleigh` | `` | x_tick_rayleigh: The x axis finest major tick in Km for the Rayleigh plot. Defaults to: 2. |
| `y_lims_calibration` | `` | y_lims_calibration: The y axis limits (lower and upper) of the gain ratios at +-45. Used for the pol. calibration subplot. Defaults to: automatic selection |
| `y_lims_rayleigh` | `` | y_lims_rayleigh: The y axis limits (lower and upper) of the volume depolarization ratio. Used for the Rayleigh subplot. Defaults to: automatic selection |
| `K` | `` | K: The K value for each channel pair. Defaults to: 1 for all channels |
| `G_R` | `` | G_R: The G value for the reflected channel of the pair. Defaults to: 1 for all channels (assuming no correction for the emission and the receiver) |
| `G_T` | `` | G_T: The G value for the transmitted channel of the pair. Defaults to: 1 for all channels (assuming no correction for the emission and the receiver) |
| `H_R` | `` | H_R: The H value for the reflected channel of the pair. Defaults to:<br>• 1 for all co-polar (p) reflected channels<br>• -1 for all  cross-polar (c) reflected channels<br>(assuming no correction for the emission and the receiver) |
| `H_T` | `` | H_T: The H value for the transmitted channel of the pair. Defaults to: 1 or -1 for all co-polar (p) and cross-polar (c) transmitted channels, respectively (no receiver optics + emitted pcb. state correction)<br>• 1 for all co-polar (p) transmitted channels<br>• -1 for all  cross-polar (c) transmitted channels<br>(assuming no correction for the emission and the receiver) |
| `R_to_T_transmission_ratio` | `` | R_to_T_transmission_ratio: The transmission ratio of the R to T channels per pair (TR/TT). TR = 1 and/or TT = 1 if no filter was applied during calibration. Defaults to: 1 for all pairs |
| `smooth` | `` | smooth: If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. Defaults to: True |
| `smoothing_exponential` | `` | smoothing_exponential: This variable is ignored if smooth = False and also if the upper and lower values of the smoothing_window variable are the same. Choose one of:<br>• True: a smoothing window that exponentially increases with height/range will be applied.<br>• False:  a smoothing window that exponentially increases with height/range will be applied<br>Defaults to: False. |
| `smoothing_range` | `` | smoothing_range: Set the first and last height/range boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 0.25, 15. |
| `smoothing_window` | `` | smoothing_window: Half smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the a single value twice to apply a constant window. Defaults to: smoothing_window = 500. |
| `dpi` | `` | dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi |
| `color_reduction` | `` | color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True |

## Full file (verbatim)

```ini
[converter]
#debug: If set to True, debugging files will be generated in the  <path_to_the_converter_folder>/debug folder. Note that the path to the converter folder Defaults to: <parent_folder>/netcdf but it can be different if manually provided by the user. These included the metadata gathered from the configuration file, the licel header, and the combination of the two. Defaults to: False
debug =

#trim_overflows: This options determines how overflow values in the raw input files will be treated if detected. Choose among:
#    • 0: the algorithm will stop and provide a diagnostic error (default)
#    • 1: the files containing at least one overflow value will be screened out
#    • 2: overflows will be interpolated (use with care and only only a few bins per profile have overflow values)
#    • 3: overflows will not be excluded, use this only for debugging purposes
trim_overflows =

#files_per_sector: The number of telecover files per sector (integer). If provided, the telecover files must be placed in a single folder (this should be <path_to_parent_folder>/tlc according to the default folder structure). An automated assignment of the telecover files in different sectors will be attempted serially assuming the following temporal sequence of sectors: north – east – south – west. Note that the telecover test can have more than 1 rounds.
files_per_sector =

#files_per_ring: The number of telecover files per ring (integer). If provided, the telecover files must be placed in a single folder (this should be <path_to_parent_folder>/tlc_rin according to the default folder structure). An automated assignment of the telecover files in different sectors will be attempted serially assuming the following temporal sequence of sectors: inner - outer. Note that the telecover test can have more than 1 rounds.
files_per_ring =

#rsonde_skip_header: Radiosonde parser option. Number of lines to skip at the beginning of the radiosonde ascii file. Defaults to: 1 (1 line reserved for header info)
rsonde_skip_header = 

#rsonde_skip_footer: Radiosonde parser option. Number of lines to skip at the end of the radiosonde file. Defaults to: 0 (no footer assumed)
rsonde_skip_footer = 

#rsonde_delimiter: Radiosonde parser option. The delimiter that separates columns in the radiosonde file choose one of:
#    • S: space
#    • C: comma
#Defaults to: S.
rsonde_delimiter = 

#rsonde_column_index: Radiosonde parser option. The index (integer) of the column that contains the Height, Pressure, Temperature, and Relative Humidity (optional) information. For example, rsonde_columns = 1, 3, 2, 6 means: 
#    • Height: 1st column in the radiosonde file 
#    • Pressure: 3rd column 
#    • Temperature: 2nd column
#    • Relative Humidity: 6th column
#The relative humidity column is OPTIONAL and can be omitted. Defaults to: 2 1 3
rsonde_column_index = 

#rsonde_column_units: Radiosonde parser option. The units of Height, Pressure, Temperature, and Relative Humidity (optional) columns in the radiosonde file. The number of values must be the same as in rsonde_column_index. Supported units for:
#    • height: m_asl (default), km_asl, m_agl, km_agl
#    • pressure: hPa (default), Pa, atm
#    • temperature: K (default), C, Cx10
#    • relative humidity: percent (default), fraction
#If the Height units are in agl (above ground level) then the station altitude must be provided in the rsonde_geodata
rsonde_column_units = 

#rsonde_latitude: The radiosonde station latitude in degrees. Defaults to None.
rsonde_latitude = 

#rsonde_longitude: The radiosonde station longitude in degrees. Defaults to None.
rsonde_longitude = 

#rsonde_altitude: The radiosonde station altitude in meters. Defaults to None. Mandatory if the radiosonde_column_units of the altitude is either m_agl or km_agl.  Defaults to None. In that case the lowermost altitude reported in the radiosonde file is used instead.
rsonde_altitude = 

#rsonde_station_name: The name of the station where the radiosonde measurement was performed. Defaults to: None.
rsonde_station_name = 

#rsonde_wmo_number: The WMO number of the radiosonde station. Defaults to: None
rsonde_wmo_number = 

#rsonde_wban_number: The WBAN number of the radiosonde station. Defaults to: None
rsonde_wban_number = 

[preprocessor]
#vertical_trimming: If set to True then bins above a certain height above the lidar (20km by default) will be removed. Can speed up computations and save RAM. Defaults to: True
vertical_trimming =

#vertical_limit: The maximum height above the lidar in km above which no calculations will be performed. Solar background calculations are performed prior to vertical signal trimming to enable background calculations up to the maximum signal range. Defaults to: 20km 
vertical_limit =

#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels =

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type =

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_channel_type =

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode =

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. 
exclude_channel_subtype =

#-------------------------------------------------------------------------------------------------
# Visualizer options following below. Please note:
#    • If new data is set to False, the converter and preprocessor sections will not be taken into
#      account as long as there are preprocessed output files. 
#    • A combination of the visualize, process, and process_qck options will determine which of
#      the following sections will be taken into account
#-------------------------------------------------------------------------------------------------

[quicklooks]
#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels =

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type =

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_channel_type =

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode =

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. 
exclude_channel_subtype =

#use_range: If set to False, the y axis units of the quicklook will correspond to the height in km above the lidar. If set to True, the y axis units corresponds to the range above the lidar. This variable determines the units of the following variables: 
#    • y_lims 
#    • z_max_zone
#    • smoothing_range
#Defaults to: True
use_range = 

#t_lims: The x axis limits in time units. Use the following format: hhmm, hhmm for both limits (not meant to be used for daylong quicklooks). Defaults to: automatic selection.
t_lims =

#t_tick: The major tick for the time axis (same as x axis but with different units). Defaults to: automatic selection
t_tick =

#y_lims: The y axis limits (height/range) in km (lower and upper). Defaults to: 0., 14.
y_lims =

#y_tick: The y axis finest major tick in km. Defaults to: 1km 
y_tick = 

#z_lims: The colorscale limits of the normalized range-corrected signal normalized with the mean value in the z_max_zone region. This setting is useful if for example the normalization takes place in a very strong aerosol layer making all other layers difficult to discern.  Defaults to: 
#    • z_lims = 0., 1. when use_log_scale = False and 
#    • z_lims = 1E-5, 1. when use_log_scale = True 
z_lims =

#z_max_zone: Provide a zone (min and max height/range) in km. The maximum range-corrected signal value encountered in the zone will be used for the normalization of the range-corrected signal for the quicklook. Particularly useful in order to avoid scaling the colors with a cloud. Defaults to: 0.1, 2.
z_max_zone = 

#smooth: If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. The smoothing_exponential, smoothing_exponential, and smoothing_window will be ignored if smooth is set to False. Defaults to: False
smooth = 

#smoothing_range: Set the first and last height/range boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 0.05, 14.
smoothing_range = 

#smoothing_window: The smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the same value twice to apply a constant window. Defaults to: smoothing_window = 50., 500.
smoothing_window =

#smoothing_exponential: This variable is ignored if the upper and lower values of smoothing_window are the same. Choose one of:
#    • True: a smoothing window that exponentially increases with height/range will be applied. 
#    • False:  a smoothing window that exponentially increases with height/range will be applied
#Defaults to: False.
smoothing_exponential =

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi =

#color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True
color_reduction = 

[rayleigh_fit]
#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels =

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type =

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_channel_type =

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode =

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subty Defaults to: b, s, a, w, c
exclude_channel_subtype =

#use_range: If set to False, the x axis units of the Rayleigh fit will correspond to the height in km above the lidar. If set to True, the x axis units corresponds to the range above the lidar. This variable determines the units of the following variables: 
#    • reference_height
#    • x_lims 
#    • smoothing_range
#Defaults to: True
use_range = 

#use_lin_scale: If set to True, a linear scale will be used for the y axis (signal axis). If set to False a logarithmic scale will be applied instead. Defaults to False:
use_lin_scale =

#cross_check_lower_limit: Lower limit applied for the cross-check criterion for the calculation of the Rayleigh fit mask. Use in case a high distance of full overlap is causing the cross-check test to fail. Defaults to 2km 
cross_check_lower_limit = 

#normalization_region: The lower and upper limits of the region used for normalizing the signal in the Rayleigh fit in km. If use_range is set to True, the limits correspond to distance. If not provided, the normalization region will be automatically identified (default by CARS) Example: normalization_region = 4., 6.
normalization_region =

#fit_mask_region:  Lower and upper thresholds for the region where the Rayleigh fit mask will be applied in km. The cross criterion might might fail if the lower limit is below the distance of full overlap. Defaults to: 3 to 34 km Example: fit_mask_region = 4., 20. 
fit_mask_region = 

#fit_mask_window: The size limits (min and max) of the window used for Rayleigh fit mask. The mask will be applied for window sizes ranging between the two provided values. Defaults to: 2., 6. km Example: fit_mask_window = 0.5, 4. 
fit_mask_window =

#fit_mask_window_step: The center of the window used for the Rayleigh fit mask will range between the fit_mask_window with a step equal to the value provided here in km. Keep in mind that smaller values will increase the resolution of the mask but will make processing slower. Defaults to: 0.2 km. Example: fit_mask_window_step = 0.1
fit_mask_window_step =

#rsem_threshold: The relative standard error of the mean threshold that is applied by the Rayleigh fit mask. Regions with relative SEM above this threshold will not be considered molecular. Defaults to 0.2 (20%). Example: rsem_threshold = 0.02 
rsem_threshold =

#first_derivative_threshold: The standard deviation threshold applied by the Rayleigh fit mask. Regions where the first derivative is less than its uncertainty multiplied by the derivative_threshold will not be considered molecular. Defaults to 2. (2 * sigma). Example: first_derivative_threshold = 1.
first_derivative_threshold =

#second_derivative_threshold: The standard deviation threshold applied by the Rayleigh fit mask. Regions where the second derivative is less than its uncertainty multiplied by the derivative_threshold will not be considered molecular. Defaults to 2. (2 * sigma). Example: second_derivative_threshold = 1.
second_derivative_threshold =

#shapiro_wilk_threshold: The p value threshold for the Shapiro-Wilk test applied by the Rayleigh fit mask. Regions where the Shapiro-Wilk test returns p values smaller than this threshold will not be considered molecular because the noise does not follow a Gaussian distribution. Defaults to 0.05. Example: shapiro_wilk_threshold = 0.10 
shapiro_wilk_threshold =

#cross_criterion_threshold: The cross criterion threshold applied by the Rayleigh fit mask when checking if negative differences of the normalized rangecorrected signal and the molecular attenuated backscatter are significant. If even one region below the normalization region is found where the difference between the normalized rangecorected signal and the molecular profile is smaller than the negative uncertainty of the normalized rangecorrected signal multiplied by the cross_criterion_threshold means that the corresponding normalization region is not molecular. Defaults to 1. (1 * sigma). Example: cross_criterion_threshold = 1.
cross_criterion_threshold = 

#durbin_watson_threshold: Lower and upper thresholds for the Durbin-Watson test. Values close to 2. indicate no auto-correlation. Defaults to [1., 3.]. Example: cross_criterion_threshold = 1.5, 2.5
durbin_watson_threshold = 

#x_lims: Set the x axis (height/range) limits in km (lower and upper). Defaults to: 0., 20.
x_lims =

#x_tick: The x axis finest major tick in km. Defaults to: 2km 
x_tick =

#y_lims: The y axis (signal) limits (lower and upper) of the normalized rangecorrected signal in m-1 sr-1 (pseudo-units). It is recommended to skip this variable and use the automatic selection because channels in deifferent wavelengths have different optimal limits. Defaults to: automatic selection
y_lims =

#smooth: Refer to the smooth option in the quicklook section. Defaults to: False
smooth = 

#smoothing_exponential: Refer to the smooth option in the quicklook section. Defaults to: False.
smoothing_exponential =

#smoothing_range: Refer to the smooth option in the quicklook section. Defaults to: 0.05, 20.
smoothing_range = 

#smoothing_window: Refer to the smooth option in the quicklook section. Defaults to: smoothing_window = 100., 100.
smoothing_window =

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi =

#color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True
color_reduction = 

[telecover]
#channels: Provided it to process only selected channels. By default, no channel is excluded if this variable is not provided. The channel name must follow the nomenclature of section 4.3.
channels =

#exclude_telescope_type: Provide it to entirely exclude all channels of a specific telescope_type. By default, no channel is excluded if this variable is not provided. 
exclude_telescope_type =

#exclude_channel_type: Provide it to entirely exclude all channels of a specific channel_type. By default:
#    • if isday = False --> no channel is excluded if this variable is not provided
#    • if isday = True --> exclude_channel_type = v, f (vibrational Raman and fluorescence channels excluded)
exclude_channel_type =

#exclude_acquisition_mode: Provide it to entirely exclude all channels of a specific acquisition_mode. By default, no channel is excluded if this variable is not provided. 
exclude_acquisition_mode =

#exclude_channel_subtype: Provide it to entirely exclude all channels of a specific channel_subtype. By default, no channel is excluded if this variable is not provided. 
exclude_channel_subtype =

#use_range: If set to False, the x axis units of the Telecover test will correspond to the height in km above the lidar. If set to True, the x axis units corresponds to the range above the lidar. This variable determines the units of the following variables: 
#    • normalization_height
#    • x_lims 
#    • smoothing_range
#Defaults to: True
use_range = 

#x_lims: Set the x axis (height/range) limits for the near range subplots of the telecover test in km (lower and upper). Defaults to: 0., 2.5 for l, m, n, and x telescope types and: 0., 5. for the rest of the types 
x_lims =

#x_tick: The x axis finest tick for the near range subplots of the telecover test in km. Defaults to: 0.5 for l, m, n, and x telescope types and: 1. for the rest of the types 
x_tick =

#y_lims: The y axis (signal) limits (lower and upper) of the normalized rangecorrected signal in m-1 sr-1 (pseudo-units). It is recommended to skip this variable and use the automatic selection because channels in deifferent wavelengths have different optimal limits. Defaults to: automatic selection
y_lims =

#use_non_rangecor: If set to True, the non range corrected signals will be used for the left subplot of the telecover test. If set to False the range corrected signals will be used instead. Defaults to: False
use_non_rangecor = 

#use_last: While set to True an additional purple line will be added in all subplots of the telecover test. In the first 3 subplots the lines corresponds to the last sector (e.g. N2), while in the 4th subplot (deviations) it is the difference between the normalized signal of the last sector and the normalized mean signal of the correspondin sector e.g. N2 – N. Set to False to not visualize this line. Defaults to True.
use_last =

#auto_fit: If set to True an automatic identification of the normalization region will be attempted. If the automatic procedure is successful, the normalization_region variable will be ignored. If the procedure is not successful or auto_fit is set to False, the manually-provided/default values will be used. Defaults to True
auto_fit = 

#normalization_region: The reference height/range in km where the signals will be normalized for the Telecover test (relevant to middle and right subplots). It should be well above the suspected full overlap height. Defaults to: normalization_region = 1.8, 2.2
normalization_region =

#smooth: Refer to the smooth option in the quicklook section. Defaults to: True
smooth = 

#smoothing_exponential: Refer to the smooth option in the quicklook section. Defaults to: False.
smoothing_exponential =

#smoothing_window: Refer to the smooth option in the quicklook section. Defaults to: smoothing_window = 100., 100.
smoothing_window =

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi =

#color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True
color_reduction = 

[polarization_calibration]
#ch_r: Provide here channel names from the available ones with the reflected channel_subtype (e.g. 0355xpar). The number of reflected channels must be the same as the number of the respective transmitted channels provided by ch_t and by the GHK parameters. By default, all available reflected channels will be matched to all available transmitted channels that share the same telescope type, detection mode, and wavelength. WARNING! The field is mandatory if non-default GHK values are applied. 
ch_r = 
 
#ch_t: Provide here channel names from the available ones with the transmitted channel_subtype (e.g. 0355xcat). The number of reflected channels must be the same as the number of the respective transmitted channels provided by ch_r and by the GHK parameters. By default, all available reflected channels will be matched to all available transmitted channels that share the same telescope type, detection mode, and wavelength. WARNING! The field is mandatory if non-default GHK values are applied. 
ch_t =

#auto_fit: If set to True an automatic identification of the calibration and Rayleigh regions will be attempted. If the automatic procedure is successful, the calibration_region and rayleigh_region variables will be ignored. If the procedure is not successful or auto_fit is set to False, the manually-provided/default values will be used. Defaults to True
auto_fit = 

#calibration_region: The lower and upper limits of the region used for Δ90 calibration. If use_range is called, the limits correspond to distance. If auto_fit is set to True and the automatic identification is successful for a specific channel combination, the calibration_region variable will be ignored even if provided. Defaults to: 2., 4. km 
calibration_region = 

#rayleigh_region: The lower and upper limits of the region used for the comparison with the Rayleigh atmosphere. If auto_fit is set to True and the automatic identification is successful for a specific channel combination, the rayleigh_region variable will be ignored even if provided. Defaults to: 8.5, 9.5
rayleigh_region = 

#x_lims_calibration: The x axis (height/range) limits in km (lower and upper) for the pol. calibration subplot. Defaults to: 0, 15
x_lims_calibration =

#x_lims_rayleigh: The x axis (height/range) limits in km (lower and upper) for the Rayleigh comparison subplot. Defaults to: 0, 15
x_lims_rayleigh =

#x_tick_calibration: The x axis finest major tick in km for the pol. calibration sunplot. Defaults to: 0.5
x_tick_calibration =

#x_tick_rayleigh: The x axis finest major tick in Km for the Rayleigh plot. Defaults to: 2.
x_tick_rayleigh =

#y_lims_calibration: The y axis limits (lower and upper) of the gain ratios at +-45. Used for the pol. calibration subplot. Defaults to: automatic selection
y_lims_calibration =

#y_lims_rayleigh: The y axis limits (lower and upper) of the volume depolarization ratio. Used for the Rayleigh subplot. Defaults to: automatic selection
y_lims_rayleigh = 

#K: The K value for each channel pair. Defaults to: 1 for all channels 
K =

#G_R: The G value for the reflected channel of the pair. Defaults to: 1 for all channels (assuming no correction for the emission and the receiver)
G_R =

#G_T: The G value for the transmitted channel of the pair. Defaults to: 1 for all channels (assuming no correction for the emission and the receiver)
G_T =

#H_R: The H value for the reflected channel of the pair. Defaults to: 
#    • 1 for all co-polar (p) reflected channels
#    • -1 for all  cross-polar (c) reflected channels
#(assuming no correction for the emission and the receiver)
H_R =

#H_T: The H value for the transmitted channel of the pair. Defaults to: 1 or -1 for all co-polar (p) and cross-polar (c) transmitted channels, respectively (no receiver optics + emitted pcb. state correction)
#    • 1 for all co-polar (p) transmitted channels
#    • -1 for all  cross-polar (c) transmitted channels
#(assuming no correction for the emission and the receiver)
H_T = 

#R_to_T_transmission_ratio: The transmission ratio of the R to T channels per pair (TR/TT). TR = 1 and/or TT = 1 if no filter was applied during calibration. Defaults to: 1 for all pairs
R_to_T_transmission_ratio = 

#smooth: If set to True, a sliding average smoothing filter will be applied on the signals across y axis for better visualization. Defaults to: True
smooth =

#smoothing_exponential: This variable is ignored if smooth = False and also if the upper and lower values of the smoothing_window variable are the same. Choose one of:
#    • True: a smoothing window that exponentially increases with height/range will be applied. 
#    • False:  a smoothing window that exponentially increases with height/range will be applied
#Defaults to: False.
smoothing_exponential =

#smoothing_range: Set the first and last height/range boundaries in km where smoothing should be applied. If they exceed the actual signal boundaries the actual boundaries will be used instead. Defaults to: 0.25, 15.
smoothing_range = 

#smoothing_window: Half smoothing window in the first and last bin of the smoothing region, in m. The widow progressively changes between from the first to the last value. Use the a single value twice to apply a constant window. Defaults to: smoothing_window = 500.
smoothing_window =

#dpi: The dots per inch (dpi) resolution of the exported figures. Defaults to: 300 dpi 
dpi =

#color_reduction: If set to True, an adaptive color reduction will be applied to save space when exporting the figures. Defaults to True
color_reduction =
```
