# ATLAS intercomparison initialization reference

!!! note "Generated reference"
    This page is generated from the intercomparison parser schemas and flavor text.

## Structure

The file contains one `[general]` section, one optional `[plotting]` section, repeated `[dataset:<id>]` sections, and repeated `[channel_group:<id>]` and `[pair_group:<id>]` sections.

Dataset sections define where and how exported data are read. Group sections do not map datasets directly anymore; instead they contain **named entries**. Each named entry points to a dataset and selects one ATLAS channel or pair. Several entries may point to the same dataset, so multiple channels or pairs from one dataset can be compared in the same plot.

Exactly one dataset must still set `reference = True`. In every active channel or pair group, exactly one named entry must also set `<entry_id>.reference = True`, and that entry must point to the global reference dataset. This identifies the exact channel/pair used as the normalization and difference reference.

### Named-entry syntax

```ini
[channel_group:355_comparison]
a_parallel.dataset = dataset_a_reference
a_parallel.atlas_channel_id = 0355xpgx
a_parallel.label = A parallel
a_parallel.reference = True

a_analog.dataset = dataset_a_reference
a_analog.atlas_channel_id = 0355xagx
a_analog.label = A analog

b_parallel.dataset = dataset_b
b_parallel.atlas_channel_id = 0355xpgx
b_parallel.label = B parallel
```

The prefix (`a_parallel`, `a_analog`, `b_parallel`) is the entry ID. Entry IDs only need to be unique within their group.

## `general`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `output_folder` | Folder where intercomparison plots, tables, and cached products are written. Relative paths are resolved against the folder containing this INI file. | `Path` | `analysis` |  | `./analysis` |
| `overwrite_output` | If True, existing intercomparison outputs with the same names may be overwritten. | `bool` | `False` |  | `False` |
| `default_qa_test` | Universal fallback QA test. A dataset-specific qa_test takes priority. | `str` | `ray` | ray, drk, pcb, tlc, tlc_rin, ray_pcb, pcb_aux, trg, dtm, cam | `ray` |
| `default_signal_source` | Universal fallback exported-stage parameter containing channel signals. A dataset-specific signal_source takes priority. | `str` | `profile` |  | `profile` |
| `default_signal_error_source` | Universal fallback exported-stage parameter containing channel uncertainties. A dataset-specific signal_error_source takes priority. | `str` | `profile_error` |  | `profile_error` |
| `default_pair_source` | Universal fallback exported-stage parameter containing channel-pair products. A dataset-specific pair_source takes priority. | `str` | `pol_cal_ratio_mean` |  | `pol_cal_ratio_mean` |
| `default_pair_error_source` | Universal fallback exported-stage parameter containing channel-pair uncertainties. A dataset-specific pair_error_source takes priority. | `str` | `pol_cal_ratio_error_mean` |  | `pol_cal_ratio_error_mean` |
| `vertical_scale` | Vertical coordinate used by vertical harmonization and later plotting. Allowed values are bins, range, height_agl, and height_asl. All four coordinates are loaded when available. Default: height_asl. | `str` | `height_asl` | bins, range, height_agl, height_asl | `height_asl` |
| `vertical_method` | Vertical harmonization method. interpolation maps every participating entry to the coarsest native grid in that group so no entry is upscaled. vertical_binning applies conservative overlap-weighted binning onto a common regular grid. Default: interpolation. | `str` | `interpolation` | interpolation, vertical_binning | `interpolation` |
| `vertical_bin_width` | Default conservative vertical-bin width used when vertical_method=vertical_binning. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Default: 0.1 km for physical scales. A group-specific vertical_bin_width overrides this value. If the requested width is smaller than the coarsest nominal native step in a group, ATLAS warns and uses that coarsest native step instead. | `float` | `0.1` | &gt; 0.0 | `0.1` |
| `first_bin_left_edge` | Left edge of the first bin in the common conservative grid. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Default: 0. Each output vertical-scale value is assigned to the center of its bin. Leading bins with no source overlap are retained and filled with NaN. | `float` | `0.0` |  | `0.0` |
| `vertical_min` | Optional lower retained limit for the common binned grid. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Empty starts from first_bin_left_edge. | `float` | `` |  | `0.5` |
| `vertical_max` | Optional upper retained limit for the common binned grid. For physical vertical scales it is specified in kilometres; for bins it is in bin units. Empty extends to the largest available vertical extent in the group. | `float` | `` |  | `15.0` |
| `plot_native_scale` | If True, vertical harmonization is skipped for plotting and every dataset is shown on its native selected vertical scale. The right-hand relative-difference panel is left without curves because the vertical samples are not aligned. Default: False. | `bool` | `False` |  | `False` |
| `slice_measurement` | Optional temporal slices as repeating start, stop pairs. Accepted formats match call_atlas.ini: HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, or yyyymmdd_HHMMSS. | `list[str]` | `` |  | `20260801_2100, 20260802_0200` |
| `exclude_measurement` | Optional temporal exclusions as repeating start, stop pairs, using the same time formats as slice_measurement. | `list[str]` | `` |  | `20260801_2330, 20260801_2345` |
| `default_channel_background_correction` | Default background-correction switch for channel groups. A value set directly in a [channel_group:&lt;id&gt;] section overrides this default. | `bool` | `False` |  | `False` |
| `default_channel_background_region` | Default channel background region in kilometres. A channel-group background_region overrides it. | `list[float]` | `18.0, 22.0` | size=2 | `18.0, 22.0` |
| `default_channel_normalisation` | Default normalization switch for channel groups. A channel-group normalisation value overrides it. | `bool` | `True` |  | `True` |
| `default_channel_normalisation_region` | Default channel normalization region in kilometres. A channel-group normalisation_region overrides it. | `list[float]` | `4.0, 6.0` | size=2 | `4.0, 6.0` |
| `default_channel_normalise_to_molecular` | Default channel normalization target. If True, channel groups are normalized to the reference entry&#x27;s molecular profile; if False, they are normalized to the reference entry&#x27;s measured signal. | `bool` | `True` |  | `True` |
| `default_channel_plot_molecular` | Default switch controlling whether the reference molecular channel profile is included in plots. | `bool` | `True` |  | `True` |
| `default_pair_background_correction` | Default background-correction switch for pair groups. A value set directly in a [pair_group:&lt;id&gt;] section overrides this default. | `bool` | `False` |  | `False` |
| `default_pair_background_region` | Default pair background region in kilometres. A pair-group background_region overrides it. | `list[float]` | `18.0, 22.0` | size=2 | `18.0, 22.0` |
| `default_pair_normalisation` | Default normalization switch for pair groups. Pair normalization is always to the reference entry&#x27;s measured pair ratio, never to the molecular ratio. | `bool` | `False` |  | `False` |
| `default_pair_normalisation_region` | Default pair normalization region in kilometres. A pair-group normalisation_region overrides it. | `list[float]` | `7.5, 9.0` | size=2 | `7.5, 9.0` |
| `default_pair_plot_molecular` | Default switch controlling whether the reference molecular ratio is included in pair plots. | `bool` | `True` |  | `True` |
| `dpi` | Resolution of exported figures in dots per inch. | `int` | `150` | &gt; 1 | `150` |
| `color_reduction` | If True, apply the ATLAS image color-reduction workflow to exported figures. | `bool` | `False` |  | `False` |

## `plotting`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `channel_x_lims` | Default horizontal limits for channel-group intercomparison plots. Empty means determine them automatically from the plotted vertical data. Units follow vertical_scale: km for physical scales, bin units for bins. | `list[float]` | `` | size=2 | `0.0, 20.0` |
| `pair_x_lims` | Default horizontal limits for pair-group intercomparison plots. Empty means determine them automatically from the plotted vertical data. Units follow vertical_scale: km for physical scales, bin units for bins. | `list[float]` | `` | size=2 | `0.0, 10.0` |
| `x_tick` | Default major horizontal-axis tick spacing for intercomparison plots. Units follow vertical_scale. | `float` | `2.0` | &gt; 0.0 | `1.0` |
| `channel_difference_y_lims` | Default right-panel y-axis limits for channel-group relative differences. Group-specific difference_y_lims can override these values. | `list[float]` | `-0.4, 0.4` | size=2 | `-0.4, 0.4` |
| `pair_difference_y_lims` | Default right-panel y-axis limits for pair-group absolute differences. Empty means determine them automatically from the SNR-filtered absolute differences. Group-specific difference_y_lims can override these values. | `list[float]` | `` | size=2 | `-0.1, 0.1` |
| `channel_y_lims` | Default left-panel y limits for channel-group plots. Empty means determine them automatically from all plotted channel signals and, when enabled, the reference molecular profile. | `list[float]` | `` | size=2 |  |
| `channel_smooth` | Default smoothing switch for channel-group plots. Smoothing/local-STD estimation is used only for data that are not conservatively binned. Harmonized vertical_binning plots use the binned signal and propagated error directly. | `bool` | `True` |  | `True` |
| `channel_smoothing_range` | Default channel smoothing range. Units follow vertical_scale: km for physical scales, bin units for bins. | `list[float]` | `0.05, 35.0` | size=2 | `0.05, 35.0` |
| `channel_smoothing_window` | Default channel smoothing window. Units follow vertical_scale: km for physical scales, bin units for bins. | `float` | `0.5` | &gt; 0.0 | `0.5` |
| `pair_y_lims` | Default left-panel y limits for pair-group plots. Empty means determine them automatically. | `list[float]` | `` | size=2 |  |
| `pair_smooth` | Default smoothing switch for pair-group plots. Smoothing/local-STD estimation is used only for data that are not conservatively binned. Harmonized vertical_binning plots use the binned values and propagated errors directly. | `bool` | `True` |  | `True` |
| `pair_smoothing_range` | Default pair smoothing range. Units follow vertical_scale: km for physical scales, bin units for bins. | `list[float]` | `0.05, 10.0` | size=2 | `0.05, 10.0` |
| `pair_smoothing_window` | Default pair smoothing window. Units follow vertical_scale: km for physical scales, bin units for bins. | `float` | `0.5` | &gt; 0.0 | `0.5` |

## `dataset:<dataset_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `stage_path` | Absolute or relative path to one exported ATLAS stage directory. Different datasets may point to different stage paths, or multiple datasets may intentionally point to the same stage path with different QA/source selections. | `Path` | `` |  | `../dataset_a/exported/preprocessing_complete` |
| `reference` | Set True for exactly one dataset. This defines the global reference dataset; every active comparison group must choose exactly one reference entry that points to it. | `bool` | `False` |  | `True` |
| `system_label` | Optional label for the physical lidar/system that produced this dataset. Multiple datasets may share the same system_label. | `str` | `` |  | `Lidar A` |
| `dataset_label` | Optional label describing this particular dataset or processing realization. When empty, later code may fall back to the dataset section ID. | `str` | `` |  | `Rayleigh processing` |
| `qa_test` | QA test used when reading this dataset. When empty, default_qa_test from [general] is used. | `str` | `` | ray, drk, pcb, tlc, tlc_rin, ray_pcb, pcb_aux, trg, dtm, cam | `ray` |
| `signal_source` | Exported-stage parameter containing channel signals for this dataset. When empty, default_signal_source is used. | `str` | `` |  | `profile` |
| `signal_error_source` | Exported-stage parameter containing channel uncertainties for this dataset. When empty, default_signal_error_source is used. | `str` | `` |  | `profile_error` |
| `pair_source` | Exported-stage parameter containing pair products for this dataset. When empty, default_pair_source is used. | `str` | `` |  | `pol_cal_ratio_mean` |
| `pair_error_source` | Exported-stage parameter containing pair uncertainties for this dataset. When empty, default_pair_error_source is used. | `str` | `` |  | `pol_cal_ratio_error_mean` |

## `channel_group:<group_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `label` | Optional channel-group display label. When empty, the reference entry&#x27;s atlas_channel_id is used. | `str` | `` |  | `355 nm parallel` |
| `background_correction` | Optional channel-group override for background correction. Empty inherits default_channel_background_correction from [general]. | `bool` | `` |  | `False` |
| `background_region` | Optional channel-group override for the background interval in kilometres. Empty inherits default_channel_background_region from [general]. | `list[float]` | `` | size=2 | `18.0, 22.0` |
| `normalisation` | Optional channel-group override for normalization. Empty inherits default_channel_normalisation from [general]. | `bool` | `` |  | `True` |
| `normalisation_region` | Optional channel-group override for the normalization interval in kilometres. Empty inherits default_channel_normalisation_region from [general]. | `list[float]` | `` | size=2 | `4.0, 6.0` |
| `normalise_to_molecular` | Optional channel-group override. True normalizes every participating channel to the reference entry&#x27;s molecular profile; False normalizes every channel to the reference entry&#x27;s measured signal. | `bool` | `` |  | `True` |
| `plot_molecular` | Optional channel-group override controlling whether the reference molecular profile is plotted. Empty inherits default_channel_plot_molecular from [general]. | `bool` | `` |  | `True` |
| `vertical_bin_width` | Optional channel-group conservative bin-width override. Empty uses [general] vertical_bin_width. If the requested width is smaller than the coarsest nominal native step in this group, ATLAS warns and uses the coarsest native step instead. | `float` | `` | &gt; 0.0 | `0.03` |
| `smooth` | Optional channel-group plotting smoothing override. Empty inherits channel_smooth from [plotting]. | `bool` | `` |  | `True` |
| `smoothing_range` | Optional channel-group smoothing-range override. Empty inherits channel_smoothing_range from [plotting]. | `list[float]` | `` | size=2 | `0.05, 35.0` |
| `smoothing_window` | Optional channel-group smoothing-window override. Empty inherits channel_smoothing_window from [plotting]. | `float` | `` | &gt; 0.0 | `0.5` |
| `x_lims` | Optional channel-group x-axis limits. Empty inherits channel_x_lims from [plotting]; an empty resolved value triggers automatic limits. | `list[float]` | `` | size=2 | `0.0, 20.0` |
| `x_tick` | Optional channel-group major x-axis tick spacing. Empty inherits x_tick from [plotting]. | `float` | `` | &gt; 0.0 | `2.0` |
| `y_lims` | Optional channel-group left-panel y limits. Empty inherits channel_y_lims from [plotting]; an empty resolved value triggers automatic limits. | `list[float]` | `` | size=2 |  |
| `difference_y_lims` | Optional channel-group right-panel y-axis limits for relative differences. Empty inherits channel_difference_y_lims from [plotting]. | `list[float]` | `` | size=2 | `-0.4, 0.4` |
| `use_log_y_scale` | Channel-group logarithmic-y switch. Default: True. Set False in this group to use a linear left-panel y axis. | `bool` | `True` |  | `True` |

### Dynamic named channel entries

- `<entry_id>.dataset = <dataset_id>` selects the source dataset.
- `<entry_id>.atlas_channel_id = <atlas_channel_id>` selects one channel from that dataset.
- `<entry_id>.label = ...` is optional and controls the legend label.
- `<entry_id>.reference = True` must appear on exactly one entry in each active group.

## `pair_group:<group_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `label` | Optional pair-group display label. When empty, the reference entry&#x27;s atlas_pair_id is used. | `str` | `` |  | `VLDR 355 nm` |
| `background_correction` | Optional pair-group override for background correction. Empty inherits default_pair_background_correction from [general]. | `bool` | `` |  | `False` |
| `background_region` | Optional pair-group override for the background interval in kilometres. Empty inherits default_pair_background_region from [general]. | `list[float]` | `` | size=2 | `18.0, 22.0` |
| `normalisation` | Optional pair-group override for normalization. Pair products are normalized to the reference entry&#x27;s measured pair ratio; molecular normalization is not used for pairs. | `bool` | `` |  | `False` |
| `normalisation_region` | Optional pair-group override for the normalization interval in kilometres. Empty inherits default_pair_normalisation_region from [general]. | `list[float]` | `` | size=2 | `7.5, 9.0` |
| `plot_molecular` | Optional pair-group override controlling whether the reference molecular ratio is plotted. Empty inherits default_pair_plot_molecular from [general]. | `bool` | `` |  | `True` |
| `vertical_bin_width` | Optional pair-group conservative bin-width override. Empty uses [general] vertical_bin_width. If the requested width is smaller than the coarsest nominal native step in this group, ATLAS warns and uses the coarsest native step instead. | `float` | `` | &gt; 0.0 | `0.03` |
| `smooth` | Optional pair-group plotting smoothing override. Empty inherits pair_smooth from [plotting]. | `bool` | `` |  | `True` |
| `smoothing_range` | Optional pair-group smoothing-range override. Empty inherits pair_smoothing_range from [plotting]. | `list[float]` | `` | size=2 | `0.05, 10.0` |
| `smoothing_window` | Optional pair-group smoothing-window override. Empty inherits pair_smoothing_window from [plotting]. | `float` | `` | &gt; 0.0 | `0.5` |
| `x_lims` | Optional pair-group x-axis limits. Empty inherits pair_x_lims from [plotting]; an empty resolved value triggers automatic limits. | `list[float]` | `` | size=2 | `0.0, 10.0` |
| `x_tick` | Optional pair-group major x-axis tick spacing. Empty inherits x_tick from [plotting]. | `float` | `` | &gt; 0.0 | `1.0` |
| `y_lims` | Optional pair-group left-panel y limits. Empty inherits pair_y_lims from [plotting]; an empty resolved value triggers automatic limits. | `list[float]` | `` | size=2 |  |
| `difference_y_lims` | Optional pair-group right-panel y-axis limits for absolute differences. Empty inherits pair_difference_y_lims from [plotting]; if both are empty, limits are calculated automatically. | `list[float]` | `` | size=2 | `-0.1, 0.1` |
| `use_log_y_scale` | Pair-group logarithmic-y switch. Default: False. Set True in this group only when a logarithmic left-panel y axis is explicitly desired. | `bool` | `False` |  | `False` |

### Dynamic named pair entries

- `<entry_id>.dataset = <dataset_id>` selects the source dataset.
- `<entry_id>.atlas_pair_id = <atlas_pair_id>` selects one pair from that dataset.
- `<entry_id>.label = ...` is optional and controls the legend label.
- `<entry_id>.reference = True` must appear on exactly one entry in each active group.

## Internal group representation

Parsed and prepared groups are keyed by `entries`, not by `datasets`. Harmonized group arrays use an `entry` dimension. Each prepared entry retains its `dataset_id`, so dataset-level metadata and stage caching remain shared even when several entries point to the same dataset.

## Deferred defaults

Entry labels may remain empty after parsing. Plotting then falls back to the dataset label/system label plus the selected ATLAS product ID. Background, normalization, and molecular-plot controls are resolved per group from `[general]`. Plotting controls are resolved per group from `[plotting]`, except `use_log_y_scale`, whose group-level defaults are True for channel groups and False for pair groups.
