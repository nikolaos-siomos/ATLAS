# ATLAS intercomparison initialization reference

!!! note "Generated reference"
    This page is generated from the intercomparison parser schemas and flavor text.

## Structure

The file contains one `[general]` section, repeated `[dataset:<id>]` sections, and repeated `[channel_group:<id>]` and `[pair_group:<id>]` sections.

Dataset sections define the exported stage, QA test, source parameters, and labels for each independently selectable comparison dataset. Different datasets may use different stages, QA tests, or source parameters, even when they originate from the same physical lidar.

Exactly one dataset must set `reference = True`. Group sections map participating datasets to one `atlas_channel_id` or `atlas_pair_id`. Empty or omitted mappings exclude that dataset. A group with no provided IDs is ignored. For an active group, the reference dataset must provide its ID.

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
| `vertical_scale` | Vertical coordinate to use later for harmonization and plotting. Allowed values are bins, range, height_agl, and height_asl. All four coordinates are loaded when available; this option only selects which one downstream processing will use. Default: height_asl. | `str` | `height_asl` | bins, range, height_agl, height_asl | `height_asl` |
| `vertical_method` | Method used to place datasets on a common physical vertical grid. interpolation uses the reference-dataset grid; vertical_binning creates common altitude intervals. | `str` | `interpolation` | interpolation, vertical_binning | `interpolation` |
| `vertical_binning` | Vertical bin width in kilometres. Required only when vertical_method is vertical_binning. | `float` | `` | &gt; 0.0 | `0.03` |
| `vertical_min` | Optional lower comparison and plotting limit in kilometres. Empty uses the common valid overlap. | `float` | `` |  | `0.5` |
| `vertical_max` | Optional upper comparison and plotting limit in kilometres. Empty uses the common valid overlap. | `float` | `` |  | `15.0` |
| `slice_measurement` | Optional temporal slices as repeating start, stop pairs. Accepted formats match call_atlas.ini: HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, or yyyymmdd_HHMMSS. | `list[str]` | `` |  | `20260801_2100, 20260802_0200` |
| `exclude_measurement` | Optional temporal exclusions as repeating start, stop pairs, using the same time formats as slice_measurement. | `list[str]` | `` |  | `20260801_2330, 20260801_2345` |
| `default_channel_background_correction` | Default background-correction switch for channel groups. A value set directly in a [channel_group:&lt;id&gt;] section overrides this default. | `bool` | `False` |  | `False` |
| `default_channel_background_region` | Default channel background region in kilometres. A channel-group background_region overrides it. | `list[float]` | `18.0, 22.0` | size=2 | `18.0, 22.0` |
| `default_channel_normalisation` | Default normalization switch for channel groups. A channel-group normalisation value overrides it. | `bool` | `True` |  | `True` |
| `default_channel_normalisation_region` | Default channel normalization region in kilometres. A channel-group normalisation_region overrides it. | `list[float]` | `7.5, 9.0` | size=2 | `7.5, 9.0` |
| `default_channel_normalise_to_molecular` | Default channel normalization target. If True, channel groups are normalized to the reference dataset&#x27;s molecular profile; if False, they are normalized to the reference measured signal. | `bool` | `True` |  | `True` |
| `default_channel_plot_molecular` | Default switch controlling whether the reference molecular channel profile is included in plots. | `bool` | `True` |  | `True` |
| `default_pair_background_correction` | Default background-correction switch for pair groups. A value set directly in a [pair_group:&lt;id&gt;] section overrides this default. | `bool` | `False` |  | `False` |
| `default_pair_background_region` | Default pair background region in kilometres. A pair-group background_region overrides it. | `list[float]` | `18.0, 22.0` | size=2 | `18.0, 22.0` |
| `default_pair_normalisation` | Default normalization switch for pair groups. Pair normalization is always to the reference measured pair ratio, never to the molecular ratio. | `bool` | `False` |  | `False` |
| `default_pair_normalisation_region` | Default pair normalization region in kilometres. A pair-group normalisation_region overrides it. | `list[float]` | `7.5, 9.0` | size=2 | `7.5, 9.0` |
| `default_pair_plot_molecular` | Default switch controlling whether the reference molecular ratio is included in pair plots. | `bool` | `True` |  | `True` |
| `dpi` | Resolution of exported figures in dots per inch. | `int` | `150` | &gt; 1 | `150` |
| `color_reduction` | If True, apply the ATLAS image color-reduction workflow to exported figures. | `bool` | `False` |  | `False` |

## `dataset:<dataset_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `stage_path` | Absolute or relative path to one exported ATLAS stage directory. Different datasets may point to different stage paths, or multiple datasets may intentionally point to the same stage path with different QA/source selections. | `Path` | `` |  | `../dataset_a/exported/preprocessing_complete` |
| `reference` | Set True for exactly one dataset. Its vertical grid, molecular profile, product IDs, and metadata defaults define the comparison reference. | `bool` | `False` |  | `True` |
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
| `label` | Optional channel-group display label. When empty, the reference dataset&#x27;s atlas_channel_id is used. | `str` | `` |  | `355 nm parallel` |
| `background_correction` | Optional channel-group override for background correction. Empty inherits default_channel_background_correction from [general]. | `bool` | `` |  | `False` |
| `background_region` | Optional channel-group override for the background interval in kilometres. Empty inherits default_channel_background_region from [general]. | `list[float]` | `` | size=2 | `18.0, 22.0` |
| `normalisation` | Optional channel-group override for normalization. Empty inherits default_channel_normalisation from [general]. | `bool` | `` |  | `True` |
| `normalisation_region` | Optional channel-group override for the normalization interval in kilometres. Empty inherits default_channel_normalisation_region from [general]. | `list[float]` | `` | size=2 | `7.5, 9.0` |
| `normalise_to_molecular` | Optional channel-group override. True normalizes every participating channel to the reference dataset&#x27;s molecular profile; False normalizes every channel to the reference measured signal. | `bool` | `` |  | `True` |
| `plot_molecular` | Optional channel-group override controlling whether the reference molecular profile is plotted. Empty inherits default_channel_plot_molecular from [general]. | `bool` | `` |  | `True` |

### Dynamic channel-group entries

- `<dataset_id>.atlas_channel_id = <atlas_channel_id>` selects one channel from that dataset. Empty or omitted mappings exclude the dataset. If no dataset provides an ID, the group is ignored. For an active group, the reference dataset must provide an ID. Example: `0355xpgx`.

## `pair_group:<group_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `label` | Optional pair-group display label. When empty, the reference dataset&#x27;s atlas_pair_id is used. | `str` | `` |  | `VLDR 355 nm` |
| `background_correction` | Optional pair-group override for background correction. Empty inherits default_pair_background_correction from [general]. | `bool` | `` |  | `False` |
| `background_region` | Optional pair-group override for the background interval in kilometres. Empty inherits default_pair_background_region from [general]. | `list[float]` | `` | size=2 | `18.0, 22.0` |
| `normalisation` | Optional pair-group override for normalization. Pair products are normalized to the reference measured pair ratio; molecular normalization is not used for pairs. | `bool` | `` |  | `False` |
| `normalisation_region` | Optional pair-group override for the normalization interval in kilometres. Empty inherits default_pair_normalisation_region from [general]. | `list[float]` | `` | size=2 | `7.5, 9.0` |
| `plot_molecular` | Optional pair-group override controlling whether the reference molecular ratio is plotted. Empty inherits default_pair_plot_molecular from [general]. | `bool` | `` |  | `True` |

### Dynamic pair-group entries

- `<dataset_id>.atlas_pair_id = <atlas_pair_id>` selects one pair from that dataset. Empty or omitted mappings exclude the dataset. If no dataset provides an ID, the group is ignored. For an active group, the reference dataset must provide an ID. Example: `0355UVAX`.

## Deferred metadata defaults

Dataset labels may remain empty after parsing. Background, normalization, and molecular-plot controls are resolved per group: an explicit group value overrides the corresponding channel/pair default from [general].
