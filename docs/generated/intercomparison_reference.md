# ATLAS intercomparison initialization reference

!!! note "Generated reference"
    This page is generated from the intercomparison parser schemas and flavor text.

## Structure

The file contains one `[general]` section, repeated `[system:<id>]` sections, and repeated `[channel:<id>]` and `[pair:<id>]` comparison sections.

Only one system may set `reference = True`. In every channel and pair section, the mapping for that reference system is mandatory. Missing non-reference system mappings inherit the reference system's data ID; `off` excludes a non-reference system.

System-specific source overrides use dotted names such as `system_b.signal_source = derived_profile`. A missing source inherits the reference system's source, which in turn falls back to the corresponding general default.

## `general`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `output_folder` | Folder where intercomparison plots, tables, and cached products are written. Relative paths are resolved against the folder containing this INI file. | `Path` | `analysis` |  | `./analysis` |
| `overwrite_output` | If True, existing intercomparison outputs with the same names may be overwritten. | `bool` | `False` |  | `False` |
| `default_qa_test` | General QA-test source used by channel and pair comparisons unless force_qa_test is provided in that comparison section. | `str` | `ray` | Allowed values: ray, drk, pcb, tlc, tlc_rin, ray_pcb, pcb_aux, trg, dtm, cam | `ray` |
| `default_signal_source` | Fallback exported-stage parameter containing channel signals. A source specified for the reference system in a channel bundle becomes the effective default for all systems in that bundle. | `str` | `profile` |  | `profile` |
| `default_signal_error_source` | Fallback exported-stage parameter containing channel uncertainties. A source specified for the reference system in a channel bundle becomes the effective default for all systems in that bundle. | `str` | `profile_error` |  | `profile_error` |
| `default_pair_source` | Fallback exported-stage parameter containing channel-pair products. A source specified for the reference system in a pair bundle becomes the effective default for all systems in that bundle. | `str` | `pol_cal_ratio_mean` |  | `pol_cal_ratio_mean` |
| `default_pair_error_source` | Fallback exported-stage parameter containing channel-pair uncertainties. A source specified for the reference system in a pair bundle becomes the effective default for all systems in that bundle. | `str` | `pol_cal_ratio_error_mean` |  | `pol_cal_ratio_error_mean` |
| `vertical_scale` | Physical vertical coordinate used for all comparisons. height_asl is strongly recommended because systems may be located at different station altitudes. | `str` | `height_asl` | Allowed values: height_asl, height_agl, range | `height_asl` |
| `vertical_method` | Method used to place systems on a common physical vertical grid. interpolation uses the reference-system grid; vertical_binning creates common altitude intervals. | `str` | `interpolation` | Allowed values: interpolation, vertical_binning | `interpolation` |
| `vertical_binning` | Vertical bin width in kilometres. Required only when vertical_method is vertical_binning. | `float` | `empty parameter --&gt; ignored or automatic selection` | Minimum: 0.0 | `0.03` |
| `vertical_min` | Optional lower comparison and plotting limit in kilometres. Empty uses the common valid overlap. | `float` | `empty parameter --&gt; ignored or automatic selection` |  | `0.5` |
| `vertical_max` | Optional upper comparison and plotting limit in kilometres. Empty uses the common valid overlap. | `float` | `empty parameter --&gt; ignored or automatic selection` |  | `15.0` |
| `slice_measurement` | Optional temporal slices as repeating start, stop pairs. Accepted formats match call_atlas.ini: HHMM, yyyymmdd, yyyymmdd_HH, yyyymmdd_HHMM, or yyyymmdd_HHMMSS. | `list[str]` | `empty parameter --&gt; ignored or automatic selection` |  | `20260801_2100, 20260802_0200` |
| `exclude_measurement` | Optional temporal exclusions as repeating start, stop pairs, using the same time formats as slice_measurement. | `list[str]` | `empty parameter --&gt; ignored or automatic selection` |  | `20260801_2330, 20260801_2345` |
| `background_correction` | General switch for subtracting the mean signal in background_region. Disabled by default and overridable per channel comparison. | `bool` | `False` |  | `False` |
| `background_region` | General fallback background interval in kilometres. A channel-specific value takes priority; when empty, later metadata loading may use the reference channel&#x27;s stored region. | `list[float]` | `empty parameter --&gt; ignored or automatic selection` | Required list size: 2 | `18.0, 22.0` |
| `normalisation` | General switch for channel normalization. It may be disabled or overridden per channel comparison. Pair products are not normalized. | `bool` | `True` |  | `True` |
| `normalisation_region` | General fallback normalization interval in kilometres. A channel-specific value takes priority; when empty, later metadata loading uses the reference channel&#x27;s stored Rayleigh-fit region. | `list[float]` | `empty parameter --&gt; ignored or automatic selection` | Required list size: 2 | `7.5, 9.0` |
| `normalise_to_molecular` | If True, normalize measured channel profiles to the reference system&#x27;s molecular profile using arithmetic means over the resolved normalization region. | `bool` | `True` |  | `True` |
| `plot_molecular` | If True, include the reference system&#x27;s molecular profile in channel-comparison plots. | `bool` | `True` |  | `True` |
| `dpi` | Resolution of exported figures in dots per inch. | `int` | `150` | Minimum: 1 | `150` |
| `color_reduction` | If True, apply the ATLAS image color-reduction workflow to exported figures. | `bool` | `False` |  | `False` |

## `system:<system_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `stage_path` | Absolute or relative path to one exported ATLAS stage directory. The path may point anywhere and no common parent-folder layout is required. | `Path` | `empty parameter --&gt; ignored or automatic selection` |  | `../reference_system/exported/preprocessing_complete` |
| `reference` | Set True for exactly one system. Its vertical grid, molecular profile, data IDs, source parameters, and channel-specific metadata defaults define the comparison reference. | `bool` | `False` |  | `True` |
| `label` | Optional display label. When empty, metadata loading will use lidar_name and finally the system section ID as fallback. | `str` | `empty parameter --&gt; ignored or automatic selection` |  | `Reference lidar` |

## `channel:<comparison_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `label` | Optional comparison label. When empty, the reference system&#x27;s atlas_channel_id is used. | `str` | `empty parameter --&gt; ignored or automatic selection` |  | `355 nm parallel` |
| `force_qa_test` | Optional QA-test source forced for all systems in this channel comparison. When empty, default_qa_test is used. | `str` | `empty parameter --&gt; ignored or automatic selection` | Allowed values: ray, drk, pcb, tlc, tlc_rin, ray_pcb, pcb_aux, trg, dtm, cam | `ray_pcb` |
| `background_correction` | Enable or disable background correction for this complete channel-comparison bundle. | `bool` | `empty parameter --&gt; ignored or automatic selection` |  | `True` |
| `background_region` | Channel-specific background interval in kilometres, applied to every system in this comparison. Empty falls back to the general value or later reference-channel metadata. | `list[float]` | `empty parameter --&gt; ignored or automatic selection` | Required list size: 2 | `18.0, 22.0` |
| `normalisation` | Enable or disable normalization for this complete channel-comparison bundle. | `bool` | `empty parameter --&gt; ignored or automatic selection` |  | `True` |
| `normalisation_region` | Channel-specific normalization interval in kilometres, applied to every system in this comparison. Empty falls back to the general value or later reference-channel metadata. | `list[float]` | `empty parameter --&gt; ignored or automatic selection` | Required list size: 2 | `7.5, 9.0` |
| `normalise_to_molecular` | Override whether this channel bundle is normalized to the reference molecular profile. | `bool` | `empty parameter --&gt; ignored or automatic selection` |  | `True` |
| `plot_molecular` | Override whether the reference molecular profile is plotted for this channel bundle. | `bool` | `empty parameter --&gt; ignored or automatic selection` |  | `True` |

### Dynamic channel entries

- `<system_id> = <atlas_channel_id>` maps a system. The reference-system mapping is mandatory; omitted non-reference IDs inherit it; `off` excludes a non-reference system.
- `<system_id>.signal_source` and `<system_id>.signal_error_source` override exported parameters. Missing values inherit the reference system and then the general defaults.

## `pair:<comparison_id>`

| Parameter | Description | Type | Default | Allowed / limits | Example |
| --- | --- | --- | --- | --- | --- |
| `label` | Optional pair-comparison label. When empty, the reference system&#x27;s atlas_pair_id is used. | `str` | `empty parameter --&gt; ignored or automatic selection` |  | `VLDR 355 nm` |
| `force_qa_test` | Optional QA-test source forced for all systems in this pair comparison. When empty, default_qa_test is used. | `str` | `empty parameter --&gt; ignored or automatic selection` | Allowed values: ray, drk, pcb, tlc, tlc_rin, ray_pcb, pcb_aux, trg, dtm, cam | `pcb` |

### Dynamic pair entries

- `<system_id> = <atlas_pair_id>` maps a system. The reference-system mapping is mandatory; omitted non-reference IDs inherit it; `off` excludes a non-reference system.
- `<system_id>.pair_source` and `<system_id>.pair_error_source` override exported parameters. Missing values inherit the reference system and then the general defaults.

## Deferred metadata defaults

Empty system labels, background regions, and normalization regions may remain unresolved after parsing. They are resolved later from the reference system metadata, separately for every channel comparison.
