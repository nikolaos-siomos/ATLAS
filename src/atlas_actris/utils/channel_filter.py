import numpy as np


CHANNEL_FIELDS = {
    "exclude_wavelength": slice(0, 4),
    "exclude_telescope_type": 4,
    "exclude_channel_type": 5,
    "exclude_acquisition_mode": 6,
    "exclude_channel_subtype": 7,
}


def _as_list(value):
    """
    Convert None, string, scalar, or list-like input to a clean list.
    """

    if value is None:
        return []

    if isinstance(value, str):
        value = value.strip()
        if not value:
            return []
        return [value]

    try:
        return list(value)
    except TypeError:
        return [value]


def _normalize_wavelengths(values):
    """
    Ensure wavelengths such as 355 become '0355'.
    """

    return [str(value).strip().zfill(4) for value in values if str(value).strip()]


def get_channels_to_keep(channel_ids, caller_info):
    """
    Build a list of channels to keep based on caller_info.
    """

    channel_ids = [str(channel_id) for channel_id in channel_ids]

    for channel_id in channel_ids:
        if len(channel_id) != 8:
            raise ValueError(
                f"Invalid channel ID {channel_id!r}. "
                "Expected a string of length 8."
            )

    # Keep only selected full channel IDs, if provided.
    selected_channels = [
        str(channel_id).strip()
        for channel_id in _as_list(caller_info.get("channels"))
        if str(channel_id).strip()
    ]

    if selected_channels:
        keep = np.array(
            [channel_id in selected_channels for channel_id in channel_ids],
            dtype=bool,
        )
    else:
        keep = np.ones(len(channel_ids), dtype=bool)

    # Exclude full channel IDs directly.
    exclude_channels = [
        str(channel_id).strip()
        for channel_id in _as_list(caller_info.get("exclude_channels"))
        if str(channel_id).strip()
    ]

    if exclude_channels:
        keep &= ~np.isin(channel_ids, exclude_channels)

    exclusions = {
        "exclude_wavelength": _normalize_wavelengths(
            _as_list(caller_info.get("exclude_wavelength"))
        ),
        "exclude_telescope_type": [
            str(value).strip()
            for value in _as_list(caller_info.get("exclude_telescope_type"))
            if str(value).strip()
        ],
        "exclude_channel_type": [
            str(value).strip()
            for value in _as_list(caller_info.get("exclude_channel_type"))
            if str(value).strip()
        ],
        "exclude_acquisition_mode": [
            str(value).strip()
            for value in _as_list(caller_info.get("exclude_acquisition_mode"))
            if str(value).strip()
        ],
        "exclude_channel_subtype": [
            str(value).strip()
            for value in _as_list(caller_info.get("exclude_channel_subtype"))
            if str(value).strip()
        ],
    }

    # Exclude based on the relevant channel-ID character fields.
    for exclusion_key, excluded_values in exclusions.items():
        if not excluded_values:
            continue

        field_position = CHANNEL_FIELDS[exclusion_key]

        for index, channel_id in enumerate(channel_ids):
            if channel_id[field_position] in excluded_values:
                keep[index] = False

    return [
        channel_id
        for channel_id, should_keep in zip(channel_ids, keep)
        if should_keep
    ]


def filter_channel_object(obj, caller_info, channel_dim="channel"):
    """
    Apply channel filtering to one xarray DataArray or Dataset.

    If the object has no dims attribute, or has no channel dimension,
    it is returned unchanged.
    """

    if not hasattr(obj, "dims"):
        return obj

    if channel_dim not in obj.dims:
        return obj

    channels_to_keep = get_channels_to_keep(
        obj[channel_dim].values,
        caller_info,
    )

    return obj.sel({channel_dim: channels_to_keep})


def filter_channel_data(data_dict, caller_info, channel_dim="channel"):
    """
    Apply channel filtering to all xarray DataArray/Dataset objects
    inside a dictionary.

    Values without a channel dimension are returned unchanged.
    """

    filtered_dict = {}

    for key, obj in data_dict.items():
        filtered_dict[key] = filter_channel_object(
            obj,
            caller_info,
            channel_dim=channel_dim,
        )

    return filtered_dict


def filter_channels(caller_info, profiles, metadata):
    """
    Apply channel filtering to profiles and metadata.

    profiles is expected to be a dictionary of xarray objects.
    metadata is expected to be a dictionary whose values are dictionaries.
    Non-xarray values are returned unchanged.
    """

    profiles = filter_channel_data(profiles, caller_info)

    for key, value in metadata.items():
        if isinstance(value, dict):
            metadata[key] = filter_channel_data(value, caller_info)
        else:
            metadata[key] = filter_channel_object(value, caller_info)

    return profiles, metadata