"""Channel availability checks shared by initialization and direct QA calls.

Checks inspect coordinates and small pair metadata only; signal data stay lazy.
"""

from utils.error_classes import CustomWarning

TELECOVER_GROUPS = {
    "tlc": ("tlc_north", "tlc_east", "tlc_south", "tlc_west"),
    "tlc_rin": ("tlc_outer", "tlc_inner"),
}
CALIBRATION_GROUPS = {
    "pcb": ("pcb_p45", "pcb_m45"),
    "pcb_aux": ("pcb_aux_p45", "pcb_aux_m45"),
}


def warn_once(caller_info, message):
    emitted = caller_info.setdefault("_channel_warnings", [])
    if message not in emitted:
        CustomWarning(message)
        emitted.append(message)


def resolve_measurement(key, store, caller_info):
    """Resolve only configured loading aliases; never substitute ray_pcb."""
    aliases = caller_info.get("loading_map", {})
    seen = set()
    while key not in seen:
        if key in store:
            return key
        seen.add(key)
        key = aliases.get(key)
        if key is None:
            break
    return None


def drop_measurement(key, profiles, metadata, caller_info):
    profiles.pop(key, None)
    for group in metadata.values():
        if isinstance(group, dict):
            group.pop(key, None)
    removed = caller_info.setdefault("_removed_channel_measurements", [])
    if key not in removed:
        removed.append(key)


def telecover_channels_match(group, profiles, caller_info):
    """Compare available sectors; a sector removed for empty channels fails QA.

    Entirely absent sectors retain the existing partial-sector policy.
    """
    sectors = TELECOVER_GROUPS[group]
    removed = caller_info.get("_removed_channel_measurements", [])
    aliases = caller_info.get("loading_map", {})
    rejected = [key for key in sectors if key in removed or aliases.get(key) in removed]
    if rejected:
        warn_once(caller_info, f"{group}: empty sector {', '.join(rejected)}; test skipped.")
        return False
    expected = caller_info.get("_channel_validation", {}).get("telecover_members", {}).get(group, [])
    lost = [key for key in expected if resolve_measurement(key, profiles, caller_info) is None]
    if lost:
        warn_once(caller_info, f"{group}: missing sector {', '.join(lost)}; test skipped.")
        return False
    available = [(key, resolve_measurement(key, profiles, caller_info)) for key in sectors]
    available = [(key, source) for key, source in available if source is not None]
    if not available:
        return False
    reference = set(profiles[available[0][1]].channel.values.tolist())
    if not reference or any(
        set(profiles[source].channel.values.tolist()) != reference
        for _, source in available[1:]
    ):
        warn_once(caller_info, f"{group}: sector channels differ; test skipped.")
        return False
    return True


def telecover_test_allowed(group, profiles, caller_info):
    previous = caller_info.get("_channel_validation", {}).get("telecover", {})
    if previous.get(group) is False:
        return False
    return telecover_channels_match(group, profiles, caller_info)


def calibration_pair_ids(profiles, pol_cal_info, caller_info):
    """Return active base-pair IDs for primary and auxiliary calibration.

    Auxiliary eligibility depends on primary eligibility, never the reverse.
    Full base metadata on ray remain available for independent MLDR calculations.
    """
    active = {"pcb": [], "pcb_aux": []}
    resolved = {
        key: resolve_measurement(key, profiles, caller_info)
        for key in ("pcb_p45", "pcb_m45", "pcb_aux_p45", "pcb_aux_m45", "ray")
    }
    if not any(resolved[key] for key in CALIBRATION_GROUPS["pcb"] + CALIBRATION_GROUPS["pcb_aux"]):
        return active
    required = ("pcb_p45", "pcb_m45", "ray")
    missing = [key for key in required if resolved[key] is None]
    if missing:
        warn_once(caller_info, f"pcb: missing {', '.join(missing)}; calibration skipped.")
        return active

    def pair_table(logical):
        source = resolved[logical]
        return pol_cal_info.get(logical, pol_cal_info.get(source))

    info = pair_table("pcb_p45")
    if info is None or "pair" not in info.dims:
        return active
    if not {"ch_r", "ch_t"}.issubset(set(info.parameters.values.tolist())):
        return active
    channel_sets = {
        key: set(profiles[source].channel.values.tolist())
        for key, source in resolved.items() if source is not None
    }
    pair_channels = {}
    for pair in info.pair.values.tolist():
        r = info.sel(parameters="ch_r", pair=pair).item()
        t = info.sel(parameters="ch_t", pair=pair).item()
        pair_channels[pair] = (r, t)
        missing = [f"{key}: {','.join(str(ch) for ch in (r, t) if ch not in channel_sets[key])}"
                   for key in required if any(ch not in channel_sets[key] for ch in (r, t))]
        if missing:
            warn_once(caller_info, f"pcb pair {r}/{t}: missing {'; '.join(missing)}; skipped.")
        else:
            active["pcb"].append(pair)

    aux_keys = CALIBRATION_GROUPS["pcb_aux"]
    if not any(resolved[key] for key in aux_keys):
        return active
    missing = [key for key in aux_keys if resolved[key] is None]
    if missing:
        warn_once(caller_info, f"pcb_aux: missing {', '.join(missing)}; evaluation skipped.")
        return active
    for pair in active["pcb"]:
        r, t = pair_channels[pair]
        missing = [f"{key}: {','.join(str(ch) for ch in (r, t) if ch not in channel_sets[key])}"
                   for key in aux_keys if any(ch not in channel_sets[key] for ch in (r, t))]
        if missing:
            warn_once(caller_info, f"pcb_aux pair {r}/{t}: missing {'; '.join(missing)}; skipped.")
        else:
            active["pcb_aux"].append(pair)
    return active


def validate_channels(caller_info, profiles, metadata):
    """Record QA eligibility without changing diagnostic plotting requests."""
    caller_info["_channel_validation"] = {
        "telecover": {
            group: telecover_channels_match(group, profiles, caller_info)
            for group in TELECOVER_GROUPS
        },
        "telecover_members": {
            group: [key for key in sectors if resolve_measurement(key, profiles, caller_info) is not None]
            for group, sectors in TELECOVER_GROUPS.items()
        },
        "calibration_pairs": calibration_pair_ids(
            profiles, metadata.get("pol_cal_info", {}), caller_info
        ),
    }
    return caller_info
