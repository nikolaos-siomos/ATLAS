# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 15:05:17 2026

@author: owner
"""

def create_glue_channels(channels, verbose=True):
    """
    Create glue channels (g) from the analog (a) and photon counting (p)

    :param channels: list of channels (with 8 characters )
    :param verbose: if True, prints the AN and PC channels
    :return: list of final channels (sum & others)
    """

    channel_map = {}
    used_channels = set()
    final_channels = []

    for ch in channels:
        if len(ch) != 8:
            continue

        key = ch[:6] + ch[7]  
        channel_map.setdefault(key, []).append(ch)

    # Create glue channels
    for group in channel_map.values():
        a_channel = next((c for c in group if c[6] == "a"), None)
        p_channel = next((c for c in group if c[6] == "p"), None)

        if a_channel and p_channel:
            sum_channel = a_channel[:6] + "g" + a_channel[7]

            if verbose:
                print(f"SUM: {a_channel} + {p_channel} -> {sum_channel}")

            final_channels.append(sum_channel)
            used_channels.update([a_channel, p_channel])

    # Add channels that didn't used
    for ch in channels:
        if ch not in used_channels:
            final_channels.append(ch)

    return final_channels

def create_sum_channels(channels, verbose=True):
    """
    Create sum channels (s) from the trabsmitted (t) and reflected (r)

    :param channels: list of channels (with 8 characters )
    :param verbose: if True, prints the Tr and Refl channels
    :return: list of final channels (sum & others)
    """

    channel_map = {}
    used_channels = set()
    final_channels = []

    for ch in channels:
        if len(ch) != 8:
            continue

        key = ch[:5] + ch[6:7]   
        channel_map.setdefault(key, []).append(ch)

    # Create sum channels
    for group in channel_map.values():
        if len(group) == 2:
            ch1, ch2 = group

            sum_channel = ch1[:5] + "t" + ch1[6] + "x"

            if verbose:
                print(f"SUM2: {ch1} + {ch2} -> {sum_channel}")

            final_channels.append(sum_channel)
            used_channels.update([ch1, ch2])

    # Add channels that didn't used
    for ch in channels:
        if ch not in used_channels:
            final_channels.append(ch)

    return final_channels

def create_telescope_glue_channels(channels, verbose=True):
    """
    Create glue channels from different telescope types.
    Channels must differ only in the 1st character (index 0).

    :param channels: list of channels (8 characters)
    :param verbose: if True, prints the telescope pairs/groups
    :return: list of final channels (glue & others)
    """

    channel_map = {}
    used_channels = set()
    final_channels = []

    for ch in channels:
        if len(ch) != 8:
            continue

        key = ch[1:]   
        channel_map.setdefault(key, []).append(ch)

    # Create telescope glue channels
    for group in channel_map.values():
        if len(group) >= 2:
            glue_channel = "a" + group[0][1:]

            if verbose:
                print(f"TEL_GLUE: {' + '.join(group)} -> {glue_channel}")

            final_channels.append(glue_channel)
            used_channels.update(group)

    # Add channels that didn't used
    for ch in channels:
        if ch not in used_channels:
            final_channels.append(ch)

    return final_channels