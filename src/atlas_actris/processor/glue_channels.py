# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 15:03:12 2026

@author: owner
"""

channels = [
    "0355xpar", "0355xppr", "0355xcat", "0355xcpt",
    "0387xvan", "0387xvpn",
    "0532xpat", "0532xppt", "0532xcar", "0532xcpr",
    "0607xvan", "0607xvpn",
    "1064xpar", "1064xcat"
]

pairs = []
glue_channels = []
used_channels = set()     
final_channels = []


channel_map = {}

for ch in channels:
    key = ch[:6] + ch[7]  
    channel_map.setdefault(key, []).append(ch)

for key, group in channel_map.items():
    a_channel = next((c for c in group if c[6] == "a"), None)
    p_channel = next((c for c in group if c[6] == "p"), None)

    if a_channel and p_channel:
        pairs.append((a_channel, p_channel))

        glue = a_channel[:6] + "g" + a_channel[7]
        glue_channels.append(glue)

        used_channels.update([a_channel, p_channel])
        final_channels.append(glue)


for ch in channels:
    if ch not in used_channels:
        final_channels.append(ch)


print("Δυάδες analog / photon counting:")
for a, p in pairs:
    print(f"{a} <-> {p}")

print("\nGlue (sum) κανάλια:")
for g in glue_channels:
    print(g)

print("\nΤελική λίστα καναλιών:")
for ch in final_channels:
    print(ch)

sum_pairs = []
sum_channels = []

used_sum_channels = set()
final_after_sum = []


sum_map = {}

for ch in final_channels:
    key = ch[:5] + ch[6:7]
    sum_map.setdefault(key, []).append(ch)

for group in sum_map.values():
    if len(group) == 2:
        ch1, ch2 = group
        sum_pairs.append((ch1, ch2))

        sum_channel = ch1[:5] + "t" + ch1[6] + "x"
        sum_channels.append(sum_channel)

        used_sum_channels.update([ch1, ch2])
        final_after_sum.append(sum_channel)

for ch in final_channels:
    if ch not in used_sum_channels:
        final_after_sum.append(ch)

print("\nGlue δυάδες για SUM:")
for c1, c2 in sum_pairs:
    print(f"{c1} <-> {c2}")

print("\nSUM κανάλια:")
for s in sum_channels:
    print(s)
    
print("\nΤελική λίστα μετά το SUM:")
for ch in final_after_sum:
    print(ch)
