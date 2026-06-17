#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 19:57:33 2026

@author: ime_rasenberg
"""

import os
import json

master_dict = {}
decoder = json.JSONDecoder()

# ==========================================================
# find all beta folders
# ==========================================================
beta_folders = [
    f for f in os.listdir(".")
    if os.path.isdir(f) and f.startswith("Data ")
]

print(f"Found {len(beta_folders)} beta folders")

for data_folder in beta_folders:

    # --------------------------------------
    # extract beta from folder name
    # --------------------------------------
    try:
        beta = float(data_folder.split(" ", 1)[1])
    except ValueError:
        print(f"Skipping folder {data_folder}")
        continue

    print(f"\nLoading beta = {beta}")

    if beta not in master_dict:
        master_dict[beta] = {}

    # --------------------------------------
    # collect json files
    # --------------------------------------
    filepaths = [
        os.path.join(data_folder, f)
        for f in os.listdir(data_folder)
        if os.path.isfile(os.path.join(data_folder, f))
        and f.endswith(".json")
    ]

    print(f"Found {len(filepaths)} JSON files")

    for filepath in filepaths:

        filename = os.path.basename(filepath)

        # --------------------------------------
        # scrape metadata
        # --------------------------------------
        name_no_ext = os.path.splitext(filename)[0]
        parts = name_no_ext.split("__")

        metadata = {}

        for part in parts:
            if "_" in part:
                key, value = part.split("_", 1)
                metadata[key] = value

        if not all(k in metadata for k in ["D", "Hz", "I"]):
            print(f"Skipping {filename}")
            continue

        try:
            D = float(metadata["D"])
            Hz = float(metadata["Hz"])
            I = int(metadata["I"])
        except ValueError:
            print(f"Malformed metadata in {filename}")
            continue

        # --------------------------------------
        # initialise nested structure
        # master_dict[beta][D][Hz][I]
        # --------------------------------------
        if D not in master_dict[beta]:
            master_dict[beta][D] = {}

        if Hz not in master_dict[beta][D]:
            master_dict[beta][D][Hz] = {}

        master_dict[beta][D][Hz][I] = {}

        # --------------------------------------
        # load concatenated json objects
        # --------------------------------------
        with open(filepath, "r") as f:
            content = f.read().strip()

        pos = 0

        while pos < len(content):

            while pos < len(content) and content[pos].isspace():
                pos += 1

            if pos >= len(content):
                break

            try:
                obj, idx = decoder.raw_decode(content[pos:])
                pos += idx

                step = obj["step"]

                master_dict[beta][D][Hz][I][step] = obj

            except json.JSONDecodeError as e:
                print(f"Error in {filename} near position {pos}: {e}")
                break

print("Finished loading.")