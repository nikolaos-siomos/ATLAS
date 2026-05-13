#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 23:11:27 2026

@author: nikos
"""

def run_linear_recipe(
    processor,
    recipe,
    initial_input,
    checkout_id=None,
):
    input_id = initial_input

    for output_id, stage_name in recipe:
        processor.run(
            output_id=output_id,
            input_id=input_id,
            stage_name=stage_name,
        )

        input_id = output_id

    if checkout_id is not None:

        processor.checkout(
            output_id=checkout_id,
            input_id=input_id,
        )

        return checkout_id
