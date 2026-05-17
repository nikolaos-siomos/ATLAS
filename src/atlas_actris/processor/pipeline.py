"""
@authors: N. Siomos, P. Paschou, Ioannis Binietoglou 
based on SULA project (https://react-gitlab.space.noa.gr/ReACT/eve/data-processing)
and also on https://gitlab.com/ioannis_binietoglou/lidar-processing/

Processing routines for signals 

=================================
Signal in 3D xarray dataset with dimensions [time, channel, bin/range]

Fucntions
 -- average_by_time: Average signals across the timeframes
 -- background_calculation: Calculates the solar background per timeframe and channel 
 -- background_correction: Performs the background correction on signals
 -- dark_correction: Removes the dark signals from the normal ananlog signals
 -- dead time correction: Performs the dead time correction onphoton channels
 -- detect_saturation: Identifies regions where signals are saturated
 -- height_calculation: Calculates the height above the lidar values per bin and channel
 -- range calculation: Calculates the range above the lidar values per bin and channel
 -- range_correction: Performs the range correction on signals
 -- smoothing: Smooths the signals (sliding average)
 -- trigger_correction: Perform the trigger correction per channel
 -- trim_vertically: Trim channels up to a maximum altitude
 -- unit_conv_counts_to_MHz: Converts raw counts to MHz for the photon channels

"""

import copy
from typing import Any, Dict
from utils.error_classes import OverwriteError
from utils.printouts import print_header

from dataclasses import dataclass
from collections import defaultdict

from processor.handle_overflows import compute_check_for_overflows
from processor.signal_flagging import compute_detect_saturation
from processor.signal_trimming import (
    compute_slice_and_exclude, 
    compute_screen_low_shots,
    )

from processor.signal_processing import (
    compute_height_and_range_calculation,
    compute_unit_conv_counts_to_MHz,
    compute_dead_time_correction,
    compute_background_calculation,
    compute_averaging_by_time_single,
    compute_averaging_by_time_low_res,
    compute_averaging_by_time_high_res,
    compute_trim_vertically,
    compute_background_correction,
    compute_range_correction,
    compute_dark_correction,
    compute_signal_noise,
    )

from processor.packaging import combine_QA_pack

from processor.signal_gluing import (
    compute_gluing_region, 
    compute_gluing,
    )

from processor.molecular_calculations import compute_molecular_calculations

from processor.selectors import (
    get_from_selector,
    set_from_selector,
    get_io_map,
    get_source_map
)

# Shared, read-only context
@dataclass(frozen=True)
class Context:
    processing_info: Dict[str, Any]
    settings_info: Dict[str, Any]
    starting_dataset: Dict[str, Any]
    
@dataclass(frozen=True)
class Stage:
    item_id: str
    stage_name: str
    input_map: str
    output_map: dict
    run_id: int

class Processor():
    
    def __init__(self, ctx: Context):
        self.processing_info = ctx.processing_info
        self.settings_info = ctx.settings_info
        self.starting_dataset = ctx.starting_dataset
        self.run_id = 0
        self.stage_info = {}
    
        self.stage_map = {
            "height_and_range_calculation":{
                "function":compute_height_and_range_calculation,
                "header":"Calculation of ranges, heights_agl, and heights_asl"
                },
            "slice_and_exclude":{
                "function":compute_slice_and_exclude,
                "header":"Slicing and excluding measurement parts"
                },        
            "screen_low_shots":{
                "function":compute_screen_low_shots,
                "header":"Removing profiles with too few shots"
                },    
            "handling_overflows":{
                "function":compute_check_for_overflows,
                "header":"Handling overflow values"
                },
            "check_saturation":{
                "function":compute_detect_saturation,
                "header":"Checking signal saturation"
                },
            "photon_units_conversion":{
                "function":compute_unit_conv_counts_to_MHz,
                "header":"Converting photon units to MHz"
                },
            "dead_time_correction":{
                "function":compute_dead_time_correction,
                "header":"Performing deadtime correction"
                },
            "background_calculation":{
                "function":compute_background_calculation,
                "header":"Calculating temporally resolved background"
                },
            "averaging_by_time":{
                "function":compute_averaging_by_time_single,
                "header":"Averaging along time dimension to a single profile"
                },
            "averaging_by_time_low_res":{
                "function":compute_averaging_by_time_low_res,
                "header":"Averaging profiles along dime dimension - low resolution"
                },
            "averaging_by_time_high_res":{
                "function":compute_averaging_by_time_high_res,
                "header":"Averaging profiles along dime dimension - high resolution"
                },
            "trim_vertically":{
                "function":compute_trim_vertically,
                "header":"Vertical trimming:\nRemoving parts below the zero bin and above the maximum provided height"
                },
            "background_correction":{
                "function":compute_background_correction,
                "header":"Background correction"
                },
            "range_correction":{
                "function":compute_range_correction,
                "header":"Range correction"
                },
            "dark_correction":{
                "function":compute_dark_correction,
                "header":"Dark correction"
                },
            "signal_noise_calculation":{
                "function":compute_signal_noise,
                "header":"Signal noise calculation"
                },
            "gluing_region":{
                "function":compute_gluing_region,
                "header":"Identify gluing region"
                },
            "gluing":{
                "function":compute_gluing,
                "header":"Gluing signals"
                },
            "molecular_calculations":{
                "function":compute_molecular_calculations,
                "header":"Molecular calculations"
                },
            # "pol_calibration_factor":{
            #     "function":,
            #     "header":"Calculating polarization calibration factor"
            #     },
            # "summing_pol":{
            #     "function":,
            #     "header":"Summing polarization signals"
            #     },
            }
        

        print_header("Initializing debugging-related data arrays")
        
        output_id = 'init'
    
        # Get I/O map and validate
        output_map = get_io_map(output_id)
                
        # Register output entries
        self.register_item(
            item_id=output_id,
            stage_name='initialization',
            output_map=output_map
        )
            
        # Save output entries
        for s in output_map.keys():
            if s in self.starting_dataset:
                set_from_selector(self, output_map[s], self.starting_dataset[s])
            else:
                set_from_selector(self, output_map[s], {})
                
        
    def register_item(
        self, item_id, stage_name, input_map=None, output_map=None
    ):
        if item_id in self.stage_info:
            raise OverwriteError(
                f"-- The provided item_id ({item_id}) already exists. Overwriting is not allowed"
            )
    
        self.run_id += 1
    
        stage_info = Stage(
            item_id=item_id,
            stage_name=stage_name,
            input_map=input_map,
            output_map=output_map,
            run_id=self.run_id
        )
    
        self.stage_info[item_id] = stage_info
               
    def export_stage(self, input_id):

        input_map = get_io_map(input_id)
        
        # Gather input data
        exported_data = {s: get_from_selector(self, input_map[s]) 
                         for s in input_map.keys()}
        
        return exported_data
    

    def export_test_from_stage(self, input_id):

        input_map = get_io_map(input_id)
        
        # Gather input data
        exported_data = {s: get_from_selector(self, input_map[s]) 
                         for s in input_map.keys()}
        
        exported_data_swapped = defaultdict(dict)
        
        for k1, inner in exported_data.items():
            for k2, value in inner.items():
                exported_data_swapped[k2][k1] = value
        
        exported_data_swapped = dict(exported_data_swapped)
        
        return exported_data_swapped
             
    def prepare_output(self, output_id, stage_name):
        output_map = get_io_map(output_id)
        
        # Register output entries
        self.register_item(
            item_id=output_id,
            stage_name=stage_name,
            output_map=output_map
        )
        
        return output_map
    
    def save_output(self, output, output_map):
        
        for s in output_map:
            set_from_selector(self, output_map[s], output[s])
        
    def delete_entry(self, entry_id):
        
        source_map = get_source_map(entry_id)
        
        db_list = [db_name['db'] for db_name in source_map.values()]
        
        for db_name in db_list:
            db = getattr(self, db_name)
            if entry_id in db:
                del db[entry_id]
        
        del self.stage_info[entry_id]
    
    def delete_registry(self, entry_id):
        
        del self.stage_info[entry_id]

    def run(self, output_id, input_id, stage_name):
    
        # Print process header
        header = self.stage_map[stage_name]['header']
        print_header(header)
            
        # Get input data and output map
        input_data = self.export_stage(input_id)
        
        # Get output map and register
        output_map = self.prepare_output(output_id, stage_name = stage_name)
    
        # Processing
        output_data = self.stage_map[stage_name]['function'](self.processing_info, input_data)
        
        # Save output entries
        self.save_output(output_data, output_map)
        
    def copy(self, output_id, input_id):
    
        # Print process header
        print_header(f"Copying stages: {input_id} to {output_id}")
            
        # Get input data and output map
        input_data = self.export_stage(input_id)
        
        # Get output map and register
        output_map = self.prepare_output(output_id, stage_name = "copy")

        # Make a copy of input data
        output_data = copy.deepcopy(input_data)

        # Save output entries
        self.save_output(output_data, output_map)
        
    def package(self, output_id, input_id):

        # Print process header
        print_header(f"Packaging stage: {input_id} to {output_id}")
            
        # Get input data and output map
        input_data = self.export_stage(input_id)
        
        # Get output map and register
        output_map = self.prepare_output(output_id, stage_name = "packaging")
    
        # Processing
        output_data = combine_QA_pack(input_data)
        
        # Save output entries
        self.save_output(output_data, output_map)
        
    def checkout(self, output_id, input_id):
    
        # Print process header
        print_header(f"Calculating stage {input_id} and save to stage {output_id}")
            
        # Get input data and output map
        input_data = self.export_stage(input_id)
        
        # Get output map and register
        output_map = self.prepare_output(output_id, stage_name = "copy")

        # Make a 2-level copy
        output_data = {
            key_1: input_data[key_1].copy()
            for key_1 in input_data.keys()
        }
    
        # Persist everything that can be persisted
        for key_1 in output_data.keys():
            for key_2 in output_data[key_1].keys():
                value = output_data[key_1][key_2]
    
                if hasattr(value, "persist"):
                    output_data[key_1][key_2] = value.persist()
    
        # Save output entries
        self.save_output(output_data, output_map)

              