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
    compute_asign_dark_blocks,
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
    compute_smoothing_dark,
    compute_dark_correction,
    compute_signal_noise,
    compute_signal_smoothing,
    compute_mean_arrays,
    )

from processor.packaging import combine_QA_pack

from processor.signal_gluing import (
    compute_gluing_region, 
    compute_gluing,
    )

from processor.pol_cal_calculations import (
    compute_gain_ratio,
    compute_calibration_factor,
    compute_calibrated_ratio,
    compute_mean_calibrated_ratio,
    compute_mldr,
    compute_mean_vldr,
    compute_vldr
    # compute_eta,
    # compute_calibrated_ratio_pol_cal,
    # compute_vldr
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
        self.starting_dataset = ctx.starting_dataset
        self.run_id = 0
        self.stage_info = {}
        self.stage_history = {}
        self.stage_history_ids = {}
    
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
            "assign_dark":{
                "function":compute_asign_dark_blocks,
                "header":"Assign dark measurement blocks"
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
            "smoothing_dark":{
                "function":compute_smoothing_dark,
                "header":"Smoothing dark profiles"
                },
            "dark_correction":{
                "function":compute_dark_correction,
                "header":"Dark correction"
                },
            "signal_noise_calculation":{
                "function":compute_signal_noise,
                "header":"Signal noise calculation"
                },
            "signal_smoothing":{
                "function":compute_signal_smoothing,
                "header":"Signal smoothing"
                },
            "gluing_region":{
                "function":compute_gluing_region,
                "header":"Identify gluing region"
                },
            "gluing":{
                "function":compute_gluing,
                "header":"Gluing signals"
                },
            "computing_mean":{
                "function":compute_mean_arrays,
                "header":"Compute mean arrays"
                },
            "molecular_calculations":{
                "function":compute_molecular_calculations,
                "header":"Molecular calculations"
                },
            "gain_ratio": {
                "function": compute_gain_ratio,
                "header": "Calculating gain ratios for pcb_x45 and pcb_aux_x45",
                },
            "calibration_factor": {
                "function": compute_calibration_factor,
                "header": "Calculating polarization calibration factor eta",
                },
            "calibrated_ratio": {
                "function": compute_calibrated_ratio,
                "header": "Calculating Rayleigh calibrated depolarization ratios",
                },
            "calibrated_ratio_mean": {
                "function": compute_mean_calibrated_ratio,
                "header": "Calculating mean Rayleigh calibrated depolarization ratios",
                },
            "mldr": {
                "function": compute_mldr,
                "header": "Calculating molecular linear depolarization ratio",
                },
            "vldr": {
                "function": compute_vldr,
                "header": "Calculating volume linear depolarization ratio",
                },
            "vldr_mean": {
                "function": compute_mean_vldr,
                "header": "Calculating mean volume linear depolarization ratio",
                },
            # "eta": {
            #     "function": compute_eta,
            #     "header": "Calculating pol. calibration factor",
            #     },
            # "pol_calibrated_ratio": {
            #     "function": compute_calibrated_ratio_pol_cal,
            #     "header": "Calculating pol. calibrated ratio",
            #     },
            # "vldr": {
            #     "function": compute_vldr,
            #     "header": "Calculating VLDR",
            #     },
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
        self.stage_history[output_id] = ()
        self.stage_history_ids[output_id] = ()
            
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


    def register_history(self, output_id, input_ids):
        """Store the names and processing IDs of stages preceding ``output_id``.

        History is kept separately from ``stage_info`` so the existing Stage
        structure and public method signatures remain unchanged. The original
        ``stage_history`` registry continues to expose only external stage names.
        """
        if isinstance(input_ids, str):
            input_ids = [input_ids]

        history = []
        history_ids = []

        for input_id in input_ids:
            history.extend(self.stage_history.get(input_id, ()))
            history_ids.extend(self.stage_history_ids.get(input_id, ()))

            history.append(input_id)
            history_ids.append(self.stage_info[input_id].stage_name)

        # Stage names are unique, so use them to remove duplicate ancestry while
        # preserving dependency order and keeping IDs aligned with their names.
        unique_history = {}
        for stage_name, stage_id in zip(history, history_ids):
            unique_history.setdefault(stage_name, stage_id)

        self.stage_history[output_id] = tuple(unique_history.keys())
        self.stage_history_ids[output_id] = tuple(unique_history.values())

    def get_stage_history(self, item_id):
        """Return the external names of all stages preceding ``item_id``."""
        if item_id not in self.stage_info:
            raise KeyError(
                f"Unknown stage '{item_id}'. "
                f"Available stages: {list(self.stage_info)}"
            )

        return self.stage_history.get(item_id, ())

    def get_stage_history_ids(self, item_id):
        """Return the processing IDs of all stages preceding ``item_id``."""
        if item_id not in self.stage_info:
            raise KeyError(
                f"Unknown stage '{item_id}'. "
                f"Available stages: {list(self.stage_info)}"
            )

        return self.stage_history_ids.get(item_id, ())

    def get_stage_history_details(self, item_id):
        """Return preceding stages as ``(name, id)`` pairs."""
        return tuple(zip(
            self.get_stage_history(item_id),
            self.get_stage_history_ids(item_id),
        ))

    def has_stage(self, item_id, stage_name, include_self=False):
        """Return whether ``stage_name`` occurs in the history of ``item_id``.

        Each history entry is exposed as ``(stage_name, stage_id)``. This
        method explicitly checks index 0, i.e. the externally assigned stage
        name, rather than the processing-function ID at index 1.
        """
        if item_id not in self.stage_info:
            raise KeyError(
                f"Unknown stage '{item_id}'. "
                f"Available stages: {list(self.stage_info)}"
            )

        if any(
            history_entry[0] == stage_name
            for history_entry in self.get_stage_history_details(item_id)
        ):
            return True

        if include_self:
            return item_id == stage_name

        return False

    def has_stage_id(self, item_id, stage_id, include_self=False):
        """Return whether ``stage_id`` occurs in the history of ``item_id``.

        Each history entry is exposed as ``(stage_name, stage_id)``. This
        method explicitly checks index 1, i.e. the processing-function ID,
        rather than the externally assigned stage name at index 0.
        """
        if item_id not in self.stage_info:
            raise KeyError(
                f"Unknown stage '{item_id}'. "
                f"Available stages: {list(self.stage_info)}"
            )

        if any(
            history_entry[1] == stage_id
            for history_entry in self.get_stage_history_details(item_id)
        ):
            return True

        if include_self:
            return self.stage_info[item_id].stage_name == stage_id

        return False
               
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
        self.stage_history.pop(entry_id, None)
        self.stage_history_ids.pop(entry_id, None)
    
    def delete_registry(self, entry_id):
        
        del self.stage_info[entry_id]
        self.stage_history.pop(entry_id, None)
        self.stage_history_ids.pop(entry_id, None)

    def run(self, output_id, input_id, stage_name):
    
        # Print process header
        header = self.stage_map[stage_name]['header']
        print_header(header)
            
        # Get input data and output map
        input_data = self.export_stage(input_id)
        
        # Get output map and register
        output_map = self.prepare_output(output_id, stage_name = stage_name)
        self.register_history(output_id, input_id)
    
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
        self.register_history(output_id, input_id)

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
        self.register_history(output_id, input_id)
    
        # Processing
        output_data = combine_QA_pack(input_data)
        
        # Save output entries
        self.save_output(output_data, output_map)
        
    def package_from_stage(self, input_id):

        print_header(f"Packaging stage directly: {input_id}")

        input_data = self.export_stage(input_id)

        exported_data = combine_QA_pack(input_data)

        exported_data_swapped = defaultdict(dict)
        
        for k1, inner in exported_data.items():
            for k2, value in inner.items():
                exported_data_swapped[k2][k1] = value
        
        exported_data_swapped = dict(exported_data_swapped)
        
        return exported_data_swapped
        
    def checkout(self, output_id, input_id):
    
        # Print process header
        print_header(f"Saving to stage {output_id}")
            
        # Get input data and output map
        input_data = self.export_stage(input_id)
        
        # Get output map and register
        output_map = self.prepare_output(output_id, stage_name = "copy")
        self.register_history(output_id, input_id)

        # Make a 2-level copy
        output_data = {
            key_1: input_data[key_1].copy()
            for key_1 in input_data.keys()
        }
    
        # # Persist everything that can be persisted
        # for key_1 in output_data.keys():
        #     for key_2 in output_data[key_1].keys():
        #         value = output_data[key_1][key_2]
    
        #         if hasattr(value, "persist"):
        #             output_data[key_1][key_2] = value.persist()
    
        # Save output entries
        self.save_output(output_data, output_map)

              