"""Noninteractive end-of-run actions selected in the GUI."""
import os


def _enabled(name):
    return os.environ.get(name, "N").strip().lower() in ("y", "yes")


def validate_gui_post_run(caller_info):
    """Validate parsed settings before any measurements are processed."""
    if os.environ.get("ATLAS_GUI_POST_RUN") != "1":
        return
    if _enabled("ATLAS_EXPORT_STAGE_ANSWER"):
        if _enabled("ATLAS_DELETE_EXPORTED_ANSWER"):
            raise ValueError("Choose either Clean stages or Export stages, not both.")
        if not caller_info.get("export_stages"):
            raise ValueError("Export stages requires a nonempty export_stages selection in call_atlas.ini.")


def run_gui_post_run(processor, caller_info, export_stages, delete_stages, clean_cache):
    """Called only after processing/report success; export before cache removal.

    An export exception propagates, preserving measurement caches for a retry.
    The existing stage utilities enforce path validation and export overwrite rules.
    """
    validate_gui_post_run(caller_info)
    if _enabled("ATLAS_EXPORT_STAGE_ANSWER"):
        export_stages(processor=processor, stage_names=caller_info["export_stages"],
                      output_folder=caller_info["output_folder"], overwrite=True,
                      ask=False, default_answer=False, print_estimated_size=True)
    elif _enabled("ATLAS_DELETE_EXPORTED_ANSWER"):
        delete_stages(output_folder=caller_info["output_folder"], ask=False)
    clean_cache(caller_info)
