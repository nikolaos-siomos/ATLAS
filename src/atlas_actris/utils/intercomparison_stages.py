#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Selective access to exported ATLAS stages for intercomparison datasets."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from utils.export_processing_stage import import_processing_entry


def _split_stage_path(stage_path: str | Path) -> tuple[Path, str]:
    path = Path(stage_path).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"Exported stage directory does not exist: {path}")
    if path.parent.name != "exported":
        raise ValueError(
            "Configured stage_path must point to an exported ATLAS stage directory "
            f"with layout '<output_folder>/exported/<stage_name>': {path}"
        )
    return path.parent.parent, path.name


class IntercomparisonStageStore:
    """On-demand, cached access to the exported stage behind each dataset."""

    def __init__(self, intercomparison_info: Mapping[str, Any]):
        datasets = intercomparison_info.get("datasets")
        if not isinstance(datasets, Mapping) or not datasets:
            raise ValueError(
                "intercomparison_info does not contain any configured datasets"
            )

        self.datasets: dict[str, dict[str, Any]] = {}
        self.cache: dict[str, dict[str, dict[str, Any]]] = {}
        self._entry_cache: dict[tuple[str, str, str, str], Any] = {}

        for dataset_id, dataset_info in datasets.items():
            stage_path = dataset_info.get("stage_path")
            if stage_path is None:
                raise ValueError(f"Dataset {dataset_id!r} does not define stage_path")

            output_folder, stage_name = _split_stage_path(stage_path)
            dataset_id = str(dataset_id)

            self.datasets[dataset_id] = {
                "system_label": dataset_info.get("system_label"),
                "dataset_label": dataset_info.get("dataset_label"),
                "reference": bool(dataset_info.get("reference", False)),
                "stage_path": Path(stage_path),
                "output_folder": output_folder,
                "stage_name": stage_name,
                "qa_test": dataset_info.get("qa_test"),
                "signal_source": dataset_info.get("signal_source"),
                "signal_error_source": dataset_info.get("signal_error_source"),
                "pair_source": dataset_info.get("pair_source"),
                "pair_error_source": dataset_info.get("pair_error_source"),
            }
            self.cache[dataset_id] = {}

    def get_entry(self, dataset_id: str, qa_test: str, parameter: str) -> Any:
        """Return one required exported entry, importing it only once per dataset."""
        return self._get_entry(dataset_id, qa_test, parameter, required=True)

    def get_optional_entry(self, dataset_id: str, qa_test: str, parameter: str) -> Any:
        """Return one optional exported entry, or None when that parameter is absent."""
        return self._get_entry(dataset_id, qa_test, parameter, required=False)

    def _get_entry(
        self,
        dataset_id: str,
        qa_test: str,
        parameter: str,
        *,
        required: bool,
    ) -> Any:
        """Import one entry with shared caching and optional-metadata handling."""

        dataset_id = str(dataset_id)
        qa_test = str(qa_test)
        parameter = str(parameter)

        if dataset_id not in self.datasets:
            raise KeyError(
                f"Unknown intercomparison dataset {dataset_id!r}. "
                f"Available datasets: {sorted(self.datasets)}"
            )

        qa_cache = self.cache[dataset_id].setdefault(qa_test, {})
        if parameter in qa_cache:
            return qa_cache[parameter]

        dataset = self.datasets[dataset_id]
        shared_key = (
            str(dataset["output_folder"]),
            dataset["stage_name"],
            qa_test,
            parameter,
        )
        if shared_key in self._entry_cache:
            value = self._entry_cache[shared_key]
            qa_cache[parameter] = value
            return value

        try:
            value = import_processing_entry(
                output_folder=str(dataset["output_folder"]),
                stage_name=dataset["stage_name"],
                qa_test=qa_test,
                parameter=parameter,
            )
        except KeyError as exc:
            if not required:
                self._entry_cache[shared_key] = None
                qa_cache[parameter] = None
                return None
            # Keep the original importer detail available to the caller, but
            # avoid wrapping KeyError in another quoted KeyError message.
            detail = exc.args[0] if exc.args else str(exc)
            raise ValueError(
                f"Could not load parameter {parameter!r} for dataset {dataset_id!r} "
                f"from QA test {qa_test!r} in exported stage "
                f"{dataset['stage_name']!r}.\n"
                f"{detail}"
            ) from None
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Could not load data for dataset {dataset_id!r} from exported "
                f"stage {dataset['stage_name']!r}.\n"
                f"{exc}"
            ) from None

        self._entry_cache[shared_key] = value
        qa_cache[parameter] = value
        return value

    def loaded_entries(self) -> dict[str, dict[str, tuple[str, ...]]]:
        """Summarize entries currently present in the cache."""

        return {
            dataset_id: {
                qa_test: tuple(parameters)
                for qa_test, parameters in qa_cache.items()
                if parameters
            }
            for dataset_id, qa_cache in self.cache.items()
        }


def collect_intercomparison_stages(
    intercomparison_info: Mapping[str, Any],
) -> IntercomparisonStageStore:
    """Create an on-demand stage store without importing stage data yet."""

    return IntercomparisonStageStore(intercomparison_info)
