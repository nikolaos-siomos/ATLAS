#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command-line entry point for exporting an ATLAS configuration from SCC."""

import sys
from .utils.get_scc_config import parse_args, export_scc_config, SCCDownloadError


def main(argv=None):
    """Run the SCC exporter using the existing command-line arguments."""
    cmd_args = parse_args(argv)
    try:
        export_scc_config(
            scc_configuration_id=cmd_args.scc_configuration_id,
            atlas_configuration_file=cmd_args.atlas_configuration_file,
            export_hoi_cfg=cmd_args.export_hoi_cfg,
            output_folder=cmd_args.output_folder,
            csv_output_folder=cmd_args.output_folder,
            scc_compatible_format=cmd_args.scc_compatible_format,
        )
    except SCCDownloadError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("SCC export cancelled.", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
