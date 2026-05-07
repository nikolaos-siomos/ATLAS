GET_CONFIG - Retrieve and Export SCC and ATLAS Configuration Files

get_config retrieves an SCC HOI configuration file, exports it as a
.txt file, and generates a corresponding ATLAS config_file.ini.
Both output locations can be controlled through command-line arguments.
The name of the config file will have the following format:
config_file_<scc_config_id>_exported_<creation_date>_<creation_time>

MINIMAL USAGE
get_config -i <scc_configuration_id> -o <hoi_output_folder>

ARGUMENTS
Required:

-i, --scc_configuration_id (string)
The SCC configuration ID to retrieve.

-o, --hoi_output_folder (string)
Folder where the exported SCC HOI file (.txt) will be written.

Optional:

-c, --atlas_configuration_folder (string)
Folder where the ATLAS .ini configuration file will be created.
If not provided, the HOI output folder (-o) is used.

-v, --verbose / --no-verbose (bool)
Enables extra debug information during processing.

NOTES

The SCC ID (-i) and HOI output folder (-o) are mandatory.

The ATLAS configuration file is always generated after the SCC file.

If -c is omitted, the ATLAS configuration file is placed in the same
folder given by -o.


EXAMPLES
get_config -i 962 -o ./output
get_config -i 1203 -o ~/Downloads -v
get_config -i 2001 -o ./SCC_HOI -c ./configurations

