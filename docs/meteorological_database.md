# Meteorological station database

ATLAS uses this table to fill empty meteorological initialization settings for the configuration file's `station_id`. See [station defaults and the distance limit](radiosonde_selection.md#station-defaults-from-the-meteorological-database) for the rules and user-override behavior.

The CSV remains at `src/atlas_actris/molecular/meteorological_database.csv`. This page reads it during each documentation build; edit that file to update both the application defaults and the published table, then rebuild the site.

`station_id` is the EARLINET code. `cloudnet_station_name` contains the Cloudnet site identifier. Radiosonde distances are in kilometres, and `rsonde_station_altitude` is in metres above sea level, shown to one decimal place. Empty cells mean that the corresponding value is not supplied. A blank radiosonde distance prevents radiosonde defaults from being used.

The database includes fixed EARLINET sites after excluding stations marked not active or without a listed PI. Radiosonde candidates include historical records; the listed site and WMO ID do not establish present operational status or Wyoming download availability. The default distance limit still applies to every row.

<!-- ATLAS_METEOROLOGICAL_DATABASE -->
