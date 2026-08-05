README for eia_{%Y%m}_{sat}{id}_ascii.txt
-------------------------------------------------------------------------------
COLUMNS
-------
~ the following columns were created in
    pyDetEPI.detection.stats.cand_info.epi_stats ~

        pass_id, time_start, time_end, mlat_start, mlat_end, glon_start,
        glon_end, glat_start, glat_end, lt, eia_type, crest1_mlat, crest1_np,
        crest2_mlat, crest2_np, crest3_mlat, crest3_np, alt

-------------------------------------------------------------------------------
Notes
------
    This output file utilizes pyValEIA to detect the EIA types.
    /Crest{#}/ is used in place of 1 through 3 depedning on the number of
    EIA crests detected
-------------------------------------------------------------------------------
Columns
-------
pass_id : format %Y%m%d_%H%M_{sat}{#} for each satellite pass meant to be used
            as an identifier for detections during the same pass
time_start : starting universal time of pass between +/-30 degrees MLat
time_end : ending universal time of pass between +/-30 degrees MLat
mlat_start : starting magnetic latitude of pass
mlat_end : ending magnetic latitude of pass
mlon_start : starting magnetic longitude of pass
mlon_end : ending magnetic longitude of pass
glon_start : starting geographic longitude of pass
glon_end : ending geographic longitude of pass
glat_start : starting geographic latitude of pass
glat_end : ending geographic latitude of pass
lt : local time at the magentic equator
eia_type : EIA categorization
    Types
        trough, flat (north, south, symmetric), peak (north, south, symmetric),
        eia classic ((north, south, symmetric),
        eia saddle (north, south, symmetric),
        eia ghost (north, south, symmetric), eia ghost peak (north, south)
crest{#}_mlat : magnetic latitude of crests
crest{#}_np : plasma density of crests
min_lat# : latitude of minimum density between EIA crests
            if less than 2 crests, min_lat1 and min_lat2 are Nan
            if 2 crests, min_lat1 is a number and min_lat2 is Nan
            if 3 crests present, min_lat1 is southern and min_lat2 is northern
alt : alitude of satelltie during pass
