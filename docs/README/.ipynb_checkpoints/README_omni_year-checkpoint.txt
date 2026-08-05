README for supplemental/solarwind/solarwind_{%Y%m}_{sat}{id}_ascii.txt
-------------------------------------------------------------------------------
COLUMNS
-------

    pyDetEPI.detection.stats.solarwind_stats.swind_year and saved by
    pyDetEPI.depi_io.write.omni_write.sw_year_save

    time_start, time_end, time_{min or max}, {param}_{min or max}

-------------------------------------------------------------------------------
Notes
------
    ~ {param} can be any OMNI solar wind parameter.
    ~ {max or min} depends on if below is True (min) or False (max) in
        swind_year
    ~ Currently, the saved parameters are the AL index (al_index)
        and SYM-H (symh), which below is True so the minimum values are saved
    ~ set minimum or maximum by user, should probably have saved that
        for SYM-H, -30 and for AL index -100
    ~ Parentheses are used in column description for the case of if below is
        False
-------------------------------------------------------------------------------
Columns
-------
time_start : time when {param} first went below (above) set minimum (maximum)
time_end : time when {param} went back above (below) set minimum (maximum)
time_{min or max} : time of deepest minimum (highest maximum) between
                    time_start and time_end
{param}_{min or max} : deepest mininum (highest maximum) {param} value