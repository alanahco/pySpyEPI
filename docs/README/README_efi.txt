README for efi_{%Y%m}_{sat}{id}_ascii.txt
-------------------------------------------------------------------------------
COLUMNS
-------
~ the following columns were created in
    pyDetEPI.detection.stats.swarm_stats.EFI_stats ~
        center_time, center_te, center_te_std, edge1_te, edge1_te_std,
        edge2_te, edge2_te_std, te_med, te_std, te_percnan, lp_flag,
        time_lpsweep, center_dte_max, edge1_dte_max, edge2_dte_max,
        center_dte_min, edge1_dte_min, edge2_dte_min, center_dte_med,
        edge1_dte_med, edge2_dte_med

-------------------------------------------------------------------------------
Notes
------
    EFI is a 2 Hz dataset
    {param}_{measure} is used for measure of max (maximum), min (minimum), and
        med (median)

    /center/ refers to deepest point in bubble based on plasma density
    /edge1/ refers to the edge of the plasma density halfwidth detected
        by the satellite first in time/chronology.
    /edge2/ refers to the edge of the plasma density halfwidth detected
        by the satellite second in time/chronology.
    When referring to center, edge1, and edge2, {loc} will be used
    temperature refers to electron temperature
    for 0702 files (update), Te_error may be empty
-------------------------------------------------------------------------------
Columns
-------
center_time : universal time at EPI center
{loc}_te : temperature at {loc} (3 point average of before, {loc}, and after)
{loc}_te_std : standard deviation of 3 point temperature at {loc}
te_med : median temperature of EPI halfwidth
te_std : standard deviation of temepratuer of EPI halfwidth
te_percnan :  percentage of NaN temperature values between edges
lp_flag : LP flag for high gain and low gain probes
            if a sweep occurred during halfwidth, lp_flag = 9
            if all high gain, lp_flag = 1
            if all low gain, lp_flag = 3
            if only original datapoint used (0602), lp flag = 7
            if none of the above are met (or mixed), lp_flag = 2
            Sweep (9) takes precidence
        When analyzing Te, best to only consider lp_flag = 1
te_jump : if te jumped from one point to the next (not continuous) more than
        500 K during epi, then te_jump is 1 otherwise it is 0 after removing
        outliers but before a median filter
te_dif_b4 : temperature difference from one point to the next before median
            filtering
te_dif_af : temperature difference from one point to the next after median
            filtering. This with te_dif_b4 could help determine if the te_jump
            is still useable.
time_lpsweep : seconds between last LP sweep and first/left edge
                will be 0 seconds if an LP sweep occurs during depletion
{loc}_dte_{measure} : measure of Te error at {loc}
fnum : file number string. If it is the old file format, this number will be
        '0602', and if it is the new format, the number will be '0701', this
        will help when dte is smaller because original error is not taken into
        account because error is empty in new files
pass_id : format %Y%m%d_%H%M_{sat}{#} for each satellite pass meant to be used
            as an identifier for detections during the same pass. Will match
            detection file.