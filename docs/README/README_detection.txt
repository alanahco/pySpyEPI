README for detection_{%Y%m}_{sat}{id}_ascii.txt
-------------------------------------------------------------------------------
COLUMNS
-------
~ the following columns were created in
    pyDetEPI.detection.stats.cand_info.epi_stats ~
        center_time, edge1_time, edge2_time, center_np, edge1_np,
        edge2_np, center_prom, center_bnp, center_mlat, edge1_mlat, edge2_mlat,
        center_mlon, edge1_mlon, edge2_mlon, center_lat, edge1_lat, edge2_lat,
        center_lon, edge1_lon, edge2_lon, center_lshell, edge1_lshell,
        edge2_lshell, center_lt, edge1_lt, edge2_lt, barrel_slope,
        barrel_slope_left, barrel_slope_right, altitude, satellite, fpeak_prom,
        barrel_high_south_ne_max, barrel_high_south_ne_min,
        barrel_high_south_mlat_max, barrel_high_south_mlat_min,
        barrel_mid_south_ne_max, barrel_mid_south_ne_min,
        barrel_mid_south_mlat_max, barrel_mid_south_mlat_min,
        barrel_eq_south_ne_max, barrel_eq_south_ne_min,
        barrel_eq_south_mlat_max, barrel_eq_south_mlat_min,
        barrel_eq_north_ne_max, barrel_eq_north_ne_min,
        barrel_eq_north_mlat_max, barrel_eq_north_mlat_min,
        barrel_mid_north_ne_max, barrel_mid_north_ne_min,
        barrel_mid_north_mlat_max, barrel_mid_north_mlat_min,
        barrel_high_north_ne_max, barrel_high_north_ne_min,
        barrel_high_north_mlat_max, barrel_high_north_mlat_min,


~ the following columns were created in
    pyDetEPI.detection.stats.test_candidates.test_peaks ~
        epi_flag, Ldnp_dif, Rdnp_dif, Hdnp_dif, Hdnp_med, skew, percent_depth,
        bub_in_bub, flagged_np, max_freq_dnp20s, high_freq

~ the following columns were created in
    pyDetEPI.detection.stats.flag_detections.mlat_flag ~
    mlat_flag

~ the following columns were created in
    pyDetEPI.detection.stats.flag_detections.density_flag ~
        den_flag

~ the following columns were created in
    pyDetEPI.detection.stats.cand_info.pass_info ~
        pass_id, pass_time1, pass_time2, pass_mlat1, pass_mlat2,
        pass_max_np, pass_min_np, pass_mlat_max_np, pass_mlat_min_np,
        pass_eqtime, eia_state

-------------------------------------------------------------------------------
Notes
------
    /center/ refers to deepest point in bubble based on plasma density
    /edge1/ refers to the edge of the plasma density halfwidth detected
        by the satellite first in time/chronology.
    /edge2/ refers to the edge of the plasma density halfwidth detected
        by the satellite second in time/chronology.
    When referring to either edge, edge{#} will be used
    When referring to center, edge1, and edge2, {loc} will be used

    detrended density is barrel density - original density
    the magnetic latitude, magnetic longitude, and L shell are calculated
        by apexpy
    background density is barrel density

    Version 1 uses Delta Ne 40s and 10 barrels ranging from largest to smallest
        (for Swarm 80 seconds to 8 seconds)
    Version 2 uses Np / Weighted Barrel and 5 barrels from largest to smallest
        (for Swarm 80 seconds to 16 seconds)
        and Weighted barrel is made of 5 small barrels rolled over and under
        (for Swarm 8 seconds to 40 seconds)
-------------------------------------------------------------------------------
Columns
------------------
"About" parameters
--------------------
{loc}_time : universal time at {loc}
{loc}_np : plasma density at {loc}
center_prom : prominence at center from detrended density
center_bnp : background plasma density from barrel trend
{loc}_mlat : magnetic latitude at {loc}
{loc}_mlon : magnetic longitude at {loc}
{loc}_lat : geographic latitude at {loc}
{loc}_lon : geographic longitude at {loc}
{loc}_lshell : L Shell at {loc}
{loc}_lt : local time at {loc}
barrel_slope : slope between the plasma density at edge1 and at edge 2
                calculated as (edge2_np - edge1_np)/(edge2_tloc - edge1_tloc)
                where tloc is time in seconds. This could help distinguish
                between blobs and bubbles using epi_flag = 9
barrel_slope_left : slope between the plasma density at peak and at edge 1
                calculated as (peak_np - edge1_np)/(peak_tloc - edge1_tloc)
                where tloc is time in seconds. This could help distinguish
                between blobs and bubbles using epi_flag = 9
barrel_slope_right : slope between the plasma density at edge 2 and at peak
                calculated as (edge2_np - peak_np)/(edge2_tloc - peak_tloc)
                where tloc is time in seconds. This could help distinguish
                between blobs and bubbles using epi_flag = 9
altitude : satellite altitude. remains the same between passes, but may change
            over time
satellite : satellite letter or number
fpeak_prom : prominence recorded bu find_peaks
------------------------
EPI detection parameters
------------------------
epi_flag : flag that indicates if a candidate is considered an EPI
            Flags overwrite each other in following order: 3, 9, 8, 4, 2, 7
            A negative flag indicates that > 30% of the flags were 2
                indicating large amounts of erroneous data
                or a bub in bub EPI detection had a bad frequency
                but the outer one did not

            Flag values and meanings:
            1 : EPI
            2 : peak flagged because frequency is > 0.2 Hz at max amplitude
                or a > 0.2 Hz freq of amplitude > 0.5 * max(amplitude)
            3 : peak flagged because of a bubble inside a in bubble
            4 : peak flagged because prominence does not meet threshold (10%)
            7 : peak flagged becuase during the halfwidth, Np is flagged
            8 : peak flagged because the prominence is < 10**4,
                will be overwritten if 4 is also true

            Version 1
            ----------
            9 : peak flagged becuase of the deltaNe40s criteria
                H_dif (halfwidth 95th - 5th percentile of dne)
                L_dif (edge1 to peak 95th - 5th percentile of dne)
                R_dif (peak to edge2 95th - 5th percentile of dne)
                if H_dif < 0.1 AND abs(skew) > 1
                if H_dif < 0.01
                if H_dif < 0.1 AND R_dif/L_dif < 0.1 or > 10 (factor of 10)
                if R_dif or L_dif == 0
            Version 2
            ----------
            9 : peak flagged becuase of the Np / W_B (Weighted Barrel)
                l20 : 20th percentile of EPI left edge (edge 1) to center
                r20 : 20th percentile of EPI center to right edge (edge 2)
                med_l : median of left edge (edge 1) to center
                med_r : median of center to right edge (edge 2)
                lr_15 : 15th percentile of pass Np / W_B
                lr_set : 0.995 (0.5 %)
                Remove if any of the following are True
                (l20 > lr_15)
                (r20 > lr_15)
                (med_l > lr_set)
                (med_r > lr_set)
                * Also recommend keeping if lr_15 is less than 0.95 and
                    l20 and r20 are below 0.95 (can retroactively add this
                    back in)
percent_depth : depth between center and correpsonding barrel point
                (abs(barrel_df["filt_det_np"].iloc[xp]
                   / barrel_df["filt_np"].iloc[xp])) * 100
                where filt_np is barrel density, filt_det_np is detrended
                density, and xp is index at center
bub_in_bub : 1 if an EPI is detected inside of another EPI
            0 means detection is alone or on the outside of another EPI
flagged_np : 1 if there are flagged points based on plasma denisty
            between the edges, otherwise it is 0
max_freq_dnp20s : maximum frequency from FFT of delta Np 20s
high_freq : 1 if a freq is  > 0.2 Hz freq of amplitude > 0.5 * max(amplitude)
            0 if not
barrel_{high, mid, eq}_{south, north}_{ne, mlat}_{max, min} :
        maximum and minimum barrel background trend density (ne) and
        corresponding magnetic latitude (mlat) in high mag latitudes
        (+/- 40 to +/-50) (high), mid mag latitudes (+/-25 to +/-40), and
        equatorial mag latitudes (0 to +/-25) all are inclusive (<= and >=)
        note: NaN if no data present for respective region, nanmax and nanmin
        used
    Version 1
    ---------
    Ldnp_dif : 95th - 5th percentile of dNp 40s from edge1 to center
    Rdnp_dif : 95th - 5th percentile of dNp 40s from center to edge2
    Hdnp_dif : halfwidth 95th - 5th percentile of dNp 40s
    Hdnp_med : median of halfwidth dNp 40s
    skew : skew of dNp 40s from edge1 to edge2
    Version 2
    ----------
    p20_wdr_{h, l, r} : 20th percentile of Np / weighted barrel for
                        halfwidth (h), left side (l), and right side (r) of EPI
    med_wdr_{h, l, r} : median of Np / weighted barrel for
                        halfwidth (h), left side (l), and right side (r) of EPI
    pass_wdr_p{15, 25} : 15th and 20th percentile of Np / weighted barrel
                        for whole pass

----------------
Additional Flags
----------------
mlat_flag : Depletions above 40 degrees need to be treated skeptically
            Additionally if the only depletions are "super", be skeptical
            0 : no flag
            5, 10, 15 : abs(peak mlat) >= 40 degrees
                5 : at least one 1 flag between +/- 25
                10 : no flags between +/- 25
                15 : flags present but no 1 flag between +/- 25
                Value can be X.0, X.1, or X.2
                X.0 means no flags between 25 and 40
                X.1 means a 1 flag between 25 -> 35 AND 35 -> 40
                X.2 means a 1 flag between 25 and 35 but not 35 to 40
                X.3 means a 1 flag between 35 and 40 but not 25 to 35
                X.4 means a non 1 flag between 25 and 40
            4, 8, 12 : abs(peak mlat) < 40 degrees and >= 35
                4 : at least one 1 flag between +/- 25
                8 : no flags between +/- 25
                12 : flags present but no 1 flag between +/- 25
                Value can be X.0, X.1, or X.2
                X.0 means no flags between 25 and 35
                X.1 means a 1 flag between 25 and 35
                X.2 means a non 1 flag between 25 and 35
            3, 6, 9 : abs(peak mlat) < 35 degrees and > 25
                3 : at least one 1 flag between +/- 25
                6 : no flags between +/- 25
                9 : flags present but no 1 flag between +/- 25
            Generally best to keep flags 0, 3, 4.1, and maybe 5.1

----------------
Pass Information
----------------
pass_id : format %Y%m%d_%H%M_{sat}{#} for each satellite pass meant to be used
            as an identifier for detections during the same pass
pass_time1 : starting universal time of pass
pass_time2 : ending universal time of pass
pass_mlat1 : magnetic latitude at pass_time1
pass_mlat2 : magnetic latitude at pass_time2
        pass_mlat1 and pass_mlat2 can be used to determine if the trough dropped
        below 50 degrees magnetic latitude
pass_max_np : maximum plasma density measured during pass
pass_min_np : minimum plasma density measured during pass
pass_mlat_max_np : magnetic latitude of maximum plasma density during pass
pass_mlat_min_np : magnetic latitude of minimum plasma density during pass
pass_eqtime : universal time at mangetic equator for pass
eia_state : EIA category/orientation based on pyValEIA
