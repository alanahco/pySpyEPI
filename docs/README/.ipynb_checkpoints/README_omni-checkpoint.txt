README for omni_{%Y%m}_{sat}{id}_ascii.txt
-------------------------------------------------------------------------------
COLUMNS
-------
~ the following columns were created in
    pyDetEPI.omni_write.sw_epi ~
        center_time, epi_by_gsm, epi_bz_gsm, epi_flow_speed,
        epi_vx, epi_vy, epi_vz, epi_proton_density,
        epi_pressure, epi_efield, epi_plasma_beta, epi_bsn_x,
        epi_bsn_y, epi_bsn_z, epi_ae_index, epi_al_index,
        epi_au_index, epi_symh, time_sym_dip_b4,
        time_symh_min_b4, symh_min_b4, time_sym_rec_b4,
        time_sym_dip_af, time_symh_min_af, symh_min_af,
        time_sym_rec_af, time_al_dip_b4, time_al_min_b4,
        al_min_b4, time_al_rec_b4, time_al_dip_af,
        time_al_min_af, al_min_af, time_al_rec_af
-------------------------------------------------------------------------------
Notes
------
    if a column starts with epi_, then the saved value is the closest time to
        the EPI center
    if a column contains b4 or af, it is related to the supplemental/solarwind/
        yearly files. For these parameters, symh and al will be referred to as
        param.
-------------------------------------------------------------------------------
Columns
-------

Closest to EPI centers
----------------------
center_time : universal time at epi center
epi_by_gsm : By (nT), GSM, determined from post-shift GSE components
epi_bz_gsm : Bz (nT), GSM, determined from post-shift GSE components
epi_flow_speed : solar wind Flow Speed (km/s), GSE
epi_vx : solar wind Vx Velocity (km/s), GSE
epi_vy : solar wind Vy Velocity (km/s), GSE
epi_vz : solar wind Vz Velocity (km/s), GSE
epi_proton_density : solar wind Proton density (n/cc)
epi_pressure : solar wind Flow pressure (nPa)
epi_efield : solar wind Electric Field (mV/m)
epi_plasma_beta : solar wind Plasma beta
epi_bsn_x : Bow Shock Nose (Re) location, X, GSE
epi_bsn_y : Bow Shock Nose (Re) location, Y, GSE
epi_bsn_z : Bow Shock Nose (Re) location, Z, GSE
epi_ae_index : AE - 5-minute AE-index, from WDC Kyoto
epi_al_index : AL - 5-minute AL-index, from WDC Kyoto
epi_au_index : AU - 5-minute AU-index,from WDC Kyoto
epi_symh : SYM/H - 5-minute SYM/H index,from WDC Kyoto

Yearly file columns
-------------------
time_{param}_dip_b4 : closest time before center EPI when param dipped below
                        -100 for AL index or -30 for SYM H
time_{param}_min_b4 : time when param minimum occurred related to dip_b4
{param}_min_b4 : minimum param at time_{param}_min_b4
time_{param}_rec_b4 : time when param recovered to above -100 for AL and -30
                        for SYM H
time_{param}_dip_af : closest time after center EPI when param dipped below
                        -100 for AL index or -30 for SYM H
time_{param}_min_af : time when param minimum occurred related to dip_af
{param}_min_af : minimum param at time_{param}_min_af
time_{param}_rec_af : time when param recovered to above -100 for AL and -30
                        for SYM H
