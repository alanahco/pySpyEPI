README for irr_{%Y%m}_{sat}{id}_ascii.txt
-------------------------------------------------------------------------------
COLUMNS
-------
~ the following columns were created in
    pyDetEPI.detection.stats.swarm_stats.IRR_stats ~
        center_time, gne100_max, gne100_min, gne100_med, gne100_std, gne50_max,
        gne50_min, gne50_med, gne50_std, gne20_max, gne20_min, gne20_med,
        gne20_std, dne10_max, dne10_min, dne10_med, dne10_std, dne20_max,
        dne20_min, dne20_med, dne20_std, dne40_max, dne40_min, dne40_med,
        dne40_std, rod_max, rod_min, rod_med, rod_std, rodi10_max, rodi10_min,
        rodi10_med, rodi10_std, rodi20_max, rodi20_min, rodi20_med, rodi20_std,
        mvtec_max, mvtec_min, mvtec_med, mvtec_std, mrot_max, mrot_min,
        mrot_med, mrot_std, mroti10_max, mroti10_min, mroti10_med, mroti10_std,
        mroti20_max, mroti20_min, mroti20_med, mroti20_std, tec_std_max,
        tec_std_min, tec_std_med, tec_std_std, ibi_max, ipir_max, irf_center,
        irf_edge1, irf_edge2, pass_id

-------------------------------------------------------------------------------
Notes
------
    IPIR/IRR is a 1 Hz dataset
    {param}_{measure} is used for measure of max (maximum), min (minimum),
        med (median), and std (standard deviation)

    /center/ refers to deepest point in bubble based on plasma density
    /edge1/ refers to the edge of the plasma density halfwidth detected
        by the satellite first in time/chronology.
    /edge2/ refers to the edge of the plasma density halfwidth detected
        by the satellite second in time/chronology.
    When referring to center, edge1, and edge2, {loc} will be used

    Parameter descriptions from Swarm Level 2 Product Defintions
    (file:///Users/alanahco/Downloads/Swarm-Level-2-Product-Definitions.html)

    Grad_Ne@100km	CDF_DOUBLE	The electron density gradient in a running
                    window calculated via linear regression over 27 data points
                    for the 2 Hz electron density data	cm^-3/m
    Grad_Ne@50km	CDF_DOUBLE	The electron density gradient in a running
                    window calculated via linear regression over 13 data points
                    for the 2 Hz electron density data	cm^-3/m
    Grad_Ne@20km	CDF_DOUBLE	The electron density gradient in a running
                    window calculated via linear regression over 5 data points
                    for the 2 Hz electron density data	cm^-3/m
    ROD	CDF_DOUBLE	Rate Of change of Density	cm^-3/s
    RODI10s	CDF_DOUBLE	Rate Of change of Density Index (RODI) is the standard
                deviation of ROD over 10 seconds	cm^-3/s
    RODI20s	CDF_DOUBLE	Rate Of Density Index (RODI) is the standard deviation
                of ROD over 20 seconds	cm^-3/s
    mVTEC	CDF_DOUBLE	Median of VTEC from all available GPS satellites above
                30 degrees	TECU
    mROT	CDF_DOUBLE	Median of Rate Of TEC (ROT) from all available GPS
            satellites above 30 degrees	TECU/s
    mROTI10s	CDF_DOUBLE	Median of Rate Of TEC Index (ROTI) from all
                available GPS satellites above 30 degrees. The ROTI of each
                satellite is the standard deviation of ROT over 10 seconds
                TECU/s
    mROTI20s	CDF_DOUBLE	Median of Rate Of TEC Index (ROTI) from all
                available GPS satellites above 30 degrees. The ROTI of each
                satellite is the standard deviation of ROT over 20 seconds
                TECU/s
    IBI_flag	CDF_INT4	Plasma Bubble Index
                0 no bubble, 1 bubble, -1 not analyzed
    Ionopshere_region_ flag	CDF_INT4	0: equator, 1: mid-latitudes;
                2: auroral oval; 3: polar cap
    IPIR_index	CDF_INT4	0-3 low, 4-5 medium, and > 6 high level of
                fluctuations in the ionospheric plasma density
    TEC_STD	CDF_DOUBLE	Standard deviation of VTEC from GPS satellites	TECU

    deltaNe was recalculated so that the absolute value was not taken
    delta_Ne10s	CDF_DOUBLE	Derived by subtracting Ne by its median filtered
                value in 10 seconds; indicates the electron density
                fluctuations smaller than 75 km	cm^-3
    delta_Ne20s	CDF_DOUBLE	Derived by subtracting Ne by its median filtered
                value in 20 seconds; indicates the electron density
                fluctuations smaller than 150 km	cm^-3
    delta_Ne40s	CDF_DOUBLE	Derived by subtracting Ne by its median filtered
                value in 40 seconds; indicates the electron density
                fluctuations smaller than 300 km	cm^-3
-------------------------------------------------------------------------------
Columns
-------
center_time : universal time at EPI center
gne{km}_{measure} : measure of Grad_Ne@ 100, 50, and 20 km
dne{s}_{measure} : measure of deltaNe 10, 20, and 40s
rod_{measure} : measure of ROD
rodi{s}_{measure} : measure of RODI 10 and 20s
mvtec_{measure} : measure of mVTEC
mrot_{measure} : measure of mROT
mroti{s}_{measure} : measure of mROTI 10 and 20s
tec_std_{measure} : measure of TEC_STD
ibi_max : maximum IBI during EPI halfwidth
ipir_max : maximum IPIR index detected during EPI halfwidth
irf_{loc} : Ionosphere Region Flag for each {loc}
pass_id : format %Y%m%d_%H%M_{sat}{#} for each satellite pass meant to be used
            as an identifier for detections during the same pass







