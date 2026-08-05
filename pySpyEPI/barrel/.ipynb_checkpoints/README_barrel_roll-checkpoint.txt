barrel_roll.py README

barrel_df columns from get_barrel_info function:
x : x axis variable such as time or lat or lon
ne : original density untouched
barrel_ne : barrel density filtered by savitzky-golay filter (optional)
dne_dd : barrel_ne - ne (dd = density difference)
dne_dd_filt : filtered dne_dd
ne_scale : scaled density
x_scale : scaled x axis

barrel_df columns from triple_barrel function:

from large upper barrel:
ne : original density untouched
barrel_ne : barrel density filtered by savitzky-golay filter (optional)
dne_dd : barrel_ne - ne (dd = density difference)
dne_dd_filt : filtered dne_dd
ne_scale : scaled density
x_scale : scaled x axis

plus:
weight_barrel : weighted barrel from narrow barrels above and below density
dne_wdr : ne / weighted_barrel for weighted density ratio
narlo_barrel : barrel_ne from lower narrow barrel
narup_barrel : barrel_ne from upper narrow barrel

