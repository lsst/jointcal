# The cfht_minimal data only has the sdss refcat.
config.connections.photometryRefCat = "sdss_dr9_fink_v5b"
config.photometryRefObjLoader.filterMap = {}

# Need a different set of selection flags, as these catalogs don't have all
# the new columns.
config.sourceSelector['science'].flags.bad = ['pixelFlags_edge',
                                              'pixelFlags_saturated',
                                              'pixelFlags_interpolatedCenter',
                                              'pixelFlags_interpolated',
                                              'pixelFlags_crCenter',
                                              'pixelFlags_bad',
                                              'apFlux_12_0_flag',
                                              ]
