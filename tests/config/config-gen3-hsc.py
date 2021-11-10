# Config for gen3 tests, which cannot use astrometrySourceSelector, and which
# do source selection on Parquet tables, which have different field names.

# TODO DM-31889: Update the field names to the new names when
# testdata_jointcal is updated to the new column names.

# TODO DM-29008: move these configs into the jointcal defaults before gen2 removal.
# make science selector look like astrometry selector (see note in tests/config/config.py)
config.sourceSelector['science'].flags.bad = ['PixelFlags_edge',
                                              'PixelFlags_saturated',
                                              'PixelFlags_interpolatedCenter',
                                              'PixelFlags_interpolated',
                                              'PixelFlags_crCenter',
                                              'PixelFlags_bad',
                                              'HsmPsfMoments_flag',
                                              'ApFlux_12_0_flag',
                                              ]
config.sourceSelector['science'].doUnresolved = False
config.sourceSelector['science'].signalToNoise.fluxField = 'ApFlux_12_0_instFlux'
config.sourceSelector['science'].signalToNoise.errField = 'ApFlux_12_0_instFluxErr'
config.sourceSelector['science'].isolated.parentName = 'parentSourceId'
config.sourceSelector['science'].isolated.nChildName = 'Deblend_nChild'
config.sourceFluxType = "ApFlux_12_0"

# TODO DM-27843: workaround for gen3 not supporting anyFilterMapsToThis
config.astrometryRefObjLoader.filterMap = {'g': 'phot_g_mean',
                                           'r': 'phot_g_mean',
                                           'i': 'phot_g_mean',
                                           'z': 'phot_g_mean',
                                           'y': 'phot_g_mean',
                                           'N921': 'phot_g_mean',
                                           }
config.astrometryRefObjLoader.anyFilterMapsToThis = None
