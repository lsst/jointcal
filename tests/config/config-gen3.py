# Config for gen3 tests, which cannot use astrometrySourceSelector, and which
# do source selection on Parquet tables, which have different field names.

# TODO DM-29008: move these configs into the jointcal defaults before gen2 removal.
# make science selector look like astrometry selector (see note in tests/config/config.py)
config.sourceSelector['science'].flags.bad = ['pixelFlags_edge',
                                              'pixelFlags_saturated',
                                              'pixelFlags_interpolatedCenter',
                                              'pixelFlags_interpolated',
                                              'pixelFlags_crCenter',
                                              'pixelFlags_bad',
                                              'hsmPsfMoments_flag',
                                              'apFlux_12_0_flag',
                                              ]
config.sourceSelector['science'].doUnresolved = False
config.sourceSelector['science'].signalToNoise.fluxField = 'apFlux_12_0_instFlux'
config.sourceSelector['science'].signalToNoise.errField = 'apFlux_12_0_instFluxErr'
config.sourceSelector['science'].isolated.parentName = 'parentSourceId'
config.sourceSelector['science'].isolated.nChildName = 'deblend_nChild'
config.sourceFluxType = 'apFlux_12_0'

# TODO DM-27843: workaround for gen3 not supporting anyFilterMapsToThis
config.astrometryRefObjLoader.filterMap = {'g': 'phot_g_mean',
                                           'r': 'phot_g_mean',
                                           'i': 'phot_g_mean',
                                           'z': 'phot_g_mean',
                                           'y': 'phot_g_mean',
                                           'N921': 'phot_g_mean',
                                           }
config.astrometryRefObjLoader.anyFilterMapsToThis = None
