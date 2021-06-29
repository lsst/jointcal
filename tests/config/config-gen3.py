# Config for gen3 tests, which cannot use astrometrySourceSelector, and which
# do source selection on Parquet tables, which have different field names.

# TODO DM-29008: move these configs into the jointcal defaults before gen2 removal.
# make science selector look like astrometry selector (see note in tests/config/config.py)
config.sourceSelector['science'].flags.bad = ['PixelFlags_edge',
                                              'PixelFlags_saturated',
                                              'PixelFlags_interpolatedCenter',
                                              'Centroid_flag',
                                              'PsfFlux_flag',
                                              'PixelFlags_suspectCenter',
                                              ]
config.sourceSelector['science'].doUnresolved = True
config.sourceSelector['science'].unresolved.name = 'extendedness'
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
