# Config for gen3 HSC tests, which have different field names because it's older data.

# TODO DM-31889: Remove this file and just use the default config when
# testdata_jointcal is updated to the new column names.

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
