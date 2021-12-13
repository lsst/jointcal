# for the cfht_minimal tests with a minSnr threshold
config.doAstrometry = False
config.photometryModel = "simpleFlux"
config.sourceSelector['science'].doSignalToNoise = True
config.sourceSelector['science'].signalToNoise.minimum = 10000
