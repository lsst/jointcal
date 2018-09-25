from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
config.astrometryRefObjLoader.retarget(LoadIndexedReferenceObjectsTask)
config.astrometryRefObjLoader.ref_dataset_name = "gaia"
config.astrometryRefObjLoader.filterMap = {'u': 'phot_g_mean_mag',
    'g': 'phot_g_mean_mag',
    'r': 'phot_g_mean_mag',
    'i': 'phot_g_mean_mag',
    'z': 'phot_g_mean_mag',
    'y': 'phot_g_mean_mag',
}
