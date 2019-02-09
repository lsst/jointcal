from lsst.pipe.tasks.colorterms import ColortermLibrary, ColortermDict, Colorterm
# this test needs the sdss r-band colorterm
config.applyColorTerms = True
# The values in this library were taken from obs_subaru as of stack version 16.0
config.colorterms = ColortermLibrary(data={
    "sdss*": ColortermDict(data={
        'r': Colorterm(primary="r", secondary="i", c0=0.00231810, c1=0.01284177, c2=-0.03068248)
})
})
