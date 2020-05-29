from lsst.pipe.tasks.colorterms import ColortermLibrary, ColortermDict, Colorterm
# this test needs the ps1 r-band colorterm
config.applyColorTerms = True
# The values in this library were taken from obs_subaru as of stack version 19.0
config.colorterms = ColortermLibrary(data={
    "ps1*": ColortermDict(data={
        'r': Colorterm(primary="r", secondary="i", c0=-0.000144, c1=0.001369, c2=-0.008380),
})
})
