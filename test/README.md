# LArPhysicsContent testing

To perform validation of the ntuple mechanics:
1. Run Pandora using the `PandoraSettings_NtupleTest.xml` settings file to produce an ntuple ROOT file; e.g. 

```PandoraInterface -i PandoraSettings_NtupleTest.xml -e [events] -g [geometry] -r [mode]```

2. Run the validation macro `ValidateNtuple.c` on the ntuple ROOT file; e.g. 

```root -l -q 'ValidateNtuple.c("PandoraNtuple.root")'```

The validation macro will run detailed tests on the ntuple to facilitate debugging.
