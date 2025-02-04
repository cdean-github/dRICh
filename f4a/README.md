# dRICh Fun4All porting

## Setup and Usage
- first, make sure you can run `fun4all`
  - if you have the `fun4all` `singularity` image, start a new container; you
    can do this with `util/startFun4all.sh /path/to/fun4all/Singularity`, where
    the path points to [Fun4all Singularity](https://github.com/eic/Singularity.git);
    this script will automatically source `eic_setup.sh`, once in the container
  - if your usage of `fun4all` is different, please note that so far this port
    has been developed using the `singularity` container; if you encounter
    issues, send an issue on `github`, or a pull request if you manage to fix
    it
- build the dRICh by running `util/buildDetector.sh` or `util/buildDetector.sh clean`
  - building with the `clean` option will delete the build directory before compiling,
    useful if you want a fresh build
  - libraries and headers will appear in `./install`, if successful
- set environment variables with `source env.sh`
  - most scripts in the `util/` directory do this automatically, so this step
    is not necessary if you use those scripts
  - `echo $LD_LIBRARY_PATH` should include `./install/lib`
- the script `util/shootPion.sh` is being used for testing and development
  - it throws a `pi+` to a specific location on the dRICh
  - run this script with no arguments for a print out of the usage
  - this script will call some macros in `./drich/macros`:
    - `shootPion.C` is responsible for running the simulation
    - `hitHistos.C` will draw some useful histograms
      - if you edit `hitHistos.C`, you can use `./util/rerunHitHistos.sh` to
        quickly rerun this script, without having to rerun the simulation
  - ROOT files that were produced will appear in `./out`
    - `dRIChTree.root` contains the output hit `TTree`
    - `dRIChHists.root` contains histograms
    - `G4HitNtuple.root` contains the DST tree (not tested)
  - if you start the GUI, use the command `/Fun4All/run 1` to shoot one
    particle and visualize

## Development Notes
- see [diagram.pdf](doc/diagram.pdf) for a schematic of the classes and
  inheritance for this `fun4all` dRICh implementation

### changelog for shared code, or ESCalate implementation
This is a summary of changes needed in code outside this
directory, in order to allow the `fun4all` implementation to
build and run and not complain. Hopefully they won't
conflict with the `ESCalate` implementation development.
- `g4dRIChOptics.hh`
  - replace `fmt` dependence with `stdio` functions
    - `fmt` is not in `fun4all` singularity container;
      this is just a quick fix
  - `setOpticalParams` are empty virtual methods with
    non-void return
    - fix: added return dummy value `-1`
  - in `g4dRIChAerogel`, mode 0 (Vorobiev), `rho` has no
    value assigned to it, but it is needed later in
    `scaledE` and `scaledS`
    - quick fix: initialize `rho` with dummy value `0`
    - mode 3 (CLAS12) is currently being used
- `ci_DRICH_Config.hh`
  - comment out `InitializationContext.hh` include directive
    (not sure what this is)
  - add `G4GenericMessenger.hh` include directive
  - should this header be moved to `share/source`?
- `drich-g4model.txt`
  - the construction by text in `fun4all` was getting confused about which
    volume is the mother; as a result, sometimes `fun4all` would take an aerogel
    sector, or a single photo-sensor, as the mother volume, therefore only that volume would be drawn
    - other times, it would correctly take the vessesl as the mother volume, and
      everything is fine
    - regardless of how it runs, several errors regarding the conflict of
      mother volumes were printed during construction, indicating the confusion
  - a simple fix is placing the gas-filled vessel, before defining other
    volumes and petal assemblies; while this practically guaranteed that the
    correct mother volume would be chosen, the conflict errors were still
    appearing
  - a more involved fix, which has been implemented, is the following:
    - `motherName` is now renamed to `vesselName`; this volume is aluminum, and
      will contain all the petals, and we want it to be the mother
      - we don't place `vesselName` directly, but we make sure the first volume
        that we do place has `vesselName` as its mother, so `vesselName` will
        be placed "indirectly" (and first)
    - the original `vesselName`, which was the gas volume, is now called
      `petalName`, and now only spans one petal in azimuth (60 degrees)
      - 6 copies of these are placed, to span the full azimuth, and their mother
        is the new `vesselName` volume
      - the 6 `petalName` copies are touching, but hopefully the boundary between
        each petal volume won't matter, since it is all the same gas
    - then we place the aerogel, filter, mirror, and photo-sensors, all with mother
      `petalName`
      - we only need to place them in one petal, since the petal was already copied;
        these volumes will therefore also be copied to each petal
      - when placing the mirror, the sign of the x coordinate needed to be flipped,
        from what it was when the `VOLU_ASSEMBLY` method was used; it would be good
        to cross check with the previous configuration, to make sure this is correct
    - after these changes, the heirarchy of volumes is more clear to `fun4all`, and
      the construction no longer throws conflict errors and seems to be stable
