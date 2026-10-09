# DISMO integration validation on Windows

The integration fixes retain the DSSAT 4.8.6.4 CROPGRO and FORAGE interfaces and functionality, register the DISMO GET/PUT variables, initialize the classic pest fallback, and write a text-only DISMO.OUT header. Windows linking no longer uses `/FORCE`, so unresolved symbols fail the build. CMake installation includes the disease parameter file and skips the Unix launcher on Windows.

## Build

Validated with Visual Studio 2022, Intel Fortran Classic 2021.11.0, CMake 3.29.0-rc3, and x64:

```powershell
.\Build-Windows.ps1 -Configuration Release
.\Build-Windows.ps1 -Configuration Debug
```

The executable is `build/bin/<Configuration>/dscsm048.exe`. The script requires the Intel Fortran integration for Visual Studio. Release compiled with zero errors and warnings. Debug compiled with zero errors and two existing warnings about the `SYSTEM` declaration in `Utilities/UTILS.for`.

## Runtime checks performed

These results use local copies of the LOND0501 and TAMA0501 soybean experiments and their associated weather, soil, and calibrated genotype data. They are examples of the validation performed, not universal expected values for different datasets. Those local datasets are not included in this repository.

| Case | Yield (kg/ha) | Maximum severity (%) | DISMO rows |
| --- | ---: | ---: | ---: |
| LOND0501, disease on, PHOTO=L | 3079 | 62.81 | 195 |
| LOND0501, disease off, PHOTO=L | 5684 | — | No DISMO output |
| LOND0501, disease on, PHOTO=C | 3319 | 61.12 | 195 |
| LOND0501, disease off, PHOTO=C | 5516 | — | No DISMO output |
| TAMA0501, disease on, PHOTO=L | 2009 | 2.39 | 166 |

Additional checks:

- `DISES=Y` without `disease_parameters.txt` used the classic pest pathway and produced 5684 kg/ha for the LOND case, matching the disease-off control. The uninitialized factor previously caused zero yield in this case.
- A nonexistent target disease produced the expected DSSAT diagnostic, ERROR.OUT, and exit code 99.
- Two identical LOND entries in one batch produced 390 DISMO rows and identical results per run, confirming season-state reset for that test.
- Debug and Release produced matching yield and severity. Spore-count columns had small floating-point differences within relative tolerance 2e-6 and absolute tolerance 0.001.
- The outputs contained finite values, severity within 0–100%, nondecreasing severity within each run, and no NUL bytes in the header. Debug reported no array-bounds or floating-point exception in these tests.
- CMake configuration, compilation, and installation completed successfully.

Run from the experiment directory with a DSSAT profile and associated input data available. For example:

```powershell
& '<repository>\build\bin\Release\dscsm048.exe' CRGRO048 A LOND0501.SBX 1
```

Use `DISES=Y` and put `disease_parameters.txt` in the current working directory to activate DISMO. Inspect both process exit status and DSSAT diagnostic files when assessing a run.

## Humidity settings and reference comparison

The supplied source uses `USE_WTH_RH=.TRUE.` and `USE_FUNGICIDE=.FALSE.`. Biological parameters were not changed by these integration fixes.

The local reference used DSSAT 4.8.5.9 and `USE_WTH_RH=.FALSE.`, with a spray buffer of 23 days and a DVIP threshold of 1. A temporary build of the corrected integration with those same settings reproduced every DISMO data column and row from the rebuilt reference for LOND, TAMA, PHOTO=C, and the repeated batch. The source was then restored to the fork's supplied settings: weather RH, a 16-day spray buffer, and threshold 6.

Changing the humidity mode can change simulated severity and yield. The buffer and threshold are separate source constants; changing `USE_WTH_RH` alone does not automatically change them.

## DSSAT data compatibility and scope

Maize FLSC8101, alfalfa AGZG1219, and brachiaria CBAN8401 also completed. Maize and alfalfa required the 4.8.6 genotype files provided under `Data/Genotype`. Older installed genotype files caused an input conversion error for maize and a missing `SURVIVING WINTER COLD AND DROUGHT` section for alfalfa; a separate profile pointing to the repository genotype directory resolved those failures.

These checks establish computational behavior for the stated cases. They do not establish scientific calibration for other locations, test every crop or rotation/seasonal mode, or validate fungicide efficacy with `USE_FUNGICIDE` enabled.
