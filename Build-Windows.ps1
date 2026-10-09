param(
    [ValidateSet('Release', 'Debug')]
    [string]$Configuration = 'Release'
)
$ErrorActionPreference = 'Stop'
$repoRoot = $PSScriptRoot
$buildRoot = Join-Path $repoRoot 'build'
& cmake -S $repoRoot -B $buildRoot -G 'Visual Studio 17 2022' -A x64 `
    -DCMAKE_BUILD_TYPE=Release `
    '-DCMAKE_Fortran_FLAGS=/nologo /fpp /libs:static /threads /fpe:0 /extend-source:132 /Qdiag-disable:10448' `
    '-DCMAKE_Fortran_FLAGS_RELEASE=/O2 /DNDEBUG /Qunroll /Ob2' `
    '-DCMAKE_Fortran_FLAGS_DEBUG=/Od /debug:full /dbglibs /warn:all /traceback /check:bounds' `
    '-DCMAKE_Fortran_FLAGS_TESTING=/O2' `
    '-DCMAKE_EXE_LINKER_FLAGS=/machine:x64'
if ($LASTEXITCODE -ne 0) { throw 'CMake configuration failed.' }
& cmake --build $buildRoot --config $Configuration --parallel 4
if ($LASTEXITCODE -ne 0) { throw 'DSSAT build failed.' }
Write-Host ('Executable: ' + (Join-Path $buildRoot "bin/$Configuration/dscsm048.exe"))
