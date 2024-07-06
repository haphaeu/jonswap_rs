# JONSWAP

Lib and CLI to calculate JONSWAP wave spectrum and its realisation in time domain.

## Use

```
Use:
  jonswap hs tp [-n nharms] [-g gamma] [-s seed] [-d duration] [-t timestep] [-h]
where:
  hs: significant wave height in meters.
  tp: spectral peak period in seconds
  -n: number of harmonics used to discretise the spectrum.
  -g: value for gamma. If ommited, DNV is used.
  -s: seed number for phase randomisation.
  -d: duration - timetrace will be shown.
  -t: time step. default is 0.1 seconds.
  -h: show this help.
```
