;Download a magnetogram or He II eq. width from solis for a given date
pro chole_getsolis, heii=heii, mag=mag, flist=flist

;mag:
fmag='ftp://solis.nso.edu/synoptic/level2/vsm/2012/07/svsm_m1100_S2_20120706_1728.fts.gz'

;heii:
fheii='ftp://solis.nso.edu/synoptic/level2/vsm/2012/07/svsm_e3100_S2_20120706_2007.fts.gz'

;ends up being wrong by seconds-hrs... 
;datevms=anytim(1342764000.+anytim('1-jan-1970'),/vms)



;sock_copy

;flist=








end