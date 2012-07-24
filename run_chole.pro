;A wrapper to run the coronal hole detection code
;Optional:
; Either supply these:
;	chole_path = the path pointing to the CHOLE root directory (where the parameter file etc is sitting)
;	chole_param = the full path file name of the parameter file to be used
; Or this:
;	param = a structure holding the parameters for running the code
pro run_chole, time=time, fheii=fheii, fmag=fmag, $
	chole_path=chole_path, chole_param=chole_param, params=inparams

;TEMP!
fheii='~/science/procedures/chole_detect/data/heii/svsm_e3100_S2_20120706_2007.fts.gz'
fmag='~/science/procedures/chole_detect/data/mag/svsm_m1100_S2_20120706_1728.fts.gz'
time=anytim('06-jul-2012 16:00',/vms)

;Define paths and files to run the code
chole_setup, chole_path=chole_path, chole_param=chole_param

;Load the parameter file- read in the parameters for running the code
if n_elements(params) eq 1 then params=inparams $
	else params=chole_loadparam()


;Run the detection code
chole_detect, time, fheii, fmag, params, _extra=_extra, $
	outstruct=outstruct, outfile=outfile, err=err


end