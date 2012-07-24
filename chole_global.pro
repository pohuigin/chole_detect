;----------------------------------------------------------------------------->

function chole_global, fparam=fparam, proot=proot

retval=''

if keyword_set(fparam) then retval=!CHOLE_PATH+!CHOLE_PARAM ;'~/science/procedures/chole_detect/chole_param.txt'

;if n_elements(pdata) gt 0 then retval='~/science/procedures/chole_detect/data/'

if keyword_set(proot) then retval=!CHOLE_PATH




return, retval

end

;----------------------------------------------------------------------------->