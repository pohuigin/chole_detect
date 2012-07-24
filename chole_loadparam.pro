;----------------------------------------------------------------------------->

;Load the parameter meta data file for the detection algorithm.
function chole_loadparam, fparam=fparamin, meta=outmeta

;Determine name of meta data file
if not keyword_set(fparamin) then fparam=chole_global(/fparam) else fparam=fparamin

;Read parameters from meta data file
readcol, fparam, param, val, type, meta, comment='#', format='A,A,A,A', delim=';'
param=strtrim(param,2)
val=strtrim(val,2)
type=strtrim(type,2)
meta=strtrim(meta,2)

;Make array of data types for each field in structure
dataspec=strjoin(type,',')

;Create empty structure
create_struct, paramstruct, '', param, dataspec

;Fill the structure
for i=0,n_elements(param)-1 do paramstruct.(i)=val[i]

;Output the description of each parameter
outmeta=[[param],[meta]]

return, paramstruct

end

;----------------------------------------------------------------------------->