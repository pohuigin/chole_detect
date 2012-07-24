;A setup file that makes the environment variables to run the code
pro chole_setup, chole_path=inchole_path, chole_param=inchole_param

;Set default variables
chole_path='~/science/procedures/chole_detect/'
chole_param='chole_param.txt'
;chole_param='chole_param2.txt'

;Allow for variable names to be set externally
if n_elements(inchole_path) eq 1 then chole_path=inchole_path
if n_elements(inchole_param) eq 1 then chole_param=inchole_param

;Set the environment variables
DEFSYSV, '!CHOLE_PATH', chole_path
DEFSYSV, '!CHOLE_PARAM', chole_param


end