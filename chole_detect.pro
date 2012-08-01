;----------------------------------------------------------------------------->
;+
; Detect coronal holes using the method from Henley and Harvey 2005 using 
; He II 1083 and magnetogram data.
; 
; USAGE:
;
; ...
;
; HISTORY:
; 	Written - Paul A. Higgins - 20120717
;
; TODO:
; 	1. Error handling
; 	2. Do equal shape projection when segmenting, then project back to HPC
;	3. Use coronal images to validate darkness of CHs?
;	4. Do validation of code by comparing stability of detections to EIT ones (charm/spoca), predicted ones (solis), over time etc.
;	5. Validate by how do the boundaries evolve?
;
;-
;----------------------------------------------------------------------------->

;----------------------------------------------------------------------------->

pro chole_getmag, time, outflist, err
err=-1

pdata=chole_global(/pdata)

;PLACE HOLDER!!
;find and download data to PDATA
;chole_getsolis,/mag,flist=flist
;chole_getsolis,/heii,flist=flist

end

;----------------------------------------------------------------------------->

pro chole_getheii, time, outflist, err
err=-1

pdata=chole_global(/pdata)

;find and download data to PDATA


end

;----------------------------------------------------------------------------->

function chole_fread, flist

mreadfits,flist,ind,dat

;Make a map with the full FITS header in addition to the map-specific keywords
mindex2map, ind, dat, map

;Gets around a bug in DROT_MAP
;TEMP??
map=add_tag(map,map.time,'rtime',/no_copy)

return,map

end

;----------------------------------------------------------------------------->

function chole_getdataparam, time, map

;Create World Coordinate System structure
wcs = fitshead2wcs(map[0])

;Get radius of Sun in arcsecs, particular to the data source
dsunmm=wcs.position.DSUN_OBS/1d6 ;put dist. to sun in Mm
rsunmm=WCS_RSUN(units='Mm') ;radius of sun in Mm
radpasec=2.*!pi/(360.*3600.)
rsunasec=atan(rsunmm/dsunmm)/radpasec

;Determine Degrees per pixel
degperpx=1./(2.*!pi*rsunasec/360./(map[0]).dx)

;Determine Mm per pixel
mmppx=degperpx*(rsunmm*2.*!pi)/360.

dparams = {mmppx:mmppx, degppx:degperpx, rsunasec:rsunasec, wcs:wcs}

;Update time in WCS and map structure to differentially rotated time
for i=0,n_elements(tag_names(dparams.wcs.time))-1 do dparams.wcs.time.(i)=time

;Make DX, DY positive in WCS...
;TEMP??
dparams.wcs.cdelt=abs(dparams.wcs.cdelt)

outdparams=dparams
return, outdparams
end

;----------------------------------------------------------------------------->

;Cut off bad bits of the data and change into equivalent width from intensity if necessary.
function chole_procheii, time, inmheii, params=params, dparams=dparams
mheii=inmheii

;Make intensity images into equivalent line width
if params.dolinwidth then begin
	
	for i=0,n_elements(mheii)-1 do begin
		thisdat=mheii[i].data
		thisdat=thisdat-median(thisdat)
		mheii[i].data=thisdat
	endfor
	
endif

;Differentially rotate images to TIME
mheii=drot_map(mheii, time=time, miss=params.heiimiss)

;Update time in map
mheii.time=time

;Average images together
if n_elements(mheii) gt 1 then begin
	datavg=average(mheii.data, 3, miss=params.heiimiss)
	mheii.data=datavg
endif

return,mheii
end

;----------------------------------------------------------------------------->

;Cut off bad bits of the data and change into equivalent width from intensity if necessary.
function chole_procmag, time, inmmag, params=params, dparams=dparams
mmag=inmmag

;Do cosine correction to magnetic field values
if params.domagcoscor then begin
	COORD = WCS_GET_COORD(dparams.WCS)
	
	;Get coordinate arrays
	xx=reform(coord[0,*,*])
	yy=reform(coord[1,*,*])
	rr=(xx^(2.)+yy^(2.))^(.5)

	coscor=rr
	coscor=1./cos(asin(coscor/dparams.rsunasec))
	coscor[where(rr gt dparams.rsunasec)]=1
	;Limit correction to maximum area a pixel can cover on edge of disk
	thetalim=asin(1.-(dparams.wcs.cdelt)[0]/(dparams.rsunasec)) ;should it be half a pixel or a full pixel? asin(1.-(wcs.cdelt)[0]/(2.*rsunasec))
	coscor=coscor < 1./cos(thetalim)
	
	for i=0,n_elements(mmag)-1 do begin
		thisdat=mmag[i].data
		thisdat=thisdat*coscor
		mmag[i].data=thisdat
	endfor
	
endif

;Differentially rotate images to TIME
mmag=drot_map(mmag, time=time, miss=params.magmiss)

;Update time in map to differentially rotated time
mmag.time=time

;Average images together
if n_elements(mmag) gt 1 then begin
	datavg=average(mmag.data, 3, miss=params.magmiss)
	mmag.data=datavg
endif

return,mmag
end

;----------------------------------------------------------------------------->

;Remap the data into a stereographic projection
function chole_project, inmap, params=params, dparams=dparams, $
	inverse=inverse, inproject=inproject, instgx=instgx, instgy=instgy, $ ;input all of the above to produce an inverse projection
	outstgx=outstgx, outstgy=outstgy, $ ;output these to  allow the inverse projection later on
	outprojlon=outprojlon, outprojlat=outprojlat, outprojdc=outprojdc ;Set these to get the projected HG lat. lon. arrays for overplotting grid lines on the projection

map=inmap
image=map.data

wcs=dparams.wcs

;stop

;Size of input image
imgsz=size(image,/dim)

;Determine the coordinate arrays for the image
COORD = WCS_GET_COORD(WCS) ;coord in sol-x and sol-y
xx=reform(coord[0,*,*])
yy=reform(coord[1,*,*])

;stop

;WTF! WHY is XX a mirror image???
;xx=-xx
;yy=-yy

if keyword_set(inverse) then begin
	;Setup input arrays
	projdat=inproject
	stgx=instgx
	stgy=instgy
	imgszp=size(projdat,/dim)

	;Determine matching tie points between projection and deprojection
	wmatch=where(finite(stgx) eq 1)
	winsx=wmatch mod imgsz[0] ;positions in the output (HPC) image
	winsy=wmatch/imgsz[1]
	wpx=stgx[wmatch] ;positions in the input (stereographic) image
	wpy=stgy[wmatch]
	
	;De-project and interpolate to get back to an HPC projection
	deproj = reform( WARP_TRI( winsx, winsy, wpx, wpy, projdat, output_size=imgsz ))
	deprojmap=map
	deprojmap.data=deproj
	
	return,deprojmap
	
endif

;Get radial array of degrees from D.C. as if image was looking down on the pole.
rrfrac=((xx^(2.)+yy^(2.))^(.5))/dparams.rsunasec

;Find minimum and maximum X and Y pixel bounds for limb in image
limbmask=fltarr(imgsz[0],imgsz[1])
limbmask[where(rrfrac le 1)]=1.
xbound=minmax(where(total(limbmask,2)<1 gt 0))
ybound=minmax(where(total(limbmask,1)<1 gt 0))

;Identify off disk pixels
woff=where(rrfrac gt 1.)
rrfrac[woff]=params.nan
rrdeg=asin(rrfrac)/!dtor*(-1.)+90.

;Get rid of pixels within X number of degrees from the edge (assuming center of image is 90)
rrdeg[where(rrdeg lt (90.-params.threshlon))]=params.nan

;Get azimuthal coordinates
theta=rect2pol(xx,yy)
thsz=size(theta,/dim)
theta=theta[thsz[0]/2:*,*]
;Shift values to positive degrees (and change direction of theta)
theta=(theta-min(theta))/!dtor
theta[woff]=params.nan

;Project to stereographic coordinates
STEREOGRAPHIC, theta, rrdeg, stgx, stgy ;output is normalised to limb of Sun

;Account for clipping of edges (stretch a bit further)
stgx=stgx*(90./params.threshlon)
stgy=stgy*(90./params.threshlon)

;Re-scale the projection coordinates -> scale values to an array a specified factor larger than the original image
stgx=(((stgx+1.)/2.)*(xbound[1]-xbound[0])+xbound[0])*params.projscl
stgy=(((stgy+1.)/2.)*(ybound[1]-ybound[0])+ybound[0])*params.projscl

;Insert data values into locations in projected array determined by stereographic X, Y arrays
wins=where(finite(stgx) eq 1)
winsx=stgx[wins]
winsy=stgy[wins]

;Project and Interpolate the image
projdat = reform( WARP_TRI( winsx, winsy, (wins mod imgsz[0]), (wins/imgsz[1]), map.data, output_size=imgsz*params.projscl ))

if keyword_set(params.doplot) then begin
	;Project the coordinate arrays
	wcs_convert_from_coord, wcs, coord, ’hg’, hglon, hglat 
	outprojlon = reform( WARP_TRI( winsx, winsy, (wins mod imgsz[0]), (wins/imgsz[1]), hglon, output_size=imgsz*params.projscl ))
	outprojlat = reform( WARP_TRI( winsx, winsy, (wins mod imgsz[0]), (wins/imgsz[1]), hglat, output_size=imgsz*params.projscl ))
	;HG Angle from disk center (latitute if looking down on pole)
	outprojdc = reform( WARP_TRI( winsx, winsy, (wins mod imgsz[0]), (wins/imgsz[1]), rrdeg, output_size=imgsz*params.projscl ))

	chole_plotproj, projheii, projmag, projlon, projlat, projdc
	
;	rrdeg=asin(rrfrac)/!dtor*(-1.)+90.
;	outprojdc = reform( WARP_TRI( winsx, winsy, (wins mod imgsz[0]), (wins/imgsz[1]), rrdeg, output_size=imgsz*params.projscl ))
;	!p.multi=[0,2,1]
;	plot_image,map.data>(-100)<100
;	contour,rrdeg,lev=reverse([89,80,70,60,50,40,30,20,10]),c_lines=0,color=0,/over
;	contour,rrdeg,lev=reverse([89,80,70,60,50,40,30,20,10]),c_lines=2,color=255,/over
;	plot_image,rot(projdat,180)>(-100)<100
;	contour,outprojdc,lev=reverse([89,80,70,60,50,40,30,20,10]),c_lines=0,color=255,/over
;	contour,outprojdc,lev=reverse([89,80,70,60,50,40,30,20,10]),c_lines=2,color=0,/over	
	
;	stop
endif

;Output the stereographic position arrays to allow the inverse projection later on
outstgx=stgx
outstgy=stgy

;Output the projections
;dparams.wcs=wcs
outproj=projdat
return,outproj
end

;----------------------------------------------------------------------------->

;Segment the He II image
function chole_segment,indheii, params=params, dparams=dparams
dheii=indheii
imgsz=size(dheii,/dim)

;filter NaNs
wnan=where(dheii eq params.nan)
if wnan[0] ne -1 then dheii[wnan]=0.

;Pull out the non-zero values into an array
wfilt=where(dheii ne 0)
heiifilt=(dheii)[wfilt]

;Determine the mode of values
heiidist = Histogram(heiifilt, MIN=Min(heiifilt))
heiimode = (Where(heiidist EQ Max(heiidist)) + Min(heiifilt))[0]

;Determine segmentation threshold for the image
segthresh=0.1*median(heiifilt[where(heiifilt gt 0)])
;segthresh=median(heiifilt)+0.1*(abs(min(heiifilt))+heiimode);median(heiifilt))
;;segthresh=median(heiifilt)+stddev(heiifilt)
;segthresh=heiimode+stddev(heiifilt)


;Segment the image
chmask=fltarr(imgsz[0],imgsz[1])
;chmask[where(dheii le segthresh)]=0.
chmask[where(dheii gt segthresh)]=1.

loadct,0
plot_image,dheii>(-100)<100
setcolors,/sys,/sil
contour,chmask,level=0.5,/fill,color=!red,/over
stop

;Close the detection
kernalclose=dist(params.kernclosrad/dparams.mmppx)
kernalclose[where(kernalclose le params.kernclosrad/dparams.mmppx)]=1.
kernalclose[where(kernalclose gt params.kernclosrad/dparams.mmppx)]=0.
chmaskgrw = MORPH_CLOSE(chmask, kernalclose)

;Number the detections
chmaskgrw = LABEL_REGION(chmaskgrw)

stop

;Necessary??
;Remove regions that are too small
yhistarea=histogram(chmaskgrw,bin=1.,locations=xhistarea)
wsmall=where(yhistarea lt params.thresharea/(dparams.mmppx)^2.)
if wsmall[0] ne -1 then begin
	match,chmaskgrw,yhistarea[wsmall],wremo
	chmaskgrw[wremo]=0
endif

stop

;Necessary??
;Mask out blobs not in grown contours
chmaskgrw[where(chmaskgrw gt 0)]=1.
chmask=LABEL_REGION(chmask)
for i=1,max(chmask) do begin
	wreg=where(chmask eq i)
	wovlap=where(chmaskgrw[wreg] eq 1)
	if wovlap[0] eq -1 then chmask[wreg]=0.
endfor
chmask=chmask < 1

stop

;Open the detection
kernalopen=dist(params.kernopenrad/dparams.mmppx)
kernalopen[where(kernalopen le params.kernopenrad/dparams.mmppx)]=1.
kernalopen[where(kernalopen gt params.kernopenrad/dparams.mmppx)]=0.
chmaskgrw = MORPH_OPEN(chmask, kernalopen)



stop

;Update numbering of the detections
chmaskgrw = LABEL_REGION(chmaskgrw)

return,outdetmask

end

;----------------------------------------------------------------------------->

;Filter the blobs by their unipolarity
function chole_validate, indetmask, inmag, params=params, dparams=dparams
detmask=indetmask
mag=inmag

for i=1,max(detmask) do begin
	
	;Determine unipolarity
	wthis=where(detmask eq i)
	thispol=mag[wthis]
	wpos=where(thispol gt 0)
	wneg=where(thispol lt 0)
	unipol=(n_elements(wpos) > n_elements(wneg))/float(n_elements(thispol))
	
	;Filter using unipolarity threshold
	if unipol lt params.threshuni then detmask[wthis]=0
	
endfor

;Update numbering
detmask=LABEL_REGION(detmask < 1)

outdetmask=detmask
return, outdetmask
end

;----------------------------------------------------------------------------->

;Plot the projected data arrays
pro chole_plotproj, projheii, projmag, projlon, projlat, projdc

window,xs=1000,ys=500
!p.multi=[0,2,1]

plot_image,projheii>(-100)<100
contour,projdc,level=[20,30,40,50,60,70,80],c_lines=2,color=0,/over
contour,projlon,level=0,c_lines=2,color=0,/over
contour,projlat,level=0,c_lines=2,color=0,/over

plot_image,projheii>(-100)<100
contour,projdc,level=[20,30,40,50,60,70,80],c_lines=2,color=0,/over
contour,projlon,level=0,c_lines=2,color=0,/over
contour,projlat,level=0,c_lines=2,color=0,/over

end

;----------------------------------------------------------------------------->

;Input the desired detection time, 1 or more days worth of He II and LOS magnetogram data
pro chole_detect, time, fheii, fmag, params, _extra=_extra, $
	outstruct=outstruct, outfile=outfile, err=err

;Download files if need be
if n_elements(fheii) eq 0 then $
	chole_getheii, time, fheii, errfheii
;if not errfheii then save error code

if n_elements(fmag) eq 0 then $
	chole_getmag, time, fmag, errfmag
;if not errfmqg then save error code

;Read in the He II files
mheii=chole_fread(fheii)

;Read in the Magnetograms
mmag=chole_fread(fmag)

;Generate coordinate and parameters from the data to be used later
dparamsheii=chole_getdataparam(time, mheii)
dparamsmag=chole_getdataparam(time, mmag)

;Process the He II files
;differentially rotate images in stack to given input time
;generate an averaged line-width and averaged image from the stack of images
mheii=chole_procheii(time, mheii, params=params, dparams=dparamsheii)

;Process the Magnetograms
;generate averaged image and los correct the B vals
mmag=chole_procmag(time, mmag, params=params, dparams=dparamsmag)

;Remap to equal shape or area projection
if params.doproject then begin
	projheii=chole_project(mheii, params=params, dparams=dparamsheii, outstgx=heiistgx, outstgy=heiistgy) ;, outprojlon=projlon, outprojlat=projlat, outprojdc=projdc)
	projmag=chole_project(mmag, params=params, dparams=dparamsmag, outstgx=magstgx, outstgy=magstgy)

stop

endif else begin
	projheii=mheii.data
	projmag=mmag.data
endelse	

;Make binary detections of CHs on line-width image
projdetmask=chole_segment(projheii, params=params, dparams=dparamsheii)

;Compare detections to magnetic fields to filter non-unipolar detections
projdetmask=chole_validate(projdetmask, projmag, params=params, dparams=dparamsmag)

;De-project from equal shape or area projection back to HPC
if params.doproject then begin
	detmask=chole_project(projdetmask, params=params, dparams=dparamsheii, /inverse)
endif

;Output detection mask, chain code of detections, averaged line-width image and averaged LOS magnetogram differentially rotated to the same time
mdetmask=mheii & mdetmask.data=detmask
outstruct = {mask:mdetmask, heii:mheii, mag:mmag}
;outstruct = {mask:'', chain:'', heii:'', mag:''}

;if keyword_set(outfile) then save,outstruct,file=''
;Save in fits format? compressed? save chain codes in ascii?

stop

end

;----------------------------------------------------------------------------->