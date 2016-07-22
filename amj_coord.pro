pro amj_coord, image_in, hdr_in, CRD_out, instr, seg_const=seg_const, display=display, info=info, correct = correct
;+
; NAME:
;       amj_coord
;
; PURPOSE:
;       Using a given magnetogram and header calculate heliographic coordinates,
;       correct for LOS components, and calculate elements of area and elements of flux
;       
; CALLING SEQUENCE:
;       amj_coord, image_in, hdr_in, CRD_out, instr
;   
;   
; INPUTS:   
;       image_in: Magnetogram image
;
;       hdr_in: header file of the input magnetogram
;
;   	instr: 	Variable that indicates which instrument is used:
;   			1. KPVT-512 files
;   			2. SPMG files
;   			3. MDI files
;   			4. HMI files
;   
; OPTIONAL INPUTS:
;       none
;   
;   
; OUTPUTS: 
;       CRD_out: a structure holding the original and corrected images, as well as the calculated
;                heliographic coordinates
;            1.  im_raw: original magnetogram
;            2.  im_crr: magnetogram corrected using LOS component and assuming al field is radial
;            3.  hdr: header file of the input magnetogram
;            4.  mgnt_ar: elements of area corresponding to each pixel
;            5.  mgnt_flx_raw: flux elements corresponding to each pixel - uncorrected field
;            5.  mgnt_flx_corr: flux elements corresponding to each pixel - corrected field
;            6.  Xar: x coordinate in a heliographic coordinate system 
;            7.  Yar: y coordinate in a heliographic coordinate system 
;            8.  Zar: z coordinate in a heliographic coordinate system
;            9.  Lath: heliographic latitude
;           10.  Lonh: heliographic longitude
;             
;       
; KEYWORDS:
;       display: plots intermediate steps
;       info: prints process information on screen
;       
; MODIFICATION HISTORY:
;   2012/07/11:  Andres Munoz-Jaramillo:   initiated.
;   2015/09/16:  Andres Munoz-Jaramillo:   Adapted to work with different instruments
;-

print, 'Updating supporting variables'

;define the usage
if n_params() lt 4 then begin
    print,'USAGE:  amj_coord, image_in, hdr_in, CRD_out, instr'
    return
endif


;get the input image and dimension
im=image_in
sz=size(im)

;Constant with parameters for plotting, make sure the values are the same as in pnr_dt.rpo
seg_const_def={k_sig:15.0, valid_range:[-20000.,20000.], deg_lim:75.0}

if not keyword_set(seg_const) then begin
    seg_const=seg_const_def
endif 

;
;set the display window size and display thresholds ----------------------------------------------------------

display_zoom=0.5
display_xsize=sz[1]*display_zoom & display_ysize=sz[2]*display_zoom  

print_zoom=1.0
print_xsize=sz[1]*print_zoom & print_ysize=sz[2]*print_zoom 

display_thresholdu = 350
display_thresholdd = - display_thresholdu       


;Reading header variables

;HMI uses structures for header values
if instr eq 4 then begin
	date = hdr_in.DATE_OBS

	;Define center and radius
	hfx = hdr_in.CRPIX1 ;  Location of the center in x pixels 
	hfy = hdr_in.CRPIX2 ;    Location of the center in y pixels
	di = hdr_in.RSUN_OBS/hdr_in.CDELT1;

	;Load Solar Coordinates
	P0 = 0.0
	RD = hdr_in.DSUN_OBS/hdr_in.RSUN_REF
	B0 = hdr_in.CRLT_OBS
	L0 = hdr_in.CRLN_OBS

	;Observer Coordinates
	X_scl = hdr_in.CDELT1/60.0
	Y_scl = hdr_in.CDELT2/60.0
	
endif else begin


	date = fxpar(hdr_in, 'DATE_OBS')

	;KPVT-512
	if instr eq 1 then begin
	
		;Define center and radius
		hfx = fxpar(hdr_in, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
		hfy = fxpar(hdr_in, 'CRPIX2A');+1.0;    Location of the center in y pixels
		di = fxpar(hdr_in,'EPH_R0');

		;Load Solar Coordinates
		P0 = 0.0
		RD = !values.f_nan
		B0 = fxpar(hdr_in, 'EPH_B0')
		L0 = fxpar(hdr_in, 'EPH_L0')

		;Observer Coordinates
		X_scl = fxpar(hdr_in, 'CDELT1')*fxpar(hdr_in, 'CRR_SCLX')/60.0
		Y_scl = fxpar(hdr_in, 'CDELT2')*fxpar(hdr_in, 'CRR_SCLY')/60.0

	endif

    ;KPVT-SPMG
    if instr eq 2 then begin
    
        ;Define center and radius
        hfx = fxpar(hdr_in, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
        hfy = fxpar(hdr_in, 'CRPIX2A');+1.0;    Location of the center in y pixels
        di = fxpar(hdr_in,'EPH_R0')/fxpar(hdr_in,'SCALE');

        ;Load Solar Coordinates
        P0 = 0.0
        RD = !values.f_nan
        B0 = fxpar(hdr_in, 'EPH_B0')
        L0 = fxpar(hdr_in, 'EPH_L0')

        ;Observer Coordinates
        X_scl = fxpar(hdr_in, 'CDELT1')*fxpar(hdr_in, 'CRR_SCLX')/60.0
        Y_scl = fxpar(hdr_in, 'CDELT2')*fxpar(hdr_in, 'CRR_SCLY')/60.0

    endif


	;MDI
	if instr eq 3 then begin
	
		;Define center and radius
		hfx = fxpar(hdr_in, 'X0');  Location of the center in x pixels 
		hfy = fxpar(hdr_in, 'Y0');  Location of the center in y pixels
		di = fxpar(hdr_in,'R_SUN');

		;Load Solar Coordinates
		P0 = 0.0
		RD = fxpar(hdr_in, 'OBS_DIST')/0.0046491
		B0 = fxpar(hdr_in, 'B0')
		L0 = fxpar(hdr_in, 'L0')

		;Observer Coordinates
		X_scl = fxpar(hdr_in, 'XSCALE')/60.0
		Y_scl = fxpar(hdr_in, 'YSCALE')/60.0	
	
	endif

endelse



print, 'Calculating heliospheric coordinates...'
;Calculation of heliospheric coordinates---------------------------------------------------------------------------


Xg = transpose(double(floor(findgen(sz[1],sz[1])/sz[1]))) - hfx 
Yg = double(floor(findgen(sz[1],sz[1])/sz[1])) - hfy 
R = sqrt( Xg^2 + Yg^2 )

Xg = Xg*X_scl
Yg = Yg*Y_scl

helio = arcmin2hel(Xg[findgen(long(sz[1]*sz[1]))], Yg[findgen(long(sz[1]*sz[1]))], date = date, p = P0, b0 = B0, l0 = L0, sphere = 1)

Lath = reform(helio[0,*],sz[1],sz[1])
Lonh = reform(helio[1,*],sz[1],sz[1])


nps_pl = 100.0
Lon_pl = findgen(nps_pl+1)/nps_pl*360.0 - 180.0

;Plotting Circles
LatN_pl = fltarr(nps_pl+1) + 70.0
LatE_pl = fltarr(nps_pl+1)
LatS_pl = fltarr(nps_pl+1) - 70.0

N_pl = hel2arcmin(LatN_pl, Lon_pl, vsblN, p = P0, b0 = B0, r0 = RD, l0 = L0 , rsun = di)*60.0
N_pl_ind = where(vsblN eq 0)
N_pl[0,N_pl_ind] = !values.f_nan
N_pl[1,N_pl_ind] = !values.f_nan

E_pl = hel2arcmin(LatE_pl, Lon_pl, vsblE, p = P0, b0 = B0, r0 = RD, l0 = L0 , rsun = di)*60.0
E_pl_ind = where(vsblE eq 0)
E_pl[0,E_pl_ind] = !values.f_nan
E_pl[1,E_pl_ind] = !values.f_nan

S_pl = hel2arcmin(LatS_pl, Lon_pl, vsblS, p = P0, b0 = B0, r0 = RD, l0 = L0 , rsun = di)*60.0
S_pl_ind = where(vsblS eq 0)
S_pl[0,S_pl_ind] = !values.f_nan
S_pl[1,S_pl_ind] = !values.f_nan

;Circle
nps_cir = 36000.0
Th_cir = findgen(nps_pl+1)/nps_pl*360.0 - 180.0

X_cir = (di)*cos(Th_cir*!dtor)+hfx
Y_cir = (di)*sin(Th_cir*!dtor)+hfy

;
;Raw image Display------------------------------------------------
;

;;if keyword_set(display) then begin
;
;    set_plot,'X'
;    loadct,0,/silent
;    window,3,xsize=display_xsize,ysize=display_ysize,retain=2 
;    plot,[1,1],/nodata,xstyle=5,ystyle=5
;    tv,bytscl(congrid(Lonh,display_xsize,display_ysize))
;    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
;    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
;    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
;;endif
;
;stop

if keyword_set(display) then begin

    set_plot,'X'
    loadct,0,/silent
    window,0,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color=255,/device,thick=3
endif

if keyword_set(prt) then begin
    set_plot,'PS'
    device, filename='Mgntgram_raw.eps'
    device,/color,bits_per_pixel=8,/portr,/inches,xsize=6.0,ysize=6.0,$
    xoff=0.5,yoff=0.5
    !p.position=[0.0,0.0,1.0,1.0]  
    !x.window=[0.0,1.0]
    !y.window=[0.0,1.0]
    px = !x.window * !d.x_vsize ;Get size of window in device units
    py = !y.window * !d.y_vsize
    
    plot,[1,1],/nodata,xstyle=5,ystyle=5    
    tv,bytscl(congrid(im,print_xsize,print_ysize),min=display_thresholdd,max=display_thresholdu)
    device, /close
    set_plot,'X'

endif
imgs0=im



print, 'Correcting LOS field...'
;Correct line of sight magnetic field assuming field is radial

;Observers' unit vector
Xobs = cos(B0*!dtor)*cos(L0*!dtor)
Yobs = cos(B0*!dtor)*sin(L0*!dtor)
Zobs = sin(B0*!dtor)

M_corr = cos(Lath*!dtor)*cos(Lonh*!dtor)*Xobs + cos(Lath*!dtor)*sin(Lonh*!dtor)*Yobs + sin(Lath*!dtor)*Zobs
im_corr = im/M_corr

im_corr_ind = where( (R gt di*sin(seg_const.deg_lim*!dtor)) and finite(im_corr), count )

if count NE 0 then im_corr[im_corr_ind] = !values.f_nan 
;
;Corrected image Display------------------------------------------------
;
if keyword_set(display) then begin

    set_plot,'X'
    loadct,0,/silent
    window,1,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(im,display_xsize,display_ysize),min=display_thresholdd,max=display_thresholdu)
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color='FFFFFF'x,/device,thick=3
    tvcircle,di*display_zoom*sin(seg_const.deg_lim*!dtor), hfx*display_zoom,hfy*display_zoom,color=255,/device,thick=3
endif

;
;Calculate the Element of Area--------------------------------------------------------
print, 'Calculating element of Area...'

;Converting to Cartesian
Xar = cos(Lath*!dtor)*cos(Lonh*!dtor)
Yar = cos(Lath*!dtor)*sin(Lonh*!dtor)
Zar = sin(Lath*!dtor)

;Finding the first diagonal vector
Parx = dblarr(sz[1],sz[2])
Parx[1:sz[1]-2,1:sz[2]-2] = Xar[0:sz[1]-3,0:sz[2]-3] - Xar[2:sz[1]-1,2:sz[2]-1]
Pary = dblarr(sz[1],sz[2])
Pary[1:sz[1]-2,1:sz[2]-2] = Yar[0:sz[1]-3,0:sz[2]-3] - Yar[2:sz[1]-1,2:sz[2]-1]
Parz = dblarr(sz[1],sz[2])
Parz[1:sz[1]-2,1:sz[2]-2] = Zar[0:sz[1]-3,0:sz[2]-3] - Zar[2:sz[1]-1,2:sz[2]-1]

P = [reform(Parx,1,sz[1]*sz[2]), reform(Pary,1,sz[1]*sz[2]), reform(Parz,1,sz[1]*sz[2])]

;Finding the second diagonal vector
Qarx = dblarr(sz[1],sz[2])
Qarx[1:sz[1]-2,1:sz[2]-2] = Xar[2:sz[1]-1,0:sz[2]-3] - Xar[0:sz[1]-3,2:sz[2]-1]
Qary = dblarr(sz[1],sz[2])
Qary[1:sz[1]-2,1:sz[2]-2] = Yar[2:sz[1]-1,0:sz[2]-3] - Yar[0:sz[1]-3,2:sz[2]-1]
Qarz = dblarr(sz[1],sz[2])
Qarz[1:sz[1]-2,1:sz[2]-2] = Zar[2:sz[1]-1,0:sz[2]-3] - Zar[0:sz[1]-3,2:sz[2]-1]

Q = [reform(Qarx,1,sz[1]*sz[2]), reform(Qary,1,sz[1]*sz[2]), reform(Qarz,1,sz[1]*sz[2])]

Atmp = crosspn(P,Q)

mgnt_ar = reform(sqrt(Atmp[0,*]^2.0 + Atmp[1,*]^2.0 + Atmp[2,*]^2.0)*6.955e10*6.955e10/8.0, sz[1], sz[2])
ind_npf = where(R gt di)
mgnt_ar(ind_npf) = !values.f_nan

if keyword_set(display) then begin
    loadct,0,/silent
    if not keyword_set(ps) then window,2,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt_ar,display_xsize,display_ysize),min=max(mgnt_ar,/nan)/100,max=max(mgnt_ar,/nan))
    loadct, 13
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color=200*3/6,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color=200*2/6,/device,thick=3
    tvcircle,di*display_zoom*sin(70.0*!dtor), hfx*display_zoom,hfy*display_zoom,color=255,/device,thick=3
    tvcircle,di*display_zoom*sin(50.0*!dtor), hfx*display_zoom,hfy*display_zoom,color=255,/device,thick=3
    tvcircle,di*display_zoom*sin(30.0*!dtor), hfx*display_zoom,hfy*display_zoom,color=255,/device,thick=3
    
    loadct, 0
endif

;print, junk

;Calculating flux---------------------------------------------
print, 'Calculating magnetic flux...'

mgnt_flux_corr = mgnt_ar*im_corr
mgnt_flux_raw = mgnt_ar*imgs0

if keyword_set(display) then begin
    loadct,0,/silent
    if not keyword_set(ps) then window,3,xsize=display_xsize,ysize=display_ysize,retain=2 
    plot,[1,1],/nodata,xstyle=5,ystyle=5
    tv,bytscl(congrid(mgnt_flux_raw, display_xsize,display_ysize),min=min(mgnt_flux_raw/10.0,/nan),max=max(mgnt_flux_raw/10.0,/nan))
    loadct, 13
    plots, (N_pl[0,*] + hfx)*display_zoom, (N_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    plots, (E_pl[0,*] + hfx)*display_zoom, (E_pl[1,*] + hfy)*display_zoom ,color=200*3/6,/device,thick=3
    plots, (S_pl[0,*] + hfx)*display_zoom, (S_pl[1,*] + hfy)*display_zoom ,color=200,/device,thick=3
    tvcircle,di*display_zoom, hfx*display_zoom,hfy*display_zoom,color=200*2/6,/device,thick=3
    loadct, 0
endif

CRD_out = {im_raw: imgs0, im_crr: im, hdr:hdr_in, mgnt_ar:mgnt_ar, mgnt_flux_raw: mgnt_flux_raw, mgnt_flux_corr: mgnt_flux_corr, Xar:Xar, Yar:Yar, Zar:Zar, Lath:Lath, Lonh:Lonh}
return
end
