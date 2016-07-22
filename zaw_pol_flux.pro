pro zaw_pol_flux, start_date, end_date, instr
;+
; NAME:
;    ZAW_POL_FLUX
;
; PURPOSE:
;   Identifies the pixels belonging to the norhtern and southern poles and calculates
;   data on them.
;
; CATEGORY:
;    DATA ANALYSIS
;
; CALLING SEQUENCE:
;    zaw_pol_flux, start_date, end_date, instr
;
; INPUTS:
;    start_date                 - starting date to parse from
;    end_date                   - ending date to stop
;    instr                      - instrument being analyzed
;
; OUTPUTS:
;    Writes out a CSV of the struct PF_data containing norther and southern hemisphere data
;
; MODIFICATION HISTORY:
;   2011/09/26: Andres Munoz-Jaramillo:  Adapted from Jie Zhang's AR detection code
;   2016/02/06: Zach Werginz: Adapted to work with amj_coord.pro and amj_file_read.pro
;   2016/03/12: Zach Wergin: Added new fields to PF_data for maximum pixel field and flux
;-


;Define the usage
if n_params() lt 2 then begin
    print,'USAGE:  zaw_pol_flux, start_date, end_date, instr'
    print, "instr: 1 = KPVT 512"
    print, "instr: 2 = KPVT SPMG"
    print, "instr: 3 = MDI"
    print, "instr: 4 = HMI"
    return
endif

; Runtime parameters
deg_lim = 65.0
inv_px_tol = .85 
stat = 1

;Set instrument parameters
;First reference day for keeping track time easily
if instr eq 1 then DayOff = julday(1,1,1970);   KPVT 512
if instr eq 2 then DayOff = julday(1,1,1990);   KPVT SPMG
if instr eq 3 then DayOff = julday(1,1,1993);   MDI
if instr eq 4 then DayOff = julday(1,1,2009);   HMI

;Setting dates
mdi_i  = julday(strmid(start_date,5,2),strmid(start_date,8,2),strmid(start_date,0,4))-DayOff
mdi_f  = julday(strmid(end_date,5,2),strmid(end_date,8,2),strmid(end_date,0,4))-DayOff

;Data structure for poles
PF_data = {pf_day, mdi_i: 0L, date: '', intf_n: !values.f_nan, intfc_n: !values.f_nan, unsflux_n: !values.f_nan, unsfluxc_n: !values.f_nan, sflux_n: !values.f_nan, sfluxc_n: !values.f_nan, posfluxc_n: !values.f_nan, $
            negfluxc_n: !values.f_nan, vnpc_pxn: !values.f_nan, visarea_n: !values.f_nan, max_pxflux_n: !values.f_nan, max_pxf_n: !values.f_nan, max_pxfc_n: !values.f_nan, n_swt: !values.f_nan, $
            intf_s: !values.f_nan, intfc_s: !values.f_nan, unsflux_s: !values.f_nan, unsfluxc_s: !values.f_nan, sflux_s: !values.f_nan, sfluxc_s: !values.f_nan, posfluxc_s: !values.f_nan, negfluxc_s: !values.f_nan, $
            vnpc_pxs: !values.f_nan, visarea_s: !values.f_nan, max_pxflux_s: !values.f_nan, max_pxf_s: !values.f_nan, max_pxfc_s: !values.f_nan, s_swt: !values.f_nan}

head = ['mdi_i','date','intf_n', 'intfc_n', 'unsflux_n', 'unsfluxc_n', 'sflux_n', 'sfluxc_n', 'posfluxc_n', 'negfluxc_n', 'vnpc_pxn', 'visarea_n', 'max_pxflux_n', 'max_pxf_n', 'max_pxfc_n', 'n_swt', $
                       'intf_s', 'intfc_s', 'unsflux_s', 'unsfluxc_s', 'sflux_s', 'sfluxc_s', 'posfluxc_s', 'negfluxc_s', 'vnpc_pxs', 'visarea_s', 'max_pxflux_s', 'max_pxf_s', 'max_pxfc_s', 's_swt']

filename = string('PF_data' + start_date + '_' + end_date + '.csv')

REPEAT begin
    
    ;Reading files 
    repeat begin
            caldat, mdi_i + DayOff, Month, Day, Year
            date = strtrim(string(Year), 2)+'-'+strtrim(string(Month,format='(I02)'), 2)+'-'+strtrim(string(Day,format='(I02)'), 2)
            print, date
            mg = amj_file_read( date, hdr, instr )
        
            s = size(mg)    
            
            mdi_i = mdi_i + 1          
    endrep until (s[2] eq 8 or mdi_i gt mdi_f)

    if (date eq end_date) then stat = 0

    ;Reading header variables

    ;HMI uses structures for header values
    if instr eq 4 then begin
        date = hdr.DATE_OBS

        ;Define center and radius
        hfx = hdr.CRPIX1 ;  Location of the center in x pixels 
        hfy = hdr.CRPIX2 ;    Location of the center in y pixels
        di = hdr.RSUN_OBS/hdr.CDELT1;

        ;Load Solar Coordinates
        P0 = 0.0
        RD = hdr.DSUN_OBS/hdr.RSUN_REF
        B0 = hdr.CRLT_OBS
        L0 = hdr.CRLN_OBS

        ;Observer Coordinates
        X_scl = hdr.CDELT1/60.0
        Y_scl = hdr.CDELT2/60.0
        
    endif else begin


        date = fxpar(hdr, 'DATE_OBS')

        ;KPVT-512
        if instr eq 1 then begin
        
            ;Define center and radius
            hfx = fxpar(hdr, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
            hfy = fxpar(hdr, 'CRPIX2A');+1.0;    Location of the center in y pixels
            di = fxpar(hdr,'EPH_R0');

            ;Load Solar Coordinates
            P0 = 0.0
            RD = !values.f_nan
            B0 = fxpar(hdr, 'EPH_B0')
            L0 = fxpar(hdr, 'EPH_L0')

            ;Observer Coordinates
            X_scl = fxpar(hdr, 'CDELT1')*fxpar(hdr, 'CRR_SCLX')/60.0
            Y_scl = fxpar(hdr, 'CDELT2')*fxpar(hdr, 'CRR_SCLY')/60.0

        endif

        ;KPVT-SPMG
        if instr eq 2 then begin
        
            ;Define center and radius
            hfx = fxpar(hdr, 'CRPIX1A');35;'CRPIX1');  Location of the center in x pixels 
            hfy = fxpar(hdr, 'CRPIX2A');+1.0;    Location of the center in y pixels
            di = fxpar(hdr,'EPH_R0')/fxpar(hdr,'SCALE');

            ;Load Solar Coordinates
            P0 = 0.0
            RD = !values.f_nan
            B0 = fxpar(hdr, 'EPH_B0')
            L0 = fxpar(hdr, 'EPH_L0')

            ;Observer Coordinates
            X_scl = fxpar(hdr, 'CDELT1')*fxpar(hdr, 'CRR_SCLX')/60.0
            Y_scl = fxpar(hdr, 'CDELT2')*fxpar(hdr, 'CRR_SCLY')/60.0

        endif        

        ;MDI
        if instr eq 3 then begin
        
            ;Define center and radius
            hfx = fxpar(hdr, 'X0');  Location of the center in x pixels 
            hfy = fxpar(hdr, 'Y0');  Location of the center in y pixels
            di = fxpar(hdr,'R_SUN');

            ;Load Solar Coordinates
            P0 = 0.0
            RD = fxpar(hdr, 'OBS_DIST')/0.0046491
            B0 = fxpar(hdr, 'B0')
            L0 = fxpar(hdr, 'L0')

            ;Observer Coordinates
            X_scl = fxpar(hdr, 'XSCALE')/60.0
            Y_scl = fxpar(hdr, 'YSCALE')/60.0    
        
        endif

    endelse

    ;Calculate heliospheric data    
    amj_coord, mg.img, hdr, CRD, instr;, /disp
    sz = size(CRD.im_raw)
    Xg = transpose(double(floor(findgen(sz[1],sz[1])/sz[1]))) - hfx 
    Yg = double(floor(findgen(sz[1],sz[1])/sz[1])) - hfy 
    R = sqrt( Xg^2 + Yg^2 )

    Xg = Xg*X_scl
    Yg = Yg*Y_scl




    ;;North
    n_swt = 0.0
    pc_pxindn = where(R le di and CRD.Lath ge deg_lim, npc_pxn);   Polar cap pixels
    vpc_pxindn = where(R le di and finite(CRD.im_raw) and CRD.Lath ge deg_lim, vnpc_pxn);  Valid pixels inside the polar cap
    pc_posn = where(R le di and CRD.Lath ge deg_lim and CRD.im_crr gt 0.0, npc_posn); Pixels with positive corrected field inside the polar cap
    pc_negn = where(R le di and CRD.Lath ge deg_lim and CRD.im_crr lt 0.0, npc_negn); Pixels with negative corrected field inside the polar cap

    if ( npc_negn eq 0 or npc_posn eq 0 ) then begin    
        print, 'Northern hemisphere no polarity mixture' 
        n_swt = 3.0
    endif    
    if ( float(vnpc_pxn)/float(npc_pxn) le inv_px_tol or vpc_pxindn[0] eq - 1 ) then begin
        print, 'Northen hemisphere has more than ' + string((1.0-inv_px_tol)*100.0) + '% invalid pixels' 
        n_swt = 2.0
    endif
    if ( n_swt ne 0.0 ) then begin 

        intf_n = !values.f_nan 
        intfc_n = !values.f_nan
        unsflux_n = !values.f_nan 
        unsfluxc_n = !values.f_nan
        sflux_n = !values.f_nan 
        sfluxc_n = !values.f_nan
        posfluxc_n = !values.f_nan 
        negfluxc_n = !values.f_nan
        vnpc_pxn = !values.f_nan 
        visarea_n = !values.f_nan
        max_pxflux_n = !values.f_nan
        max_pxf_n = !values.f_nan
        max_pxfc_n = !values.f_nan                 
    endif else begin   
            
        intf_n = mean( CRD.im_raw(vpc_pxindn),/double, /nan)                    ;Uncorrected Mean polar cap field
        intfc_n = mean( CRD.im_crr(vpc_pxindn),/double, /nan)                   ;Corrected Mean polar cap field
        unsflux_n = total( abs(CRD.mgnt_flux_raw(vpc_pxindn)),/double, /nan)    ;Total polar unsigned flux - uncorrected
        unsfluxc_n = total( abs(CRD.mgnt_flux_corr(vpc_pxindn)), /double, /nan) ;Total polar unsigned flux - corrected
        sflux_n = total( CRD.mgnt_flux_raw(vpc_pxindn), /double, /nan)          ;Total polar signed flux   - uncorrected
        sfluxc_n = total( CRD.mgnt_flux_corr(vpc_pxindn), /double, /nan)        ;Total polar signed flux   - corrected
        posfluxc_n = total( CRD.mgnt_flux_corr(pc_posn),/double, /nan)           ;Total positive polar flux - corrected
        negfluxc_n = total( CRD.mgnt_flux_corr(pc_negn),/double, /nan)           ;Total negative polar flux - corrected
        visarea_n = total( CRD.mgnt_ar(vpc_pxindn), /double, /nan)              ;Total visible area
        max_pxflux_n = max(abs( CRD.mgnt_flux_corr(vpc_pxindn) ), /nan)         ;Maximum pixel flux - corrected
        max_pxf_n = max( CRD.im_raw(vpc_pxindn) ,/nan)                          ;Maximum uncorrected pixel field
        max_pxfc_n =  max( CRD.im_crr(vpc_pxindn) ,/nan)                        ;Maximum corrected pixel field
        if (unsfluxc_n eq 0.0) then begin
            unsfluxc_n = !values.f_nan
            sfluxc_n = !values.f_nan
            posfluxc_n = !values.f_nan
            negfluxc_n = !values.f_nan
            max_pxflux_n = !values.f_nan
            n_swt = 4                   ;for corrected field zero identifier
        endif
    endelse


    ;South
    s_swt = 0.0
    pc_pxinds = where(R le di and CRD.Lath le -deg_lim, npc_pxs);   Polar cap pixels
    vpc_pxinds = where(R le di and finite(CRD.im_raw) and CRD.Lath le -deg_lim, vnpc_pxs);  Valid pixels inside the polar cap
    pc_poss = where(R le di and CRD.Lath le -deg_lim and CRD.im_crr gt 0.0, npc_poss); Pixels with possitive corrected field inside the polar cap
    pc_negs = where(R le di and CRD.Lath le -deg_lim and CRD.im_crr lt 0.0, npc_negs); Pixels with negative corrected field inside the polar cap

    if ( npc_negs eq 0 or npc_poss eq 0 ) then begin
        print, 'Southern hemisphere no polarity mixture' 
        s_swt = 3.0
    endif
    if ( float(vnpc_pxs)/float(npc_pxs) le inv_px_tol or vpc_pxinds[0] eq -1 ) then begin
        print, 'Southern hemisphere has more than ' + string((1.0-inv_px_tol)*100.0) + '% invalid pixels' 
        s_swt = 2.0
    endif

        
    if ( s_swt ne 0.0 ) then begin 

        intf_s = !values.f_nan
        intfc_s = !values.f_nan
        unsflux_s = !values.f_nan
        unsfluxc_s = !values.f_nan
        sflux_s = !values.f_nan
        sfluxc_s = !values.f_nan
        posfluxc_s = !values.f_nan
        negfluxc_s = !values.f_nan
        vnpc_pxs = !values.f_nan
        visarea_s = !values.f_nan
        max_pxflux_s = !values.f_nan
        max_pxf_s = !values.f_nan  
        max_pxfc_s = !values.f_nan 
        
    endif else begin   
            
        intf_s = mean( CRD.im_raw(vpc_pxinds),/double, /nan)                    ;Uncorrected Mean polar cap field
        intfc_s = mean( CRD.im_crr(vpc_pxinds),/double, /nan)                   ;Corrected Mean polar cap field
        unsflux_s = total( abs(CRD.mgnt_flux_raw(vpc_pxinds)),/double, /nan)    ;Total polar unsigned flux - uncorrected
        unsfluxc_s = total( abs(CRD.mgnt_flux_corr(vpc_pxinds)),/double, /nan)  ;Total polar unsigned flux - corrected
        sflux_s = total( CRD.mgnt_flux_raw(vpc_pxinds),/double, /nan)           ;Total polar signed flux   - corrected
        sfluxc_s = total( CRD.mgnt_flux_corr(vpc_pxinds),/double, /nan)         ;Total polar signed flux   - corrected
        posfluxc_s = total( CRD.mgnt_flux_corr(pc_poss),/double, /nan)          ;Total positive polar flux - corrected
        negfluxc_s = total( CRD.mgnt_flux_corr(pc_negs),/double, /nan)          ;Total negative polar flux - corrected
        visarea_s = total( CRD.mgnt_ar(vpc_pxinds), /double, /nan)              ;Total visible area
        max_pxflux_s = max(abs( CRD.mgnt_flux_corr(vpc_pxinds) ), /nan)         ;Maximum pixel flux        - corrected
        max_pxf_s = max( CRD.im_raw(vpc_pxinds) ,/nan)                          ;Maximum uncorrected pixel field
        max_pxfc_s =  max( CRD.im_crr(vpc_pxinds) ,/nan)                        ;Maximum corrected pixel field 
        if (unsfluxc_s eq 0.0) then begin
            unsfluxc_s = !values.f_nan
            sfluxc_s =   !values.f_nan
            posfluxc_s = !values.f_nan
            negfluxc_s = !values.f_nan
            max_pxflux_s = !values.f_nan
            s_swt = 4                   ;for corrected field zero identifier
        endif
                
    endelse
    tmp_PF = {pf_day, mdi_i, date, intf_n, intfc_n, unsflux_n, unsfluxc_n, sflux_n, sfluxc_n, posfluxc_n, negfluxc_n, vnpc_pxn, visarea_n, max_pxflux_n, max_pxf_n, max_pxfc_n, n_swt, intf_s, intfc_s, unsflux_s, unsfluxc_s, sflux_s, sfluxc_s, posfluxc_s, negfluxc_s, vnpc_pxs, visarea_s, max_pxflux_s, max_pxf_s, max_pxfc_s, s_swt}
    PF_data = [PF_data, tmp_PF] 
    write_csv, filename, PF_data, HEADER = head
ENDREP UNTIL (mdi_i gt mdi_f)


filename = string('PF_data' + start_date + '_' + end_date + '.csv')
head = ['mdi_i','date','intf_n', 'intfc_n', 'unsflux_n', 'unsfluxc_n', 'sflux_n', 'sfluxc_n', 'posfluxc_n', 'negfluxc_n', 'vnpc_pxn', 'visarea_n', 'max_pxflux_n', 'max_pxf_n', 'max_pxfc_n', 'n_swt', $
                       'intf_s', 'intfc_s', 'unsflux_s', 'unsfluxc_s', 'sflux_s', 'sfluxc_s', 'posfluxc_s', 'negfluxc_s', 'vnpc_pxs', 'visarea_s', 'max_pxflux_s', 'max_pxf_s', 'max_pxfc_s', 's_swt']
write_csv, filename, PF_data, HEADER = head

END