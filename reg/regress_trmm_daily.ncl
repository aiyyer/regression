;--------------------------------------------------------------------------------------
; This program will read a reference time series from specified input file
; and regress trmm data on the reference time series
;
;

; Aiyyer  09/2015
; updated 01/2015 now for multiple levels
; updated 01/2015 used the append function to speed up data read
; updated 09/2016 multiple reference time series can be used in one shot
; updated 01/2018 to remove the reference timeseries creation and instead
;                 now code will read existing time series
;
;
;--------------------------------------------------------------------------------------
; Notes
; The program reads band-passed filtered data created earlier. These data have time dimension
; that is YYYYMMDDHH. 
;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../..//ncl/aiyyer_string_routines.ncl"

;
begin

; 
; Define the timeSeries name and read parameters
  fileName = "ref_ts.nc"
  inFile = addfile( fileName, "r" )
  latRef=inFile->latRef
  lonRef=inFile->lonRef
  levRef=inFile->levRef
  ystart=inFile->ystart
  nyears=inFile->nyears
  sTime=inFile->sTime
  eTime=inFile->eTime
  TS=inFile->TS



; Since TRMM data begin in 1998, we need to remove all date prior to it
; if they exist in the data
  if ( ystart .lt. 1998 ) then
    print ( "ystart = " + ystart)
    print ("ystart prior to TRMM data. Adjusting the years")
    nyearsCut = 1998-ystart     
    
    dims = dimsizes(TS)
    ntimes = dims(0)
; number of days in each year in the reference timeseries is  
    ndaysperyear = ntimes/nyears
    N = nyearsCut*ndaysperyear
    TS2 = TS(N:,:)
    TS  :=TS2
    delete(TS2)
    print (nyears + " " + nyearsCut + " " + ndaysperyear + " " + N)
  end if
  printVarSummary(TS)
  
; now reset ystart and nyears

  nyears = ystart + nyears - 1998
  ystart = 1998
  print ("Number of years; ystart =  " + nyears + " " + ystart)
  
;
  dh       = 24  ; data is every 6 hours
  lagHours = 72 ; lag +/- lagHours
  nlag     = lagHours/dh
  lagArray = ispan(-lagHours,lagHours,dh)
  
  print ( sTime + "  " + eTime )

;
; Now read the data to be regressed against the reference time series

  print (nyears)

  varName = "raindaily" ;"rain"

  ; the trmm data file
  
  trmmFile = "/share01/trmm/carl_3b42/trmm3b42_daily.nc"
  in   := addfile (trmmFile, "r")
                                ; read the units of time in this file

  trmmUnits = in->time@units

  trmmDates = in->time

  
  print (trmmUnits)

; Define the output file name

; fileName = "trmm" +  "_Vfilt_levRef"+ levRef+ "latRef"+latRef+ ".nc"
  fileName = "trmm" +  "_trmm_" + "latRef"+latRef+ ".nc"

  system( "rm " + fileName )        


;---------------------------------------------------------------------------------------------
    firstTime = True
    do iy = ystart,ystart+nyears-1
      sDate = (iy*1000000 + sTime)/100
      eDate = (iy*1000000 + eTime)/100

      ; convert the starting and ending dates to
      ; units of the trmm data

      ; call custom aiyyer routines
      dateBegin = ymddate(sDate,trmmUnits)
      dateEnd   = ymddate(eDate,trmmUnits)
;      print (cd_calendar(dateBegin,-2) + " " + cd_calendar(dateEnd,-2))
;      print (sDate + " " + eDate)

; now read the data
;      print (dateBegin + " " + dateEnd)
      dat := lonFlip(in->$varName$({dateBegin:dateEnd},{-45.:45.},:))
      ;printVarSummary(dat)

      datIn := dat(:,:,{-90.:90.})
      if ( firstTime ) then
        print ("First Time")
        datAppend := datIn
        firstTime = False
      else        
        datTemp   := array_append_record(datAppend, datIn,0)
        datAppend := datTemp
        delete(datTemp)
      end if
    end do
;---------------------------------------------------------------------------------------------
       
    delete(dat)
    delete(datIn)
    datAppend!0    = "time"
    datAppend&time = TS&time
    DATA = datAppend(lat|:,lon|:,time|:)
    delete(datAppend)
    dimsTS = dimsizes(TS)
    nT = dimsTS(0)/nyears
    nL = dimsTS(1)
    delete(dimsTS)

    printVarSummary(DATA)
;
; Now that we have the reference timeseries and the data to regress on to it
; We proceed to lag-regressions
;
    ilag = 0
    do ih = -nlag,nlag,1
      
; define the new arrays to hold the lagged data
; Since we have seasonal data, we will lose ih points for each year
      ntnew  = (nT - abs(ih))  ; number of times per year
;
; ntAll defines the total number of temporal points over nyears
;
      ntAll = nyears*ntnew     ; number of times over all years
;
; create new arrays to hold the time series and data to be regressed
      dimData = dimsizes(DATA)
     
      ny =  dimData(0)
      nx =  dimData(1)
      A  =  new ( (/ntAll,nL/)      ,   typeof(TS) )
      B  =  new ( (/ny,nx,ntAll/), typeof(DATA) )
      A = -999.
      B = -999.
;
      do iy = ystart,ystart+nyears-1
;      
; the bounds of the data for this season
        i1 = (iy-ystart)*nT
        i2 = i1 + nT - 1
;           
        i1new = (iy-ystart)*ntnew
        i2new = i1new + ntnew -1 
;      
; the new bounds accounting for the lag
;      
        if ( ih .lt. 0 ) then
; For negative lag, the reference time series is Ahead in time relative to the data to be regressed
; So we shift the reference TS forward (-ih is positve!)
          j1 = i1+abs(ih)
          j2 = j1 + ntnew-1
          k1 = i1
          k2 = k1 + ntnew -1
        end if
        if ( ih .ge. 0 ) then
; For positive lag, the reference timeseries is Behind in time relative to the data to be regressed
; So we shift the data forward
          k1 = i1+abs(ih)
          k2 = k1 + ntnew-1
          j1 = i1
          j2 = j1 + ntnew -1           
        end if      
        A(i1new:i2new,:)   =     TS(j1:j2,:)
        B(:,:,i1new:i2new) =   DATA(:,:,k1:k2)         
        print ("  " )     
        print ( TS&time(j1) + " " +  TS&time(j2) )
        print ( DATA&time(k1) + " " +  DATA&time(k2) )
      end do
      
      
; at this point we have the correctly phased timeseries for all years and the data to be regressed for all years
      
      tval       =  new ( (/ny,nx/),  typeof ( A ) )
      nptsxy     =  new ( (/ny,nx/),       integer )
      
      do iL = 0,nL-1
        rc = regcoef(A(:,iL),B,tval,nptsxy)
;
; compute the std dev of A and multiply the regression coeff by it
; to get anomalies
        stdA = dim_stddev_Wrap( A(:,iL) )
        rc = rc*stdA
        print ( "Lag = " + ih + " stdA =" + stdA )
        df    = nptsxy-2   ; degrees of freedom
        bb     = tval            ; b must be same size as tval (and df)
        bb     = 0.5
        prob  = betainc(df/(df+tval^2),df/2.0,bb)
;
; Create the arrays to hold the regressed data      
        if ( ilag .eq. 0 .and. iL .eq. 0)then
          rcAll   = new ( (/ny,nx,2*nlag+1,nL/), "float" )
          prAll   = new ( (/ny,nx,2*nlag+1,nL/), "float")
          
          rcAll!0 = "lat"
          rcAll!1 = "lon"
          rcAll!2 = "lag"
          rcAll!3 = "lonRef"
          
          rcAll&lat   = DATA&lat
          rcAll&lon   = DATA&lon
          rcAll&lag  = lagArray
          copy_VarMeta(rcAll,prAll)
        end if
;        
        prAll(:,:,ilag,iL) = (/prob/)
        rcAll(:,:,ilag,iL) = (/rc/)
        delete(prob)
        delete(rc)
      end do
      ilag = ilag + 1      
;
      delete (A)
      delete (B)
      delete (tval)
      delete (nptsxy)
      delete (df)
      delete (bb)   
    end do
    
    printVarSummary( rcAll )
    printVarSummary( prAll )
    
    outFile = addfile( fileName, "c" )
    outFile->rcAll=rcAll
    outFile->prAll=prAll
    

    delete(rcAll)
    delete(prAll)
    delete (DATA)


end 
