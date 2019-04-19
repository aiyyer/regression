;--------------------------------------------------------------------------------------
; This program will read specified data files and create an output file with
; timeseries of data for selected locations for a fixed level.
;
; Aiyyer 11/2018
;--------------------------------------------------------------------------------------
; Notes:
; as configured currently, reads trmm filtered fields
; reference timeseries to be used for regressions (Nov 2018)
;
;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
begin
;
;
; Set the level for the reference time series
  levRef = 0 ; not needed for trmm but keep for consistency with era data
;
  
; as of now we go to 2014 as filtered data goes up to sep 1 in 2015
; since we used nwght=12 in the bp filter program. Aiyyer 11/2018
; set the starting year
  ystart = 1998
; Number of years to read
  nyears = 15
  yend   = ystart + nyears - 1
  
; for each year define the date range to read
; Define the JAS period with some padding on either side for lag regressions
; make sure that the lagHours are not more than the padding hours
; this is set within the year loop
;  sTime = 062400
;  eTime = 100700
;
; set the reference time series location 
; latRef defines the reference lat
  nRef = 9 ; this is the number of locations for the timeseries
  lonRef = new(nRef,"float")
  latRef = 10. ;18. ; only one lat allowed
  lonRef = (/-40.,-30.,-20.,-10.,0.,10.,20.,30.,40./)

  fileName = "reference_timeseries_trmm_daily" + latRef + "N.nc"
  system( "rm " + fileName )
  outFile = addfile( fileName, "c" )

  print (fileName)

  dir  = "/share01/trmm/carl_3b42/"
  fNam = "trmm3b42_daily_2to12day_filtered.nc"
  varNam  = "rain2to12"
  fil  = dir+fNam
  in   = addfile(fil, "r")
  ;time = in->time
  timeUnits =  in->time@units
  print(timeUnits)

  sc=0
  mn=0
  hh=0


;
;----------------------------------------------------------------------------------
  firstTime = 1

  do iy = ystart,yend
    mm=06
    dd=24
    sTime = mm*10000 + dd*100
    timeStart   = cd_inv_calendar(iy,mm,dd,hh,mn,sc,timeUnits, 0)
    mm=10
    dd=07  
    timeEnd   = cd_inv_calendar(iy,mm,dd,hh,mn,sc,timeUnits, 0)
    eTime = mm*10000 + dd*100
    
    ndays = timeEnd-timeStart+1
    print(timeStart + " " + timeEnd + " " + ndays)
    print(sTime + " " + eTime)
    print (cd_calendar(timeStart,-2) + "  " + cd_calendar(timeEnd,-2))
                                ; read the data for the given days

; read the data for the reference point (or the nearest location)
; refTS is the reference time series against which we will regress
; unfiltered fields


    data=in->$varNam$({timeStart:timeEnd},{latRef},:)
  
    printVarSummary(data)
    refTS = data(:,{lonRef})
    if ( iy .eq. ystart )
      TS = refTS
    else
      refTemp = array_append_record (TS,refTS,0)
      delete(TS)
      TS = refTemp
      delete(refTemp)
    end if    
  end do

; write the timeseries to output file
  outFile->TS = TS
  outFile->latRef=latRef
  outFile->lonRef=lonRef
  outFile->levRef=levRef
  outFile->ystart=ystart
  outFile->nyears=nyears
  outFile->sTime=sTime
  outFile->eTime=eTime

  end 
