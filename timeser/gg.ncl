;--------------------------------------------------------------------------------------
; This program will read specified data files and create an output file with
; timeseries of data for selected locations for a fixed level.
;
; Aiyyer 11/2018
;--------------------------------------------------------------------------------------
; Notes:
; as configured currently, reads gpcp filtered fields
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
  
; set the starting year
  ystart = 1998
; Number of years to read
  nyears = 17
  yend   = ystart + nyears - 1

  fileName = "gpcp_daily.nc"
  system( "rm " + fileName )
  outFile = addfile( fileName, "c" )
  print (fileName)

  dir  = "/tamrind/gpcp/"

; read one of the files to get info on units
  in   = addfile(fil, "r")


  timeUnits =  "days since 1970-01-01 00:00:00 0:00"
  print(timeUnits)
  sc=0
  mn=0
  hh=0
;
  varNam = "rainBP"
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


; read the data for the reference point (or the nearest location)
; refTS is the reference time series against which we will regress
; unfiltered fields

    fNam = "gpcp_2to10day_" + iy + ".nc"
    fil  = dir+fNam
    in   = addfile(fil, "r")
    data=lonFlip(in->$varNam$({timeStart:timeEnd},{latRef},:))
  
   ; printVarSummary(data)
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
