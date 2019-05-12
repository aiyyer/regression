;--------------------------------------------------------------------------------------
; This program will read specified data files and create an output file with
; timeseries of data for selected locations for a fixed level.
;
; Aiyyer 01/2018
;--------------------------------------------------------------------------------------
; Notes:
; as configured currently, reads filetered V wind from yearly data to create
; reference timeseries to be used for regressions (Jan 2018)
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
  levRef = 650
;
; set the starting year
  ystart = 1998
; Number of years to read
  nyears = 13; 36
;
; for each year define the date range to read
; Define the JAS period with some padding on either side for lag regressions
; make sure that the lagHours are not more than the padding hours
;
  sTime = 062400
  eTime = 100700
;
; set the reference time series location 
; latRef defines the reference lat
  nRef = 6 ; this is the number of locations for the timeseries
  lonRef = new(nRef,"float")
  latRef = 10. ;18. ; only one lat allowed
  lonRef = (/-40.,-20.,0.,20.,40.,60./)

; 
; Define the output file name and write some data
;  fileName = "reference_timeseries_" + latRef + "N.nc"

  fileName = "reference_timeseries_daily_" + levRef +"hPa_" + latRef + "N.nc"


  system( "rm " + fileName )
  outFile = addfile( fileName, "c" )

  print (fileName)

;  dirFiltData = "/typhoon/projects/data/filtered/"

; for daily data
  dirFiltData = "/typhoon/projects/data/filtered/2to12/
;
;----------------------------------------------------------------------------------
  firstTime = 1
;
  do iy = ystart,ystart+nyears-1
    sDate = iy*1000000 + sTime
    eDate = iy*1000000 + eTime
    print ( sDate + "  " + eDate )

; set the input file
; here it is defined as the filtered meridional wind data for each year
    
    inFile = dirFiltData + "V2T10Day" + iy + ".nc"
    inFile = dirFiltData + "V10DayHP" + iy + ".nc"
    inFile = dirFiltData + "2-12DayBP_daily"+ iy + ".nc"


    in = addfile (inFile,"r")

;    
; read the data for the reference point (or the nearest location)
; refTS is the reference time series against which we will regress
; unfiltered fields
    

;;;; EDITS IN PROGRESS need to fix dates for the subsetting to work

    dat   = lonFlip(in->VBP({sDate:eDate},{levRef},{latRef},:))
    refTS = dat(:,{lonRef})
    delete(in)
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
