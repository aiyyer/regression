;--------------------------------------------------------------------------------------
; This program will read a reference time series from specified input file
; and regress a given field on to the reference time series
;
; This version is for processing 1xdaily fields (i.e. daily averages) 
; currently configured for daily erai files on our system
;
; This version is for 2D fields such as surface fluxes
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
;NCSU Tropical Dynamics
;A. Aiyyer
;
;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------
;
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
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

;
  dh       = 24  ; data is every 6 hours
  lagHours = 72 ; lag +/- lagHours
  nlag     = lagHours/dh
  lagArray = ispan(-lagHours,lagHours,dh)
  
  print ( sTime + "  " + eTime )
  printVarSummary(TS)
  dimsTS = dimsizes(TS)
  nTimeRef= dimsTS(0) 

  nTimeRefPerYear = nTimeRef/nyears
  print (nTimeRef + " " + nyears + " " + nTimeRefPerYear)

;
; Now read the data to be regressed against the reference time series
; this is a quick way to create an array containing the time stamp of each 6 hourly data 
; we will use the array to create another array of filenames to do a mass read of data

  print (nyears)


; list of variables that we need
 derived = False
 variables = (/"SSR","TSR","TTR","STR","SSHF","SLHF","TP"/)


;  derived = True
; variables = (/"ColAvgQ","ColAvgMse","ColAvgDse","ColAvgzonalAdvMse","ColAvgmeridAdvMse","ColAvgvertAdvMse","ColAvgzonalAdvDse","ColAvgmeridAdvDse","ColAvgvertAdvDse","ColAvgzonalAdvQ","ColAvgmeridAdvQ","ColAvgvertAdvQ","TendColAvgMse","TendColAvgDse","TendColAvgQ"/)

print(variables)
;
;
; the names of variables in the datafiles
  varName = variables
;
; loop over the variables
  nvar = dimsizes(variables)
  do ivar = 0,nvar-1

; Define the output file name
    fileName = "./" + "regressed_" + variables(ivar)+"_latRef"+latRef+ ".nc"

    system( "rm " + fileName )        
    print (fileName)

    outFile = addfile( fileName, "c" )
    outFile->latRef=(/latRef/)
    outFile->lonRef=(/lonRef/)
    outFile->levRef=(/levRef/)
    outFile->ystart=(/ystart/)
    outFile->nyears=(/nyears/)
    outFile->sTime=(/sTime/)
    outFile->eTime=(/eTime/)


;---------------------------------------------------------------------------------------------
    firstTime = True
    do iy = ystart,ystart+nyears-1
      sDate = (iy*1000000 + sTime)/100
      eDate = (iy*1000000 + eTime)/100
      allTimes := yyyymmdd_time(iy,iy, "integer")
      TIME     := allTimes({sDate:eDate})


; now define all files 
; are the fields derived variables or original variables?      	
      if ( derived ) then
        fils := "/tamrind/eradaily/2D_derived/" + TIME + ".nc" 
      else
        fils := "/tamrind/eradaily/2D/" + TIME + ".nc" 
      end if



; TP is in 3D files from 1998-2010 
      if ( ivar .eq. 6) 
        if ( iy .lt. 2011 ) then
          fils := "/tamrind/eradaily/3D/" + TIME + ".nc" 
        else
          fils := "/tamrind/eradaily/2D/" + TIME + ".nc" 
        end if
      end if
      


      print ( "year = " + iy)
      print(fils)
      in   := addfiles (fils, "r")
      ListSetType (in, "join")             ; concatenate or "merge" (default)
      delete(allTimes)
      delete(TIME) 


; now read the data subset and lonflip to center on 0
      dat := lonFlip(in[:]->$varName(ivar)$(:,{-30.:36.},:))      
      printVarSummary(dat)
      delete(in)
      datIn := dat(:,:,{-90.:90.})
      
      dimsDatIn = dimsizes(datIn)
      nTimesDat = dimsDatIn(0)
      print (nTimesDat + " " + nTimeRefPerYear)
      if ( nTimeRefPerYear .ne. nTimesDat) then
        print ("Missing data! Stop!")
        exit
      end if
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
    exit
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
; This yields anomalies of the same units as the data being regresssed

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
          rcAll@lonRef = lonRef
          rcAll@long_name=DATA@long_name
          rcAll@units=DATA@units
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
    
    outFile->rcAll=rcAll
    outFile->prAll=prAll
    

    delete(rcAll)
    delete(prAll)
    delete (DATA)
; end loop over all variables
    end do

end 
