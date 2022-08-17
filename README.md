# ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System User Manual
## Alan Chou 09/08/2022 alan_chou@ymail.com, ac6522@ic.ac.uk

## Table of Contents
- [Flow diagram of the slope stability modelling system](#Flow-diagram-of-the-slope-stability-modelling-system)
- [1. ERA5 Data extraction and data pre-processing](#ERA5-Data-extraction-and-pre-processing)
   - [1.1 ERA5 data downloading by CDS API](#ERA5-data-downloading-by-CDS-API)
   - [1.2 Grib file data extraction and stored as variables](#Grib-file-data-extraction-and-stored-as-variables)
   - [1.3 Merging all data into 3 list](#Merging-all-data-into-3-lists)
   - [1.4 Data pre-processing](#Data-pre-processing)
   - [1.5 Saving data as csv file](#Saving-data-as-csv-file)
- [2. Hydrological Processing](#Hydrological-Processing)
   - [2.1 DAYMOD soil moisture model](#DAYMOD-Soil-moisture-model)
   - [2.2 Reading ERA5 data in R](#Reading-ERA5-data-in-R)
   - [2.3 Processing ERA5 soil moisture data](#Processing-ERA5-soil-moisture-data)
   - [2.4 Saving ERA5 soil moisture data as csv file](#Saving-ERA5-soil-moisture-data-as-csv-file)
- [3. Calibrating all parameters in the system](#Calibrating-all-parameters-in-the-system)
   - [3.1 Extracting the data and determine the FoS in the period of interest](#Extracting-the-data-and-determine-the-FoS-in-the-period-of-interest)
- [4. Sensitivity testing on all parameters](#Sensitivity-testing-on-all-parameters)
   - [4.1 Extracting the control data and all testing data](#Extracting-the-control-data-and-all-testing-data)
   - [4.2 Processing data in geotehnical model and stored in individual lists](#Processing-data-in-geotechnical-model-and-stored-in-individual-lists)
   - [4.3 Organising and plotting difference in FoS caused by varied parameter values](#Organising-and-plotting-difference-in-FoS-caused-by-varied-parameter-values)
- [5. Additional sensitivity testing on parameters](#Additional-sensitivity-testing-on-parameters)
   - [5.1 Testing significance of Alpha value](#Alpha-value-testing)
   - [5.2 Testing significance of drying to wetting reduce percentage](#Drying-to-wetting-reduce-percentage-testing)
   - [5.3 Testing significance of Field capacity](#Testing-significance-of-Field-capacity)
   - [5.4 Testing significance of Effective friction angle](#Testing-significance-of-Effective-friction-angle)
   - [5.5 Testing significance of slope angle](#Testing-significance-of-slope-angle)
   - [5.6 Testing significance of smax value](#Testing-significance-of-smax_value)
- [6. Initial approach to integrate Monte-Carlo simulation](#Initial-approach-to-integrate-Monte-Carlo-simulation)


<a name="Flow-diagram-of-the-slope-stability-modelling-system"></a>
## Flow diagram of the slope stability modelling system

[![Screenshot-2022-08-15-at-17-56-12.png](https://i.postimg.cc/CL1PLb48/Screenshot-2022-08-15-at-17-56-12.png)](https://postimg.cc/yWMXLgt6)

<a name="ERA5-Data-extraction-and-pre-processing"></a>
## 1. ERA5 Data extraction and data pre-processing
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/1. ERSS data extraction and pre-processing.ipynb

The aim of this process is to download the entire range of data that is given 
from ERA5, i.e. 1959-present. Which would allow the system to pick up 
a wider range of data and spot the periods in the past which the slope 
were unstable when the FoS is at similar range as the case studies.

<a name="ERA5-data-downloading-by-CDS-API"></a>
### 1.1 ERA5 data downloading by CDS API 
Import all packages required for the tasks
```python
import cdsapi
import xarray as xr
import cfgrib
import cf2cdm
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.dates as mdates
from itertools import chain
from datetime import datetime
from datetime import time
import datetime
import math
from plotnine import *
from matplotlib.dates import DateFormatter
import pyreadr
import rpy2
import cdsapi
```

Downloading data from the Climate Data Store API, this section of the code is currently mark as markdown as an alternative 
method was used for the case studies, the data were downloaded from the Copernicus website directly as it offers a faster 
and clearer method of downloading data, where once a request was sent, the downloading process can be 
done offline until the request is completed for download.

Input
- `product_type` = the type of data required, in the example, ERA5 reanalysis data is required
- `variable` = parameters required for the model, in the example, potential evaporation and total precipitation is needed
- `year`, `month`, `day`, `time` = specify the dates and the time that are required, resolution of the data given is 1hour
- `area` = Specifying the location of interest, resolution of the ERA5 model is 0.1degrees by 0.1 degrees
- `format` = Data from ECMWF on Seasonal forecast is given as grib file
- `name.grib` = Enter the `name` for the grib file

```python
c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
             'potential_evaporation', 'total_precipitation',
         ],
        'year': [
            '2009', '2010', '2011',
            '2012',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            56.8, -4.7, 56.7,
            -4.6,
        ],
        'format': 'grib',
    },
    'name.grib')
```
* Please note that there are a size limit to the downloading of data for 
both the API and the website, when all hourly data are required the maximum size of 
each data file can only be a maximum of 5 years period.

<a name="Grib-file-data-extraction-and-stored-as-variables"></a>
### 1.2 Grib file data extraction and stored as variables 

The `combine_data` function is designed to read the given dataset and locate all 
data with the longitude and latitude set and output the 3 data lists. 
This function is then repeated for all dataset taken for a case-study and stored in individual variables.

The grib files of the examples can be found in 
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/Data file/1. CDS Downloaded Data

Input
- `data_set` = grib data file obtained from CDS
- `-4.84` = Longitude of the point of interest, resolution of 0.1 degrees, Example - Rest and be Thankful Pass
- `56.23` = Latitude of the point of interest, resolution of 0.1 degrees, Example - Rest and be Thankful Pass
```python
def combine_data(data_set):
    ds = xr.open_dataset(data_set, engine='cfgrib',
                      backend_kwargs={'read_keys': ['experimentVersionNumber']})
    df = ds.to_dataframe()
    df_index=df.index
    test = df_index.to_frame()
    df = test.join(df)
    df.reset_index(drop=True, inplace=True)
    df = df.loc[df['longitude'] == -4.84]
    df = df.loc[df['latitude'] == 56.23]
    df.dropna(subset = ['tp'], inplace=True)
    dp = df['tp']
    pev = df['pev']
    valid_time = df['valid_time']
    return dp, pev, valid_time

dp_1959_1963 = combine_data('Rest_and_Be_Thankful_Pass_1959_1963.grib')[0]
pev_1959_1963 = combine_data('Rest_and_Be_Thankful_Pass_1959_1963.grib')[1]
valid_time_1959_1963 = combine_data('Rest_and_Be_Thankful_Pass_1959_1963.grib')[2]

dp_1964_1968 = combine_data('Rest_and_Be_Thankful_Pass_1964_1968.grib')[0]
pev_1964_1968 = combine_data('Rest_and_Be_Thankful_Pass_1964_1968.grib')[1]
valid_time_1964_1968 = combine_data('Rest_and_Be_Thankful_Pass_1964_1968.grib')[2]

dp_1969_1973 = combine_data('Rest_and_Be_Thankful_Pass_1969_1973.grib')[0]
pev_1969_1973 = combine_data('Rest_and_Be_Thankful_Pass_1969_1973.grib')[1]
valid_time_1969_1973 = combine_data('Rest_and_Be_Thankful_Pass_1969_1973.grib')[2]

dp_1974_1978 = combine_data('Rest_and_Be_Thankful_Pass_1974_1978.grib')[0]
pev_1974_1978 = combine_data('Rest_and_Be_Thankful_Pass_1974_1978.grib')[1]
valid_time_1974_1978 = combine_data('Rest_and_Be_Thankful_Pass_1974_1978.grib')[2]

dp_1979_1983 = combine_data('Rest_and_Be_Thankful_Pass_1979_1983.grib')[0]
pev_1979_1983 = combine_data('Rest_and_Be_Thankful_Pass_1979_1983.grib')[1]
valid_time_1979_1983 = combine_data('Rest_and_Be_Thankful_Pass_1979_1983.grib')[2]

dp_1984_1988 = combine_data('Rest_and_Be_Thankful_Pass_1984_1988.grib')[0]
pev_1984_1988 = combine_data('Rest_and_Be_Thankful_Pass_1984_1988.grib')[1]
valid_time_1984_1988 = combine_data('Rest_and_Be_Thankful_Pass_1984_1988.grib')[2]

dp_1989_1993 = combine_data('Rest_and_Be_Thankful_Pass_1989_1993.grib')[0]
pev_1989_1993 = combine_data('Rest_and_Be_Thankful_Pass_1989_1993.grib')[1]
valid_time_1989_1993 = combine_data('Rest_and_Be_Thankful_Pass_1989_1993.grib')[2]

dp_1994_1998 = combine_data('Rest_and_Be_Thankful_Pass_1994_1998.grib')[0]
pev_1994_1998 = combine_data('Rest_and_Be_Thankful_Pass_1994_1998.grib')[1]
valid_time_1994_1998 = combine_data('Rest_and_Be_Thankful_Pass_1994_1998.grib')[2]

dp_1999_2003 = combine_data('Rest_and_Be_Thankful_Pass_1999_2003.grib')[0]
pev_1999_2003 = combine_data('Rest_and_Be_Thankful_Pass_1999_2003.grib')[1]
valid_time_1999_2003 = combine_data('Rest_and_Be_Thankful_Pass_1999_2003.grib')[2]

dp_2004_2008 = combine_data('Rest_and_Be_Thankful_Pass_2004_2008.grib')[0]
pev_2004_2008 = combine_data('Rest_and_Be_Thankful_Pass_2004_2008.grib')[1]
valid_time_2004_2008 = combine_data('Rest_and_Be_Thankful_Pass_2004_2008.grib')[2]

dp_2009_2013 = combine_data('Rest_and_Be_Thankful_Pass_2009_2013.grib')[0]
pev_2009_2013 = combine_data('Rest_and_Be_Thankful_Pass_2009_2013.grib')[1]
valid_time_2009_2013 = combine_data('Rest_and_Be_Thankful_Pass_2009_2013.grib')[2]

dp_2014_2018 = combine_data('Rest_and_Be_Thankful_Pass_2014_2018.grib')[0]
pev_2014_2018 = combine_data('Rest_and_Be_Thankful_Pass_2014_2018.grib')[1]
valid_time_2014_2018 = combine_data('Rest_and_Be_Thankful_Pass_2014_2018.grib')[2]

dp_2019_2021 = combine_data('Rest_and_Be_Thankful_Pass_2019_2021.grib')[0]
pev_2019_2021 = combine_data('Rest_and_Be_Thankful_Pass_2019_2021.grib')[1]
valid_time_2019_2021 = combine_data('Rest_and_Be_Thankful_Pass_2019_2021.grib')[2]
```

<a name="Merging-all-data-into-3-lists"></a>
### 1.3 Merging all data into 3 lists

Combining all data from different years range into big lists, so that all data is listed from 1959-2021
* All data files to be downloaded in 5 years length from 1959 to 2021
```python
dp = []
pev = []
valid_time = []
dp.extend(dp_1959_1963)
dp.extend(dp_1964_1968)
dp.extend(dp_1969_1973)
dp.extend(dp_1974_1978)
dp.extend(dp_1979_1983)
dp.extend(dp_1984_1988)
dp.extend(dp_1989_1993)
dp.extend(dp_1994_1998)
dp.extend(dp_1999_2003)
dp.extend(dp_2004_2008)
dp.extend(dp_2009_2013)
dp.extend(dp_2014_2018)
dp.extend(dp_2019_2021)

pev.extend(pev_1959_1963)
pev.extend(pev_1964_1968)
pev.extend(pev_1969_1973)
pev.extend(pev_1974_1978)
pev.extend(pev_1979_1983)
pev.extend(pev_1984_1988)
pev.extend(pev_1989_1993)
pev.extend(pev_1994_1998)
pev.extend(pev_1999_2003)
pev.extend(pev_2004_2008)
pev.extend(pev_2009_2013)
pev.extend(pev_2014_2018)
pev.extend(pev_2019_2021)

valid_time.extend(valid_time_1959_1963)
valid_time.extend(valid_time_1964_1968)
valid_time.extend(valid_time_1969_1973)
valid_time.extend(valid_time_1974_1978)
valid_time.extend(valid_time_1979_1983)
valid_time.extend(valid_time_1984_1988)
valid_time.extend(valid_time_1989_1993)
valid_time.extend(valid_time_1994_1998)
valid_time.extend(valid_time_1999_2003)
valid_time.extend(valid_time_2004_2008)
valid_time.extend(valid_time_2009_2013)
valid_time.extend(valid_time_2014_2018)
valid_time.extend(valid_time_2019_2021)
```
<a name="Data-pre-processing"></a>
### 1.4 Data pre-processing 

This part of the codes convert the data into defined parameters
and units used in the hydrological model, where precipitation is converted from m to mm,
and potential evaporation is converted from neagative values in m as evaporation to positive values in mm

```python
def m_to_mm(m_data):
    mm = []
    for j in m_data:
        mm.append(j*1000)
    return mm
dp_mm = m_to_mm(dp)

def pev_to_mm(m_data): 
    mm = []
    for j in m_data:
        mm.append(j*-1000)
    return mm
pev_mm = pev_to_mm(pev)
```

<a name="Saving-data-as-csv-file"></a>
### 1.5 Saving data as csv file

This part of the codes set data lists as dictionary to be organised
in a dataframe which is then gathered as a csv file

Input
- `df.to_csv` = name of the csv file
```python
dict = {'Observed Date': valid_time, 'Hourly Precipitation': dp_mm, 'Potential Evaporation Rate': pev_mm}
df = pd.DataFrame(dict)
df.to_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_r_read_data.csv')
```

The example file is saved in 
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/2. Hydrological Model Read Data

<a name="Hydrological-Processing"></a>
## 2. Hydrological Processing

File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/2. ERSS ERA5 R script.R

Following converting the data into csv files, 
the data can then be processed in the hydrological model to model
the change in soil moisture across the period of interest. 
The model was written by Dr Thomas Kjeldsen and was modified to target
on the need of interest of the research.

<a name="DAYMOD-Soil-moisture-model"></a>
### 2.1 DAYMOD Soil moisture model

The model shown is an extension to the existing DAYMOD model which is
capable of estimating the soil moisture of the soil by taking in 
parameter such as:
- `smax` = The maximum amount of water that the soil can hold
- `FC` = Field capacity of the soil, the upper limit which when soil 
mositure exceed the limit, drainage to deeper soil layer will take place
- `RC` = Root Capacity of the soil, the lower limit which when 
soil mositure fall below the limit, evaporation does not occur as the
suction in the soil is extremely strong which prevent evaporation
- `dkt` = Drainage coefficient, the parameter that determine the amount
drainage when soil moisture is higher than the FC
- `m` = Initial soil moisture
- `rrain` = Hourly precipitation
- `EP` = Hourly potential evaporation
- `d.t` = time step of the model, in this example the model runs
with a resolution of 1 hour

```R
pdsoil<-function(smaxi,FC,RC,dkt,m,rrain,Ep,d.t) 
{
  
  # define jump matrix
  SM<-1/smaxi
  m.limit<-c((RC*SM),(FC*SM),SM) #capcity in % * smaxi
  jump<-matrix(c(0,1,2,3,0,4,5,6,0),nrow=3,ncol=3)
  konst<-matrix(c(0,0,Ep/(RC*SM), Ep,0,0, Ep,FC*SM,dkt),ncol=3,nrow=3,byrow=TRUE) # p,fc,k
  
  # initialise
  stmean<-0
  erf<-0
  evap<-0
  rechg<-0
  
  # find starting zone, and set ea, k and fc accordingly
  zone0<-zone.check(m,RC,FC,SM)
  
  # Calculate soil moisture
  m0<-m
  m1<-m.t(d.t,m0,rrain,konst,SM,zone0)
  if(m1<0) {m1<-0}
  if(m1>SM) {m1<-SM}
  
  # check if soil moisture has crossed boundary
  zone1<-zone.check(m1,RC,FC,SM)
  
  # if zone has changed (i.e. jump!=0) then
  
  # First, find the fraction of d.t spend in zone0, tt
  # Then calculate mean soil moisture in zone0
  # Finally calculate soil moisture for (1-tt)*d.t time in zone 1
  
  if(jump[zone0,zone1]!=0)  
  {
    # find fraction tt that makes difference between m.limit and m0 as small as possible
    # in other words, how long does it take to get from m0 to fill-up zone0
    # Find value of tt that gives (m.t-m.limit)==0
    
    # These two lines are necessary as the reference limits (m.limit) are conditional on
    # increasing or decreasing soil moisture
    if(m1>m0) {zone<-zone0}
    if(m1<m0) {zone<-zone1}
    
    # Pass zone0, as this function is finding the time (tt) spend in zone0 before crossing the limit
    #root.tt<-uniroot(fun.cross.m,interval=c(0,1),lower=0,upper=1,m0,rrain,konst,SM,zone0,m.limit[zone])
    root.tt<-uniroot(fun.cross.m,interval=c(0,1),lower=0,upper=1,m0,rrain,konst,SM,zone0,m.limit[zone])
    tt<-root.tt$root
    
    d.tt<-d.t*tt #time spent in a certain zone
    stmean<-(m0+m.limit[zone])/2 
    smean<-stmean*d.tt 
    
    if(m1<RC*SM) { #Zone 1 Evaporation
      evap<-d.tt*stmean/(RC*SM)
      rechg<-0 }
    if(m1>=RC*SM & m1<FC) { #Zone 2 Evaporation
      evap<-d.tt
      rechg<-0}
    if(m1>FC*SM) {#Zone 3 Evaporation+Drainage
      evap<-d.tt
      rechg<-d.tt*dkt*(stmean-FC*SM) }
    
    # Now calculate soil moisture after (d.t-d.tt) above/below the threshold
    m0<-m.limit[zone]
    m1<-m.t((d.t-d.tt),m0,rrain,konst,SM,zone1)  #zone1, as this is after the threshold is crossed
    if(m1<0) {m1<-0}
    if(m1>SM) {m1<-SM}
    
    stmean<-(m1+m.limit[zone])/2
    if(m1<RC*SM) {
      evap<-evap+d.tt*stmean/(RC*SM)
      rechg<-rechg+0 }
    if(m1>=RC*SM & m1<FC*SM) {
      evap<-evap+d.tt
      rechg<-rechg+0}
    if(m1>FC*SM) {
      evap<-evap+d.tt
      rechg<-rechg+d.tt*dkt*(stmean-FC*SM) }
  }
  
  # If no change of zone is observed
  if(jump[zone0,zone1]==0)
  {
    stmean<-(m0+m1)/2
    smean<-stmean*d.t
    
    if(m1<RC*SM) {
      evap<-d.t*stmean/(RC*SM)
      rechg<-0 }
    if(m1>=RC*SM & m1<FC*SM) {
      evap<-d.t
      rechg<-0 }
    if(m1>=FC*SM) {
      evap<-d.t
      rechg<-d.t*dkt*(stmean-FC*SM) }
  }
  
  # Form ouput variables
  
  stmean<-stmean*smaxi
  if(stmean>1) stmean<-1
  erf<-(1-sqrt(1-stmean))*rrain
  evap<-evap*Ep
  
  pdsoil<-c(m1,erf,evap,rechg)
  pdsoil
  
}

m.t<-function(d.t,m0,rrain,konst,SM,zone)
{
  
  p<-konst[zone,1]
  fc<-konst[zone,2]
  k<-konst[zone,3]
  
  # define convenience variable
  i.star<-rrain*d.t/(2*SM)
  E<-(p*d.t-k*d.t*fc)/SM
  k.star<-1+k*d.t/2
  
  
  m.ratio<-(m0/SM)*(2/k.star-1)-E/k.star-(i.star/k.star)^2+(i.star/k.star)*sqrt((i.star/k.star)^2+4*(k.star-m0/SM+E/2)/k.star)
  
  m.t<-m.ratio*SM
  if(m.t<=0) m.t<-0
  m.t
}

zone.check<-function(m,RC,FC,SM)
{
  if(m<RC*SM) zone.check<-1
  if(m>=RC*SM & m<FC*SM) zone.check<-2
  if(m>=FC*SM) zone.check<-3
  zone.check
}-m.t(d.t,m0,rrain,konst,SM,zone0)
  fun.cross.m<-(m.limit.s-m1)
  fun.cross.m
}
```

Setting out the environment of the model
```R
smax<-100
smaxi<-1/smax
FC<-0.8 
RC<-0.214 
dkt<-2.0/24
m<-79.49
```

<a name="Reading-ERA5-data-in-R"></a>
### 2.2 Reading ERA5 data in R

File - Documents/General/1. Existing Model/1.1 Existing ERA5 Model/ERA5 R script.R

```R
y <- read.csv("Rest_and_Be_Thankful_Pass_1959_to_2021_r_read_data.csv", TRUE, ",")


p<-y$Hourly.Precipitation #Extracting Hourly Precipitation data  
Ep<-y$Potential.Evaporation.Rate #Extracting Hourly Evaporation data
Od<-y$Observed.Date #Extracting valid datetime
d.t<-1 # Time step, one unit of time step of ERA5 data (hours)
np<-length(p) #Length of the data
```

<a name="Processing-ERA5-soil-moisture-data"></a>
### 2.3 Processing ERA5 soil moisture data 

As ERA5 data only contain a set of reanalyzed data, therefore only 1 loop is required.

```R
for(i in 1:np)
{
  out<-pdsoil(smaxi,FC,RC,dkt,soil.b4,p[i],Ep[i],d.t)
  soil.b4<-out[1] #soil moisture in previous time step
  soil.m[i]<-out[1]#soil moisture in current time step
  qs[i]<-out[2] #surface runoff
  evapa[i]<-out[3] #evaporation
  rch[i]<-out[4] #recharge, draiange
  }
}
```

<a name="Saving-ERA5-soil-moisture-data-as-csv-file"></a>
### 2.4 Saving ERA5 soil moisture data as csv file

The data will then be saved as csv file to allow the Geotechnical 
model accessing the soil moisture data

Input
- `data.frame(Od,soil.m)` = all varaiables that is of interest
of the geotechnical model, if required other variables such as 
surface runoff and etc can be obtained
```R
df<- data.frame(Od,soil.m)
write.csv(df,"C:/Users/Alan/Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv", row.names = FALSE)
```

The control dataset is saved in 
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/3. Control Dataset
The 10%-20% varied dataset is saved in 
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/4. Varied parameter datasets
Files for additional testing in saved in 
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.3 Varied FC datasets
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.6 Varied smax datasets

<a name="Calibrating-all-parameters-in-the-system"></a>
## 3. Calibrating all parameters in the system

File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/3. Calibration model.ipynb

This is a model developed to test and calibrate the parameters to achieve FoS
values in a range that are realistic, this would require a more detailed sensitivity testing
to achieve and therefore this model should be used with the `sensitivity_test_overall.ipynb`

<a name="Extracting-the-data-and-determine-the-FoS-in-the-period-of-interest"></a>
### 3.1 Extracting the data and determine the FoS in the period of interest

Import all packages required for the tasks
```python
import pandas as pd
%matplotlib inline
from matplotlib import pyplot as plt
import numpy as np
import csv
import matplotlib
import matplotlib.dates as mdates
import pyreadr
from itertools import chain
import math
from datetime import datetime
from statistics import NormalDist
import tempfile
import statistics
```

Reading the calibrated/control data and store the data as variables,
the example shown is a calibrated data with smax = 500, FC = 0.8, RC = 0.214, dkt = 1/12.
The period of interest is chosen to be from 2000-2021 and therefore relative 
constraint is given to the soil data and the datatime

- `smax` = The maximum amount of water that the soil can hold
- `FC` = Field capacity of the soil, the upper limit which when soil 
mositure exceed the limit, drainage to deeper soil layer will take place
- `RC` = Root Capacity of the soil, the lower limit which when 
soil mositure fall below the limit, evaporation does not occur as the
suction in the soil is extremely strong which prevent evaporation
- `dkt` = Drainage coefficient, the parameter that determine the amount
drainage when soil moisture is higher than the FC
- `soil_moisture_list_2000_2021` = Soil moisture of the typical soil in the slope,
- where the maximum soil moisture of the soil is smax

Input
- `pd.read_csv` = name of the csv file
- `359393:552257` = constraint given to the data to get the period of interest
```python
df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv')

df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')
df['date'] = pd.to_datetime(df['Od']).dt.date
df['time'] = pd.to_datetime(df['Od']).dt.time

datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]
 
od = df['date']
od_l = od.tolist()

ot = df['time']
ot_l = ot.tolist()

soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]
```


Converting hourly soil moisture data and datetime data to daily data to reduce
the overall size of the dataset to enable plotting with matplotlib
```python
def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily = hourly_to_daily_data(soil_moisture_list_2000_2021)
date_daily = hourly_to_daily_date(datetime_list_2000_2021)
```

Estimated value for the parameters required in the model, these parameters can be 
updated once calibrated for different case-studies

Input
- `e_cohesion` = Effective Cohesion in kPa
(Assumption from literature)
- `e_friction_ang` = Degrees Effective angle of internal friction 
(Assumption from literature)
- `rate_S_to_suction` = Degrees Angle indicating the rate of increase
of shear strength relative to increase matric suction (Assumption from literature)
- `soil_weight` = Soil Unit Weight in kN/m3 (Estimated)
- `slope_angle` = Degrees Slope angle (Assumed)
- `n` = Van Genuchten parameter n (Estimated for glacial till materials)
- `m` = Van Genuchten parameter m (Derived from n)
- `alpha` = Van Genuchten parameter alpha 
(Estimated for glacial till materials)
- `r_wc` = Residual water content (Estimated for glacial till materials)
- `s_wc` = Saturated water content (Estimated for glacial till materials)
- `porosity` = porosity of the soil (Estimated for glacial till materials)
```python
e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25 
n = 2.3 
m = 1-1/n 
alpha = 0.1 #0.0054 
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75
smax = 500
```

The failure depth of the slope was found to be a function of the soil moisture and porosity of the typical soil.
This function is designed to estimate the failure depth across the period of interest.

```python
def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d

fail_plane_d = finding_failure_depth(soil_moisture_list_2000_2021, porosity, smax)
fail_plane_d_daily = finding_failure_depth(soil_m_daily, porosity, smax)
```

The `soil_m_to_wetting_drying_suction`function determine the relative soil suction with the soil moisture. 
As the soil suction behaviour varies depending on whether the soil is wetting or drying. 
This model accounts for that effect by taking an average of the soil moisture in the `accounted_time` 
(set to 5 days in advance of the day of interest in example) and compared to the day of interest, 
if the average soil moisture is higher than the soil moisture on that day, the soil locates in the drying curve, 
if the average soil moisture is lower than the soil moisture on that day, the soil locates in the wetting curve and a 
reduced percentage is given to estimate the difference between the wetting and drying curve.
An intermediate percentage is also given when the soil moisture on that day is higher than the previous day,
the intermediate percentage will apply, acting as an intermediate step. 

- `accounted_time` = days accounted in advance of the day of interest to be included
- `reduce_percentage` = the reduced percentage given to the wetting curve compared against the drying curve
- `intermidate_percentage` = the percentage given as a step of the soil before reaching from drying to wetting curve
```python
def soil_m_to_wetting_drying_suction(soil_m, r_wc, s_wc, n, m, alpha, accounted_time,
                                     reduce_percentage, intermediate_percentage, smax):
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/m)
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                if soil_m[i] > soil_m[i-1]:
                    actual_suction = suction[i]*intermediate_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction

matric_suction_updated = soil_m_to_wetting_drying_suction(soil_moisture_list_2000_2021, r_wc, s_wc, n, m, alpha, 5,0.5,0.8, smax)
matric_suction_updated_1 = soil_m_to_wetting_drying_suction(soil_m_daily, r_wc, s_wc, n, m, alpha, 5,0.5,0.8, smax)
```

The `FoS_cal` function determine the relative FoS for a given slope given the soil suction value is known.

```python
def FoS_cal(e_cohesion, matric_suction, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d, e_friction_ang):
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS

FoS = FoS_cal(e_cohesion, matric_suction_updated, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d, e_friction_ang)
FoS_1 = FoS_cal(e_cohesion, matric_suction_updated_1, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d_daily, e_friction_ang)
```

Plotting FoS against time graph and indicating the landslide events in the period of interest
```python
fig, ax = plt.subplots()

feb_obs_plot = ax.plot(date_daily[1096:1461], FoS_1[1096:1461],
                      color = 'b',
                      label = 'Historic observed FoS' )

plt.plot(date_daily[1120], FoS_1[1120],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red", label = 'Landslide event')

plt.plot(date_daily[1428], FoS_1[1428],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale
ax.set_ylim(1.7,3)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval = 1))
plt.xlabel('Obeserved Period')
plt.ylabel('FoS')
plt.title('Variation of FoS in Rest and Be Thankful Pass over Observed Period')
leg = plt.legend(loc = 1)
fig.set_size_inches((20, 10))

fig, ax = plt.subplots()

feb_obs_plot = ax.plot(datetime_list_2000_2021, FoS,
                      color = 'b',
                      label = 'Historic observed FoS' )

plt.plot(datetime_list_2000_2021[26879], FoS[26879],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red", label = 'Landslide event')

plt.plot(datetime_list_2000_2021[34271], FoS[34271],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[35495], FoS[35495],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[68567], FoS[68567],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[72335], FoS[72335],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[77231], FoS[77231],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[84911], FoS[84911],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[112943], FoS[112943],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[139607], FoS[139607],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[140207], FoS[140207],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

plt.plot(datetime_list_2000_2021[185471], FoS[185471],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale 
ax.xaxis.set_major_locator(mdates.YearLocator(1))
#ax.set(xlim = [datetime_list_2000_2021[0], datetime_list_2000_2021[192862]])
plt.xlabel('Obeserved Period')
plt.ylabel('FoS')
plt.title('Variation of FoS in Rest and Be Thankful Pass over Observed Period')
leg = plt.legend(loc = 1)
fig.set_size_inches((20, 10))

fig, ax = plt.subplots()

feb_obs_plot = ax.plot(datetime_list_2000_2021[26279:35063], FoS[26279:35063],
                      color = 'b',
                      label = 'Historic observed FoS' )

plt.plot(datetime_list_2000_2021[26879], FoS[26879],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red", label = 'Landslide event')

plt.plot(datetime_list_2000_2021[34271], FoS[34271],
         marker="o", markersize=10, markeredgecolor="red", 
         markerfacecolor="red")

# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale 
ax.xaxis.set_major_locator(mdates.MonthLocator(interval = 1))
ax.set(xlim = [datetime_list_2000_2021[26279], datetime_list_2000_2021[35063]])
plt.xlabel('Obeserved Period')
plt.ylabel('FoS')
plt.title('Variation of FoS in Rest and Be Thankful Pass over Observed Period')
leg = plt.legend(loc = 1)
fig.set_size_inches((20, 10))
```


<a name="Sensitivity-testing-on-all-parameters"></a>
## 4. Sensitivity testing on all parameters

File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/Sensitivity test model.ipynb

<a name="Extracting-the-control-data-and-all-testing-data"></a>
### 4.1 Extracting the control data and all testing data

Import all packages required for the tasks
```python
import pandas as pd
%matplotlib inline
from matplotlib import pyplot as plt
import numpy as np
import csv
import matplotlib
import matplotlib.dates as mdates
from itertools import chain
import math
from datetime import datetime
from statistics import NormalDist
import tempfile
import statistics
```

As some parameters are specific to the hydrological model only,
sensitivity tests were performed to the data with +/- 10% and 20% to the parameter,
e.g control value of smax = 500, +/- 10% = 450, 550, +/- 20% = 400, 600 and etc.

The control model, df are performed using smax = 500, FC = 0.8, RC = 0.214, dkt = 2/24

Input
- `pd.read_csv()` = csv file name to be read
```python
df_smax_0_8 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_400.csv') 
df_smax_0_9 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_450.csv') 
df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') 
df_smax_1_1 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_550.csv')
df_smax_1_2 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_600.csv')

df_dkt_0_8 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_dkt_80%.csv') 
df_dkt_0_9 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_dkt_90%.csv') 
df_dkt_1_1 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_dkt_110%.csv') 
df_dkt_1_2 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_dkt_120%.csv') 

df_FC_0_8 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_80%.csv')
df_FC_0_9 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_90%.csv')
df_FC_1_1 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_110%.csv')
df_FC_1_2 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_120%.csv')
```

Obtain the control soil moisture data and its relative date and time reference in form of a datetime variable,
the data is then constraint to between 2000 - 2021 and stored as list in `datetime_list_2000_2021`, `soil_moisture_list_2000_2021`.
```python
df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')

datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]

soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]
```

Similar process to the control dataset, all other datasets based on varied hydrological parameters are constraint and stored in individual list
```python
soil_moisture_0_8 = df_smax_0_8['soil.m']
soil_moisture_list_0_8 = soil_moisture_0_8.tolist()
sm_list_2000_2021_smax_0_8 = soil_moisture_list_0_8[359393:552257]

soil_moisture_0_9 = df_smax_0_9['soil.m']
soil_moisture_list_0_9 = soil_moisture_0_9.tolist()
sm_list_2000_2021_smax_0_9 = soil_moisture_list_0_9[359393:552257]

soil_moisture_1_1 = df_smax_1_1['soil.m']
soil_moisture_list_1_1 = soil_moisture_1_1.tolist()
sm_list_2000_2021_smax_1_1 = soil_moisture_list_1_1[359393:552257]

soil_moisture_1_2 = df_smax_1_2['soil.m']
soil_moisture_list_1_2 = soil_moisture_1_2.tolist()
sm_list_2000_2021_smax_1_2 = soil_moisture_list_1_2[359393:552257]



soil_moisture_0_8 = df_dkt_0_8['soil.m']
soil_moisture_list_0_8 = soil_moisture_0_8.tolist()
sm_list_2000_2021_dkt_0_8 = soil_moisture_list_0_8[359393:552257]

soil_moisture_0_9 = df_dkt_0_9['soil.m']
soil_moisture_list_0_9 = soil_moisture_0_9.tolist()
sm_list_2000_2021_dkt_0_9 = soil_moisture_list_0_9[359393:552257]

soil_moisture_1_1 = df_dkt_1_1['soil.m']
soil_moisture_list_1_1 = soil_moisture_1_1.tolist()
sm_list_2000_2021_dkt_1_1 = soil_moisture_list_1_1[359393:552257]

soil_moisture_1_2 = df_dkt_1_2['soil.m']
soil_moisture_list_1_2 = soil_moisture_1_2.tolist()
sm_list_2000_2021_dkt_1_2 = soil_moisture_list_1_2[359393:552257]



soil_moisture_0_8 = df_FC_0_8['soil.m']
soil_moisture_list_0_8 = soil_moisture_0_8.tolist()
sm_list_2000_2021_FC_0_8 = soil_moisture_list_0_8[359393:552257]

soil_moisture_0_9 = df_FC_0_9['soil.m']
soil_moisture_list_0_9 = soil_moisture_0_9.tolist()
sm_list_2000_2021_FC_0_9 = soil_moisture_list_0_9[359393:552257]

soil_moisture_1_1 = df_FC_1_1['soil.m']
soil_moisture_list_1_1 = soil_moisture_1_1.tolist()
sm_list_2000_2021_FC_1_1 = soil_moisture_list_1_1[359393:552257]

soil_moisture_1_2 = df_FC_1_2['soil.m']
soil_moisture_list_1_2 = soil_moisture_1_2.tolist()
sm_list_2000_2021_FC_1_2 = soil_moisture_list_1_2[359393:552257]
```

<a name="Processing-data-in-geotechnical-model-and-stored-in-individual-lists"></a>
### 4.2 Processing data in geotechnical model and stored in individual lists

With the aim to reduce the size of the dataset to enable clearer presentation in 
graphs and better processing performance, all hourly soil moisture data from all datasets are converted into daily data
```python
def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily = hourly_to_daily_data(soil_moisture_list_2000_2021)
date_daily = hourly_to_daily_date(datetime_list_2000_2021)



soil_m_daily_smax_0_8 = hourly_to_daily_data(sm_list_2000_2021_smax_0_8)
soil_m_daily_smax_0_9 = hourly_to_daily_data(sm_list_2000_2021_smax_0_9)
soil_m_daily_smax_1_1 = hourly_to_daily_data(sm_list_2000_2021_smax_1_1)
soil_m_daily_smax_1_2 = hourly_to_daily_data(sm_list_2000_2021_smax_1_2)


soil_m_daily_dkt_0_8 = hourly_to_daily_data(sm_list_2000_2021_dkt_0_8)
soil_m_daily_dkt_0_9 = hourly_to_daily_data(sm_list_2000_2021_dkt_0_9)
soil_m_daily_dkt_1_1 = hourly_to_daily_data(sm_list_2000_2021_dkt_1_1)
soil_m_daily_dkt_1_2 = hourly_to_daily_data(sm_list_2000_2021_dkt_1_2)


soil_m_daily_FC_0_8 = hourly_to_daily_data(sm_list_2000_2021_FC_0_8)
soil_m_daily_FC_0_9 = hourly_to_daily_data(sm_list_2000_2021_FC_0_9)
soil_m_daily_FC_1_1 = hourly_to_daily_data(sm_list_2000_2021_FC_1_1)
soil_m_daily_FC_1_2 = hourly_to_daily_data(sm_list_2000_2021_FC_1_2)
```

Setting the base environment for the model, where the values shown in the scripts are the control values for the variables.

Input
- `e_cohesion` = Effective Cohesion in kPa
(Assumption from literature)
- `e_friction_ang` = Degrees Effective angle of internal friction 
(Assumption from literature)
- `rate_S_to_suction` = Degrees Angle indicating the rate of increase
of shear strength relative to increase matric suction (Assumption from literature)
- `soil_weight` = Soil Unit Weight in kN/m3 (Estimated)
- `slope_angle` = Degrees Slope angle (Assumed)
- `fail_plane_d` = Depth of failure plane in meters (assumed)
- `n` = Van Genuchten parameter n (Estimated for glacial till materials)
- `m` = Van Genuchten parameter m (Derived from n)
- `alpha` = Van Genuchten parameter alpha (Estimated from self calibration)
- `r_wc` = Residual water content (Estimated for glacial till materials)
- `s_wc` = Saturated water content (Estimated for glacial till materials)
- `porosity` = porosity of the soil (Estimated for glacial till materials)
- `local_smax` = maximum soil moisture content in m
```python
e_cohesion = 0 
e_friction_ang = 29 
rate_S_to_suction = 29 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
alpha = 0.1
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75 
local_smax = 500
```

Failure depth of the slope is estimated to be a function of the soil moisture and porosity, critical depth = soil moisture/porosity. 
This was set as a function and was performed with the control dataset and other datasets that has varied parameters which affect the soil moisture value.
The function was also repeated for another 4 times with the control soil moisture dataset with varying porosity to
investigate the significance of porosity to the result.
```python
def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d


fail_plane_d = finding_failure_depth(soil_moisture_list_2000_2021, porosity, local_smax)
fail_plane_d_daily = finding_failure_depth(soil_m_daily, porosity, local_smax)


fail_plane_d_smax_0_8 = finding_failure_depth(soil_m_daily_smax_0_8, porosity, local_smax*0.8)
fail_plane_d_smax_0_9 = finding_failure_depth(soil_m_daily_smax_0_9, porosity, local_smax*0.9)
fail_plane_d_smax_1_1 = finding_failure_depth(soil_m_daily_smax_1_1, porosity, local_smax*1.1)
fail_plane_d_smax_1_2 = finding_failure_depth(soil_m_daily_smax_1_2, porosity, local_smax*1.2)


fail_plane_d_porosity_0_8 = finding_failure_depth(soil_m_daily, porosity*0.8, local_smax)
fail_plane_d_porosity_0_9 = finding_failure_depth(soil_m_daily, porosity*0.9, local_smax)
fail_plane_d_porosity_1_1 = finding_failure_depth(soil_m_daily, porosity*1.1, local_smax)
fail_plane_d_porosity_1_2 = finding_failure_depth(soil_m_daily, porosity*1.2, local_smax)


fail_plane_d_dkt_0_8 = finding_failure_depth(soil_m_daily_dkt_0_8, porosity, local_smax)
fail_plane_d_dkt_0_9 = finding_failure_depth(soil_m_daily_dkt_0_9, porosity, local_smax)
fail_plane_d_dkt_1_1 = finding_failure_depth(soil_m_daily_dkt_1_1, porosity, local_smax)
fail_plane_d_dkt_1_2 = finding_failure_depth(soil_m_daily_dkt_1_2, porosity, local_smax)


fail_plane_d_FC_0_8 = finding_failure_depth(soil_m_daily_FC_0_8, porosity, local_smax)
fail_plane_d_FC_0_9 = finding_failure_depth(soil_m_daily_FC_0_9, porosity, local_smax)
fail_plane_d_FC_1_1 = finding_failure_depth(soil_m_daily_FC_1_1, porosity, local_smax)
fail_plane_d_FC_1_2 = finding_failure_depth(soil_m_daily_FC_1_2, porosity, local_smax)
```

The `soil_m_to_wetting_drying_suction`function determine the relative soil suction with the soil moisture. 
As the soil suction behaviour varies depending on whether the soil is wetting or drying. 
This model accounts for that effect by taking an average of the soil moisture in the `accounted_time` 
(set to 5 days in advance of the day of interest in example) and compared to the day of interest, 
if the average soil moisture is higher than the soil moisture on that day, the soil locates in the drying curve, 
if the average soil moisture is lower than the soil moisture on that day, the soil locates in the wetting curve and a 
reduced percentage is given to estimate the difference between the wetting and drying curve.

- `accounted_time` = days accounted in advance of the day of interest to be included
- `reduce_percentage` = the reduced percentage given to the wetting curve compared against the drying curve

The function was then used to process all soil moisture datasets and the control dataset. 
In order to perform sensitivity testing on the geotechnical model parameters, 
the control (soil moisture) dataset was repeated numerous times with varied geotechnical parameters
```python
def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = 0.5


matric_suction = soil_m_to_wetting_drying_suction(
    soil_moisture_list_2000_2021, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, local_smax)
matric_suction_daily = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha,
    accounted_time, reduce_percentage, local_smax)

ms_daily_n_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n*0.8, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_n_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n*0.9, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_n_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n*1.1, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_n_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n*1.2, alpha, accounted_time, reduce_percentage, local_smax)


ms_daily_alpha_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha*0.8, accounted_time, reduce_percentage, local_smax)
ms_daily_alpha_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha*0.9, accounted_time, reduce_percentage, local_smax)
ms_daily_alpha_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha*1.1, accounted_time, reduce_percentage, local_smax)
ms_daily_alpha_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha*1.2, accounted_time, reduce_percentage, local_smax)


ms_daily_reduce_perc_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage*0.8, local_smax)
ms_daily_reduce_perc_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage*0.9, local_smax)
ms_daily_reduce_perc_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage*1.1, local_smax)
ms_daily_reduce_perc_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage*1.2, local_smax)


ms_daily_r_wc_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc*0.8, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_r_wc_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc*0.9, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_r_wc_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc*1.1, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_r_wc_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc*1.2, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)


ms_daily_s_wc_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc*0.8, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_s_wc_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc*0.9, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_s_wc_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc*1.1, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_s_wc_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc*1.2, n, alpha, accounted_time, reduce_percentage, local_smax)


ms_daily_dkt_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily_dkt_0_8, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_dkt_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily_dkt_0_9, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_dkt_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily_dkt_1_1, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_dkt_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily_dkt_1_2, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)


ms_daily_FC_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_8, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_FC_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_9, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_FC_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_1_1, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)
ms_daily_FC_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_1_2, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax)

ms_daily_smax_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_0_8, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax*0.8)
ms_daily_smax_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_0_9, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax*0.9)
ms_daily_smax_1_1 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_1_1, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax*1.1)
ms_daily_smax_1_2 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_1_2, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, local_smax*1.2)
```

Taking a failure point of the slope and a safe point, i.e a low FoS (landslide occuring event), and a high FoS point
The 2 points is taken in all of the soil suction datasets and is appended into their relative list,
i.e 80% friction angle, slope angle and etc failure points is appended in `FoS_test_data_0_8_lp`. 
The order of which the varied parameters are appended into the list are given in `FoS_test_list`.
Therefore, varied friction angle results are located at place [0] and slope angle is placed at [1], etc.

Input
- `low point` = failure point data location in the list
- `high point` = safe point data location in the list

The high point can be identified by taking the control FoS list and pair it with `date_daily` with function `zip()` 
and rank the lists by smallest to the largest by function`sorted()`.
```python
FoS_test_list = ['friction angle', 'slope angle', 
            'Phi b','soil weight', 'n', 'alpha', 'porosity',
            'drying to wetting reduced percentage', 
            'residual water content','saturated water content', 
            'drainage coefficient', 'field capacity', 'smax']
FoS_test_data_0_8_lp = []
FoS_test_data_0_9_lp = []
FoS_test_data_lp = []
FoS_test_data_1_1_lp = []
FoS_test_data_1_2_lp = []


FoS_test_data_0_8_hp = []
FoS_test_data_0_9_hp = []
FoS_test_data_hp = []
FoS_test_data_1_1_hp = []
FoS_test_data_1_2_hp = []


def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

FoS = FoS_cal(e_cohesion, matric_suction_daily, rate_S_to_suction, 
              soil_weight, slope_angle, fail_plane_d_daily, e_friction_ang)
FoS_test_data_lp.append(FoS[low_point])
FoS_test_data_hp.append(FoS[high_point])

FoS_fric_ang_0_8 = FoS_cal(e_cohesion, matric_suction_daily, 
                           rate_S_to_suction, soil_weight, slope_angle,
                           fail_plane_d_daily, e_friction_ang*0.8)
FoS_fric_ang_0_9 = FoS_cal(e_cohesion, matric_suction_daily, 
                           rate_S_to_suction, soil_weight, slope_angle, 
                           fail_plane_d_daily, e_friction_ang*0.9)
FoS_fric_ang_1_1 = FoS_cal(e_cohesion, matric_suction_daily, 
                           rate_S_to_suction, soil_weight, slope_angle, 
                           fail_plane_d_daily, e_friction_ang*1.1)
FoS_fric_ang_1_2 = FoS_cal(e_cohesion, matric_suction_daily, 
                           rate_S_to_suction, soil_weight, slope_angle, 
                           fail_plane_d_daily, e_friction_ang*1.2)

FoS_test_data_0_8_lp.append(FoS_fric_ang_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_fric_ang_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_fric_ang_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_fric_ang_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_fric_ang_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_fric_ang_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_fric_ang_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_fric_ang_1_2[high_point])


FoS_slope_ang_0_8 = FoS_cal(e_cohesion, matric_suction_daily, 
                            rate_S_to_suction, soil_weight, slope_angle*0.8,
                            fail_plane_d_daily, e_friction_ang)
FoS_slope_ang_0_9 = FoS_cal(e_cohesion, matric_suction_daily, 
                            rate_S_to_suction, soil_weight, slope_angle*0.9, 
                            fail_plane_d_daily, e_friction_ang)
FoS_slope_ang_1_1 = FoS_cal(e_cohesion, matric_suction_daily, 
                            rate_S_to_suction, soil_weight, slope_angle*1.1,
                            fail_plane_d_daily, e_friction_ang)
FoS_slope_ang_1_2 = FoS_cal(e_cohesion, matric_suction_daily, 
                            rate_S_to_suction, soil_weight, slope_angle*1.2, 
                            fail_plane_d_daily, e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_slope_ang_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_slope_ang_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_slope_ang_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_slope_ang_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_slope_ang_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_slope_ang_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_slope_ang_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_slope_ang_1_2[high_point])


FoS_rate_S_to_suction_0_8 = FoS_cal(e_cohesion, matric_suction_daily, 
                                    rate_S_to_suction*0.8, soil_weight, 
                                    slope_angle, fail_plane_d_daily, 
                                    e_friction_ang)
FoS_rate_S_to_suction_0_9 = FoS_cal(e_cohesion, matric_suction_daily, 
                                    rate_S_to_suction*0.9, soil_weight, 
                                    slope_angle, fail_plane_d_daily, 
                                    e_friction_ang)
FoS_rate_S_to_suction_1_1 = FoS_cal(e_cohesion, matric_suction_daily, 
                                    rate_S_to_suction*1.1, soil_weight, 
                                    slope_angle, fail_plane_d_daily, 
                                    e_friction_ang)
FoS_rate_S_to_suction_1_2 = FoS_cal(e_cohesion, matric_suction_daily, 
                                    rate_S_to_suction*1.2, soil_weight, 
                                    slope_angle, fail_plane_d_daily, 
                                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_rate_S_to_suction_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_rate_S_to_suction_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_rate_S_to_suction_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_rate_S_to_suction_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_rate_S_to_suction_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_rate_S_to_suction_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_rate_S_to_suction_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_rate_S_to_suction_1_2[high_point])

FoS_soil_weight_0_8 = FoS_cal(e_cohesion, matric_suction_daily, 
                              rate_S_to_suction, soil_weight*0.8, slope_angle,
                              fail_plane_d_daily, e_friction_ang)
FoS_soil_weight_0_9 = FoS_cal(e_cohesion, matric_suction_daily, 
                              rate_S_to_suction, soil_weight*0.9, slope_angle,
                              fail_plane_d_daily, e_friction_ang)
FoS_soil_weight_1_1 = FoS_cal(e_cohesion, matric_suction_daily, 
                              rate_S_to_suction, soil_weight*1.1, slope_angle,
                              fail_plane_d_daily, e_friction_ang)
FoS_soil_weight_1_2 = FoS_cal(e_cohesion, matric_suction_daily, 
                              rate_S_to_suction, soil_weight*1.2, slope_angle,
                              fail_plane_d_daily, e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_soil_weight_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_soil_weight_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_soil_weight_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_soil_weight_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_soil_weight_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_soil_weight_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_soil_weight_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_soil_weight_1_2[high_point])


FoS_n_0_8 = FoS_cal(e_cohesion, ms_daily_n_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_n_0_9 = FoS_cal(e_cohesion, ms_daily_n_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_n_1_1 = FoS_cal(e_cohesion, ms_daily_n_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_n_1_2 = FoS_cal(e_cohesion, ms_daily_n_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_soil_weight_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_soil_weight_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_soil_weight_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_soil_weight_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_soil_weight_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_soil_weight_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_soil_weight_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_soil_weight_1_2[high_point])


FoS_alpha_0_8 = FoS_cal(e_cohesion, ms_daily_alpha_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_alpha_0_9 = FoS_cal(e_cohesion, ms_daily_alpha_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_alpha_1_1 = FoS_cal(e_cohesion, ms_daily_alpha_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_alpha_1_2 = FoS_cal(e_cohesion, ms_daily_alpha_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_alpha_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_alpha_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_alpha_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_alpha_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_alpha_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_alpha_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_alpha_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_alpha_1_2[high_point])


FoS_porosity_0_8 = FoS_cal(e_cohesion, matric_suction_daily, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_porosity_0_8, 
                    e_friction_ang)
FoS_porosity_0_9 = FoS_cal(e_cohesion, matric_suction_daily, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_porosity_0_9, 
                    e_friction_ang)
FoS_porosity_1_1 = FoS_cal(e_cohesion, matric_suction_daily, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_porosity_1_1, 
                    e_friction_ang)
FoS_porosity_1_2 = FoS_cal(e_cohesion, matric_suction_daily, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_porosity_1_2, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_porosity_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_porosity_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_porosity_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_porosity_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_porosity_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_porosity_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_porosity_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_porosity_1_2[high_point])


FoS_reduce_perc_0_8 = FoS_cal(e_cohesion, ms_daily_reduce_perc_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_reduce_perc_0_9 = FoS_cal(e_cohesion, ms_daily_reduce_perc_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_reduce_perc_1_1 = FoS_cal(e_cohesion, ms_daily_reduce_perc_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_reduce_perc_1_2 = FoS_cal(e_cohesion, ms_daily_reduce_perc_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_reduce_perc_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_reduce_perc_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_reduce_perc_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_reduce_perc_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_reduce_perc_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_reduce_perc_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_reduce_perc_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_reduce_perc_1_2[high_point])


FoS_r_wc_0_8 = FoS_cal(e_cohesion, ms_daily_r_wc_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_r_wc_0_9 = FoS_cal(e_cohesion, ms_daily_r_wc_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_r_wc_1_1 = FoS_cal(e_cohesion, ms_daily_r_wc_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_r_wc_1_2 = FoS_cal(e_cohesion, ms_daily_r_wc_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_r_wc_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_r_wc_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_r_wc_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_r_wc_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_r_wc_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_r_wc_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_r_wc_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_r_wc_1_2[high_point])


FoS_s_wc_0_8 = FoS_cal(e_cohesion, ms_daily_s_wc_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_s_wc_0_9 = FoS_cal(e_cohesion, ms_daily_s_wc_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_s_wc_1_1 = FoS_cal(e_cohesion, ms_daily_s_wc_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_s_wc_1_2 = FoS_cal(e_cohesion, ms_daily_s_wc_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_s_wc_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_s_wc_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_s_wc_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_s_wc_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_s_wc_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_s_wc_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_s_wc_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_s_wc_1_2[high_point])


FoS_dkt_0_8 = FoS_cal(e_cohesion, ms_daily_dkt_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_dkt_0_9 = FoS_cal(e_cohesion, ms_daily_dkt_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_dkt_1_1 = FoS_cal(e_cohesion, ms_daily_dkt_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
FoS_dkt_1_2 = FoS_cal(e_cohesion, ms_daily_dkt_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_dkt_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_dkt_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_dkt_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_dkt_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_dkt_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_dkt_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_dkt_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_dkt_1_2[high_point])


FoS_FC_0_8 = FoS_cal(e_cohesion, ms_daily_FC_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_8, 
                    e_friction_ang)
FoS_FC_0_9 = FoS_cal(e_cohesion, ms_daily_FC_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_9, 
                    e_friction_ang)
FoS_FC_1_1 = FoS_cal(e_cohesion, ms_daily_FC_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_1_1, 
                    e_friction_ang)
FoS_FC_1_2 = FoS_cal(e_cohesion, ms_daily_FC_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_1_2, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_FC_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_FC_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_FC_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_FC_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_FC_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_FC_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_FC_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_FC_1_2[high_point])


FoS_smax_0_8 = FoS_cal(e_cohesion, ms_daily_smax_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_0_8, 
                    e_friction_ang)
FoS_smax_0_9 = FoS_cal(e_cohesion, ms_daily_smax_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_0_9, 
                    e_friction_ang)
FoS_smax_1_1 = FoS_cal(e_cohesion, ms_daily_smax_1_1, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_1_1, 
                    e_friction_ang)
FoS_smax_1_2 = FoS_cal(e_cohesion, ms_daily_smax_1_2, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_1_2, 
                    e_friction_ang)

FoS_test_data_0_8_lp.append(FoS_smax_0_8[low_point])
FoS_test_data_0_8_hp.append(FoS_smax_0_8[high_point])
FoS_test_data_0_9_lp.append(FoS_smax_0_9[low_point])
FoS_test_data_0_9_hp.append(FoS_smax_0_9[high_point])
FoS_test_data_1_1_lp.append(FoS_smax_1_1[low_point])
FoS_test_data_1_1_hp.append(FoS_smax_1_1[high_point])
FoS_test_data_1_2_lp.append(FoS_smax_1_2[low_point])
FoS_test_data_1_2_hp.append(FoS_smax_1_2[high_point])
```

<a name="Organising-and-plotting-difference-in-FoS-caused-by-varied-parameter-values"></a>
### 4.3 Organising and plotting difference in FoS caused by varied parameter values

The code below first find the difference of FoS between the control dataset and each individual 120% varied parameter datasets.
The difference between the datasets is then paired with the 80%, 90%, 110% and 120% varied parameter datasets and ranked from largest difference to smallest.
Which prepares the data to be presented in a tornado plot where the significance of each parameter can be shown.
```python
FoS_diff_lp = []
for i in FoS_test_data_1_2_lp:
    diff = i - FoS_test_data_lp[0]
    FoS_diff_lp.append(abs(diff))

FoS_diff_lar_to_smal = [i for i in sorted(FoS_diff_lp, reverse=True)]
FoS_0_8_lp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_lp, FoS_test_data_0_8_lp))]
FoS_0_9_lp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_lp, FoS_test_data_0_9_lp))]
FoS_1_1_lp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_lp, FoS_test_data_1_1_lp))]
FoS_1_2_lp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_lp, FoS_test_data_1_2_lp))]
list_sorted = [i for j, i in sorted(zip(FoS_diff_lp, FoS_test_list))]


FoS_diff_hp = []
for i in FoS_test_data_0_8_hp:
    diff = i - FoS_test_data_hp[0]
    FoS_diff_hp.append(abs(diff))

FoS_diff_hp_lar_to_smal = [i for i in sorted(FoS_diff_hp, reverse=True)]
FoS_0_8_hp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_hp, FoS_test_data_0_8_hp))]
FoS_0_9_hp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_hp, FoS_test_data_0_9_hp))]
FoS_1_1_hp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_hp, FoS_test_data_1_1_hp))]
FoS_1_2_hp_smal_to_lar = [i for j, i in sorted(zip(FoS_diff_hp, FoS_test_data_1_2_hp))]
list_sorted = [i for j, i in sorted(zip(FoS_diff_hp, FoS_test_list))]

pos = np.arange(0, num_list,1)
```

Tornado plot was identified to be a valid method to give a graphical presentation on how each parameter affect the FoS results,
where the difference of FoS caused by a varied parameter with all other parameter kept constant can be visually observed. 
It is also worth keeping in mind that the significance of each parameter can differ when looking at a different range of FoS.
Therefore, studies have been carried out to investigate the significance of individual parameter in its entire possible range shown in section 5.

- `zip()` = the zip function pairs the inputted lists, e.g. let a and b to be 2 lists, list(zip(a,b)) = [[a[0],b[0]], [a[1], b[1]]]
- `plt.broken_barh` = in the example the broken barh function is used in plotting both the upper and lower bound hence divided into 2 section
`[(i, FoS_test_data_lp[0] - i), (FoS_test_data_lp[0], j - FoS_test_data_lp[0])]`, the first input `i` represent the base point of the plot and 
`FoS_test_data_lp[0] - i` is the difference between the control FoS value to the 80% varied FoS, which is the width that the plot is required to bring the plot to the control axis.
similar approach was then taken for the plot for 120% varied FoS
- `k-0.4` = the height the bar to be plotted in the graph
- `0.8` = the height of each bar
- `facecolors` = color of the bar
- `edgecolors` = color of the bar boundaries
- `linewidth` = width of the boundaries
- `set_label()` = setting label for the 2 hatches
- `title` = setting title of the graph
- `axvline()` = plotting a line in the graph, indicating the controlled datasets' FoS
- `yticks()` = matching the name of the varied parameters to their relative bars
- `legend()` = showing the label of the plot as legend 
```python
axes = plt.gca()  # (gca = get current axes)


for i, j, k in zip(FoS_0_8_lp_smal_to_lar, FoS_1_2_lp_smal_to_lar, pos):
     first = plt.broken_barh([(i, FoS_test_data_lp[0] - i), 
                     (FoS_test_data_lp[0], j - FoS_test_data_lp[0])],
                     (k-0.4, 0.8),
                     facecolors=['red', 'red'], 
                     edgecolors=['black', 'black'],
                     linewidth=1,)

for i, j, k in zip(FoS_0_9_lp_smal_to_lar, FoS_1_1_lp_smal_to_lar, pos):
    second = plt.broken_barh([(i, FoS_test_data_lp[0] - i), 
                     (FoS_test_data_lp[0], j - FoS_test_data_lp[0])],
                     (k-0.4, 0.8),
                     facecolors=['white', 'white'],  
                     edgecolors=['black', 'black'],
                     linewidth=1,)

first.set_label('FoS with +/- 20% in parameters')
second.set_label('FoS with +/- 10% in parameters')
    
axes.spines['left'].set_visible(False)
axes.spines['right'].set_visible(False)
axes.spines['bottom'].set_visible(False)
axes.xaxis.set_ticks_position('top')
axes.set_title('FoS')
plt.title('Sensitivity Check on all input parameters in low FoS')
#axes.set_xlim(2, 14)
plt.axvline(FoS_test_data_lp, color='blue', label = 'FoS with default parameters')
plt.yticks(pos, list_sorted)
leg = plt.legend(loc = 4, fontsize = 'small')


axes = plt.gca()  # (gca = get current axes)


for i, j, k in zip(FoS_0_8_hp_smal_to_lar, FoS_1_2_hp_smal_to_lar, pos):
     first = plt.broken_barh([(i, FoS_test_data_hp[0] - i), 
                     (FoS_test_data_hp[0], j - FoS_test_data_hp[0])],
                     (k-0.4, 0.8),
                     facecolors=['red', 'red'],  
                     edgecolors=['black', 'black'],
                     linewidth=1,)

for i, j, k in zip(FoS_0_9_hp_smal_to_lar, FoS_1_1_hp_smal_to_lar, pos):
    second = plt.broken_barh([(i, FoS_test_data_hp[0] - i), 
                     (FoS_test_data_hp[0], j - FoS_test_data_hp[0])],
                     (k-0.4, 0.8),
                     facecolors=['white', 'white'],  
                     edgecolors=['black', 'black'],
                     linewidth=1,)

first.set_label('FoS with +/- 20% in parameters')
second.set_label('FoS with +/- 10% in parameters')
    
axes.spines['left'].set_visible(False)
axes.spines['right'].set_visible(False)
axes.spines['bottom'].set_visible(False)
axes.xaxis.set_ticks_position('top')
axes.set_title('FoS')
plt.title('Sensitivity Check on all input parameters in high FoS')
#axes.set_xlim(2, 14)
plt.axvline(FoS_test_data_hp, color='blue', label = 'FoS with default parameters')
plt.yticks(pos, list_sorted)
leg = plt.legend(loc = 4, fontsize = 'small')
```

<a name="Additional-sensitivity-testing-on-parameters"></a>
## 5. Additional sensitivity testing on parameters

In order to investigate on the significance of parameters in a range of FoS that is more realistic, 
i.e FoS between 1-1.5 when slope fails and 2-3 when the slope is safe, 
an additional sensitivity testing was performed on a few parameters.

<a name="Alpha-value-testing"></a>
### 5.1 Testing significance of alpha value
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.1 Additional testing alpha.ipynb

Since alpha is a geotechnical parameter, the same control soil moisture dataset will be used.
The code shown below is a similar approach as in the sensitivity test,
where the data between 2000-2021 is extracted and stored in variables and converted into daily data.
In addition, the soil moisture is also converted from a value between 0 and smax, to a percentage to allow for easier interpretation.

```python
import pandas as pd
%matplotlib inline
from matplotlib import pyplot as plt
import numpy as np
import math
import statistics

df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') #extract soil moisture csv data to variable



df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')


datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]


soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]




def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily = hourly_to_daily_data(soil_moisture_list_2000_2021)
date_daily = hourly_to_daily_date(datetime_list_2000_2021)


soil_m_percentage = []
for i in range(len(soil_m_daily)):
    a = soil_m_daily[i]/500
    b = a*100
    soil_m_percentage.append(b)
```

Using Controlled parameters for the model, except for alpha where a range between 0.005 to 1 with
steps of 0.005 is given as shown in `alpha = [i for i in np.arange(0.005,1,0.005)]`

```python
e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63
porosity = 0.75 
alpha = [i for i in np.arange(0.005,1,0.005)]
```

Similar to the sensitivity testing, the failure plane depth is found by the `finding_failure_depth` function, 
and the same soil suction and FoS applies.
However, as there are 199 dataset from the varied alpha value,
both function is required to repeat for 199 times in a for loop, 
which will then extract the safe point and failure point FoS results as `FoS_lp` and `FoS_hp`.
- `alpha_change` = a list of list that contains the soil suction values of any alpha value in sublists
- `a` , `c` = a list that contains both the failure and safe point FoS and its relative alpha value paired
and presented in a tuple format
```python
def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d
 
fail_plane_d = finding_failure_depth(soil_moisture_list_2000_2021, porosity, 500)
fail_plane_d_daily = finding_failure_depth(soil_m_daily, porosity, 500)

def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = 0.5

def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

alpha_change = []
FoS_lp = []
FoS_hp = []

for i in alpha:
    ms_daily_alpha = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, i, accounted_time, reduce_percentage, 500)
    alpha_change.append(ms_daily_alpha)

for i in range(len(alpha_change)):
    FoS_alpha = FoS_cal(e_cohesion, alpha_change[i], rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
    FoS_lp.append(FoS_alpha[low_point])
    FoS_hp.append(FoS_alpha[high_point])

a = []
for i in range(len(FoS_lp)):
    b = tuple([alpha[i], FoS_lp[i]])
    a.append(b)
    
c = []
for i in range(len(FoS_hp)):
    d = tuple([alpha[i], FoS_hp[i]])
    c.append(d)
```

Plotting the FoS changes with the change in alpha value
```python
fig, ax = plt.subplots()

feb_obs_plot = ax.plot(alpha, FoS_lp,
                      color = 'black' )

plt.xlabel('Alpha')
plt.ylabel('FoS')
plt.title('FoS against Alpha variation')
fig.set_size_inches((20, 10))
```

As requested by Dr Cormac Reale, the alpha value of 0.005, 0.05, 0.1, 0.2, 0.4 are taken for the plot.
The `soil_moisture_list` is a list that contains soil moisture content between residual and saturated water content with 50 steps.
The soil suction value for each alpha value in the range of soil moisture is then calculated from `soil_m_to_ss`.
Finally, the soil suction of the 4 alpha values are plotted against the mositure content list.


```python
soil_moisture_list = []
interval = 49
alpha_1 = [0.005, 0.05, 0.1, 0.2, 0.4]
suction = []
for i in range(0, interval+1):
    a = r_wc + (s_wc-r_wc)/(interval)*i
    soil_moisture_list.append(a)

def soil_m_to_ss(soil_m, r_wc, s_wc, n, m, alpha):
    matric_suction_kpa = []
    for i in soil_m:
        a = i
        if a <= r_wc:
            a = r_wc + 0.001
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e = (d-1)**(1/n)
        f = e/alpha
        matric_suction_kpa.append(f)
    return matric_suction_kpa
    
for j in alpha_1:
    ms_daily_alpha = soil_m_to_ss(soil_moisture_list, r_wc, s_wc, n, m, j)
    suction.append(ms_daily_alpha)

fig, ax = plt.subplots()

plot_0_05 = ax.plot(suction[0], soil_moisture_list,
                      linestyle = 'solid',color = 'black',
                   label = 'alpha = 0.005')

plot_0_1 = ax.plot(suction[1], soil_moisture_list,
                      linestyle = 'dashed',color = 'black',
                   label = 'alpha = 0.05')

plot_0_2 = ax.plot(suction[2], soil_moisture_list,
                      linestyle = 'dotted',color = 'black',
                   label = 'alpha = 0.1')

plot_0_4 = ax.plot(suction[3], soil_moisture_list,
                      linestyle = 'dashdot',color = 'black',
                   label = 'alpha = 0.2')

plot_0_4 = ax.plot(suction[4], soil_moisture_list,
                      linestyle = (0,(1,10)),color = 'black',
                   label = 'alpha = 0.4')

ax.set_xscale('log')
plt.xlabel('Matric Suction (kPa)')
plt.ylabel('Soil Moisture (%)')
plt.title('Soil Moisture against Matric Suction')
leg = plt.legend(loc = 1)
fig.set_size_inches((20, 10))
```

<a name="Drying-to-wetting-reduce-percentage-testing"></a>
### 5.2 Testing significance of drying to wetting reduce percentage
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.2 Additional testing drying to wetting percentage.ipynb

The same processed was done on the drying to wetting reduced percentage value to the alpha value to enable creating the graph that shows
the significance of the parameter with its entire range of values.
- `reduce_percentage` =  the drying to wetting reduce percentage, it was set to vary between 0 to 1 with steps of 0.05 
as shown in `reduce_percentage = [i for i in np.arange(0, 1, 0.05)]`, all percentage is represented as decimals
```python
df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') #extract soil moisture csv data to variable

#Extract observed datetime, date and time
df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')

#Set observed datetime data to list 
datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]

#Set observed soil moisture data to list 
soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]




def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily = hourly_to_daily_data(soil_moisture_list_2000_2021)
date_daily = hourly_to_daily_date(datetime_list_2000_2021)

soil_m_percentage = []
for i in range(len(soil_m_daily)):
    a = soil_m_daily[i]/500
    b = a*100
    soil_m_percentage.append(b)



e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75 
alpha = 0.1



def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d

fail_plane_d = finding_failure_depth(soil_moisture_list_2000_2021, porosity, 500)
fail_plane_d_daily = finding_failure_depth(soil_m_daily, porosity, 500)

def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = [i for i in np.arange(0, 1, 0.05)]



def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

soil_moisture_list = []
FoS_lp = []
FoS_hp = []

for i in range(len(reduce_percentage)):
    ms_daily = soil_m_to_wetting_drying_suction(
        soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage[i], 500)
    soil_moisture_list.append(ms_daily)

for j in range(len(soil_moisture_list)):
    FoS_alpha = FoS_cal(e_cohesion, soil_moisture_list[j], rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang)
    FoS_lp.append(FoS_alpha[low_point])
    FoS_hp.append(FoS_alpha[high_point])

a = []
for i in range(len(FoS_lp)):
    b = tuple([reduce_percentage[i], FoS_lp[i]])
    a.append(b)
    
c = []
for i in range(len(FoS_hp)):
    d = tuple([reduce_percentage[i], FoS_hp[i]])
    c.append(d)

fig, ax = plt.subplots()

feb_obs_plot = ax.plot(reduce_percentage, FoS_lp,
                      color = 'b' )


# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale 
plt.xlabel('Drying to Wetting Reduced Percentage')
plt.ylabel('FoS')
plt.title('FoS against variation of Drying to Wetting Reduced Percentage')
fig.set_size_inches((20, 10))

```

<a name="Testing-significance-of-Field_capacity"></a>
### 5.3 Testing significance of Field Capacity
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.3 Additional testing FC.ipynb

As Field capacity is a hydrological model parameter, to investigate on its significnace to the overall system,
numerous soil moisture modelling is required which result in multiple soil moisture datafile
to be read by the codes as shown below. Similar process as investigation on other geotechnical parameters was used,
however as there are multiple datasets, the data pre-processing,
soil suction modelling and FoS modelling are repeated in multiple for loops.

The realistic range of Field capacity was identified as between 0.5 to 0.95, 
the testing was done with steps of 0.05 as shown in the imported csv datasets.
```python
df_FC_0_5 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.5.csv') #extract soil moisture csv data to variable
df_FC_0_55 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.55.csv') #extract soil moisture csv data to variable
df_FC_0_6 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.6.csv') #extract soil moisture csv data to variable
df_FC_0_65 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.65.csv') #extract soil moisture csv data to variable
df_FC_0_7 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.7.csv') #extract soil moisture csv data to variable
df_FC_0_75 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.75.csv') #extract soil moisture csv data to variable
df_FC_0_8 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.8.csv') #extract soil moisture csv data to variable
df_FC_0_85 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.85.csv') #extract soil moisture csv data to variable
df_FC_0_9 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.9.csv') #extract soil moisture csv data to variable
df_FC_0_95 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_FC_value_0.95.csv') #extract soil moisture csv data to variable


#Set observed soil moisture data to list 
soil_moisture_0_5 = df_FC_0_5['soil.m']
soil_moisture_list_0_5 = soil_moisture_0_5.tolist()
sm_list_2000_2021_FC_0_5 = soil_moisture_list_0_5[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_55 = df_FC_0_55['soil.m']
soil_moisture_list_0_55 = soil_moisture_0_55.tolist()
sm_list_2000_2021_FC_0_55 = soil_moisture_list_0_55[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_6 = df_FC_0_6['soil.m']
soil_moisture_list_0_6 = soil_moisture_0_6.tolist()
sm_list_2000_2021_FC_0_6 = soil_moisture_list_0_6[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_65 = df_FC_0_65['soil.m']
soil_moisture_list_0_65 = soil_moisture_0_65.tolist()
sm_list_2000_2021_FC_0_65 = soil_moisture_list_0_65[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_7 = df_FC_0_7['soil.m']
soil_moisture_list_0_7 = soil_moisture_0_7.tolist()
sm_list_2000_2021_FC_0_7 = soil_moisture_list_0_7[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_75 = df_FC_0_75['soil.m']
soil_moisture_list_0_75 = soil_moisture_0_75.tolist()
sm_list_2000_2021_FC_0_75 = soil_moisture_list_0_75[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_8 = df_FC_0_8['soil.m']
soil_moisture_list_0_8 = soil_moisture_0_8.tolist()
sm_list_2000_2021_FC_0_8 = soil_moisture_list_0_8[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_85 = df_FC_0_85['soil.m']
soil_moisture_list_0_85 = soil_moisture_0_85.tolist()
sm_list_2000_2021_FC_0_85 = soil_moisture_list_0_85[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_9 = df_FC_0_9['soil.m']
soil_moisture_list_0_9 = soil_moisture_0_9.tolist()
sm_list_2000_2021_FC_0_9 = soil_moisture_list_0_9[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_0_95 = df_FC_0_95['soil.m']
soil_moisture_list_0_95 = soil_moisture_0_95.tolist()
sm_list_2000_2021_FC_0_95 = soil_moisture_list_0_95[359393:552257]


len(soil_moisture_list_0_8)


def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt


soil_m_daily_FC_0_5 = hourly_to_daily_data(sm_list_2000_2021_FC_0_5)
soil_m_daily_FC_0_55 = hourly_to_daily_data(sm_list_2000_2021_FC_0_55)
soil_m_daily_FC_0_6 = hourly_to_daily_data(sm_list_2000_2021_FC_0_6)
soil_m_daily_FC_0_65 = hourly_to_daily_data(sm_list_2000_2021_FC_0_65)
soil_m_daily_FC_0_7 = hourly_to_daily_data(sm_list_2000_2021_FC_0_7)
soil_m_daily_FC_0_75 = hourly_to_daily_data(sm_list_2000_2021_FC_0_75)
soil_m_daily_FC_0_8 = hourly_to_daily_data(sm_list_2000_2021_FC_0_8)
soil_m_daily_FC_0_85 = hourly_to_daily_data(sm_list_2000_2021_FC_0_85)
soil_m_daily_FC_0_9 = hourly_to_daily_data(sm_list_2000_2021_FC_0_9)
soil_m_daily_FC_0_95 = hourly_to_daily_data(sm_list_2000_2021_FC_0_95)

e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75 
alpha = 0.1

def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d
 
fail_plane_d_FC_0_5 = finding_failure_depth(soil_m_daily_FC_0_5, porosity, 500)
fail_plane_d_FC_0_55 = finding_failure_depth(soil_m_daily_FC_0_55, porosity, 500)
fail_plane_d_FC_0_6 = finding_failure_depth(soil_m_daily_FC_0_6, porosity, 500)
fail_plane_d_FC_0_65 = finding_failure_depth(soil_m_daily_FC_0_65, porosity, 500)
fail_plane_d_FC_0_7 = finding_failure_depth(soil_m_daily_FC_0_7, porosity, 500)
fail_plane_d_FC_0_75 = finding_failure_depth(soil_m_daily_FC_0_75, porosity, 500)
fail_plane_d_FC_0_8 = finding_failure_depth(soil_m_daily_FC_0_8, porosity, 500)
fail_plane_d_FC_0_85 = finding_failure_depth(soil_m_daily_FC_0_85, porosity, 500)
fail_plane_d_FC_0_9 = finding_failure_depth(soil_m_daily_FC_0_9, porosity, 500)
fail_plane_d_FC_0_95 = finding_failure_depth(soil_m_daily_FC_0_95, porosity, 500)
 
def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = 0.5
 
ms_daily_FC_0_5 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_5, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_55 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_55, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_6 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_6, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_65 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_65, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_7 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_7, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_75 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_75, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_8 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_8, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_85 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_85, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_9 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_9, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
ms_daily_FC_0_95 = soil_m_to_wetting_drying_suction(
    soil_m_daily_FC_0_95, r_wc, s_wc, n, alpha, accounted_time,
    reduce_percentage, 500)
 
def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

FoS_test_list = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
FoS_test_data_lp = []
FoS_test_data_hp = []

FoS_FC_0_5 = FoS_cal(e_cohesion, ms_daily_FC_0_5, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_5, 
                    e_friction_ang)
FoS_FC_0_55 = FoS_cal(e_cohesion, ms_daily_FC_0_55, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_55, 
                    e_friction_ang)
FoS_FC_0_6 = FoS_cal(e_cohesion, ms_daily_FC_0_6, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_6, 
                    e_friction_ang)
FoS_FC_0_65 = FoS_cal(e_cohesion, ms_daily_FC_0_65, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_65, 
                    e_friction_ang)
FoS_FC_0_7 = FoS_cal(e_cohesion, ms_daily_FC_0_7, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_7, 
                    e_friction_ang)
FoS_FC_0_75 = FoS_cal(e_cohesion, ms_daily_FC_0_75, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_75, 
                    e_friction_ang)
FoS_FC_0_8 = FoS_cal(e_cohesion, ms_daily_FC_0_8, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_8, 
                    e_friction_ang)
FoS_FC_0_85 = FoS_cal(e_cohesion, ms_daily_FC_0_85, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_85, 
                    e_friction_ang)
FoS_FC_0_9 = FoS_cal(e_cohesion, ms_daily_FC_0_9, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_9, 
                    e_friction_ang)
FoS_FC_0_95 = FoS_cal(e_cohesion, ms_daily_FC_0_95, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_FC_0_95, 
                    e_friction_ang)

FoS_test_data_lp.append(FoS_FC_0_5[low_point])
FoS_test_data_hp.append(FoS_FC_0_5[high_point])
FoS_test_data_lp.append(FoS_FC_0_55[low_point])
FoS_test_data_hp.append(FoS_FC_0_55[high_point])
FoS_test_data_lp.append(FoS_FC_0_6[low_point])
FoS_test_data_hp.append(FoS_FC_0_6[high_point])
FoS_test_data_lp.append(FoS_FC_0_65[low_point])
FoS_test_data_hp.append(FoS_FC_0_65[high_point])
FoS_test_data_lp.append(FoS_FC_0_7[low_point])
FoS_test_data_hp.append(FoS_FC_0_7[high_point])
FoS_test_data_lp.append(FoS_FC_0_75[low_point])
FoS_test_data_hp.append(FoS_FC_0_75[high_point])
FoS_test_data_lp.append(FoS_FC_0_8[low_point])
FoS_test_data_hp.append(FoS_FC_0_8[high_point])
FoS_test_data_lp.append(FoS_FC_0_85[low_point])
FoS_test_data_hp.append(FoS_FC_0_85[high_point])
FoS_test_data_lp.append(FoS_FC_0_9[low_point])
FoS_test_data_hp.append(FoS_FC_0_9[high_point])
FoS_test_data_lp.append(FoS_FC_0_95[low_point])
FoS_test_data_hp.append(FoS_FC_0_95[high_point])

fig, ax = plt.subplots()

feb_obs_plot = ax.plot(FoS_test_list, FoS_test_data_lp,
                      color = 'b' )


# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale 
plt.xlabel('FC')
plt.ylabel('FoS')
plt.title('FoS against FC variation')
fig.set_size_inches((20, 10))
```

<a name="Testing-significance-of-Effective-friction-angle"></a>
### 5.4 Testing significance of Effective friction angle
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.4 Additional testing friction angle.ipynb

Same approach was taken for the investigation on effective friction angle to
other geotechnical model parameters such as the alpha value.

The realistic range of the effective friction angle was identified to be betweeen 20-40 degrees,
the testing was done in steps of 2 degrees as shown in `e_friction_ang = [i for i in np.arange(20,42,2)] `
```python
df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') #extract soil moisture csv data to variable

#Extract observed datetime, date and time
df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')

#Set observed datetime data to list 
datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]

#Set observed soil moisture data to list 
soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]




def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily = hourly_to_daily_data(soil_moisture_list_2000_2021)
date_daily = hourly_to_daily_date(datetime_list_2000_2021)

soil_m_percentage = []
for i in range(len(soil_m_daily)):
    a = soil_m_daily[i]/500
    b = a*100
    soil_m_percentage.append(b)



e_cohesion = 0 
e_friction_ang = [i for i in np.arange(20,42,2)] 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75 
alpha = 0.1

def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d

fail_plane_d = finding_failure_depth(soil_moisture_list_2000_2021, porosity, 500)
fail_plane_d_daily = finding_failure_depth(soil_m_daily, porosity, 500)

def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = 0.5

def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

FoS_lp = []
FoS_hp = []


ms_daily = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 500)

for i in range(len(e_friction_ang)):
    FoS_alpha = FoS_cal(e_cohesion, ms_daily, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_daily, 
                    e_friction_ang[i])
    FoS_lp.append(FoS_alpha[low_point])
    FoS_hp.append(FoS_alpha[high_point])

a = []
for i in range(len(FoS_lp)):
    b = tuple([e_friction_ang[i], FoS_lp[i]])
    a.append(b)
    
c = []
for i in range(len(FoS_hp)):
    d = tuple([e_friction_ang[i], FoS_hp[i]])
    c.append(d)

fig, ax = plt.subplots()

feb_obs_plot = ax.plot(e_friction_ang, FoS_lp,
                      color = 'b' )


# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale 
plt.xlabel('Effective Friction Angle')
plt.ylabel('FoS')
plt.title('FoS against Effective Friction Angle variation')
fig.set_size_inches((20, 10))
```

<a name="Testing-significance-of-slope-angle"></a>
### 5.5 Testing significance of Slope angle
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.5 Additional testing slope angle.ipynb

Same approach was taken for the investigation on slope angle to
other geotechnical model parameters such as the effective friction angle.

The slope angles that were identified to be the possible range was between 10-85 degrees, 
the testing was done in steps of 5 degrees, as shown in `[slope_angle = i for i in np.arange(10,85, 5)]`
```python
df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') #extract soil moisture csv data to variable

#Extract observed datetime, date and time
df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')

#Set observed datetime data to list 
datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]

#Set observed soil moisture data to list 
soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]


def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily = hourly_to_daily_data(soil_moisture_list_2000_2021)
date_daily = hourly_to_daily_date(datetime_list_2000_2021)

soil_m_percentage = []
for i in range(len(soil_m_daily)):
    a = soil_m_daily[i]/500
    b = a*100
    soil_m_percentage.append(b)

e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = [i for i in np.arange(10,85, 5)]  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75 
alpha = 0.1

def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d

fail_plane_d = finding_failure_depth(soil_moisture_list_2000_2021, porosity, 500)
fail_plane_d_daily = finding_failure_depth(soil_m_daily, porosity, 500)

def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = 0.5

def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

FoS_lp = []
FoS_hp = []


ms_daily = soil_m_to_wetting_drying_suction(
    soil_m_daily, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 500)

for i in range(len(slope_angle)):
    FoS_alpha = FoS_cal(e_cohesion, ms_daily, rate_S_to_suction, 
                    soil_weight, slope_angle[i], fail_plane_d_daily, 
                    e_friction_ang)
    FoS_lp.append(FoS_alpha[low_point])
    FoS_hp.append(FoS_alpha[high_point])

a = []
for i in range(len(FoS_lp)):
    b = tuple([slope_angle[i], FoS_lp[i]])
    a.append(b)
    
c = []
for i in range(len(FoS_hp)):
    d = tuple([slope_angle[i], FoS_hp[i]])
    c.append(d)



fig, ax = plt.subplots()

feb_obs_plot = ax.plot(slope_angle, FoS_lp,
                      color = 'b' )


# Setting x and y axis limit and labels, titles and graph size
#ax.set_yscale('log') #y-axis log scale 
plt.xlabel('Slope Angle')
plt.ylabel('FoS')
plt.title('FoS against Slope Angle variation')
fig.set_size_inches((20, 10))
```

<a name="Testing-significance-of-smax_value"></a>
### 5.6 Testing significance of smax value
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/5.6 Additional testing smax.ipynb

Similar to the testing done for Field capacity, as smax is also a hydrological model parameter,
multiple data files are imported due to the parameter affecting the soil moisture model. 
The realistic range of smax identified was between 100 to 1000,
the testing was done with steps of 50 as shown in the imported csv files.

```python
df_smax_100 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_100.csv') #extract soil moisture csv data to variable
df_smax_150 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_150.csv') #extract soil moisture csv data to variable
df_smax_200 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_200.csv') #extract soil moisture csv data to variable
df_smax_250 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_250.csv') #extract soil moisture csv data to variable
df_smax_300 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_300.csv') #extract soil moisture csv data to variable
df_smax_350 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_350.csv') #extract soil moisture csv data to variable
df_smax_400 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_400.csv') #extract soil moisture csv data to variable
df_smax_450 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_450.csv') #extract soil moisture csv data to variable
df_smax_500 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') #extract soil moisture csv data to variable
df_smax_550 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_550.csv') #extract soil moisture csv data to variable
df_smax_600 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_600.csv') #extract soil moisture csv data to variable
df_smax_650 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_650.csv') #extract soil moisture csv data to variable
df_smax_700 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_700.csv') #extract soil moisture csv data to variable
df_smax_750 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_750.csv') #extract soil moisture csv data to variable
df_smax_800 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_800.csv') #extract soil moisture csv data to variable
df_smax_850 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_850.csv') #extract soil moisture csv data to variable
df_smax_900 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_900.csv') #extract soil moisture csv data to variable
df_smax_950 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_950.csv') #extract soil moisture csv data to variable
df_smax_1000 = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_1000.csv') #extract soil moisture csv data to variable

#Set observed soil moisture data to list 
soil_moisture_100 = df_smax_100['soil.m']
soil_moisture_list_100 = soil_moisture_100.tolist()
sm_list_2000_2021_smax_100 = soil_moisture_list_100[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_150 = df_smax_150['soil.m']
soil_moisture_list_150 = soil_moisture_150.tolist()
sm_list_2000_2021_smax_150 = soil_moisture_list_150[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_200 = df_smax_200['soil.m']
soil_moisture_list_200 = soil_moisture_200.tolist()
sm_list_2000_2021_smax_200 = soil_moisture_list_200[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_250 = df_smax_250['soil.m']
soil_moisture_list_250 = soil_moisture_250.tolist()
sm_list_2000_2021_smax_250 = soil_moisture_list_250[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_300 = df_smax_300['soil.m']
soil_moisture_list_300 = soil_moisture_300.tolist()
sm_list_2000_2021_smax_300 = soil_moisture_list_300[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_350 = df_smax_350['soil.m']
soil_moisture_list_350 = soil_moisture_350.tolist()
sm_list_2000_2021_smax_350 = soil_moisture_list_350[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_400 = df_smax_400['soil.m']
soil_moisture_list_400 = soil_moisture_400.tolist()
sm_list_2000_2021_smax_400 = soil_moisture_list_400[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_450 = df_smax_450['soil.m']
soil_moisture_list_450 = soil_moisture_450.tolist()
sm_list_2000_2021_smax_450 = soil_moisture_list_450[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_500 = df_smax_500['soil.m']
soil_moisture_list_500 = soil_moisture_500.tolist()
sm_list_2000_2021_smax_500 = soil_moisture_list_500[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_550 = df_smax_550['soil.m']
soil_moisture_list_550 = soil_moisture_550.tolist()
sm_list_2000_2021_smax_550 = soil_moisture_list_550[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_600 = df_smax_600['soil.m']
soil_moisture_list_600 = soil_moisture_600.tolist()
sm_list_2000_2021_smax_600 = soil_moisture_list_600[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_650 = df_smax_650['soil.m']
soil_moisture_list_650 = soil_moisture_650.tolist()
sm_list_2000_2021_smax_650 = soil_moisture_list_650[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_700 = df_smax_700['soil.m']
soil_moisture_list_700 = soil_moisture_700.tolist()
sm_list_2000_2021_smax_700 = soil_moisture_list_700[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_750 = df_smax_750['soil.m']
soil_moisture_list_750 = soil_moisture_750.tolist()
sm_list_2000_2021_smax_750 = soil_moisture_list_750[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_800 = df_smax_800['soil.m']
soil_moisture_list_800 = soil_moisture_800.tolist()
sm_list_2000_2021_smax_800 = soil_moisture_list_800[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_850 = df_smax_850['soil.m']
soil_moisture_list_850 = soil_moisture_850.tolist()
sm_list_2000_2021_smax_850 = soil_moisture_list_850[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_900 = df_smax_900['soil.m']
soil_moisture_list_900 = soil_moisture_900.tolist()
sm_list_2000_2021_smax_900 = soil_moisture_list_900[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_950 = df_smax_950['soil.m']
soil_moisture_list_950 = soil_moisture_950.tolist()
sm_list_2000_2021_smax_950 = soil_moisture_list_950[359393:552257]

#Set observed soil moisture data to list 
soil_moisture_1000 = df_smax_1000['soil.m']
soil_moisture_list_1000 = soil_moisture_1000.tolist()
sm_list_2000_2021_smax_1000 = soil_moisture_list_1000[359393:552257]

def hourly_to_daily_data(soil_m):
    daily_data = []
    for i in range(0, len(soil_m), 24):
        temp = statistics.mean(soil_m[i:i+24])
        daily_data.append(temp)
    return daily_data

def hourly_to_daily_date(insert_data):
    chunked = []
    dt = []
    for i in range(0, len(insert_data), 24):
        chunked.append(insert_data[i:i + 24])
    for j in chunked:
        dt.append(j[0])
    return dt

soil_m_daily_smax_100 = hourly_to_daily_data(sm_list_2000_2021_smax_100)
soil_m_daily_smax_150 = hourly_to_daily_data(sm_list_2000_2021_smax_150)
soil_m_daily_smax_200 = hourly_to_daily_data(sm_list_2000_2021_smax_200)
soil_m_daily_smax_250 = hourly_to_daily_data(sm_list_2000_2021_smax_250)
soil_m_daily_smax_300 = hourly_to_daily_data(sm_list_2000_2021_smax_300)
soil_m_daily_smax_350 = hourly_to_daily_data(sm_list_2000_2021_smax_350)
soil_m_daily_smax_400 = hourly_to_daily_data(sm_list_2000_2021_smax_400)
soil_m_daily_smax_450 = hourly_to_daily_data(sm_list_2000_2021_smax_450)
soil_m_daily_smax_500 = hourly_to_daily_data(sm_list_2000_2021_smax_500)
soil_m_daily_smax_550 = hourly_to_daily_data(sm_list_2000_2021_smax_550)
soil_m_daily_smax_600 = hourly_to_daily_data(sm_list_2000_2021_smax_600)
soil_m_daily_smax_650 = hourly_to_daily_data(sm_list_2000_2021_smax_650)
soil_m_daily_smax_700 = hourly_to_daily_data(sm_list_2000_2021_smax_700)
soil_m_daily_smax_750 = hourly_to_daily_data(sm_list_2000_2021_smax_750)
soil_m_daily_smax_800 = hourly_to_daily_data(sm_list_2000_2021_smax_800)
soil_m_daily_smax_850 = hourly_to_daily_data(sm_list_2000_2021_smax_850)
soil_m_daily_smax_900 = hourly_to_daily_data(sm_list_2000_2021_smax_900)
soil_m_daily_smax_950 = hourly_to_daily_data(sm_list_2000_2021_smax_950)
soil_m_daily_smax_1000 = hourly_to_daily_data(sm_list_2000_2021_smax_1000)


e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63
porosity = 0.75 
alpha = 0.1 #0.0054

def finding_failure_depth(soil_m, porosity, smax):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/smax)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d

fail_plane_d_smax_100 = finding_failure_depth(soil_m_daily_smax_100, porosity, 100)
fail_plane_d_smax_150 = finding_failure_depth(soil_m_daily_smax_150, porosity, 150)
fail_plane_d_smax_200 = finding_failure_depth(soil_m_daily_smax_200, porosity, 200)
fail_plane_d_smax_250 = finding_failure_depth(soil_m_daily_smax_250, porosity, 250)
fail_plane_d_smax_300 = finding_failure_depth(soil_m_daily_smax_300, porosity, 300)
fail_plane_d_smax_350 = finding_failure_depth(soil_m_daily_smax_350, porosity, 350)
fail_plane_d_smax_400 = finding_failure_depth(soil_m_daily_smax_400, porosity, 400)
fail_plane_d_smax_450 = finding_failure_depth(soil_m_daily_smax_450, porosity, 450)
fail_plane_d_smax_500 = finding_failure_depth(soil_m_daily_smax_500, porosity, 500)
fail_plane_d_smax_550 = finding_failure_depth(soil_m_daily_smax_550, porosity, 550)
fail_plane_d_smax_600 = finding_failure_depth(soil_m_daily_smax_600, porosity, 600)
fail_plane_d_smax_650 = finding_failure_depth(soil_m_daily_smax_650, porosity, 650)
fail_plane_d_smax_700 = finding_failure_depth(soil_m_daily_smax_700, porosity, 700)
fail_plane_d_smax_750 = finding_failure_depth(soil_m_daily_smax_750, porosity, 750)
fail_plane_d_smax_800 = finding_failure_depth(soil_m_daily_smax_800, porosity, 800)
fail_plane_d_smax_850 = finding_failure_depth(soil_m_daily_smax_850, porosity, 850)
fail_plane_d_smax_900 = finding_failure_depth(soil_m_daily_smax_900, porosity, 900)
fail_plane_d_smax_950 = finding_failure_depth(soil_m_daily_smax_950, porosity, 950)
fail_plane_d_smax_1000 = finding_failure_depth(soil_m_daily_smax_1000, porosity, 1000)

def soil_m_to_wetting_drying_suction(
    soil_m, r_wc, s_wc, n, alpha, 
    accounted_time, reduce_percentage, smax):
    
    suction = []
    w_or_d_suction = []
    avg_soil_m = []
    for i in range(len(soil_m)):
        a = (soil_m[i]/smax)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/(1-1/n))
        d = 1/c
        e= (d-1)**(1/n)
        f = (e/alpha)
        suction.append(f)
        if i <= accounted_time-1:
            w_or_d_suction.append(suction[i])  
        else:
            avg = statistics.mean(soil_m[i-accounted_time:i])
            #avg_soil_m.append(avg)
            if soil_m[i] <= avg:
                w_or_d_suction.append(suction[i])
            if soil_m[i] > avg:
                actual_suction = suction[i]*reduce_percentage
                w_or_d_suction.append(actual_suction) 
    return w_or_d_suction


accounted_time = 5
reduce_percentage = 0.5

ms_daily_smax_100 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_100, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 100)
ms_daily_smax_150 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_150, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 150)
ms_daily_smax_200 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_200, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 200)
ms_daily_smax_250 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_250, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 250)
ms_daily_smax_300 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_300, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 300)
ms_daily_smax_350 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_350, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 350)
ms_daily_smax_400 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_400, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 400)
ms_daily_smax_450 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_450, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 450)
ms_daily_smax_500 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_500, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 500)
ms_daily_smax_550 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_550, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 550)
ms_daily_smax_600 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_600, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 600)
ms_daily_smax_650 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_650, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 650)
ms_daily_smax_700 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_700, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 700)
ms_daily_smax_750 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_750, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 750)
ms_daily_smax_800 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_800, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 800)
ms_daily_smax_850 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_850, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 850)
ms_daily_smax_900 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_900, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 900)
ms_daily_smax_950 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_950, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 950)
ms_daily_smax_1000 = soil_m_to_wetting_drying_suction(
    soil_m_daily_smax_1000, r_wc, s_wc, n, alpha, accounted_time, reduce_percentage, 1000)

def FoS_cal(
    e_cohesion, matric_suction, rate_S_to_suction, soil_weight, 
    slope_angle, fail_plane_d, e_friction_ang):
    
    FoS = []
    for i in range(len(matric_suction)):
        a = e_cohesion
        b = matric_suction[i]*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d[i]*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d[i]*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return FoS


low_point = 1120 #2538
high_point = 1345

FoS_test_list = [i for i in range(100, 1050, 50)]
FoS_test_data_lp = []
FoS_test_data_hp = []

FoS_smax_100 = FoS_cal(e_cohesion, ms_daily_smax_100, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_100, 
                    e_friction_ang)
FoS_smax_150 = FoS_cal(e_cohesion, ms_daily_smax_150, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_150, 
                    e_friction_ang)
FoS_smax_200 = FoS_cal(e_cohesion, ms_daily_smax_200, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_200, 
                    e_friction_ang)
FoS_smax_250 = FoS_cal(e_cohesion, ms_daily_smax_250, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_250, 
                    e_friction_ang)
FoS_smax_300 = FoS_cal(e_cohesion, ms_daily_smax_300, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_300, 
                    e_friction_ang)
FoS_smax_350 = FoS_cal(e_cohesion, ms_daily_smax_350, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_350, 
                    e_friction_ang)
FoS_smax_400 = FoS_cal(e_cohesion, ms_daily_smax_400, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_400, 
                    e_friction_ang)
FoS_smax_450 = FoS_cal(e_cohesion, ms_daily_smax_450, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_450, 
                    e_friction_ang)
FoS_smax_500 = FoS_cal(e_cohesion, ms_daily_smax_500, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_500, 
                    e_friction_ang)
FoS_smax_550 = FoS_cal(e_cohesion, ms_daily_smax_550, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_550, 
                    e_friction_ang)
FoS_smax_600 = FoS_cal(e_cohesion, ms_daily_smax_600, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_600, 
                    e_friction_ang)
FoS_smax_650 = FoS_cal(e_cohesion, ms_daily_smax_650, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_650, 
                    e_friction_ang)
FoS_smax_700 = FoS_cal(e_cohesion, ms_daily_smax_700, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_700, 
                    e_friction_ang)
FoS_smax_750 = FoS_cal(e_cohesion, ms_daily_smax_750, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_750, 
                    e_friction_ang)
FoS_smax_800 = FoS_cal(e_cohesion, ms_daily_smax_800, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_800, 
                    e_friction_ang)
FoS_smax_850 = FoS_cal(e_cohesion, ms_daily_smax_850, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_850, 
                    e_friction_ang)
FoS_smax_900 = FoS_cal(e_cohesion, ms_daily_smax_900, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_900, 
                    e_friction_ang)
FoS_smax_950 = FoS_cal(e_cohesion, ms_daily_smax_950, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_950, 
                    e_friction_ang)
FoS_smax_1000 = FoS_cal(e_cohesion, ms_daily_smax_1000, rate_S_to_suction, 
                    soil_weight, slope_angle, fail_plane_d_smax_1000, 
                    e_friction_ang)

FoS_test_data_lp.append(FoS_smax_100[low_point])
FoS_test_data_hp.append(FoS_smax_100[high_point])
FoS_test_data_lp.append(FoS_smax_150[low_point])
FoS_test_data_hp.append(FoS_smax_150[high_point])
FoS_test_data_lp.append(FoS_smax_200[low_point])
FoS_test_data_hp.append(FoS_smax_200[high_point])
FoS_test_data_lp.append(FoS_smax_250[low_point])
FoS_test_data_hp.append(FoS_smax_250[high_point])
FoS_test_data_lp.append(FoS_smax_300[low_point])
FoS_test_data_hp.append(FoS_smax_300[high_point])
FoS_test_data_lp.append(FoS_smax_350[low_point])
FoS_test_data_hp.append(FoS_smax_350[high_point])
FoS_test_data_lp.append(FoS_smax_400[low_point])
FoS_test_data_hp.append(FoS_smax_400[high_point])
FoS_test_data_lp.append(FoS_smax_450[low_point])
FoS_test_data_hp.append(FoS_smax_450[high_point])
FoS_test_data_lp.append(FoS_smax_500[low_point])
FoS_test_data_hp.append(FoS_smax_500[high_point])
FoS_test_data_lp.append(FoS_smax_550[low_point])
FoS_test_data_hp.append(FoS_smax_550[high_point])
FoS_test_data_lp.append(FoS_smax_600[low_point])
FoS_test_data_hp.append(FoS_smax_600[high_point])
FoS_test_data_lp.append(FoS_smax_650[low_point])
FoS_test_data_hp.append(FoS_smax_650[high_point])
FoS_test_data_lp.append(FoS_smax_700[low_point])
FoS_test_data_hp.append(FoS_smax_700[high_point])
FoS_test_data_lp.append(FoS_smax_750[low_point])
FoS_test_data_hp.append(FoS_smax_750[high_point])
FoS_test_data_lp.append(FoS_smax_800[low_point])
FoS_test_data_hp.append(FoS_smax_800[high_point])
FoS_test_data_lp.append(FoS_smax_850[low_point])
FoS_test_data_hp.append(FoS_smax_850[high_point])
FoS_test_data_lp.append(FoS_smax_900[low_point])
FoS_test_data_hp.append(FoS_smax_900[high_point])
FoS_test_data_lp.append(FoS_smax_950[low_point])
FoS_test_data_hp.append(FoS_smax_950[high_point])
FoS_test_data_lp.append(FoS_smax_1000[low_point])
FoS_test_data_hp.append(FoS_smax_1000[high_point])
```

<a name="Initial-approach-to-integrate-Monte-Carlo-simulation"></a>
## 6. Initial approach to integrate Monte-Carlo simulation
File - Documents/General/2. ERSS ERA5 Reanalyzed Data Based Slope Stability Modelling System/ERSS System/6. Monte-Carlo sim development.ipynb

This section of the codes are the first attempt to develop the system to incorporate monte-carlo simulation,
however as the dependency of the parameters have yet been investigated and lack the needed research to develop the model any further
the development of this model is not completed.

Importing all required package to the script.
```python
import pandas as pd
%matplotlib inline
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import math
from statistics import NormalDist
```

Reading the soil moisture dataset and storing soil moisture data and datetime in lists. 
In order to keep the effeciency of the model, only data between 2000-2021 was processed to mitigate working with unessessary data.
```python
df = pd.read_csv('Rest_and_Be_Thankful_Pass_1959_to_2021_geotechnical_model_run_ERA5_data_smax_500.csv') 

df['Od'] = pd.to_datetime(df['Od'], format = '%Y/%m/%d %H:%M:%S')
df['date'] = pd.to_datetime(df['Od']).dt.date
df['time'] = pd.to_datetime(df['Od']).dt.time

#Set observed datetime data to list 
datetime_df = df['Od']
datetime_l = datetime_df.tolist()
datetime_list_2000_2021 = datetime_l[359393:552257]

#Set observed date data to list 
od = df['date']
od_l = od.tolist()

#Set observed time data to list 
ot = df['time']
ot_l = ot.tolist()

#Set observed soil moisture data to list 
soil_moisture = df['soil.m']
soil_moisture_list = soil_moisture.tolist()
soil_moisture_list_2000_2021 = soil_moisture_list[359393:552257]
```

`dsoil_m_over_dt` is a function to find the first order of difference of soil moisture with time.
```python
def dsoil_m_over_dt(soil_m):
    dt = 1 #1 hour timestep
    dsoil_m = []
    dsoil_m_over_dt = [0]

    for i in range(1, len(soil_m)):
        delta = soil_m[i] - soil_m[i-1]
        dsoil_m.append(delta)

    for i in range(len(dsoil_m)):
        dsoil_m_dt = dsoil_m[i]/dt
        dsoil_m_over_dt.append(dsoil_m_dt)

    return dsoil_m_over_dt

change_in_soil_moisture_with_time = dsoil_m_over_dt(soil_moisture_list)
```

Unit weight and friction angle were the 2 parameters that had know distribution given to general soil and had information regarding the distribution moments.
Further investigation and research will be required for identifying ditributions for other parameters and to validate the assumptions in this model.
Both the unit weight of soil and friction angle is identified as normally distributed and the standard deviation (std) is found by taking the max value as the 97.5% percentile
and the average (avg) is found by taking the mid point between the max and min value.
```python
unit_weight_max = 23.9
unit_weight_min = 20.0

friction_ang_max = 41
friction_ang_min = 31

unit_weight_avg = (unit_weight_max+unit_weight_min)/2
unit_weight_std = (unit_weight_max-unit_weight_avg)/NormalDist(mu=0, sigma=1).inv_cdf(0.975)

friction_ang_avg = (friction_ang_min+friction_ang_max)/2
friction_ang_std = (friction_ang_max-friction_ang_avg)/NormalDist(mu=0, sigma=1).inv_cdf(0.975)
```

The `monte_carlo_sim` function is a combined function that takes the soil moisture data,
the moments of the distribution of all identified uncertained parameters and other certained parameters.
The model then generates random numbers between 0-1 via `np.random.rand()` for each uncertain parameter and match it with its relative parameter value,
the values are then stored in lists to be further process such as `unit_weight` and `friction_angle`.
The `soil_m_to_ss`, `FoS_cal` function is the exact function as given in previous models. As both uncertained parameters does not affect the matric suction calculations, only the FoS is repeated in a for loop.
As shown in `for k in range(len(matric_suction)):
        for l in range(nr):
            FoS = FoS_cal(e_cohesion, matric_suction[k], rate_S_to_suction, unit_weight[l], slope_angle, fail_plane_d, friction_angle[l])
            FoS_list.append(FoS)`
Where k is values in the matric suction list which each value represent the suction value of the soil in a given day and l is the number of simulation that the user set for the data in any given day to run.
The `chunk` function then put the results of each simulation in a day into a sublist such that `[[day1_sim1, day1_sim2,day1_simx], [day2_sim1, day2_sim2, day2_simx],...] is achieved.
```python

def monte_carlo_sim(unit_weight_avg, unit_weight_std, friction_ang_avg, 
                    friction_ang_std, nr, soil_m, r_wc, s_wc, n, m, alpha,
                   e_cohesion, rate_S_to_suction, slope_angle, fail_plane_d):
    
    unit_weight_rand = np.random.rand(nr)
    friction_ang_rand = np.random.rand(nr)
    unit_weight = []
    friction_angle = []
    
    for i in range(len(unit_weight_rand)):
        unit_weight_sim = (NormalDist(mu=0, sigma=1).inv_cdf(unit_weight_rand[i])
                           *unit_weight_std)+unit_weight_avg
        friction_ang_sim = (NormalDist(mu=0, sigma=1).inv_cdf(friction_ang_rand[i])
                           *friction_ang_std)+friction_ang_avg
        unit_weight.append(unit_weight_sim)
        friction_angle.append(friction_ang_sim)
    
    def soil_m_to_ss(soil_m, r_wc, s_wc, n, m, alpha):
        matric_suction_kpa = []
        for i in soil_m:
            a = (i/100)*s_wc
            if a <= r_wc:
                a = r_wc + 0.01
            b = (a-r_wc)/(s_wc-r_wc)
            c = b**(1/m)
            d = 1/c
            e = (d-1)**(1/n)
            f = e/alpha
            matric_suction_kpa.append(f)
        return matric_suction_kpa

    def FoS_cal(e_cohesion, matric_suction, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d, e_friction_ang):
        FoS = []
        a = e_cohesion
        b = matric_suction*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
        return FoS

    def chunk(lst, nr): #splitting long list into multiple of sublist
        new_list = []
        for i in range(0, len(lst), nr):
            new_list.append(lst[i:i + nr])
        return new_list
    
    FoS_list = []
    matric_suction = soil_m_to_ss(soil_moisture_list, r_wc, s_wc, n, m, alpha)
    suction_bracked_list = chunk(matric_suction, nr)
    
    for k in range(len(matric_suction)):
        for l in range(nr):
            FoS = FoS_cal(e_cohesion, matric_suction[k], rate_S_to_suction, unit_weight[l], slope_angle, fail_plane_d, friction_angle[l])
            FoS_list.append(FoS)
            
    FoS_bracked_list = chunk(FoS_list, nr)  
            
    return FoS_bracked_list
```

Setting out the set parameter values for the geotechnical model.
(Most parameters in the list is known to varied even in different section of a landslide event which would contain uncertainties,
however as little research have been done in this area, the parameters are taken as constant)
```python
e_cohesion = 0 
e_friction_ang = 36 
rate_S_to_suction = 36 
soil_weight = 22 
slope_angle = 25  
n = 2.3 
m = 1-1/n 
r_wc = 0.1 
s_wc = 0.63 
porosity = 0.75 
alpha = 0.1
```

Utilising the `finding_failure_depth` function to obtain the failure depth corresponding to the soil moisture 
and obtaining the results from the monte-carlo simulation with `monte_carlo_sim` function. The format to the results given will be in form of a 
list of lists where each sublist contains the simulation result regarding to one data point, i.e FoS of a single day.
- `unit_weight_avg` = average value of the unit weight of the soil
- `unit_weight_std` = standard deviation of the unit weight of the soil
- `friction_ang_avg` = average vale of the effective friction angle
- `friction_ang_std` = standard deviation of the effective friction angle
- `10` = number of simulations for a single data point, ideally each data point should be simulated for at least 1000 times.
```python
def finding_failure_depth(soil_m, porosity):
    fail_plane_d = []
    for i in soil_m:
        critical_depth = (i/100)/porosity
        fail_plane_d.append(critical_depth)
    return fail_plane_d

fail_plane_d = finding_failure_depth(soil_moisture_list, porosity)

sim_1 = monte_carlo_sim(unit_weight_avg, unit_weight_std, friction_ang_avg, 
                        friction_ang_std, 10, soil_moisture_list_2000_2021, r_wc , s_wc, n, m, alpha,
                        e_cohesion, rate_S_to_suction, slope_angle, fail_plane_d)
```

Taking maximum and minimum, 10th, 25th, 50th, 75th, and 90th percentile of each simulation date and store data of each percentile in a list 
```python
def percentile_conv(subset_data):
    percentile_list = []
    for i in subset_data:
        sub_list = []
        percen_0 = np.percentile(i, 0)
        percen_10 = np.percentile(i, 10)
        percen_25 = np.percentile(i, 25)
        percen_50 = np.percentile(i, 50)
        percen_75 = np.percentile(i, 75)
        percen_90 = np.percentile(i, 90)
        percen_100 = np.percentile(i, 100)
        sub_list.append(percen_0)
        sub_list.append(percen_10)
        sub_list.append(percen_25)
        sub_list.append(percen_50)
        sub_list.append(percen_75)
        sub_list.append(percen_90)
        sub_list.append(percen_100)
        percentile_list.append(sub_list)
    return percentile_list

def extract_percen(sublist_data):
    l_0 = []
    l_10 = []
    l_25 = []
    l_50 = []
    l_75 = []
    l_90 = []
    l_100 = []
    for i in sublist_data:
        l_0.append(i[0])
        l_10.append(i[1])
        l_25.append(i[2])
        l_50.append(i[3])
        l_75.append(i[4])
        l_90.append(i[5])
        l_100.append(i[6])
    return l_0, l_10, l_25, l_50, l_75, l_90, l_100

FOS_percentile = percentile_conv(sim_1)

percen_0 = extract_percen(FOS_percentile)[0]
percen_10 = extract_percen(FOS_percentile)[1]
percen_25 = extract_percen(FOS_percentile)[2]
percen_50 = extract_percen(FOS_percentile)[3]
percen_75 = extract_percen(FOS_percentile)[4]
percen_90 = extract_percen(FOS_percentile)[5]
percen_100 = extract_percen(FOS_percentile)[6]
```
