# RSL Risk-Based Shallow Landslide Forecasting System User Manual
## Alan Chou 09/08/2022 alan_chou@ymail.com, ac6522@ic.ac.uk

## Table of Contents

- [Flow diagram of the landslide forecasting system](#Flow-diagram-of-the-landslide-forecasting-system)
- [1. SEAS5 API data extraction and pre-processing](#SEAS5-API-data-extraction-and-pre-processing)
   - [1.1 Import all packages required to the file](#Import-all-packages-required-to-the-file)
   - [1.2 Data extraction from climate data store](#Data-extraction-from-climate-data-store)
   - [1.3 Extract grib file as Excel](#Extract-grib-file-as-Excel)
   - [1.4 SEAS5 Data pre-processing](#SEAS5-Data-pre-processing)
- [2. ERA5 API data extraction and pre-processing](#ERA5-API-data-extraction-and-pre-processing)
   - [2.1 Import all packages required to the file](#Import-all-packages-required-to-the-file)
   - [2.2 ERA5 Data pre-processing](#ERA5-Data-pre-processing)
- [3. Hydrological Processing](#Hydrological-Processing)
   - [3.1 DAYMOD Soil moisture model](#DAYMOD-Soil-moisture-model)
   - [3.2 Reading SEAS5 data in R](#Reading-SEAS5-data-in-R)
   - [3.3 Reading ERA5 data in R](#Reading-ERA5-data-in-R)
   - [3.4 Processing SEAS5 soil moisture data](#Processing-SEAS5-soil-moisture-data)
   - [3.5 Processing ERA5 soil moisture data](#Processing-ERA5-soil-moisture-data)
   - [3.6 Saving SEAS5 soil moisture data as csv file](#Saving-SEAS5-soil-moisture-data-as-csv-file)
   - [3.7 Saving ERA55 soil moisture data as csv file](#Saving-ERA5-soil-moisture-data-as-csv-file)
- [4. Geotechnical Model and Result Presentation](#Geotechnical-Model-and-Result-Presentation)
   - [4.1 Obtaining the Factor of Safety for a typical slope in the region](#Obtaining-the-Factor-of-Safety-for-a-typical-slope-in-the-region)

<a name="Flow-diagram-of-the-landslide-forecasting-system"></a>
## Flow diagram of the landslide forecasting system

[![Flow-Diagram-1.png](https://i.postimg.cc/CxqQxLTh/Flow-Diagram-1.png)](https://postimg.cc/ts95rb4w)

Both the SEAS5 forecasting system and the ERA5 system followed the same process as shown in the flow diagram. 
First, the system requires the user to download the required SEAS5 or ERA5 precipitation and 
potential evaporation data from the Copernicus CDS either via API or the webpage. 
The data will then be pre-process in the python script to match the data with the input needed for the R hydrological model.
The hydrological model will then output the soil moisture data as csv file which allows the python geotechnical model to process.
The geotechnical model can be seperated as 2 parts, 1. the soil suction modelling, 2. the slope stability modelling.
The geotechnical model accounts for the effect of soil moisture in unsaturated soil and investigate the likelihood of 
a shallow landslide with a particular state of the topsoil.

<a name="SEAS5-API-data-extraction-and-pre-processing"></a>
## 1. SEAS5 API data extraction and pre-processing

File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.2 RSL SEAS5 Hydrological Model/Extract SEAS5 Data From API and Data Pre-processing.ipynb


The aim of this process is to extract the required SEAS5 seasonal 
forecasting data from the Climate data store and present the data in an 
easily readable data type in python which will then allow pre-processing
to be performed to match the units and expected input 
for the hydrological model

<a name="Import-all-packages-required-to-the-file"></a>
### 1.1 Import all packages required to the file

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
```

<a name="Data-extraction-from-climate-data-store"></a>
### 1.2 Data extraction from climate data store

Input
- 'originating_centre' = Specifying ECMWF as the 
organisation publishing the data
- 'system' = Specifying the version of the system, 
which in our case its SEAS5
- 'variable' = Parameters required for the anaylsis, 
in this example, its total precipitation and evaporation 
which is taken as potential evaporation 
- 'year', 'month', 'day' = Specifying the forecast issued date
- 'leadtime_hour' = Specifying the forecast hours required, data is given
in 24 hour interval
- 'area' = Specifying the location of interest, 
resolution of the SEAS5 model is 1degrees by 1 degrees
- 'format' = Data from ECMWF on Seasonal forecast is given as grib file
- 'name.grib' = Enter the `name` for the grib file 
```python
c = cdsapi.Client()

c.retrieve(
    'seasonal-original-single-levels',
    {
        'originating_centre': 'ecmwf',
        'system': '5',
        'variable': [
            'evaporation', 'total_precipitation',
        ],
        'year': '2020',
        'month': '02',
        'day': '01',
        'leadtime_hour': [
            '24', '48', '72',
            '96', '120', '144',
            '168', '192', '216',
            '240', '264', '288',
            '312', '336', '360',
            '384', '408', '432',
            '456', '480', '504',
            '528', '552', '576',
            '600', '624', '648',
            '672', '696',
        ],
        'area': [
            51.7, -3.4, 51.6,
            -3.3,
        ],
        'format': 'grib',
    },
    'name.grib')
```
As an alternative, the user can also access the CDS web page to request
data where the system will be able to record the request and download
the data in the background

<a name="Extract-grib-file-as-Excel"></a>
### 1.3 Extract grib file as Excel


First, unzip the grib file to a python readable dataframe

Input
- 'name.grib' = 'name' given to the data file
```python
ds = xr.open_dataset('name.grib', engine='cfgrib',
                      backend_kwargs={'read_keys': ['experimentVersionNumber']})
df = ds.to_dataframe()
print(df)
```

Next, Bring the index into the dataframe
```python
df_index=df.index
test = df_index.to_frame()
```

Next, join index and data

Input
- df = dataframe 'name' given to the dataset
```python
df = test.join(df)
```

Next, drop the index row 
```python
df.reset_index(drop=True, inplace=True)
```

If needed to spectate the data activate the line

Input
- 'name.csv' = 'name' given to the csv data file
```python
df.to_csv("name.csv")
```

<a name="SEAS5-Data-pre-processing"></a>
### 1.4 SEAS5 Data pre-processing 


Lock data with the longitude and latitude of interest, in this example,
the location locked is the tylorstown landslide location 

Input
- -3.4 = longitude of the point of interest
- 51.6 = latitude of the point of interest
```python
df = df.loc[df['longitude'] == -3.4] 
df = df.loc[df['latitude'] == 51.6]
```

Set parameters to series
```python
tp = df['tp']
e = df['e']
valid_time = df['valid_time']
```

Transforming precipitation data series in meters to 
data list in millimeters
```python
def m_to_mm(m_data):
    mm = []
    for j in m_data:
        mm.append(j*1000)
    return mm

tp_mm_list = m_to_mm(tp)
```

Transforming the potential evaporation data series from negative value
in meters for evaporation to in positive values data list in millimeters
```python
def e_to_mm(m_data):
    mm = []
    for j in m_data:
        mm.append(j*-1000)
    return mm

e_mm_list = e_to_mm(e)
```

Converting the corresponding forecast date and time data series to list
```python
valid_time_list = valid_time.tolist()
```

Create sublist within a list, where each sublist contain all data given 
in a single forecast

Input 
- 1st = big list that contains all data of 51 set of forecasts in
the period of interest
- n = length of each forecast, in this example, the daily data taken
is the forecast for February 2020  
```python
def chunk(lst, n):
    new_list = []
    for i in range(0, len(lst), n):
        new_list.append(lst[i:i + n])
    return new_list

forecast_time = chunk(valid_time_list, 29)
tp_mm = chunk(tp_mm_list, 29)
e_mm = chunk(e_mm_list, 29)
```

Converting accumulated precipitation and potential evaporation data to daily data

non_abs = simply subtracting the data from the day of interest to the
previous data to obtain the increase in precipitation or evaporation
abs = same process to non_abs but remove daily data that gives negative
precipitation and potential evaporation, and fill in those gaps with 0 
```python
def sub_data_by_previous_day_non_abs(data): 
    outer_list = [] 
    
    for i in data: 
        inner_list = [i[0]] 
        for j in range(1,len(i)): 
            inner_list.append(i[j]-i[j-1])
        
        outer_list.append(inner_list)
    return outer_list

dp_mm = sub_data_by_previous_day_non_abs(tp_mm)
de_mm = sub_data_by_previous_day_non_abs(e_mm)

def sub_data_by_previous_day_abs(data):
    outer_list = []
    
    for i in data:
        inner_list = [i[0]]
        for j in range(1,len(i)):
            result = i[j]-i[j-1]
            if result<0:
                result = 0
            inner_list.append(result)
        
        outer_list.append(inner_list)
    return outer_list
        
dp_mm_1 = sub_data_by_previous_day_abs(tp_mm)
de_mm_1 = sub_data_by_previous_day_abs(e_mm)
```

Plotting the precipitation over the forecasting period with intervals of 
2 days in x-axis

Input
- plt.xlabel = x-axis title 
- plt.ylabel = y-axis title
- plt.title = title of the graph
- plt.savefig = save graph in the folder, activate when needed 
```python
matplotlib.rcParams.update({'font.size': 10}) 
fig, ax = plt.subplots()
for i in range(51): 
    plt.plot(forecast_time[1],tp_mm[i]) 
ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m'))
plt.xlabel('Forecast Period(Feb 2020)')
plt.ylabel('Total Precipitation (mm)')
plt.title('Total Precipitation over Forecasting Period')
#plt.savefig('total_precipitation_over_forecasting_period.png')
fig.set_size_inches((20, 10))
plt.show()
```

Converting the daily data into hourly data by dividing the data by 24 and
repeated 24 times in the list to represent data for each individual hour
```python
dp_mm_list = list(chain(*dp_mm_1))
de_mm_list = list(chain(*de_mm_1))

def hourly_data(data):
    h_data = []
    for i in data:
        hd = i/24
        for j in range(24):
            h_data.append(hd)
    return h_data

dp_mm_per_hour_list = hourly_data(dp_mm_list)
de_mm_per_hour_list = hourly_data(de_mm_list)
```

Converting the daily forecast date into hourly interval by
adding time measurement to the date
```python
def hourly_forecast_transformation1(data):
    hf_data = []
    time = []
    for i in data:
        for j in range(24):
            hf_data.append(i)
            t = datetime.time( j, 00)
            str_t = str(t)
            time.append(str_t)
    return hf_data, time

valid_time_list_1 = hourly_forecast_transformation1(valid_time_list)[0]
valid_time_list_2 = hourly_forecast_transformation1(valid_time_list)[1]
```

Storing the data and the forecast date and time in a dictionary 
to be called as a dataframe and saved as a csv file for data transfer
to the R hydrological model 
```python
dict = {'Forecast Date': valid_time_list_1,'Forecast Time': valid_time_list_2 , 'Hourly Precipitation': dp_mm_per_hour_list, 'Evaporation Rate': de_mm_per_hour_list_abs}
df = pd.DataFrame(dict)
df.to_csv('feb_r_read_SEAS5_data.csv')
```


<a name="ERA5-API-data-extraction-and-pre-processing"></a>
## 2. ERA5 API data extraction and pre-processing

File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.1 RSL ERA5 Hydrological Model/Extract ERA5 Data From API and Data Pre-processing.ipynb
Runned similar to the data extraction and pre-processing as SEAS5 model

<a name="Import-all-packages-required-to-the-file"></a>
### 2.1 Import all packages required to the file

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
```

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
        'year': '2020',
        'month': '02',
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
            '28', '29',
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
            51.7, -3.4, 51.6,
            -3.3,
        ],
        'format': 'grib',
    },
    'name.grib')
```

As an alternative, the user can also access the CDS web page to request 
data where the system will be able to record the request and download
the data in the background

* Please note that there are a size limit to the downloading of data for 
both the API and the website, when all hourly data are required the maximum size of 
each data file can only be a maximum of 5 years period

<a name="ERA5-Data-pre-processing"></a>
### 2.2 ERA5 Data pre-processing

Similar data pre-processing process as the SEAS5 model. 
Please refer back to the description given in the SEAS5 model 

```python
ds = xr.open_dataset('name.grib', engine='cfgrib',
                      backend_kwargs={'read_keys': ['experimentVersionNumber']})
df = ds.to_dataframe()

df_index=df.index
test = df_index.to_frame()
df = test.join(df)
df.reset_index(drop=True, inplace=True)
#df.to_csv("data.csv")
df = df.loc[df['longitude'] == -3.4]
df = df.loc[df['latitude'] == 51.6]

dp = df['tp']
pev = df['pev']
valid_time = df['valid_time']

def m_to_mm(m_data):
    mm = []
    for j in m_data:
        mm.append(j*1000)
    return mm
dp_mm_series = m_to_mm(dp)

def pev_to_mm(m_data): 
    mm = []
    for j in m_data:
        mm.append(j*-1000)
    return mm
pev_mm_series = pev_to_mm(pev)

dict = {'Observed Date': valid_time, 'Hourly Precipitation': dp_mm_series, 'Potential Evaporation Rate': pev_mm_series}
df = pd.DataFrame(dict)
df.to_csv('Feb_r_read_ERA5_data.csv')
```
<a name="Hydrological-Processing"></a>
## 3. Hydrological Processing

File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.2 RSL SEAS5 Hydrological Model/SEAS5 R script.R
File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.1 RSL ERA5 Hydrological Model/ERA5 R script.R

Following converting the data into csv files, 
the data can then be processed in the hydrological model to model
the change in soil moisture across the period of interest. 
The model was written by Dr Thomas Kjeldsen and was modified to target
on the need of interest of the dissertation and research.

<a name="DAYMOD-Soil-moisture-model"></a>
### 3.1 DAYMOD Soil moisture model

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
}


fun.cross.m<-function(d.t,m0,rrain,konst,SM,zone0,m.limit.s)
{
  # Calculate difference between soil moist. at end of time step d.t and zone level
  m1<-m.t(d.t,m0,rrain,konst,SM,zone0)
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

Reading the data from the csv files and assign variables to hold the data


<a name="Reading-SEAS5-data-in-R"></a>
### 3.2 Reading SEAS5 data in R

File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.2 RSL SEAS5 Hydrological Model/SEAS5 R script.R

```R
x <- read.csv("feb_r_read_SEAS5_data.csv", TRUE, ",")


p<-x$Hourly.Precipitation
Ep<-x$Evaporation.Rate
Fd<-x$Forecast.Date
Ft <-x$Forecast.Time 
d.t<-1 
np<-length(p)
```

<a name="Reading-ERA5-data-in-R"></a>
### 3.3 Reading ERA5 data in R

File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.1 RSL ERA5 Hydrological Model/ERA5 R script.R

```R
y <- read.csv("dec_r_read_ERA5_data.csv", TRUE, ",")


p<-y$Hourly.Precipitation #Extracting Hourly Precipitation data  
Ep<-y$Potential.Evaporation.Rate #Extracting Hourly Evaporation data
Od<-y$Observed.Date #Extracting valid datetime
d.t<-1 # Time step, one unit of time step of ERA5 data (hours)
np<-length(p) #Length of the data
```

<a name="Processing-SEAS5-soil-moisture-data"></a>
### 3.4 Processing SEAS5 soil moisture data 

The SEAS5 seasonal forecasts provides 51 sets of forecast for any period
of interest and therefore the soil moisture model will be required to
be repeated for 51 times

Input
- `696` = length of the data in each forecasting period i.e 24x29=696
- `695` = length of the data in each forecasting period - 1
```R
for (i in seq(1,length(p), 696))
{soil.b4<-m
  for (j in seq(i, i + 695))
  {
  out<-pdsoil(smaxi,FC,RC,dkt,soil.b4,p[j],Ep[j],d.t)
  soil.b4<-out[1]
  soil.m[j]<-out[1]
  qs[j]<-out[2]
  evapa[j]<-out[3] 
  rch[j]<-out[4] 
  }
}
```
<a name="Processing-ERA5-soil-moisture-data"></a>
### 3.5 Processing ERA5 soil moisture data 

As ERA5 data only contain a set of reanlyzed data, therefore only 1 loop is required.

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


Converting the lists of data back into a dataframe to be contained in a csv for the geotechnical model to process

<a name="Saving-SEAS5-soil-moisture-data-as-csv-file"></a>
### 3.6 Saving SEAS5 soil moisture data as csv file

Input
- `data.frame(Fd, Ft,soil.m)` = all varaiables that is of interest
of the geotechnical model, if required other variables such as 
surface runoff and etc can be obtained
```R
df<- data.frame(Fd, Ft,soil.m)
write.csv(df,"C:/Users/Alan/dec_geotechnical_model_run_SEAS5_data.csv", row.names = FALSE)
```
<a name="Saving-ERA5-soil-moisture-data-as-csv-file"></a>
### 3.7 Saving ERA5 soil moisture data as csv file

Input
- `data.frame(Od,soil.m)` = all varaiables that is of interest
of the geotechnical model, if required other variables such as 
surface runoff and etc can be obtained
```R
df<- data.frame(Od,soil.m)
write.csv(df,"C:/Users/Alan/dec_geotechnical_model_run_ERA5_data.csv", row.names = FALSE)
```

<a name="Geotechnical-Model-and-Result-Presentation"></a>
## 4. Geotechnical Model and Result Presentation

File - Documents/General/1. RSL Risk-based Shallow Landslide Forecasting System/1.3 RSL Geotechnical Model/Geotechnical_model_for_SEAS5 and ERA5.ipynb

The Factor of Safety of a typical infinite slope can be found combining
the soil suction modelling method developed by Van Genuchten and the 
infinite slope stability modelling method developed by
Fredlund and Rahardjo.


<a name="Obtaining-the-Factor-of-Safety-for-a-typical-slope-in-the-region"></a>
### 4.1 Obtaining the Factor of Safety for a typical slope in the region

Import all required packages into the geotechnical model
```python
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import csv
import matplotlib.dates as mdates
import pyreadr
from itertools import chain
import datetime
import math
from plotnine import *
from matplotlib.dates import DateFormatter
```

Read the soil moisture data in the csv file which was generated by
the soil moisture model
```python
df_f = pd.read_csv('feb_geotechnical_model_run_SEAS5_data_exp.csv')
df_o = pd.read_csv('feb_geotechnical_model_run_ERA5_data.csv')
```

Set SEAS5 forecasting data into variables as lists
- `fd` = date of the forecast data
- `ft` = time of the forecast data
- `date_time` = date and time of the forecast data
- `forec_soil_moisture` = forecast soil moisture data
```python

fd = df_f['Fd']
fd_l = fd.tolist()

ft = df_f['Ft']
ft_l = ft.tolist()

date_time = pd.to_datetime(df_f['Fd'],dayfirst = True) + pd.TimedeltaIndex(df_f['Ft'])
date_time_l = date_time.tolist()

forec_soil_moisture = df_f['soil.m']
forec_soil_moisture_list = forec_soil_moisture.tolist()
```

Set ERA5 reanalyzed data into variables as lists
- `obs_datetime_l` = date and time of the observed data
- `od_l` = date of the observed data
- `ot_l` = time of the observed data
- `obs_soil_moisture_list` = observed soil moisture data
```python

df_o['Od'] = pd.to_datetime(df_o['Od'], format = '%d/%m/%Y %H:%M')
df_o['date'] = pd.to_datetime(df_o['Od']).dt.date
df_o['time'] = pd.to_datetime(df_o['Od']).dt.time

df_o_locked = df_o

#Set observed datetime data to list 
obs_datetime_df = df_o['Od']
obs_datetime_l = obs_datetime_df.tolist()

#Set observed date data to list 
od = df_o['date']
od_l = od.tolist()

#Set observed time data to list 
ot = df_o['time']
ot_l = ot.tolist()

#Set observed soil moisture data to list 
obs_soil_moisture = df_o['soil.m']
obs_soil_moisture_list = obs_soil_moisture.tolist()
```

Set base environment of the geotechnical model

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
- `alpha` = Van Genuchten parameter alpha 
(Estimated for glacial till materials)
- `r_wc` = Residual water content (Estimated for glacial till materials)
- `s_wc` = Saturated water content (Estimated for glacial till materials)
- `porosity` = porosity of the soil (Estimated for glacial till materials)
```python
e_cohesion = 5 
e_friction_ang = 28 
rate_S_to_suction = 20 
soil_weight = 18 
slope_angle = 50 
fail_plane_d = 1 
n = 1.25 
m = 1-1/n 
alpha = 1/0.654 
r_wc = 0.098 
s_wc = 0.459
porosity = 0.385
```

Determine the soil suction value from the soil moisture data and
store the results in variables 'forec_matric_suction' for the forecast
result and 'obs_matric_suction' for observed reanalysis result

Input
- `soil_m` = forecast or observed hourly soil moisture data list
- `r_wc` = residual water content
- `s_wc` = Saturated water content
- `n` = Van Genuchten parameter n
- `m` = Van Genuchten parameter m
- `alpha` = Van Genuchten parameter alpha
```python
def soil_m_to_ss(soil_m, r_wc, s_wc, n, m, alpha):
    matric_suction_kpa = []
    for i in soil_m:
        a = (i/100)*s_wc
        b = (a-r_wc)/(s_wc-r_wc)
        c = b**(1/m)
        d = 1/c
        e= (d-1)**(1/n)
        f = e/alpha
        matric_suction_kpa.append(f)
    return matric_suction_kpa

forec_matric_suction = soil_m_to_ss(forec_soil_moisture_list, r_wc, s_wc, n, m, alpha)
obs_matric_suction = soil_m_to_ss(obs_soil_moisture_list, r_wc, s_wc, n, m, alpha)
```


Determine the FoS of the slope across the forecast period and store 
results in variables 'forec_FoS' for FoS corresponding to the forecast
data and 'obs_FoS' for FoS corresponding to the observed reanalyzed data

Input
- `e_cohesion` = Effective Cohesion in kPa
- `matric_suction` = forecast or observed soil suction list
- `rate_S_to_suction` = Degrees Angle indicating the rate of increase of shear strength relative to increase matric suction
- `soil_weight` = Soil Unit Weight in kN/m3
- `slope_angle` = Degrees Slope angle
- `fail_plane_d` = Depth of failure plane in meters
- `e_friction_ang` = Degrees Effective angle of internal friction
```python
def FoS_cal(e_cohesion, matric_suction, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d, e_friction_ang):
    FoS = []
    for i in matric_suction:
        a = e_cohesion
        b = i*math.tan(math.radians(rate_S_to_suction))
        c = soil_weight*fail_plane_d*(math.cos(math.radians(slope_angle))**2)*math.tan(math.radians(e_friction_ang))
        d = soil_weight*fail_plane_d*math.sin(math.radians(slope_angle))*math.cos(math.radians(slope_angle))
        final = (a+b+c)/d
        FoS.append(final)
    return 

forec_FoS = FoS_cal(e_cohesion, forec_matric_suction, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d, e_friction_ang)
obs_FoS = FoS_cal(e_cohesion, obs_matric_suction, rate_S_to_suction, soil_weight, slope_angle, fail_plane_d, e_friction_ang)
```

Separate the 51 forecasts and determine the 10th, 25th, median, 75th,
and 90th percentile in addition to the min and max values for each date
and time set in the forecast. 
The chunk function work as a separator which create sub-list within the
list of the forecast i.e 
[[fore_1_day1, fore_2_day1,fore_3_day1], [fore_1_day2, fore_2_day2,fore_3_day2]]
The data_set function work as a sampler which takes the data from the
same column in the matric, allowing comparison between forecasts
The data_conv function combines the 2 functions and therefore when
processing the results only 1 function is called
The percentile_conv function determine and generate output
of the interested percentiles which can then be saved in variables

Input
- `n` = length of a single forecast or observation i.e 24x28=672
```python

def chunk(lst, n): 
    new_list = []
    for i in range(0, len(lst), n):
        new_list.append(lst[i:i + n])
    return new_list

def data_set(insert_data, n):
    outer_list = []
    for i in range(n):
        inner_list = []
        for j in insert_data:
             inner_list.append(j[i])
        outer_list.append(inner_list)
    return outer_list

def data_conv(insert_data, n):
    new_data = chunk(insert_data, n)
    return data_set(new_data, n)

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

date_forec_FOS = data_conv(forec_FoS, 672)
FOS_percentile = percentile_conv(date_forec_FOS)

date_soil_moisture = data_conv(forec_soil_moisture_list, 672)
soil_moisture_percentile = percentile_conv(date_soil_moisture)
```
As forecast data on a specific date is the accumulative data from
the past date therefore all datetime will require to bring 1 day forward

Input
- `insert_data` = datetime data
- `n` = length of a single forecast or observation i.e 24x28=672
```python
def plot_date_conv(insert_data, n):
    date = []
    for i in insert_data:
        i += datetime.timedelta(days = -1)
        date.append(i)
    new_data = chunk(date, n)
    dt = []
    for j in new_data[0]:
        dt.append(j)
    return dt

date_time_plot = plot_date_conv(date_time_l, 672)
```

Line plot of the FoS of all 51 set of forecasts against
the ERA5 reanalyzed data 

Input
- `label` = label of each individual line
- `color` = color code of the line
- `xlim`, `ylim` = x and y axis upper and lower limit
- `mdates.DayLocator(interval=1)` = x axis date interval
- `DateFormatter` = datetime presented in the graph
- `xlabel`, `ylabel` = axis title
- `title` = graph title
- `plt.legend(loc = 1)` = designed legend location
```python
fig, ax = plt.subplots()
Feb_forec_plot = ax.plot(date_time_plot, chunk(forec_FoS,672)[0],
                        color = 'b',
                        label = 'February 2020 51 forecasts')

Feb_forec_plot = ax.plot(date_time_plot, date_forec_FOS,
                        color = 'b',)

feb_obs_plot = ax.plot(obs_datetime_l, obs_FoS,
                      color = 'r',
                      label = 'February 2020 observed' )

ax.set(xlim=(0, 29),
       ylim=(1, 1.5))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
date_form = DateFormatter("%d-%m")
ax.xaxis.set_major_formatter(date_form)
ax.set(xlim = [obs_datetime_l[24], obs_datetime_l[671]])
plt.xlabel('Forecast Period(Feb 2020)')
plt.ylabel('Factor of Safety')
plt.title('Variation of Factor of safety over Forecasting Period')
leg = plt.legend(loc = 1)
fig.set_size_inches((20, 10))
```

Extract individual percentile from the FoS and 
soil moisture list containing sub-list of all percentiles 

```python
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

percen_0 = extract_percen(FOS_percentile)[0]
percen_10 = extract_percen(FOS_percentile)[1]
percen_25 = extract_percen(FOS_percentile)[2]
percen_50 = extract_percen(FOS_percentile)[3]
percen_75 = extract_percen(FOS_percentile)[4]
percen_90 = extract_percen(FOS_percentile)[5]
percen_100 = extract_percen(FOS_percentile)[6]

percen_0_1 = extract_percen(soil_moisture_percentile)[0]
percen_10_1 = extract_percen(soil_moisture_percentile)[1]
percen_25_1 = extract_percen(soil_moisture_percentile)[2]
percen_50_1 = extract_percen(soil_moisture_percentile)[3]
percen_75_1 = extract_percen(soil_moisture_percentile)[4]
percen_90_1 = extract_percen(soil_moisture_percentile)[5]
percen_100_1 = extract_percen(soil_moisture_percentile)[6]
```

Plotting the FoS percentiles generated from the seasonal forecast 
compared against the observed reanalyzed data

Input
- line 2 and 7, Activate when the max and min value is needed
- `label` = label of the line or hatch colour shown in the graph 
corresponding to the data
- `DateFormatter` = datetime presented in the graph
- `xlabel`, `ylabel` = axis title
- `title` = graph title
- `plt.legend(loc = 1)` = designed legend location
- plt.savefig, Activate when storing the graph is needed
```python
fig, ax = plt.subplots()
#ax.fill_between(date_time_plot, percen_0, percen_10, alpha=.05, linewidth=0, color = 'b', label = 'data coverage')
ax.fill_between(date_time_plot, percen_10, percen_25, alpha=.1, linewidth=0, color = 'b', label = '10 and 90 percentile')
ax.fill_between(date_time_plot, percen_25, percen_50, alpha=.25, linewidth=0, color = 'b',label = '25 and 75 percentile')
ax.fill_between(date_time_plot, percen_50, percen_75, alpha=.25, linewidth=0, color = 'b')
ax.fill_between(date_time_plot, percen_75, percen_90, alpha=.1, linewidth=0, color = 'b')
#ax.fill_between(date_time_plot, percen_90, percen_100, alpha=.05, linewidth=0, color = 'b')

feb_forec_plot = ax.plot(date_time_plot, percen_50,
                        color = 'g',alpha=1,
                        label = 'February 2020 forecast median')

feb_obs_plot = ax.plot(obs_datetime_l, obs_FoS,
                      color = 'r',
                      label = 'February 2020 observed' )

ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
date_form = DateFormatter("%d-%m")
ax.xaxis.set_major_formatter(date_form)
ax.set(xlim = [obs_datetime_l[0], obs_datetime_l[671]])
plt.xlabel('Forecast Period(Feb 2020)')
plt.ylabel('Factor of Safety')
plt.title('Variation of Factor of safety over Forecasting Period title to be writing below figure in d')
leg = plt.legend(loc = 1)
fig.set_size_inches((20, 10))
#plt.savefig('feb_forecast_against_observed.png')
```
