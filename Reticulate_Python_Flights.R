## This program modified for windos
## In windows I have to set QT paltform plugi. In mac it was Ok.
setwd("F:/OneDrive - University of Nebraska Medical Center/Python")
## https://www.kaggle.com/miquar/explore-flights-csv-airports-csv-airlines-csv#
library(reticulate)

repl_python() ### Start python
##Python start
import pandas as pd
flights = pd.read_csv("flights.csv", low_memory=False)
flights = flights[flights['DESTINATION_AIRPORT'] == "ORD"]
flights = flights[['AIRLINE', 'DEPARTURE_DELAY', 'ARRIVAL_DELAY']]
flights = flights.dropna()
flights.head()
exit 
### Exit from Python


## Using pandas daraFrame in R
library(ggplot2)
ggplot(py$flights, aes(AIRLINE, ARRIVAL_DELAY)) + geom_point() + geom_jitter()
flights <- py$flights ## Save python dataframe in R data.frame

### Using R data.frame in Python
## start of Python
repl_python()
r.flights.head()
import matplotlib.pyplot as plt
import numpy as np
import os
os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = 'C:/Users/nitish.mishra/AppData/Local/Continuum/anaconda3/Library/plugins/platforms'
plt.scatter(r.flights['DEPARTURE_DELAY'], r.flights['ARRIVAL_DELAY'])
plt.show() # 
exit
### end of Python #########
