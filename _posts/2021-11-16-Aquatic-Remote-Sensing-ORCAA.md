---
layout: post
title: Google Earth Engine and the ORCAA Tool for Aquatic Remote Sensing of Environmental Parameters
date: '2022-02-10'
categories: E5 Coral
tags: RemoteSensing
---

### Using the Google Earth Engine platform to load, process, analyze, and visualize remote sensing data

*TL;DR Notice: Start at the _"Products and Interface"_ section if you aren't interest in background*

Prior to joining the Roberts Lab this Fall, I worked for NASA Earth Applied Sciences’ [DEVELOP]( https://develop.larc.nasa.gov/) and [ARSET]( https://appliedsciences.nasa.gov/what-we-do/capacity-building/arset) programs. During this time, I became familiar a variety of environmental remote sensing techniques and platforms. One of the easiest ways to work with remote sensing data is Google Earth Engine (GEE). GEE is a cloud-based raster computing platforms that uses a JavaScript API to manipulate remote sensing data. A variety of satellite images and data products are hosted through GEE (including Landsat and Sentinel imagery). This negates data storage issues and local machine computing limits. If you would like to learn more about GEE in general, please refer to [this training]( https://appliedsciences.nasa.gov/join-mission/training/english/arset-using-google-earth-engine-land-monitoring-applications) I gave in June.

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Overall_GEE_Interface.png)
_The Google Earth Engine interface. Access to GEE is free, so you can sign up for an account at no cost. Free access is guaranteed to scientists, students, and natural resource managers for the foreseeable future._

### Aquatic environmental parameter assessment and using the Optical Reef and Coastal Area Assessment (ORCAA) Tool

For the purposes of coral reef ecology, satellite remote sensing has the potential to expand spatial and temporal coverage of parameters such as chlorophyll, turbidity, SST, and CDOM. Public satellite technology made available by government agencies (namely NASA and ESA) can provide data at a 10m-30m resolution. Temporal resolution varies between sensors, but using data from more than one satellite sensor can provide multiple images per month.

During my time at NASA DEVELOP a GEE tool was developed for remote assessment of various parameters related to reef and coastal water quality using imagery from the NASA Landsat 8 OLI sensor and the ESA Sentinel-2 MSI sensor. The tool implements a graphical user interface (GUI) giving the user the ability to call a time series of imagery over a specified period and spatial extent, providing a ready made platform to alter and specialize for the application of water quality algorithms according to the needs of research goals.

### Products and Interface

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/NDCI_LayerDisplayedFor2019Aug.png)
*ORCAA Interface displaying NDCI*
Technical information in this section was originally prepared by the creators of the ORCAA tool: Hayley Pippin, Arbyn Olarte, Roxana Pilot, and Vanessa Valenti

Data inputs from common NASA and ESA sensors were used via the [GEE internal data catalog](https://developers.google.com/earth-engine/datasets):
![iamge](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/ORCAA_RemoteSensingDataInputs.png)

Parameters include [Turbidity](https://doi.org/10.1117/12.830700), [Normalized Difference Chlorophyll Index (NDCI) and Chlorophyll-a](https://doi.org/10.1016/j.rse.2011.10.016), [Dissolved Organic Matter (CDOM)](https://doi.org/10.1117/1.JRS.11.036007), and [Sea Surface Temperature (SST)](https://modis.gsfc.nasa.gov/data/dataprod/mod28.php). Links are to publications which established the following equations for calculation of these parameters from surface reflectance bands:
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/cb37cacfbcae39d6bc6e986540f88ea684e604dd/assets/img/ORCAA_ParameterCalculationEquations.png)

Interface allows you to select dates of interest and upload your own spatial boundaries:
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/SetDateParametersGUI.png)

I chose an aggregate of data from Jan - Nov in 2019 in this example over a polygon drawn in GEE (below).

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/MooreaTestAreaPolygon.png)

Check parameters you want to retrieve, calculate, and display. Relevant colorscales will populate:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/SelectLayer_ColorscaleLegendsGUI.png)

After you run, layers are made available in the viewer. Including the highest quality true color satellite image available in your date range:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/TrueColorSatelliteLayer.png)

Chlorophyll-a layer:
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Chlorophyll-a_Layer.png)
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Chl_Legend.png)
Zoom for pixel by pixel view:
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Chlorophyll-a_Zoom.png)

Time series chart generator allows you to display data over the course of your data range, values are an average of all pixels within your spatial boundary:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/TimeSeriesChartGeneratorGUI.png)

Imagery, data layers, time series graphs, and point data are all exportable in TIFF, JPEG, and CSV format. Updates under way by the NASA DEVELOP program are implementing new cloud filtering data and improved surface reflectance values. This format of data retrieval, manipulation, and display is particularly helpful in deciding which data products may be most useful in research and allow for the quick processing of satellite imagery.

Take a look at the [ORCAA repository](https://github.com/NASA-DEVELOP/ORCAA) for further information or use the link [here](https://code.earthengine.google.com/1005c8e1910d3ad104fe8bb2a927af95) (you'll need a GEE account) to test it out yourself. 







### _More to come, post in progress_
