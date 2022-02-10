---
layout: post
title: Google Earth Engine and the ORCAA Tool for Aquatic Remote Sensing of Environmental Parameters
date: '2021-12-10'
categories: E5 Coral
tags: RemoteSensing
---

*TLDR Notice: Start at the Products and Interface section if you aren't interest in background*

### Using the Google Earth Engine platform to load, process, analyze, and visualize remote sensing data

Prior to joining the Roberts Lab this Fall, I worked for NASA Earth Applied Sciences’ [DEVELOP]( https://develop.larc.nasa.gov/) and [ARSET]( https://appliedsciences.nasa.gov/what-we-do/capacity-building/arset) programs. During this time, I became familiar a variety of environmental remote sensing techniques and platforms. One of the easiest ways to work with remote sensing data is Google Earth Engine (GEE). GEE is a cloud-based raster computing platforms that uses a JavaScript API to manipulate remote sensing data. A variety of satellite images and data products are hosted through GEE (including Landsat and Sentinel imagery). This negates data storage issues and local machine computing limits. If you would like to learn more about GEE in general, please refer to [this training]( https://appliedsciences.nasa.gov/join-mission/training/english/arset-using-google-earth-engine-land-monitoring-applications) I gave in June.

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Overall_GEE_Interface.png)
_The Google Earth Engine interface. Access to GEE is free, so you can sign up for an account at no cost. Free access is guaranteed to scientists, students, and natural resource managers for the foreseeable future._

### Aquatic environmental parameter assessment and using the Optical Reef and Coastal Area Assessment (ORCAA) Tool

For the purposes of coral reef ecology, satellite remote sensing has the potential to expand spatial and temporal coverage of parameters such as chlorophyll, turbidity, SST, and CDOM. Public satellite technology made available by government agencies (namely NASA and ESA) can provide data at a 10m-30m resolution. Temporal resolution varies between sensors, but using data from more than one satellite sensor can provide multiple images per month.

During my time at NASA DEVELOP a GEE tool was developed for remote assessment of various parameters related to reef and coastal water quality using imagery from the NASA Landsat 8 OLI sensor and the ESA Sentinel-2 MSI sensor. The tool implements a graphical user interface (GUI) giving the user the ability to call a time series of imagery over a specified period and spatial extent, providing a ready made platform to alter and specialize for the application of water quality algorithms according to the needs of research goals.

### Products and Interface



Pre-processed products include Turbidity, Chlorophyll-a, Normalized Difference Chlorophyll Index (NDCI), Dissolved Organic Matter (CDOM), and Sea Surface Temperature (SST). Links are to publications where the following equations were developed for calculation of these parameters from surface reflectance bands:
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/cb37cacfbcae39d6bc6e986540f88ea684e604dd/assets/img/ORCAA_ParameterCalculationEquations.png)



### _More to come, post in progress_
