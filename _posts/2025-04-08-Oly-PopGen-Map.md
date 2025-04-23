---
title:  "04.08.2025 Oly Pop Gen Manuscript Map"
author: "Zachary Bengtsson"
date:   2025-04-23
layout: post
tags:   [oly, mapping, pop gen]
---

> *A fully-commented walk-through of an R script that builds a publication-quality map of northwest Washington, ready for reuse or adaptation.*

------------------------------------------------------------------------

## 1 · Load the libraries

``` r
library(sf)
library(ggplot2)
library(ggmap)
library(ggspatial)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)  # for alpha()
```

**Why?**\
\* `sf`, `rnaturalearth`, and `dplyr` handle vector data and simple-feature operations.\
\* `ggplot2` and `ggspatial` layer data, add scale-bars & north arrows.\
\* `ggmap` retrieves tiled basemaps.\
\* `scales::alpha()` lets us tweak fill transparency later.

------------------------------------------------------------------------

## 2 · Authenticate with Stadia Maps

``` r
register_stadiamaps(key = "7da33b85-7352-45de-8189-217884e7376a")
```

Registers your **Stadia Maps** API key so `ggmap` can download high-resolution tiles.\
*Replace the key with an environment variable for production work.*

------------------------------------------------------------------------

## 3 · Define the map extent

``` r
nw_wash_bbox <- c(left = -124.5, bottom = 47.0, right = -122.0, top = 49.0)
```

A named numeric vector that sets the *bounding box* for northwest Washington.\
Top edge fixed at **49.0 °N** to match international border.

------------------------------------------------------------------------

## 4 · Read the basin polygons

``` r
basins_sf <- st_read("/Users/zacharybengtsson/Desktop/puget-sound-basins")
```

Loads the Puget-Sound basin shapefile into an **sf** object.

------------------------------------------------------------------------

## 5 · Remove basins we don’t need

``` r
basins_sf <- basins_sf %>% 
  filter(!regn_nm %in% c("Grays Harbor", "Willapa Bay"))
```

Drops **Grays Harbor** and **Willapa Bay** for a cleaner map.

------------------------------------------------------------------------

## 6 · Define key locations

``` r
locations <- data.frame(
  Name       = c("South Sound", "Central Sound", "North Sound"),
  Latitude   = c(47.398583,    47.591472,     48.481944),
  Longitude  = c(-122.819722,  -122.672722,   -122.574583)
)
```

Defines labels and coordinates for regional markers.

------------------------------------------------------------------------

## 7 · Convert to spatial features

``` r
locations_sf <- st_as_sf(locations,
                         coords = c("Longitude", "Latitude"),
                         crs    = 4326)
```

Creates spatial points using WGS84 (EPSG:4326).

------------------------------------------------------------------------

## 8 · Download a basemap

``` r
basemap <- get_stadiamap(bbox     = nw_wash_bbox,
                         zoom     = 8,
                         maptype  = "stamen_terrain_background")
```

Fetches a terrain-style background from Stadia Maps.

------------------------------------------------------------------------

## 9 · Create the map

``` r
map_plot <- ggmap(basemap) +
  geom_sf(data = basins_sf, aes(fill = regn_nm), inherit.aes = FALSE,
          color = "black", linewidth = 0.5, alpha = 0.6) +
  geom_sf(data = locations_sf, inherit.aes = FALSE,
          color = "black", fill = "white", shape = 21, size = 2, stroke = 1.5) +
  geom_label(data = locations, aes(x = Longitude, y = Latitude, label = Name),
             hjust = 1, vjust = -1, size = 4, color = "black", fontface = "bold",
             fill = "white", alpha = 0.8, inherit.aes = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3, style = "ticks") +
  annotation_north_arrow(location = "tl", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         height = unit(0.8, "cm"), width = unit(0.8, "cm"),
                         style = north_arrow_orienteering(
                           line_width = 0.5,
                           line_col = "black",
                           fill = c("black", "black"),
                           text_col = "black",
                           text_size = 6)) +
  ggtitle("Sound Locations and Basins in Northwest Washington") +
  theme_minimal() +
  theme(
    legend.position = c(0.2, 0.24),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.15, "cm")
  )
```

------------------------------------------------------------------------

## 10 · Export the map

``` r
ggsave("nw_washington_final_map.png", plot = map_plot,
       width = 8, height = 6, dpi = 300)
```

Saves a high-resolution PNG of the final figure.

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/nw_washington_final_map.png?raw=true)

------------------------------------------------------------------------
