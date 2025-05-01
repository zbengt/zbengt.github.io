---
layout: post
title: 04.25.2025 Plate Reader Installation
date: '2025-04-25'
categories: TA
tags: assays
---

## Overview

This post documents how we verified compatibility between our lab's Agilent BioTek **Synergy HTX multimode plate reader** and the **Cayman Chemical metabolic assays** we plan to use (Glucose, Lactate, BCA Protein, and Triglycerides).

## Plate Reader Specifications

- **Model**: Agilent BioTek Synergy HTX
- **Detection Modes**:
  - Absorbance: Monochromator (200–999 nm, 1 nm steps)
  - Fluorescence: Filter-based
- **Available Fluorescence Filters**:
  - Excitation: 360/40, 485/20, 528/20
  - Emission: 460/40, 528/20, 590/20

## Assay Requirements vs Reader Capabilities

| Assay          | Requirement                  | Reader Capability        | Compatible? |
|----------------|-------------------------------|---------------------------|-------------|
| **Glucose**    | Absorbance at 514 nm          | Monochromator (514 nm)    | ✅           |
| **Lactate**    | Excitation 530–540 nm / Emission 585–595 nm | 528/20 Ex + 590/20 Em | ✅           |
| **BCA Protein**| Absorbance at 562 nm           | Monochromator (562 nm)    | ✅           |
| **Triglycerides** | Absorbance at 540 nm        | Monochromator (540 nm)    | ✅           |

## Summary

- **Lactate assay (fluorescence)** is fully compatible with our installed filters.
- **Glucose**, **BCA protein**, and **triglyceride assays (absorbance)** are covered by the monochromator.

## Notes for Running Assays

- **Fluorescence assays** (Lactate):
  - Excitation: 528/20 nm filter
  - Emission: 590/20 nm filter

- **Absorbance assays** (Glucose, BCA, Triglycerides):
  - Select target wavelengths manually in Gen5 software.
  - Use clear, flat-bottom plates for best accuracy.

- **Validation Step**:
  - Run a test plate with standards only to confirm proper calibration and signal linearity before beginning full experiments.

## Conclusion

With our current configuration, we can move forward confidently using the Cayman Chemical assays on the Synergy HTX. No additional hardware upgrades are needed!