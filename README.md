# Nuclei_ORF1P

* **Developed for:** Rania
* **Team:** Fuchs
* **Date:** October 2022
* **Software:** Fiji

### Images description

3D images taken with a x60 objective

2 channels:
  1. *CSU_405:* DAPI nuclei
  2. *CSU_488:* ORF1p cells

### Plugin description

* Detect DAPI nuclei and ORF1p cells with Cellpose
* Keep ORF1p cells with a nucleus only
* Detect ORF1p dots with Stardist
* Give various measurements in the nucleus, inner nucleus, inner ring, outer ring and cytoplasm

### Dependencies

* **3DImageSuite** Fiji plugin
* **Cellpose** conda environment + *cyto2* model
* **StarDist** conda environment + *pmls2.zip* (homemade) model

### Version history

Version 1 released on October 14, 2022.
