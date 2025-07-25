# Protocol for integrating external information with Xenium data at single-cell resolution
<div align="center">
  <img src="https://github.com/user-attachments/assets/c0a01923-6ba4-4e3a-a2f1-6cd7ceec8d35" alt="STAR_AI_Template-01" width="500"/>
</div>

<br>

# This repository includes all the code used in the manuscript **"Protocol for integrating external information with Xenium data at single-cell resolution"**, available [here]:(#) and [here] (#).

The dataset used in the protocol can be found [here](https://zenodo.org/records/15367950).

## 1. Installation of Qupath, Fiji, and Warpy.
- **[QuPath](https://qupath.github.io/)**: Open-source software for digital pathology image analysis.
- **[Fiji](https://fiji.sc/)**: A distribution of ImageJ focused on biological-image analysis.
- **[Warpy](https://imagej.net/plugins/bdv/warpy/warpy)**: Tool for image registration and warping.
 
## 2. Setup of python environment
- **environment.yml**: file for setting up Python environment with conda
- **Install conda and add conda to path**
- **Setup environment with either method**
      
  **Method 1: Automatic installation**
   
  In Bash:
  ```bash
  git clone https://github.com/ba-sketcvh/Mask2SC.git
  cd Mask2SC
  conda env create -n mask2sc --file environment.yml
  conda activate mask2sc
  ```

  **Method 2: Manual installation**

  In Bash:   
  ```bash
  conda create -n mask2sc
  conda activate mask2sc
  pip install numpy pandas geopandas imageio shapely rasterio affine matplotlib scanpy
  conda install -c conda-forge squidpy
  ```
  
## 3. Image registration with Warpy
- Steps are availlable in our [website](#).
- More information on Warpy [website](https://imagej.net/plugins/bdv/warpy/warpy).

## 4. Exportation of registered image
- **IHC Channel Separation.groovy**: Extracts IHC channels as individual .tiff files in Qupath
- **Binary Mask Creation.groovy**: Creates a binary mask from annotations in Qupath

## 5. Integration
- **complete_python_code.ipynb**: Complete Python code for integration and visualization of IHC images and H&E binary masks
   

