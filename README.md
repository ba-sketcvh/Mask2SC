# Protocol for integrating external information with Xenium data at single-cell resolution
![STAR_AI_Template-01](https://github.com/user-attachments/assets/c0a01923-6ba4-4e3a-a2f1-6cd7ceec8d35)

This includes all the code used in the manuscript "Protocol for integrating external information with Xenium data at single-cell resolution", available here: .

The dataset used in the protocol can be found here: https://zenodo.org/records/15367950

1. Environment Setup
    - **environment.yml**: YAML file for setting up Python packages
    - **Install conda and add conda to path**
    - **Setup environment with either method**
      
      **Method 1: Automatic installation**
      ```bash
      git clone https://github.com/ba-sketcvh/Mask2SC.git
      cd Mask2SC
      conda env create -n mask2sc --file environment.yml
      conda activate mask2sc
      ```

      **Method 2: Manual installation**
      ```bash
      conda create -n mask2sc
      conda activate mask2sc
      pip install numpy pandas geopandas imageio shapely rasterio affine matplotlib scanpy squidpy
      # Squidpy installation may need to use "conda install -c conda-forge squidpy"
      ```
3. Image Registration 
    - **IHC Channel Separation.groovy**: Extracts IHC channels as individual .tiff files 
    - **Binary Mask Creation.groovy**: Creates a binary mask from annotations    
    
4. Integration
     - **complete_python_code.ipynb**: Complete Python code for integration and visualization of IHC images and H&E binary masks
   

