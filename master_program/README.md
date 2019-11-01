Version 1.0
Release date: Oct 19, 2019

# General Info:

This initial release includes Matlab programs used to process raw and gridded COVIS data in various modes and formats. Those programs are divided into `lower level` and `upper level` groups. The lower-level programs read the raw data recorded in Imaging, Diffuse-flow, and Bathymetry modes and generate gridded data after processing following the procedures described in Xu et al., (2019). The upper-level programs read those gridded data to produce various data products such as 3D plume images, maps of bathymetry and diffuse-flow distributions, etc. In the current release, the lower-level programs include generic code that processes raw COVIS data in a compressed format (e.g., `.gz` or `.7z`). There is also specialized code for processing decompressed raw data recorded in specific modes (e.g., Imaging, Diffuse-flow, Bathymetry). The upper-level programs include code for making bathymetry maps and 3D plume images.

# Instructions for User:

After checking out this [Github repo](https://www.github.com/) to a local directory, the first step before any processing is to use Matlab’s ‘Set Path’ function (under ‘Environment’) to add `master_program` and its subfolders to Matlab’s search path list.

To process a raw compressed dataset, the user should use the function `covis_raw_sweep`, which decompresses the raw dataset, identifies the data type (e.g., Imaging), and processes the raw data using a designated function (e.g., `covis_imaging_sweep.m`). Alternatively, the user can decompress the raw data in a separate step and choose a specialized function (e.g., `covis_imaging_sweep`) based on the data type to process the decompressed raw data. Once finished, the user can choose a program in the `upper_level` folder to read the gridded data and make the corresponding data product (e.g., plume images). More information on the purpose and use of each program can be found in its header.  
