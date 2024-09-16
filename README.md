# MgCST classifier v2

This application provides visualizations and tools related to [Metagenomic Community State Types v2](). It also offers a user-friendly interface to run the classifier on [VIRGO2]() output. The application is run through a Conda environment, which can be easily created using a provided YAML file. Follow the steps below to set it up:

# Run the application

1. Navigate to the folder where you want to save the application
   ```bash
    cd <path_to_your_folder>
   ```

2. Clone the repository:
    ```bash
    git clone https://github.com/JHolm-Lab/mgCST-classifier-v2.git
    ```
3. Navigate to the application directory:
    ```bash
    cd mgCST-classifier-v2
    ```
4. Download and unzip the necessary data and media into the application folder:
    ```bash
    # Download data related to the application
    wget https://figshare.com/ndownloader/files/49149751
    unzip 49149751
    
    # Download medias related to the application
    wget https://figshare.com/ndownloader/files/49149853
    unzip 49149853
    
    # Download VIRGO2 annotation files related to the application and the classifier alone
    wget https://figshare.com/ndownloader/files/49149880
    unzip 49149880
    
    # Download mgCST-classifier-master directory related to the application and the classifier alone
    wget https://figshare.com/ndownloader/files/49178470
    unzip 49178470

    # Remove zip files and hidden folder
    rm -rf 49149751 49149853 49149880 49178470 __MACOSX 
    ```
   The command ```ls -l``` should return this :
   
    ```bash
   (base) mgCST-classifier-v2 % ls -l
   total 29664
   -rw-r--r--   1 amaros  staff       955 Sep 13 12:51 0_🏠_Home.py
   -rw-r--r--   1 amaros  staff      2397 Sep 13 12:51 README.md
   drwxr-xr-x   6 amaros  staff       192 Sep 10 09:10 VIRGO2
   -rw-r--r--   1 amaros  staff  15152457 Sep 13 12:51 VIRGO2_Compiled_example.summary.NR.txt
   drwxr-xr-x  14 amaros  staff       448 Sep 11 17:00 data
   -rw-r--r--   1 amaros  staff       212 Sep 13 12:51 env.yaml
   drwxr-xr-x   7 amaros  staff       224 Sep 10 09:26 medias
   drwxr-xr-x   4 amaros  staff       128 Sep 13 12:42 mgCST-classifier-master
   -rw-r--r--   1 amaros  staff     17418 Sep 13 12:51 mgCST_classifier_v2.R
   drwxr-xr-x   5 amaros  staff       160 Sep 13 12:51 pages
    ```
    
5. Create the conda environment:
    ```bash
    conda env create -f env.yaml
    ```
6. Run the app:
    ```bash
    conda activate mgcst_app
    streamlit run 0_🏠_Home.py
    ```

    You can test it with the example file provided : ```VIRGO2_Compiled_example.summary.NR.txt```

## Notes

1. The classifier can be run without the application.
   
   Rscript mgCST_classifier_v2.R path/to/VIRGO2_Compiled_Output.txt path/to/VIRGO2 path/to/mgCST-classifier-master n_cores
   
   Example:
   ```bash
   Rscript mgCST_classifier_v2.R VIRGO2_Compiled_example.summary.NR.txt ./VIRGO2 ./mgCST-classifier-master 4
   ```
   Output is written to current directory.

   Each output file is dated.

3. Adjust the maximum size of uploaded file (if needed, currently, the limit is set to 30GB):

   a. Access the Streamlit config file (```.streamlit/config.toml```)

   b. Modify the ```config.toml``` file
   ```toml
   [server]

    maxUploadSize = 30000   # modify this value if needed
   ```
