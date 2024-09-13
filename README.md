# MgCST classifier v2

This application provides visualizations and tools related to [Metagenomic Community State Types v2] (refer to the paper). It also offers a user-friendly interface to run the classifier on [VIRGO2](reference VIRGO2 paper) output. The application is run through a Conda environment, which can be easily created using a provided YAML file. Follow the steps below to set it up:

# Run the application

1. Navigate to the folder where you want to save the application
   ```bash
    cd path/to/folder
   ```

2. Clone the repository:
    ```bash
    git clone https://github.com/JHolm-Lab/mgCST-classifier-v2.git
    ```
3. Navigate to the application directory:
    ```bash
    cd mgCST-classifier-v2
    ```

4. Download and unzip the [necessary data and media]() into the application folder:
    ```bash
    # Download data related to the application
    wget <download_url>
    unzip data.zip
    
    # Download medias related to the application
    wget <download_url>
    unzip medias.zip
    
    # Download VIRGO2 annotation files related to the application and the classifier alone
    wget <download_url>
    unzip VIRGO2.zip
    
    # Download mgCST-classifier-master directory related to the application and the classifier alone
    wget <download_url>
    unzip mgCST-classifier-master.zip
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

## Notes

1. The classifier can be run without the application.
   
   Rscript mgCST_classifier_v2.R path/to/VIRGO2_Compiled_Output.txt path/to/VIRGO2_annotation_files path/to/mgCST-classifier-master n_cores
   
   Example:
   ```bash
   Rscript mgCST_classifier_v2.R VIRGO2_Compiled.summary.NR.txt ./VIRGO2 ./mgCST-classifier-master 4
   ```
   Output is written to current directory. Each output file is dated.

3. Adjust the maximum size of uploaded file (if needed, currently, the limit is set to 30GB):

   a. Access the Streamlit config file (```.streamlit/config.toml```)
    ```bash
    cd path/to/mgCST-classifier-v2/.streamlit
    ```
   b. Modify the ```config.toml``` file
   ```toml
   [server]

    maxUploadSize = 30000
   ```
