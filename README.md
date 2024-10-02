# MgCST classifier v2

This application provides visualizations and tools related to [Metagenomic Community State Types v2](). It also offers a user-friendly interface to run the classifier on [VIRGO2]() output. The application is run through a Conda environment, which can be easily created using a provided YAML file. Follow the steps below to set it up:

# Run the application locally

1. **Download** and **unzip** the [image of the application](https://figshare.com/ndownloader/files/49573392)
   ```bash
   # Download the image via the provided link or by running the following:
   wget https://figshare.com/ndownloader/files/49573392

   # Unzip the archive
   unzip 49573392
   ```
   
2. Navigate to the application directory:
    ```bash
    cd mgCST-classifier-v2
    ```
    
3. Create the conda environment:
    ```bash
    conda env create -f env.yaml
    ```
4. Run the app:
    ```bash
    conda activate mgcst_app
    streamlit run 0_🏠_Home.py
    ```

    You can test the classifier through the application using the provided example file: ```VIRGO2_Compiled_example.summary.NR.txt```

## Notes

1. The classifier can be run without the application.
      
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
