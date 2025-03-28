# MgCST classifier v2

This application provides visualizations and tools related to Metagenomic Community State Types v2. 

It also offers a user-friendly interface to run the classifier on VIRGO2 output.

# Requirements
The application is run through Streamlit with the following requirements:

<div align="center">

| Python (v3.8+)         | R (v4.3+)          |
|------------------------|--------------------|
| streamlit              | randomForestSRC    |
| pandas                 | gplots             |
| numpy                  | data.table         |
| plotly                 | dplyr              |
| streamlit-pdf-viewer   | parallel           |

</div>

# Run the application locally

1. **Download** and **unzip** the [image of the application](https://figshare.com/ndownloader/files/53286029)
   
3. Open a terminal and navigate to the application directory
   ```bash
    cd path/to/app/directory
    ```

4. Run the app:
    ```bash
    streamlit run 0_Home.py
    ```

    You can test the classifier through the application using the provided example file.
   
## Notes

1. The classifier can be run without the application:
   
   ```bash
   Rscript path/to/mgCST_classifier_v2.R path/to/VIRGO2_output.txt.gz path/to/VIRGO2 path/to/mgCST-classifier-master n_cores

   # Example after downloading the image:
   Rscript mgCST_classifier_v2.R VIRGO2_mgCST_example.txt volume/VIRGO2 volume/mgCST-classifier-master 4
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
