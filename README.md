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

Python packages can be installed using:
```bash
pip install -r requirements.txt
```

Required R packages are incorporated in the script ```mgCST_classifier_v2.R```


# Run the application locally

1. **Download** and **unzip** the [image of the application](https://figshare.com/ndownloader/files//########)
   ```bash
   #Download the image via the provided link or by running the following:
   wget https://figshare.com/ndownloader/files/########
  
   #Unzip the archive
   tar -xvzf mgCST-classifier-v2.tar.gz
   ```
   
2. Navigate to the application directory:
    ```bash
    cd mgCST-classifier-streamlit
    ```

3. Run the app:
    ```bash
    streamlit run 0_Home.py
    ```

    You can test the classifier through the application using the provided example file: ```VIRGO2_mgCST_example.txt.gz```

## Notes

1. The classifier can be run without the application. After downloading and unzipping, run the following:
      
   Example:
   ```bash
   Rscript mgCST_classifier_v2.R VIRGO2_mgCST_example.txt.gz VIRGO2 mgCST-classifier-master 4
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
