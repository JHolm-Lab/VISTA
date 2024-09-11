# MgCST classifier v2

This application provides visualizations and tools related to [Metagenomic Community State Types v2] (refer to the paper). It also offers a user-friendly interface to run the classifier on [VIRGO2](reference VIRGO2 paper) output. The application is run through a Conda environment, which can be easily created using a provided YAML file. Follow the steps below to set it up:

## Run the classifier alone

1. Clone the repository:
    ```bash
    git clone <url> 
    ```
2. Navigate to the repo directory:
    ```bash
    cd /path/to/mgCST_classifier_v2
    ```
3. Run the classifier (choose n)
   ```bash
   Rscript mgCST_classifier_v2.R ./VIRGO2 ./mgCST-classifier-master 4
   ```
   
## Run the application

1. Clone the repository:
    ```bash
    git clone <url> 
    ```

2. Navigate to the application directory:
    ```bash
    cd /path/to/mgCST_classifier_v2
    ```

3. Download the necessary data and media:
    ```bash
    wget <download_url>
    ```

4. Create and activate the conda environment:
    ```bash
    conda env create -f env.yaml
    ```
    ```bash
    conda activate mgcst_app
    ```

5. Run the app:
    ```bash
    streamlit run 0_🏠_Home.py
    ```

## Notes

1. Adjust the maximum size of uploaded file (if needed, currently, the limit is set to 30GB):

   a. Access the Streamlit config file (```.streamlit/config.toml```)
    ```bash
    cd path/to/mgCST_classifier_v2/.streamlit
    ```
   b. Modify the ```config.toml``` file
   ```toml
   [server]

    maxUploadSize = 30000
   ```
