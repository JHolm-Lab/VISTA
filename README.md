<img src="assets/VISTA_logo.jpg" alt="VISTA Logo" width="300"/>

# VISTA: Vaginal Metagenome Community State Type Assignment

This application provides visualizations and tools for assigning **Metagenomic Community State Types (mgCSTs)** to vaginal metagenomes. It includes a user-friendly interface for running the classifier on **VIRGO2** output files.

---

## 🛠 Requirements

### Dependencies

This application runs via **Streamlit** and requires both **Python** and **R** environments:

<div align="center">

| Python (v3.8+)         | R (v4.3+)          |
|------------------------|--------------------|
| streamlit              | randomForestSRC    |
| pandas                 | pheatmap           |
| numpy                  | data.table         |
| plotly                 | dplyr              |
| seaborn                | parallel           |
| scipy                  |                    |
| streamlit-pdf-viewer   |                    |

</div>

To install the required Python libraries:

```bash
pip install streamlit pandas numpy plotly seaborn scipy streamlit-pdf-viewer
```

Required R packages are installed automatically via the classifier's R script.

---

## ⚙️ Classifier Command

The Python script calls the R classifier using the `Rscript` command. Default path (Mac/Linux):

```python
command = [
    "/usr/local/bin/Rscript",  # Update this path if needed
    r_script_path,
    temp_file_path,            # VIRGO2 output file
    VISTA_path
]
```

To verify your system’s Rscript path:

```bash
which Rscript
```

---

To set up the VISTA app and models required for `run_VISTA.R`, follow these steps:

### 1. 📥 Download VISTA Resources

Download the bundled archive from Figshare:  
**🔗 [VISTA_data.tar.gz](https://figshare.com/ndownloader/files/57689476)**  
This includes:  
- VISTA Streamlit app  
- mgSs classification models  
- mgCST reference centroids  

### 2. 📦 Move the Archive to Your VISTA Directory

Place the downloaded file into your cloned VISTA project folder:

    mv VISTA_data.tar.gz /path/to/VISTA/

### 3. 🔓 Unpack the Archive

Extract the contents of the archive:

    tar -xzvf VISTA_data.tar.gz

---

## ⌨️ Classifying mgCSTs with VISTA via Command Line

Run the mgCST classifier without launching the Streamlit app.
- The VISTA input file should be a VIRGO2 output summary, provided either as a plain text file or compressed with .gz.
```bash
# Usage:
#   Rscript mgCST_classifier_v2.R <VIRGO2_summary> <VISTA_data_dir>

Rscript path/to/mgCST_classifier_v2.R \
  path/to/VIRGO2_output_Compiled.summary.NR.txt \
  /path/to/VISTA

Rscript run_VISTA.R \
  VISTA_example.txt.gz \
  ~/bin/VISTA
```

- Output files are saved to the current working directory.
- Each output file is timestamped.

---
## 💻 Running the VISTA App
The VISTA app allows for exploration of VISTA mgCSTs and mgSs and classification of your own data. 

Launch the VISTA Streamlit interface by navigating to the app directory and running:

```bash
cd path/to/VISTA/VISTA_data
streamlit run 0_Home.py
```

By default, the app allows uploads up to 30GB. To change this limit:

Open the Streamlit configuration file located at:
```bash
.streamlit/config.toml
```
Update or add the following setting:
```bash
[server]
maxUploadSize = 30000  # Set your desired limit in megabytes
```

