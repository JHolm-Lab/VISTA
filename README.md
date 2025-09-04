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
| pandas                 | gplots             |
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
    virgo2_path,
    mgCST_classifier_master_path,
    str(num_cores)
]
```

To verify your system’s Rscript path:

```bash
which Rscript
```

---

## 🚀 Running VISTA Locally

1. **Download** and **unzip** the [application image](https://figshare.com/ndownloader/files/:

```bash
cd path/to/app/directory
```

3. Launch the app:

```bash
streamlit run 0_Home.py
```

You can test the classifier using the provided example file.

---

## 🧪 Running the Classifier Without the App

You can run the classifier directly via R:

```bash
Rscript path/to/mgCST_classifier_v2.R /path/to/VIRGO2_output_Compiled.summary.NR.txt path/to/VIRGO2 path/to/mgCST-classifier-master n_cores
```

**Example:**

```bash
Rscript mgCST_classifier_v2.R VIRGO2_mgCST_example.txt volume/VIRGO2 volume/mgCST-classifier-master 4
```

- Output files are saved to the current directory.
- Each output file is timestamped.

---

## 📁 File Upload Limit

To adjust the maximum upload size (default: 30GB):

1. Open the Streamlit config file: `.streamlit/config.toml`

2. Modify the following line:

```toml
[server]
maxUploadSize = 30000  # Change this value as needed
```

