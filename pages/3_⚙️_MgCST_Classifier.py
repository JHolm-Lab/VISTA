import streamlit as st
import subprocess
import os
import tempfile
import multiprocessing
import pandas as pd
from streamlit_pdf_viewer import pdf_viewer

st.set_page_config(layout="wide")
st.sidebar.subheader("Contact")
st.sidebar.write("jholm@som.umaryland.edu")


# Streamlit app title
st.title("MgCSTs classifier")

st.container()
col1, col2 = st.columns(2)

with col1 :
    st.subheader("Select number of CPUs:")

    # Display the number of available cores
    available_cores = multiprocessing.cpu_count()
    st.write(f"Available CPUs: {available_cores}")

    # Slider for selecting the number of cores to use
    n_cores = st.slider("CPU:", 2, available_cores, 4)

    # # Display the selected number of cores
    # st.write(f"The classifier will use {n_cores} cores")

with col2:
    st.subheader("Import VIRGO2 output")

    # File uploader widget
    uploaded_file = st.file_uploader("Select a .txt file from VIRGO2 after mapping and compiled steps", type="txt")

# Check if a file has been uploaded
if uploaded_file is not None:
    # Save the uploaded file to a temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as temp_file:
        temp_file.write(uploaded_file.read())
        temp_file_path = temp_file.name

    # Read the uploaded file
    with open(temp_file_path, "r") as file:
        file_content = file.read()

    # Execute the R script
    r_script_path = "mgCST_classifier_v2.R"
    virgo2_path = "VIRGO2/"
    mgCST_classifier_master_path = "mgCST-classifier-master"
    num_cores = n_cores

    # Build the command to run the R script
    command = [
        "Rscript",
        r_script_path,
        temp_file_path,
        virgo2_path,
        mgCST_classifier_master_path,
        str(num_cores)
    ]

    with st.spinner('Running the classifier...'):
        # Run the R script
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Display output or error
        if process.returncode == 0:
            st.success("mgCSTs classifier executed successfully!")
            st.title("Outputs")
            # st.text(stdout.decode("utf-8"))

        else:
            st.error("Error executing R script")
            # st.text(stderr.decode("utf-8"))

        # Delete temporary file after execution
        os.remove(temp_file_path)

        st.container()
        # Display CSV files from classifier
        col1, col2 = st.columns(2)
        csv_files = [file for file in os.listdir() if file.endswith(".csv")]

        for idx, file in enumerate(csv_files):
            if file.startswith("norm_"):
                with col1 :
                    st.subheader(f"{file}")
                    st.dataframe(pd.read_csv(file))
            else:
                with col2 :
                    st.subheader(f"{file}")
                    st.dataframe(pd.read_csv(file))

        # Display PDF file from classifier
        pdf_files = [file for file in os.listdir() if file.endswith(".pdf")]

        if pdf_files:
            for file in pdf_files:
                pdf_viewer(file)

                
