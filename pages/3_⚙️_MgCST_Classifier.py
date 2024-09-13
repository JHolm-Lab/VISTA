import streamlit as st
import subprocess
import os
import tempfile
import multiprocessing
import pandas as pd
from streamlit_pdf_viewer import pdf_viewer
import io
import zipfile

# Clean folder before running the classifier (remove former outputs)
for file in os.listdir():
    if file.startswith("mgCSTs_") or file.startswith("norm_counts_") or file.startswith("relabund_w_") or file.startswith("mgCST_heatmap"):
        os.remove(file)

# Set page layout
st.set_page_config(layout="wide")
st.sidebar.subheader("Contact")
st.sidebar.write("jholm@som.umaryland.edu")

st.title("MgCSTs classifier")

# Create a form for user input
with st.form(key='classifier_form'):
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("CPU")
        available_cores = multiprocessing.cpu_count()
        n_cores = st.slider("Choose a number of CPU to use", 2, available_cores, 4)
    with col2:
        st.subheader("VIRGO2 output")
        uploaded_file = st.file_uploader("Select a .txt file from VIRGO2 after mapping and compiled steps", type="txt")
    
    submit_button = st.form_submit_button(label='Submit')

if submit_button and uploaded_file is None:
    st.warning("Select an input file")

if submit_button and uploaded_file is not None:
    # Save the uploaded file to a temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as temp_file:
        temp_file.write(uploaded_file.read())
        temp_file_path = temp_file.name

    # Execute the R script
    r_script_path = "mgCST_classifier_v2.R"
    virgo2_path = "VIRGO2/"
    mgCST_classifier_master_path = "mgCST-classifier-master"
    num_cores = n_cores

    command = [
        "Rscript",
        r_script_path,
        temp_file_path,
        virgo2_path,
        mgCST_classifier_master_path,
        str(num_cores)
    ]

    # Initialize progress bar
    with st.spinner('Running the classifier...'):
    
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            os.remove(temp_file_path)

            # Collect CSV and PDF files
            file_list = []
            csv_files = [file for file in os.listdir() if file.endswith(".csv")]
            pdf_files = [file for file in os.listdir() if file.endswith(".pdf")]

            # Update progress for CSV files
            for file in csv_files:
                file_list.append({'file_name': file, 'data': pd.read_csv(file).to_csv(index=False)})

            # Update progress for PDF files
            for pdf_file in pdf_files:
                with open(pdf_file, 'rb') as f:
                    file_list.append({'file_name': pdf_file, 'data': f.read()})

            st.success("mgCSTs classifier executed successfully!")

            container = st.container()
            container.info("Don't forget to download your files", icon="⚠️")

            # Create a ZIP file in memory
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                for file_info in file_list:
                    zip_file.writestr(file_info['file_name'], file_info['data'])

            zip_buffer.seek(0)

            container.download_button(label="Download all files", data=zip_buffer, file_name="mgCST_classifier_v2_outputs.zip", mime="application/zip")

            st.container()
            st.title("Outputs")

            # Display CSV and PDF files in columns
            col1, col2 = st.columns(2)
            col3, col4 = st.columns(2)

            for file in csv_files:
                if file.startswith("norm_"):
                    with col1:
                        st.subheader(f"{file}")
                        st.dataframe(pd.read_csv(file).head(10))
                elif file.startswith("relabund_"):
                    with col2:
                        st.subheader(f"{file}")
                        st.dataframe(pd.read_csv(file).head(10))
                elif file.startswith("mgCST"):
                    with col3:
                        st.subheader(f"{file}")
                        st.dataframe(pd.read_csv(file).head(10))

            for file in pdf_files:
                with col4:
                    pdf_viewer(file)

        else:
            st.error("Error executing R script")
            st.error(stderr)