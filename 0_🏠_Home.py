import streamlit as st

st.set_page_config(page_title='Home', layout="centered")

st.sidebar.subheader("Contact")
st.sidebar.write("jholm@som.umaryland.edu")

st.image("medias/Holm_Lab_Logo.png")

st.title("Metagenomic Community State Types version 2")
# st.header("Version 2")

st.markdown('''<div style="text-align: justify;">Metagenomic CSTs (mgCSTs) were developed to explore functionally distinct vaginal microbiomes using the genetic underpinnings of the constituent vaginal bacteria. \
            Recently, the genetic repertoire of the vaginal microbiome has been substantially expanded through VIRGO2.0, necessitating an overhaul of the original metagenomic subspecies and mgCSTs. \
            Together, VIRGO2.0 with mgSs and mgCSTs assignment provides a streamlined tool for quantifying the relationships between the vaginal microbiome and urogenital health outcomes</div>''', unsafe_allow_html = True)

st.write("https://www.jbholmlab.org/")