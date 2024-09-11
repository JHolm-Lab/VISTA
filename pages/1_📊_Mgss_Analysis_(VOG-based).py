import streamlit as st
import pandas as pd
import pickle
import plotly.graph_objects as go

st.set_page_config(layout="wide")
st.sidebar.subheader("Contact")
st.sidebar.write("jholm@som.umaryland.edu")

# Function to lazy load CSV file
@st.cache_data()
def load_csv(file_path):
    return pd.read_csv(file_path)

# Function to lazy load PKL file
@st.cache_data()
def load_pkl(file_path):
    with open(file_path, 'rb') as f:
        pkl = pickle.load(f)
    return pkl

# Main function
# def main():
        
######################################## ----- DATA IMPORTATION ----- ########################################

vog_mgss_coverage = load_csv("data/vog.mgSs.coverage.stats.csv").rename(columns={"Unnamed: 0" : "sub_species", "as.factor(sample_cluster)" : "sample_cluster"}) # mgss coverage stats
VOG_GeneProduct = load_csv("data/VOG_gene_product.csv")
reorder_dataframe = load_pkl('data/reorder_dataframe.pkl')
hover = load_pkl('data/hover_dict.pkl')
vog_clusters = load_pkl('data/vog_clusters.pkl')
gene_count = load_pkl('data/gene_pa_count.pkl')

######################################## ----- DATA PROCESSING ----- ########################################

not_to_cluster = ["Alterileibacterium", "Anaerococcus", "Bacteroides", "Campylobacter",
                "Corynebacterium", "Gardnerella", "Gulosibacter", "Lactobacillus", 
                "Limosilactobacillus", "MultiGenera", "Porphyromonas", "Prevotella", "Streptococcus"]

vog_species = vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]).unique()
species = [s for s in vog_species if s not in not_to_cluster]

######################################## ----- PAGE CONTENT ----- ########################################

st.title("Metagenomic Subspecies (VOG-based)")

st.subheader("Visualizations")

option = st.selectbox("Species", species)

gene_count_df = gene_count[option]
Gene = pd.DataFrame(gene_count_df.set_index('Gene').loc[:,gene_count_df.columns[1]:].sum(axis = 1)).rename(columns = {0:"Number_of_samples"}).reset_index()
Gene['%_of_samples'] = Gene['Number_of_samples'].apply(lambda x : round((100 * x / (len(gene_count_df.columns) - 1)),2))

tab1, tab2, tab3 = st.tabs(["Species coverage", "Subspecies stats", "Presence Absence heatmap"])

with tab1:

        col1, col2 = st.columns(2)

        with col1 :
                st.image("medias/vog_mgss_coverage_png/" + option + "_subspecies_coverage_boxplot.png")

        with col2 :
                st.image("medias/vog_mgss_coverage_png/" + option + "_subspecies_coverage_by_NoVOG.png")

with tab2 :

        st.subheader("Subspecies stats")
        st.dataframe(vog_mgss_coverage[vog_mgss_coverage['sub_species'].apply(lambda x : x.split(".")[0]) == option])

with tab3 :

        col1, col2 = st.columns(2)

        with col1:
                st.image("medias/vog_heatmap_presence_absence/_" + option + "_heatmap_presence_absence.png")
        with col2:
               df_reorder = reorder_dataframe[option]
               
               heatmap = go.Figure(data=go.Heatmap(
                        z=df_reorder,
                        x=df_reorder.columns,
                        y=df_reorder.index,
                        colorscale=[[0, 'antiquewhite'], [1, 'mediumblue']],
                        showscale=False, 
                        hovertext=hover[option],
                        hovertemplate="VOG: %{y}<br>SampleID: %{x}<br>GeneProduct: %{hovertext}<extra></extra>"
                ))
               heatmap.update_layout(
                # title=f'{option} Presence-Absence Heatmap',
                # xaxis=dict(title='Samples'),
                # yaxis=dict(title='VOG'),
                width=3600,  # Set width in pixels
                height=900
                )
               # Update x axis
               heatmap.update_xaxes(
                showticklabels=False
                )
               # Update y axis
               heatmap.update_yaxes(
                showticklabels=False
                )
               # Display heatmap
               st.plotly_chart(heatmap, use_container_width=True)

