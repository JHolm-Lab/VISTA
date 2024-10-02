import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np
from streamlit_pdf_viewer import pdf_viewer

st.set_page_config(layout="wide")
st.sidebar.subheader("Contact")
st.sidebar.write("jholm@som.umaryland.edu")

######################################## ----- DATA IMPORTATION ----- ########################################

mgcsts_samples = pd.read_csv("data/samples_w_mgCSTs.csv")
mgcsts = pd.read_csv("data/mgCSTs.csv")
projects = pd.read_csv("data/VIRGO2_projects_anonymous.csv")

######################################## ----- DATA PROCESSING ----- #############################################

mgcsts_samples_project = mgcsts_samples.merge(projects, on = "sampleID", how = "left")
data_to_replace = {
    "Chinese":"China",
    "FIJI":"Fiji",
    "Irish":"Ireland",
    "PreSSMat":"Bangladesh",
    "ZAPPS":"Zambia",
    "LSVF":"USA - AL",
    "Gates-Stanford":"USA - CA, AL",
    "Gates-UMD":"USA - MD, AL",
    "MMOTH":"USA - MD",
    "STING":"USA - MD",
    "Symptoms":"USA - AL",
    "VMRC":"USA - AL",
    "XPENN":"USA - PA"
}

study_order = ["mgCST", "Bangladesh", "China", "Fiji", "Ireland", "Zambia", "USA - AL", "USA - MD", "USA - PA", "USA - MD, AL", "USA - CA, AL"]

color_pal = {
  "China":"#FDAE61",
  "Fiji":"#D73027",
  "Bangladesh":"#F46D43",
  "Ireland":"#4D9221", 
  "Zambia":"#FEE090",
  "USA - MD":"#C6DBEF",
  "USA - AL":"#6BAED6",
  "USA - MD, AL":"#4575B4",
  "USA - CA, AL":"#08519C",
  "USA - PA":"#ABD9E9"
  }

projects = projects.replace(data_to_replace)
mgcsts_samples_project = mgcsts_samples.merge(projects, on = "sampleID", how = "left")
df2 = mgcsts_samples_project.groupby(["Project", "mgCST"]).size().reset_index().pivot(columns='Project', index = 'mgCST', values =0).reset_index()
df2 = df2.reindex(columns=study_order[::-1])

# Calculate count_sample for each mgCST
count_sample = mgcsts_samples['mgCST'].value_counts().reset_index()
count_sample.columns = ['mgCST', 'count_sample']

# Merge with mgCSTs to get the color and domTaxa columns
count_sample = count_sample.merge(mgcsts[['mgCST', 'color', 'domTaxa']], on='mgCST', how='left')

# Ensure mgCST is treated as a categorical variable with sorted levels
count_sample['mgCST'] = pd.Categorical(count_sample['mgCST'], categories=sorted(count_sample['mgCST'].unique()), ordered=True)

# Create a color map for mgCST to color
color_map = dict(zip(count_sample['mgCST'], count_sample['color']))

######################################## ----- PAGE CONTENT ----- ########################################

st.container()

col1, col2, col3 = st.columns(3)

# Barplot mgcsts - number of samples colored by dominant Taxa
with col1 :

    fig = px.bar(
        count_sample,
        x='mgCST',
        y='count_sample',
        color='mgCST',
        hover_data=['domTaxa'],
        color_discrete_map=color_map, 
        labels={"count_sample" : "Number of samples"},
        title = "Distribution and prevalence of mgCSTs"
        )

    # Customize layout
    fig.update_layout(
        xaxis_title='mgCST',
        yaxis_title='Number of samples',
        template="plotly_white"
    )

    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])

    fig.update_layout(
        xaxis=dict(
            tickangle=0,  # Change the angle of the tick labels
            tickfont=dict(size=10)  # Change the size of the tick labels
        ))
    
    st.plotly_chart(fig, use_container_width=True)

# Stacked barplot mgcsts - number of samples colored by project
with col2 :

    fig = px.bar(
        df2, 
        x='mgCST', 
        y=df2.drop('mgCST', axis=1).columns.values, 
        labels={'value' : 'Number of samples'},
        color_discrete_map=color_pal,
        title = f"Distribution and prevalence of mgCSTs by region")
       
    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])

    fig.update_layout(
    xaxis=dict(
        tickangle=0,  # Change the angle of the tick labels
        tickfont=dict(size=10)  # Change the size of the tick labels
        ))
    fig.update_layout(
    legend=dict(
        title = 'Region',
        traceorder='reversed'
        )
    )
    
    st.plotly_chart(fig, use_container_width=True)

with col3 :
    new_vs_old = pd.read_csv("data/new_vs_old_mgcsts.csv")
    # old_mgCST = File_S6[['mapID', 'mgCST']]
    # old_mgCST = old_mgCST.rename(columns={'mapID':'sampleID'})
    # new_mgCST = mgcsts_samples
    # new_vs_old = pd.merge(new_mgCST, old_mgCST, on='sampleID', how='inner')
    # st.write(new_vs_old.shape)
    bubble_data = new_vs_old.groupby(['mgCST_x', 'mgCST_y']).size().reset_index(name='count')
    bubble_data = bubble_data.rename(columns={"mgCST_x":"mgCST", "mgCST_y":"old_mgCST"})

    bubble_color = []
    mgcsts['mgCST'] = mgcsts['mgCST'].astype(int)
    for i in sorted(bubble_data['mgCST'].unique()):
        bubble_color.append(mgcsts[mgcsts['mgCST'] == i]['color'].values[0])

    bubble_data['mgCST'] = bubble_data['mgCST'].astype(str)
    fig = px.scatter(
        bubble_data,
        x='mgCST',
        y = 'old_mgCST',
        color='mgCST',
        color_discrete_sequence=list(bubble_color),
        size='count',
        title = "Distribution and prevalence of mgCSTs - Comparison between mgCST v1 and mgCST v2")
    
    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])
    fig.update_yaxes(tickmode='array', tickvals=np.arange(1,28))

    fig.update_layout(
        
        xaxis=dict(
            tickangle=0,  # Change the angle of the tick labels
            tickfont=dict(size=10)  # Change the size of the tick labels
        ),

        yaxis=dict(
            tickangle=0,  # Change the angle of the tick labels
            tickfont=dict(size=10)  # Change the size of the tick labels
        ),
    )

    st.plotly_chart(fig, use_container_width=True)


st.subheader("Most abund species per mgCSTs")
with st.expander("Show table"):
    st.dataframe(mgcsts[['mgCST','domTaxa','meanRelabund','color']].set_index(['mgCST']))


# MgCSTS HEATMAP construction - button option

st.subheader("MgCSTs heatmap")

pdf_viewer("medias/mgCST_VOG_heatmap.pdf")
