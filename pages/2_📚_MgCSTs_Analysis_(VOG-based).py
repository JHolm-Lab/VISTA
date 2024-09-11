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
projects = pd.read_csv("data/VIRGO2_projects.csv")

######################################## ----- DATA PROCESSING ----- #############################################

mgcsts_samples_project = mgcsts_samples.merge(projects, on = "sampleID", how = "left")

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
        title = "Number of samples on each MgCSTs - colored by mgCST"
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

    df2 = mgcsts_samples_project.groupby(["Project", "mgCST"]).size().reset_index().pivot(columns='Project', index = 'mgCST', values =0).reset_index()

    fig = px.bar(
        df2, 
        x='mgCST', 
        y=df2.drop('mgCST', axis=1).columns.values, 
        labels={'value' : 'Number of samples'},
        title = f"Number of samples on each MgCSTs - colored by projects")
       
    fig.update_xaxes(tickmode='array', tickvals=mgcsts['mgCST'])

    fig.update_layout(
    xaxis=dict(
        tickangle=0,  # Change the angle of the tick labels
        tickfont=dict(size=10)  # Change the size of the tick labels
        ))
    
    st.plotly_chart(fig, use_container_width=True)

with col3 :

    File_S6 = pd.read_excel('data/File_S6_clean.xlsx')
    old_mgCST = File_S6[['mapID', 'mgCST']]
    old_mgCST = old_mgCST.rename(columns={'mapID':'sampleID'})
    new_mgCST = mgcsts_samples
    new_vs_old = pd.merge(new_mgCST, old_mgCST, on='sampleID', how='inner')
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
        title = "Number of samples in previous MgCST vs New MgCST")
    
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

pdf_viewer("medias/mgCST_heatmap.pdf")