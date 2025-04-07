import logging

import streamlit as st

logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)

# Set up Streamlit page configuration
st.set_page_config(page_title="Multi-page Streamlit App", page_icon="👋", layout="wide")

# Sidebar for page navigation
# st.sidebar.title("Navigation")
# page = st.sidebar.selectbox("Choose a page", ["Image Search", "Query to Database"])

# # Show the corresponding page content
# if page == "Mol Viewer":
#     # Import and run the 'image_search' page
#     import pages.mol_viewer

# elif page == "Database Explorer":
#     # Import and run the 'query_to_db' page
#     import pages.db_viewer
