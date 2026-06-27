from django.urls import path
from . import views

app_name = 'api_proxies'

urlpatterns = [
    # ── Ribosome.xyz ──────────────────────────────────────────────────
    path('ban-name/', views.proxy_ban_name, name='proxy_ban_name'),
    path('ribosome/polymers/', views.proxy_ribosome_polymers, name='proxy_ribosome_polymers'),
    path('ribosome/nomenclature/', views.proxy_ribosome_nom, name='proxy_ribosome_nom'),

    # ── PDBe ──────────────────────────────────────────────────────────
    path('pdbe/static-entry/', views.proxy_pdbe_static_entry, name='proxy_pdbe_static_entry'),
    path('pdbe/model-server/', views.proxy_pdbe_model_server, name='proxy_pdbe_model_server'),
    path('pdbe/entry-summary/', views.proxy_pdbe_entry_summary, name='proxy_pdbe_entry_summary'),
    path('pdbe/molecules/', views.proxy_pdbe_molecules, name='proxy_pdbe_molecules'),
    path('pdbe/polymer-coverage/', views.proxy_pdbe_polymer_coverage, name='proxy_pdbe_polymer_coverage'),
    path('pdbe/entry-files/', views.proxy_pdbe_entry_files, name='proxy_pdbe_entry_files'),
    path('pdbe/mappings/', views.proxy_pdbe_mappings, name='proxy_pdbe_mappings'),
    path('pdbe/superposition/', views.proxy_pdbe_superposition, name='proxy_pdbe_superposition'),
    path('pdbe/carbohydrate/', views.proxy_pdbe_carbohydrate, name='proxy_pdbe_carbohydrate'),

    # ── FR3D / RNA3DHub ───────────────────────────────────────────────
    path('fr3d/data/', views.proxy_fr3d_data, name='proxy_fr3d_data'),

    # ── AlphaFold ─────────────────────────────────────────────────────
    path('alphafold/', views.proxy_alphafold, name='proxy_alphafold'),

    # ── NCBI BLAST ────────────────────────────────────────────────────
    path('ncbi-blast/', views.proxy_ncbi_blast, name='proxy_ncbi_blast'),
    path('ncbi-blast/status/', views.proxy_ncbi_blast_status, name='proxy_ncbi_blast_status'),
    path('ncbi-blast/result/', views.proxy_ncbi_blast_result, name='proxy_ncbi_blast_result'),

    # ── RCSB ──────────────────────────────────────────────────────────
    path('rcsb-graphql/', views.proxy_rcsb_graphql, name='proxy_rcsb_graphql'),
    path('rcsb/model-server/', views.proxy_rcsb_model_server, name='proxy_rcsb_model_server'),

    # ── PDBe Volume Server ─────────────────────────────────────────────
    path('pdbe/volume-server/', views.proxy_pdbe_volume_server, name='proxy_pdbe_volume_server'),

    # ── LiteMol / Coords ─────────────────────────────────────────────
    path('litemol/coords/', views.proxy_litemol_coords, name='proxy_litemol_coords'),
]
