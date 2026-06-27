"""
API Proxy Views
===============
All outgoing external API calls are routed through these backend proxies.
Frontend calls these endpoints instead of hitting external services directly.
This avoids CORS/firewall issues at universities and restricted networks.
"""
import logging
import urllib.parse

import requests
from django.http import HttpResponse, JsonResponse, HttpResponseBadRequest
from django.views.decorators.http import require_http_methods

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Allowed external hosts (origin allowlist for CORS-style backend protection)
# Only these domains may be contacted by the proxy layer.
# ---------------------------------------------------------------------------
ALLOWED_HOSTS = {
    'api.ribosome.xyz',
    'www.ebi.ac.uk',
    'alphafold.ebi.ac.uk',
    'rna.bgsu.edu',
    'data.rcsb.org',
    'models.rcsb.org',
    'coords.litemol.org',
}

# ---------------------------------------------------------------------------
# Centralised external base URLs
# ---------------------------------------------------------------------------
EXTERNAL = {
    'ribosome_xyz_ban':       'https://api.ribosome.xyz/neo4j/get_banclass_for_chain/',
    'ribosome_xyz_polymers':  'https://api.ribosome.xyz/polymers/polynucleotide',
    'ribosome_xyz_nom':       'https://api.ribosome.xyz/neo4j/gmo_nom_class/',
    'pdbe_static_entry':      'https://www.ebi.ac.uk/pdbe/static/entry/',
    'pdbe_model_server':      'https://www.ebi.ac.uk/pdbe/model-server/v1',
    'pdbe_api':               'https://www.ebi.ac.uk/pdbe/api',
    'pdbe_graph_api':         'https://www.ebi.ac.uk/pdbe/graph-api',
    'pdbe_static':            'https://www.ebi.ac.uk/pdbe/static',
    'pdbe_entry_files':       'https://www.ebi.ac.uk/pdbe/entry-files/download/',
    'pdbe_volume_server':     'https://www.ebi.ac.uk/pdbe/volume-server',
    'rna3dhub':               'http://rna.bgsu.edu/rna3dhub/rest',
    'alphafold':              'https://alphafold.ebi.ac.uk/api/prediction/',
    'ncbi_blast_base':        'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast',
    'rcsb_graphql':           'https://data.rcsb.org/graphql',
    'rcsb_model_server':      'https://models.rcsb.org/v1',
    'litemol_coords':         'https://coords.litemol.org',
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _safe_get(url, timeout=15):
    """GET with host allowlist enforcement."""
    host = urllib.parse.urlparse(url).hostname
    if host not in ALLOWED_HOSTS:
        raise ValueError(f"Host {host!r} is not in the proxy allowlist")
    return requests.get(url, timeout=timeout)


def _safe_post(url, **kwargs):
    """POST with host allowlist enforcement."""
    host = urllib.parse.urlparse(url).hostname
    if host not in ALLOWED_HOSTS:
        raise ValueError(f"Host {host!r} is not in the proxy allowlist")
    kwargs.setdefault('timeout', 30)
    return requests.post(url, **kwargs)


def _json_ok(response):
    """Return a JsonResponse from an upstream JSON response."""
    response.raise_for_status()
    data = response.json()
    if isinstance(data, list):
        return JsonResponse(data, safe=False)
    return JsonResponse(data)


def _sanitize_path(path):
    """Strip path-traversal segments from a URL path component."""
    import posixpath
    cleaned = posixpath.normpath(path).lstrip('/')
    if '..' in cleaned.split('/'):
        raise ValueError(f"Path traversal detected in {path!r}")
    return cleaned


# Allowed mapping types for PDBe mappings endpoint
_VALID_MAPPING_TYPES = {'uniprot', 'pfam', 'cath', 'scop', 'interpro', 'go', 'ec',
                        'ensembl', 'hmmer', 'sequence_domains'}

# Allowed query prefixes for PDBe model-server endpoint
_VALID_QUERY_PREFIXES = {'residueSurroundings', 'atoms', 'full', 'residueRange',
                         'assembly', 'chains', 'ligand'}


def _binary_ok(response, content_type='application/octet-stream'):
    return HttpResponse(response.content, content_type=content_type)


def _text_ok(response):
    return HttpResponse(response.content, content_type='text/plain')


def _error(msg, status=500):
    return JsonResponse({'error': msg}, status=status)


# ═══════════════════════════════════════════════════════════════════════════
#  1. Ribosome.xyz
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["GET"])
def proxy_ban_name(request):
    """BAN class for a chain  →  api.ribosome.xyz"""
    pdb_id = request.GET.get('pdbid')
    chain_id = request.GET.get('auth_asym_id')
    if not pdb_id or not chain_id:
        return HttpResponseBadRequest("Missing pdbid or auth_asym_id")
    try:
        url = f"{EXTERNAL['ribosome_xyz_ban']}?pdbid={pdb_id}&auth_asym_id={chain_id}&format=json"
        return _json_ok(_safe_get(url, timeout=30))
    except requests.exceptions.Timeout:
        logger.warning("proxy_ban_name timed out for %s/%s", pdb_id, chain_id)
        return JsonResponse([], safe=False)
    except Exception as e:
        logger.error("proxy_ban_name: %s", e)
        return _error("Failed to fetch BAN name data")


@require_http_methods(["GET"])
def proxy_ribosome_polymers(request):
    """Polynucleotide list by RNA class  →  api.ribosome.xyz"""
    rna_class = request.GET.get('rna_class')
    if not rna_class:
        return HttpResponseBadRequest("Missing rna_class")
    try:
        url = f"{EXTERNAL['ribosome_xyz_polymers']}?rna_class={rna_class}&format=json"
        return _json_ok(_safe_get(url, timeout=30))
    except requests.exceptions.Timeout:
        logger.warning("proxy_ribosome_polymers timed out for rna_class=%s", rna_class)
        return JsonResponse([], safe=False)
    except Exception as e:
        logger.error("proxy_ribosome_polymers: %s", e)
        return _error("Failed to fetch ribosome polymer data")


@require_http_methods(["GET"])
def proxy_ribosome_nom(request):
    """Nomenclature class by BAN name  →  api.ribosome.xyz"""
    ban_name = request.GET.get('banName')
    if not ban_name:
        return HttpResponseBadRequest("Missing banName")
    try:
        url = f"{EXTERNAL['ribosome_xyz_nom']}?banName={ban_name}&format=json"
        return _json_ok(_safe_get(url, timeout=30))
    except requests.exceptions.Timeout:
        logger.warning("proxy_ribosome_nom timed out for banName=%s", ban_name)
        return JsonResponse([], safe=False)
    except Exception as e:
        logger.error("proxy_ribosome_nom: %s", e)
        return _error("Failed to fetch ribosome nomenclature data")


# ═══════════════════════════════════════════════════════════════════════════
#  2. PDBe  (European Bioinformatics Institute)
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["GET"])
def proxy_pdbe_static_entry(request):
    """RNA 2D topology JSON  →  PDBe static/entry"""
    pdb_id = request.GET.get('pdbid')
    entity_id = request.GET.get('entity_id')
    chain_id = request.GET.get('chain_id')
    if not all([pdb_id, entity_id, chain_id]):
        return HttpResponseBadRequest("Missing pdbid, entity_id or chain_id")
    try:
        url = f"{EXTERNAL['pdbe_static_entry']}{pdb_id.lower()}_{entity_id}_{chain_id}.json"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_static_entry: %s", e)
        return _error("Failed to fetch PDBe entry data")


@require_http_methods(["GET"])
def proxy_pdbe_model_server(request):
    """Atomic coordinates  →  PDBe model-server

    Supports two modes:
      1. Simple: ?pdbid=X&entity_id=Y&encoding=bcif  → /v1/{pdb}/atoms?...
      2. Full query: ?pdbid=X&query=full&encoding=bcif → /v1/{pdb}/{query}
         (used by Molstar ligand view, residueSurroundings, etc.)
    """
    pdb_id = request.GET.get('pdbid')
    encoding = request.GET.get('encoding', 'bcif')
    query_path = request.GET.get('query', '')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        if query_path:
            # Validate query prefix against allowlist
            query_prefix = query_path.split('?')[0].split('/')[0]
            if query_prefix not in _VALID_QUERY_PREFIXES:
                return HttpResponseBadRequest(f"Invalid query type: {query_prefix}")
            # Full query mode: /v1/{pdb}/{query_path}&encoding=...
            sep = '&' if '?' in query_path else '?'
            url = f"{EXTERNAL['pdbe_model_server']}/{pdb_id.lower()}/{query_path}{sep}encoding={encoding}"
            low_prec = request.GET.get('lowPrecisionCoords')
            if low_prec:
                url += f"&lowPrecisionCoords={low_prec}"
        else:
            # Simple atoms mode
            params = [f"encoding={encoding}"]
            entity_id = request.GET.get('entity_id')
            auth_asym_id = request.GET.get('auth_asym_id')
            if entity_id:
                params.append(f"label_entity_id={entity_id}")
            if auth_asym_id:
                params.append(f"auth_asym_id={auth_asym_id}")
            url = f"{EXTERNAL['pdbe_model_server']}/{pdb_id.lower()}/atoms?{'&'.join(params)}"
        resp = _safe_get(url, timeout=30)
        resp.raise_for_status()
        if encoding == 'bcif':
            return _binary_ok(resp)
        return _json_ok(resp)
    except Exception as e:
        logger.error("proxy_pdbe_model_server: %s", e)
        return _error("Failed to fetch PDBe model data")


@require_http_methods(["GET"])
def proxy_pdbe_entry_summary(request):
    """Entry summary  →  PDBe API"""
    pdb_id = request.GET.get('pdbid')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        url = f"{EXTERNAL['pdbe_api']}/pdb/entry/summary/{pdb_id.lower()}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_entry_summary: %s", e)
        return _error("Failed to fetch PDBe entry summary")


@require_http_methods(["GET"])
def proxy_pdbe_molecules(request):
    """Molecule / entity list  →  PDBe API"""
    pdb_id = request.GET.get('pdbid')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        url = f"{EXTERNAL['pdbe_api']}/pdb/entry/molecules/{pdb_id.lower()}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_molecules: %s", e)
        return _error("Failed to fetch PDBe molecules data")


@require_http_methods(["GET"])
def proxy_pdbe_polymer_coverage(request):
    """Polymer coverage  →  PDBe API"""
    pdb_id = request.GET.get('pdbid')
    chain_id = request.GET.get('chain_id')
    if not pdb_id or not chain_id:
        return HttpResponseBadRequest("Missing pdbid or chain_id")
    try:
        url = f"{EXTERNAL['pdbe_api']}/pdb/entry/polymer_coverage/{pdb_id.lower()}/chain/{chain_id}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_polymer_coverage: %s", e)
        return _error("Failed to fetch polymer coverage")


@require_http_methods(["GET"])
def proxy_pdbe_entry_files(request):
    """Binary CIF / other file download  →  PDBe"""
    pdb_id = request.GET.get('pdbid')
    fmt = request.GET.get('format', 'bcif')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        url = f"{EXTERNAL['pdbe_entry_files']}{pdb_id.lower()}.{fmt}"
        resp = _safe_get(url, timeout=30)
        resp.raise_for_status()
        return _binary_ok(resp)
    except Exception as e:
        logger.error("proxy_pdbe_entry_files: %s", e)
        return _error("Failed to download PDBe entry file")


@require_http_methods(["GET"])
def proxy_pdbe_mappings(request):
    """Domain / UniProt mappings  →  PDBe API"""
    pdb_id = request.GET.get('pdbid')
    mapping_type = request.GET.get('type', '')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        if mapping_type and mapping_type not in _VALID_MAPPING_TYPES:
            return HttpResponseBadRequest(f"Invalid mapping type: {mapping_type}")
        url = f"{EXTERNAL['pdbe_api']}/mappings/{pdb_id.lower()}"
        if mapping_type:
            url = f"{EXTERNAL['pdbe_api']}/mappings/{mapping_type}/{pdb_id.lower()}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_mappings: %s", e)
        return _error("Failed to fetch PDBe mappings")


@require_http_methods(["GET"])
def proxy_pdbe_superposition(request):
    """Superposition segments / matrices  →  PDBe graph-api / static"""
    molecule_id = request.GET.get('molecule_id')
    data_type = request.GET.get('type', 'segments')
    if not molecule_id:
        return HttpResponseBadRequest("Missing molecule_id")
    try:
        if data_type == 'segments':
            url = f"{EXTERNAL['pdbe_graph_api']}/uniprot/superposition/{molecule_id}"
        elif data_type == 'matrices':
            acc = request.GET.get('matrix_accession', molecule_id)
            url = f"{EXTERNAL['pdbe_static']}/superpose/matrices/{acc}"
        else:
            return HttpResponseBadRequest("Invalid type parameter")
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_superposition: %s", e)
        return _error("Failed to fetch superposition data")


@require_http_methods(["GET"])
def proxy_pdbe_carbohydrate(request):
    """Carbohydrate polymer info  →  PDBe API"""
    pdb_id = request.GET.get('pdbid')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        url = f"{EXTERNAL['pdbe_api']}/pdb/entry/carbohydrate_polymer/{pdb_id.lower()}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_pdbe_carbohydrate: %s", e)
        return _error("Failed to fetch carbohydrate data")


# ═══════════════════════════════════════════════════════════════════════════
#  3. FR3D / RNA3DHub  (via backend – no CORS proxy needed anymore)
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["GET"])
def proxy_fr3d_data(request):
    """Base-pair annotations  →  rna.bgsu.edu (direct from server)"""
    pdb_id = request.GET.get('pdb_id')
    chain = request.GET.get('chain')
    only_nested = request.GET.get('only_nested', 'False')
    if not pdb_id or not chain:
        return HttpResponseBadRequest("Missing pdb_id or chain")
    try:
        url = f"{EXTERNAL['rna3dhub']}/getChainSequenceBasePairs?pdb_id={pdb_id.lower()}&chain={chain}&only_nested={only_nested}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_fr3d_data: %s", e)
        return _error("Failed to fetch FR3D data")


# ═══════════════════════════════════════════════════════════════════════════
#  4. AlphaFold
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["GET"])
def proxy_alphafold(request):
    """Prediction metadata  →  alphafold.ebi.ac.uk"""
    accession = request.GET.get('accession')
    if not accession:
        return HttpResponseBadRequest("Missing accession")
    try:
        url = f"{EXTERNAL['alphafold']}{accession}"
        return _json_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_alphafold: %s", e)
        return _error("Failed to fetch AlphaFold data")


# ═══════════════════════════════════════════════════════════════════════════
#  5. NCBI BLAST  (via EBI)
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["POST"])
def proxy_ncbi_blast(request):
    """Submit BLAST job  →  EBI ncbiblast/run"""
    try:
        url = f"{EXTERNAL['ncbi_blast_base']}/run"
        resp = _safe_post(url, data=request.body, headers={
            'Content-Type': 'application/x-www-form-urlencoded',
        }, timeout=60)
        resp.raise_for_status()
        return _text_ok(resp)
    except Exception as e:
        logger.error("proxy_ncbi_blast: %s", e)
        return _error("Failed to submit BLAST job")


@require_http_methods(["GET"])
def proxy_ncbi_blast_status(request):
    """Check BLAST job status"""
    job_id = request.GET.get('job_id')
    if not job_id:
        return HttpResponseBadRequest("Missing job_id")
    try:
        url = f"{EXTERNAL['ncbi_blast_base']}/status/{job_id}"
        return _text_ok(_safe_get(url))
    except Exception as e:
        logger.error("proxy_ncbi_blast_status: %s", e)
        return _error("Failed to check BLAST status")


@require_http_methods(["GET"])
def proxy_ncbi_blast_result(request):
    """Retrieve BLAST results"""
    job_id = request.GET.get('job_id')
    result_type = request.GET.get('type', 'out')
    fmt = request.GET.get('format', '10')
    if not job_id:
        return HttpResponseBadRequest("Missing job_id")
    try:
        url = f"{EXTERNAL['ncbi_blast_base']}/result/{job_id}/{result_type}?format={fmt}"
        return _text_ok(_safe_get(url, timeout=30))
    except Exception as e:
        logger.error("proxy_ncbi_blast_result: %s", e)
        return _error("Failed to fetch BLAST results")


# ═══════════════════════════════════════════════════════════════════════════
#  6. RCSB
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["GET", "POST"])
def proxy_rcsb_graphql(request):
    """GraphQL queries  →  data.rcsb.org"""
    query = request.GET.get('query') if request.method == 'GET' else request.POST.get('query')
    if not query:
        return HttpResponseBadRequest("Missing query")
    try:
        url = EXTERNAL['rcsb_graphql']
        if request.method == 'GET':
            resp = _safe_get(f"{url}?query={urllib.parse.quote(query)}")
        else:
            resp = _safe_post(url, json={'query': query})
        resp.raise_for_status()
        return _json_ok(resp)
    except Exception as e:
        logger.error("proxy_rcsb_graphql: %s", e)
        return _error("Failed to query RCSB database")


@require_http_methods(["GET"])
def proxy_rcsb_model_server(request):
    """Atomic model download  →  models.rcsb.org"""
    pdb_id = request.GET.get('pdbid')
    encoding = request.GET.get('encoding', 'bcif')
    if not pdb_id:
        return HttpResponseBadRequest("Missing pdbid")
    try:
        url = f"{EXTERNAL['rcsb_model_server']}/{pdb_id.lower()}/atoms?encoding={encoding}"
        resp = _safe_get(url, timeout=30)
        resp.raise_for_status()
        return _binary_ok(resp)
    except Exception as e:
        logger.error("proxy_rcsb_model_server: %s", e)
        return _error("Failed to fetch RCSB model data")


# ═══════════════════════════════════════════════════════════════════════════
#  7. LiteMol / Coords
# ═══════════════════════════════════════════════════════════════════════════

@require_http_methods(["GET"])
def proxy_pdbe_volume_server(request):
    """Volume streaming data  →  PDBe volume-server (passthrough path)"""
    path = request.GET.get('path')
    if not path:
        return HttpResponseBadRequest("Missing path")
    try:
        path = _sanitize_path(path)
        url = f"{EXTERNAL['pdbe_volume_server']}/{path}"
        resp = _safe_get(url, timeout=30)
        resp.raise_for_status()
        ct = resp.headers.get('content-type', 'application/octet-stream')
        return HttpResponse(resp.content, content_type=ct)
    except Exception as e:
        logger.error("proxy_pdbe_volume_server: %s", e)
        return _error("Failed to fetch volume server data")


@require_http_methods(["GET"])
def proxy_litemol_coords(request):
    """Coordinate data  →  coords.litemol.org  (passthrough path)"""
    path = request.GET.get('path')
    if not path:
        return HttpResponseBadRequest("Missing path")
    try:
        path = _sanitize_path(path)
        url = f"{EXTERNAL['litemol_coords']}/{path}"
        resp = _safe_get(url, timeout=30)
        resp.raise_for_status()
        return _binary_ok(resp)
    except Exception as e:
        logger.error("proxy_litemol_coords: %s", e)
        return _error("Failed to fetch LiteMol coordinate data")
