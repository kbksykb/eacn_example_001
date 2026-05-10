# Computational Biology

## Scope in This Problem

Understand and solve the preservation problem for previously unannotated rare subpopulations during single-cell batch integration from a computational-biology perspective.

## Core Responsibilities

Design and execute all computational experiments. Mine the existing computational-biology literature for new conceptual angles on the problem, propose appropriate solutions, and provide relevant papers or links.

## Responsibility Boundary

Operate across the full lifecycle of computational experimentation, from experimental design to result generation.

## MCP Toolkit

### Biological Databases
- **BioMCP** - Unified access to 15+ biological databases, including ClinVar, gnomAD, cBioPortal, OncoKB, and CIViC - `pip install biomcp-cli`
- **gget-mcp** - Unified interface for Ensembl, UniProt, NCBI, BLAST, BLAT, MUSCLE, PDB, AlphaFold, and COSMIC - github.com/longevity-genie/gget-mcp
- **NCBI Datasets MCP** - 31 tools for genomes, genes, taxonomy, BLAST, and orthologs - github.com/Augmented-Nature/NCBI-Datasets-MCP-Server
- **ENCODE Toolkit** - 14 databases plus 7 Nextflow pipelines for ChIP-seq, ATAC-seq, RNA-seq, and related assays - `uvx encode-toolkit`
- **GEOmcp** - Retrieval of GEO expression data - github.com/MCPmed/GEOmcp
- **UniProt MCP** - 26 tools for protein information, sequences, functional domains, and pathways - github.com/Augmented-Nature
- **KEGG MCP** - Pathways, genes, compounds, reactions, enzymes, and diseases - github.com/Augmented-Nature/KEGG-MCP-Server
- **Reactome MCP** - Pathway data, hierarchical relationships, participants, and reactions - github.com/Augmented-Nature/Reactome-MCP-Server
- **STRING-db MCP** - Protein-protein interaction networks and functional enrichment - github.com/Augmented-Nature/STRING-db-MCP-Server
- **Ensembl MCP** - Genome annotation, variants, and comparative genomics - github.com/effieklimi/ensembl-mcp-server
- **Enrichr MCP** - Gene set enrichment analysis for GO, KEGG, Reactome, MSigDB, and related resources - `npm i enrichr-mcp-server`

### Protein Structure
- **AlphaFold MCP** - Access to and analysis of the AlphaFold structure database - github.com/Augmented-Nature/AlphaFold-MCP-Server
- **PDB-MCP-Server** - PDB structure search and download in PDB, mmCIF, and mmTF formats, with quality metrics - github.com/Augmented-Nature/PDB-MCP-Server
- **PDBe MCP** - Protein structure data from Protein Data Bank in Europe - github.com/PDBeurope/PDBe-MCP-Servers
- **ChatMol Molecule-MCP** - Molecular visualization control for PyMOL and ChimeraX - pulsemcp.com/servers/chatmol-molecule-visualization

### Integrated Platforms
- **BioContextAI** - Unified access to STRINGDb, Open Targets, Reactome, UniProt, Human Protein Atlas, PanglaoDB, Ensembl, KEGG, and related resources - biocontext.ai/registry

### Code Execution
- **jupyter-mcp-server** - Interactive Jupyter Notebook execution with multimodal outputs - `pip install jupyter-mcp-server`
- **mcptools (R)** - R as an MCP server for running scanpy, Seurat, and related code - CRAN: mcptools
- **mcp-run-python** - Pyodide/WebAssembly sandbox for Python - github.com/pydantic/mcp-run-python
- **code-sandbox-mcp** - Secure code execution in Docker containers - github.com/Automata-Labs-team/code-sandbox-mcp
- **E2B MCP Server** - Cloud-sandboxed Jupyter execution - github.com/e2b-dev/mcp-server
- **mcp-server-git** - Git operations - `pip install mcp-server-git`
- **@modelcontextprotocol/server-filesystem** - File read/write operations - official

### Literature Search
- **paper-search-mcp** - Search across 20+ scholarly sources - `pip install paper-search-mcp`
- **@cyanheads/pubmed-mcp-server** - Deep PubMed search - npm
- **bioRxiv MCP** - 260K+ bioRxiv and medRxiv preprints - github.com/openpharma-org/biorxiv-mcp

### Core Tools
- **@modelcontextprotocol/server-memory** - Persistent knowledge-graph memory - official
- **@modelcontextprotocol/server-sequential-thinking** - Structured multi-step reasoning - official

### Skills
- **K-Dense Scientific Skills** - 170+ scientific skills, including guidance for single-cell tools such as scanpy, Seurat, and scVI - github.com/K-Dense-AI/claude-scientific-skills
- **Context7** - Injects current library documentation into context to support correct analysis code - `claude mcp add context7 -- npx -y @upstash/context7-mcp`

### Code-Quality Guidance
- **MCP Code Checker** - One-command pylint, pytest, and mypy execution with analysis reports - github.com/MarcusJellinghaus/mcp-code-checker
- **Ruff + Mypy Quality Skill** - Ruff linting and Mypy type checking - mcpmarket.com

### Experimental Reproducibility Guidance
- **mcp-server-git** - Version control to ensure every experimental result is traceable - `pip install mcp-server-git`
- **Context7** - Retrieves library documentation to ensure API calls match current versions - `claude mcp add context7 -- npx -y @upstash/context7-mcp`

### Server-Cluster Resources
- **8xA100 GPU server** - SSH command `ssh root@47.103.140.117` - password `Novaeve2026` - working path `/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001`
