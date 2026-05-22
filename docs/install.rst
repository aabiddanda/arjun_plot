Installation
============

From source
-----------

Clone the repository and install with ``pip``::

  git clone https://github.com/aabiddanda/arjun_plot.git
  cd arjun_plot/
  pip install .

Directly from GitHub
--------------------

Install the latest stable version from the ``main`` branch without cloning::

  pip install git+https://github.com/aabiddanda/arjun_plot.git@main

Development install
-------------------

To install in editable mode (changes to the source are reflected immediately)::

  git clone https://github.com/aabiddanda/arjun_plot.git
  cd arjun_plot/
  pip install -e ".[dev]"

Only the ``main`` branch is recommended for general use.

Optional dependencies
---------------------

Some functions require packages that are not installed by default:

* ``statgen.plot_gene_region_worker`` / ``gene_plot`` — requires ``requests`` to query the UCSC Genome Browser API.
* ``statgen.extract_ld_matrix`` — requires ``cyvcf2`` and ``tqdm`` to read VCF files.

