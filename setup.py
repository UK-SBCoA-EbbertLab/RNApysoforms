from setuptools import setup, find_packages

setup(
    name="rna_pysoforms",
    version="0.1.0-dev",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[

    "plotly>=5.0,<6.0",
    "polars[excel]>=1.0,<2.0",
    "pyarrow>=17.0.0,<18.0.0",

    ],
    python_requires='>=3.8',
    author="Bernardo Aguzzoli Heberle",
    author_email="bernardo.aguzzoli@gmail.com",
    description="A Plotly-and-Polars-based Python implementation of ggtranscript-like functionality for vizualizing RNA isoform structure and expression.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/UK-SBCoA-EbbertLab/rna_pysoforms",
)
