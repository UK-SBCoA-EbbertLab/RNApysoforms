from setuptools import setup, find_packages

setup(
    name="RNApysoforms",
    version="1.0.1",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[

    "plotly>=5.0,<6.0",
    "polars[excel]>=1.0,<2.0",
    "pyarrow>=17.0,<18.0",
    "pandas>=1.3,<3.0"
    ],
    python_requires='>=3.8',
    author="Bernardo Aguzzoli Heberle",
    author_email="bernardo.aguzzoli@gmail.com",
    description="A Python package designed for visualizing RNA isoform structures and expression levels by leveraging Plotly for interactive plotting and Polars for efficient data manipulation, enabling the creation of fast-rendering, interactive plots.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/UK-SBCoA-EbbertLab/RNApysoforms",
)
