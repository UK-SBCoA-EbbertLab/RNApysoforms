from setuptools import setup, find_packages

setup(
    name="pytranscript",
    version="0.1",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=[
        'pandas',
        'numpy',
        'plotly',
    ],
    author="Bernardo Aguzzoli Heberle",
    author_email="bernardo.aguzzoli@gmail.com",
    description="A Plotly-based Python implementation of ggtranscript-like functionality",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/UK-SBCoA-EbbertLab/pytranscript",
)
