from setuptools import find_packages, setup

setup(
    name="decaf",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
    ],
    author="Auguste GARDETTE",
    author_email="auguste.gardette@ird.fr",
    description="DECAF: DEcontamination and Classification of Amplicon Fragment",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Aaramis/DECAF",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "decaf=decaf.cli:main",
        ],
    },
)
