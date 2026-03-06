import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="taxontabletools2",
    version="2.3.0",
    author="Till-Hendrik Macher",
    author_email="macher@uni-trier.de",
    description="Taxontabletools2 - A comprehensive, platform-independent graphical user interface software to explore and visualise DNA metabarcoding data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/taxontabletools2",
    packages=setuptools.find_packages(),
    license = 'MIT',
    install_requires=[
        "streamlit",
        "stqdm",
        "scikit-learn",
        "xlsxwriter",
        "pycirclize",
        "update_checker",
        "pandas",
        "numpy",
        "plotly",
        "matplotlib",
        "matplotlib-venn",
        "scipy",
        "biom-format",
        "statsmodels",
        "openpyxl",
        "psutil",
        "pymannkendall",
        "tqdm",
        "natsort",
        "kaleido",
        "GitPython",
        "watchdog",
    ],
    include_package_data = True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    entry_points = {
        "console_scripts" : [
            "taxontabletools2 = taxontabletools2.__main__:main",
        ]
    },
)
