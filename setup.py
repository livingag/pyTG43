from setuptools import setup

setup(
    name = "pyTG43",
    version = "0.2.0",
    author = "Alex Livingstone",
    author_email = "livingstone.alex@gmail.com",
    description = "Calculate dose to a point for a brachytherapy plan.",
    long_description = """
    Calculate dose to a point using provided treatment plan information in the form of DICOM files.

    Calculations are performed using 2D calculation framework outlined in AAPM TG-43.

    Data specifying the source parameters is required, examples and instructions are available at:

    http://www.github.com/livingag/pyTG43
    """,
    keywords = ["radiotherapy", "brachytherapy", "TG43", "medical physics"],
    url = "https://github.com/livingag/pyTG43",
    packages = ["pyTG43"],
    python_requires = '>=3.5.0',
    install_requires = [
        'numpy==1.14.2',
        'pydicom==1.0.2',
        'xlrd==1.1.0',
        'scipy==1.0.0',
        'terminaltables==3.1.0',
        'matplotlib==2.2.0+601.g197bb59'
    ],
    license = '',
    classifiers = [],
)