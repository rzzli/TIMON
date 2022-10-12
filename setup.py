from setuptools import setup, find_packages


setup(
    name = "MONet",
    version = "0.1.0",
    author = "Rick Z Li",
    author_email = "zhl022@eng.ucsd.edu",
    include_package_data=True,
    packages=find_packages(include=['MONet', 'MONet.*']),
    description = ("MONet is a tool to build transcription factor co-occurrence networks from epigenetic data."),
    install_requires=[
        'biopython ==1.77',
        'numpy ==1.21.5',
        'scipy ==1.6.1',
        'pandas ==1.3.5',
        'networkx ==2.6.3',
        'matplotlib ==3.3.2',
        'pygraphviz ==1.3',
        'jupyter'
    ],
    setup_requires=[  'flake8'],
    
    license = "MIT",
    keywords = "motif co-occurrence, transcription factors, microglia",
    url = "https://github.com/rzzli/MONet",
)
