from setuptools import setup, find_packages
import os


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='pathogen-embed',
    version='1.0.0',
    description='Reduced dimension embeddings for pathogen sequences',
    url='https://github.com/blab/pathogen-embed/',
    author='Sravani Nanduri <nandsra@cs.washington.edu> , John Huddleston <huddlej@gmail.com>',
    author_email='nandsra@cs.washington.edu',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT License',
    project_urls = {
        "Documentation": "https://blab.github.io/pathogen-embed/",
        "Bug Reports": "https://github.com/blab/pathogen-embed/issues",
        "Source Code": "https://github.com/blab/pathogen-embed/tree/main",
	"Change Log": "https://github.com/blab/pathogen-embed/tree/master/CHANGES.md",
    },
    package_dir={"": "src"},
	packages = find_packages(where="src", exclude=['test']),
    #packages=['seaborn', 'scikit-learn', 'umap-learn', 'matplotlib', 'pandas', 'numpy', 'hdbscan'],
    install_requires=['numpy',
                      'pandas',
                      "biopython",
                      'seaborn',
                      'scikit-learn',
                      'umap-learn',
                      'matplotlib',
                      'hdbscan'
                      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points = {
        "console_scripts": [
            "pathogen-embed = pathogen_embed.__main__:main",
            "pathogen-distance = pathogen_distance.__main__:main",
            "pathogen-cluster = pathogen_cluster.__main__:main"
        ]
    }
)
