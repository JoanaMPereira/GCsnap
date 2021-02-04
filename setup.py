from setuptools import setup
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='gcsnap',
      version='1.0.9',
      description='GCsnap: Interactive snapshots for the comparison of protein-coding genomic contexts',
      long_description=long_description,  # Optional
      long_description_content_type='text/markdown',  # Optional (see note above)
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='genomic contexts homology',
      url='https://github.com/JoanaMPereira/GCsnap',
      author='Joana Pereira',
      author_email='pereira.joanam@gmail.com',
      packages=['gcsnap'],
      install_requires=[
          'biopython', 
          'bokeh == 1.3.4',
          'networkx',
          'numpy',
          'pandas',
          'requests_cache',
          'scipy',
          'matplotlib'
      ],
      include_package_data=True,
      zip_safe=False,

      entry_points={  # Optional
        'console_scripts': 'GCsnap=gcsnap.GCsnap:main',}
)
