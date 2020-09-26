from setuptools import setup

setup(name='gcsnap',
      version='1.0.2',
      description='Interactive snapshots for the comparison of protein-coding genomic contexts',
      long_description='GCsnap: interactive snapshots for the comparison of protein-coding genomic contexts.',
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
          'bokeh >= 1.3.4, <= 2.1.1',
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
        'console_scripts': 'gcsnap=gcsnap.gcsnap:main',}
)
