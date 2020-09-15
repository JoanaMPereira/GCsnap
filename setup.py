from setuptools import setup

setup(name='gcsnap',
      version='0.1',
      description='Interactive snapshots for the comparison of protein-coding genomic contexts',
      long_description='GCsnap: interactive snapshots for the comparison of protein-coding genomic contexts.',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3.0',
        'Programming Language :: Python :: 3.1',
        'Topic :: Bioinformatics :: Genomic contexts',
      ],
      keywords='genomic contexts homology',
      url='http:https://github.com/JoanaMPereira/GCsnap',
      author='Joana Pereira',
      author_email='pereira.joanam@gmail.com',
      license='GNU',
      packages=['gcsnap'],
      install_requires=[
          'biopython >= 1.74', 
          'bokeh >= 1.3.4, <= 2.1.1',
          'networkx >= 2.3',
          'numpy >= 1.17.2',
          'pandas >= 0.25.1'
      ],
      include_package_data=True,
      zip_safe=False)