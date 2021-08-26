from setuptools import setup, find_packages

# def readme():
#     with open('README.md') as f:
#         return f.read()

data_files_to_include = [('', ['README.md', 'LICENSE', 'example_expanded_amino_acid_alphabet.txt'])]

setup(name='olga',
      version='1.2.4',
      description='Compute generation probability of CDR3 sequences',
      long_description='text/markdown',
      url='https://github.com/zsethna/OLGA',
      author='Zachary Sethna',
      author_email='zachary.sethna@gmail.com',
      license='GPLv3',
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Intended Audience :: Healthcare Industry',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Physics',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.6',
            ],
      packages=find_packages(),
      package_data = {
            'default_models': [],
            'default_models/human_T_alpha/': ['default_models/human_T_alpha/*'],
            'default_models/human_T_beta/': ['default_models/human_T_beta/*'],
            'default_models/mouse_T_beta/': ['default_models/mouse_T_beta/*'],
            'default_models/human_B_heavy/': ['default_models/human_B_heavy/*'],
            'default_models/human_B_kappa/': ['default_models/human_B_kappa/*'],
            'default_models/human_B_lambda/': ['default_models/human_B_lambda/*']
            },
      data_files = data_files_to_include,
      include_package_data=True,
      entry_points = {'console_scripts': [
            'olga-compute_pgen=olga.compute_pgen:main',
            'olga-generate_sequences=olga.generate_sequences:main'
            ], },
      zip_safe=False)
