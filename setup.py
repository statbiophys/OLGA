from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='olga',
      version='0.1',
      description='Compute generation probability of CDR3 sequences',
      long_description=readme(),
      author='Zachary Sethna',
      author_email='sethna@princeton.edu',
      license='GPLv3',
      packages='olga',
      scripts=['./run_pgen', './generate_synthetic_sequences', './compute_single_sequence_pgen'],
      include_package_data=True,
      zip_safe=False)
