try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(
    name='marburg_biobank',
    version='0.1',
    packages=['marburg_biobank',],
    license='BSD',
    url='/martha/imt/e/20160331_AG_Mueller_Biobank/revisions/5_201709',
    author='Florian Finkernagel',
    description = "An interface to our biobank",
    author_email='finkernagel@imt.uni-marburg.de',
    long_description='',
    install_requires=[
        'numpy>=1.3',
        'pandas>=0.16',
        ]
)
