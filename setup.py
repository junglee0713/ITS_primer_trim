from setuptools import setup

setup(
    name='primertrim',
    version='0.0.1',
    description='Trim primer sequences from FASTQ files',
    author='Kyle Bittinger',
    author_email='kylebittinger@gmail.com',
    url='https://github.com/junglee0713/ITS_primer_trim',
    packages=['primertrim'],
    entry_points = {
        'console_scripts': [
            'remove_primers.py=primertrim.remove_primers:main',
        ],
    }
)
