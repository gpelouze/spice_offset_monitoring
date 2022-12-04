import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()
with open('requirements.txt', 'r') as f:
    requirements = f.read().strip('\n').split('\n')

entry_points = {
    'console_scripts': [
        'spice_offset_monitoring=spice_offset_monitoring.cli:cli',
        ]
    }

setuptools.setup(
    name='spice_offset_monitoring',
    version='2022.12.02',
    author='Gabriel Pelouze',
    author_email='gabriel.pelouze@universite-paris-saclay.fr',
    description='Monitor the evolution of the pointing offset of SPICE',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/gpelouze/spice_offset_monitoring',
    entry_points=entry_points,
    packages=setuptools.find_packages(),
    python_requires='>=3.8',
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
    )
