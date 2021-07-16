from setuptools import setup, find_packages

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='SINATRA Pro Test',
    version='0.0.23',
    author="Wai Shing Tang",
    description = 'Python3 package for SINATRA Pro.',
    long_description = long_description,
    url = "https://github.com/lcrawlab/SINATRA-Pro",
    project_urls = {
        "Bug Tracker" : "https://github.com/lcrawlab/SINATRA-Pro/issues",
    },
    license = 'GNU General Public License v3.0',
    license_files = 'LICENSE',
    python_requires = '>=3.6',
    install_requires = [ 
        'numpy>=1.18.0',
        'scipy>=1.5.0',
        'mdanalysis>=0.20.0',
        'fast-histogram>=0.9',
        'joblib>=0.16.0',
        ],
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
)


