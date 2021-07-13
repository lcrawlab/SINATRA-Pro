from setuptools import setup, find_packages

setup(
    name='SINATRA Pro Test',    # This is the name of your PyPI-package.
    version='0.0.10',                          # Update the version number for new releases
    description = 'Python3 package for SINATRA Pro',
    author='Wai Shing Tang',
    license = 'GNU General Public License v3.0',
    license_files = 'LICENSE',
    #packages=['sinatra_pro'],
    python_requires = '>=3.6',
    install_requires = [ 
        'numpy>=1.18.0',
        'scipy>=1.5.0',
        'mdanalysis>=0.20.0',
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


