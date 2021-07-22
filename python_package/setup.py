from setuptools import setup, find_packages

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='SINATRA Pro',
    version='0.0.1',
    author="Wai Shing Tang",
    description = 'Python3 package for SINATRA Pro.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = "https://github.com/lcrawlab/SINATRA-Pro",
    project_urls = {
        "Bug Tracker" : "https://github.com/lcrawlab/SINATRA-Pro/issues",
    },
    license = 'MIT License',
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
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
)


