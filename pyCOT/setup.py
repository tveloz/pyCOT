from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

setup(
    name='pyCOT',

    version='0.1',

    description='Reaction Network library to study Chemical Organization Theory',
    long_description='Reaction Network library to study Chemical Organization Theory',

    url='https://github.com/tveloz/pyCOT',

    author='Tomas Veloz',
    author_email='tveloz@gmail.com',

    license='GNU v.3',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science',
        'Intended Audience :: Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3) ',
        'Programming Language :: Python :: 3.9',
    ],

    keywords='chemical organization theory, reaction networks, pathway analysis, systems biology',

    packages=find_packages(exclude=['tests']),

    # run-time dependencies that will be installed by pip
    # install_requires=['numpy','pandas','beautifulsoup4', 'bitarray','scipy','networkx','libroadrunner','lxml','pyvis','pypoman','matplotlib', 'joblib']
    install_requires=['numpy','pandas','beautifulsoup4', 'bitarray','scipy','networkx','libroadrunner','lxml','pyvis','matplotlib','joblib','pulp','seaborn']
)
