import io
import os
import re

from setuptools import find_packages
from setuptools import setup


def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())


setup(
    name="EFGs",
    version="0.8.4",
    url="https://github.com/HelloJocelynLu/EFGs",
    license='MIT',

    author="Jocelyn Lu",
    author_email="jl8570@nyu.edu",

    description="Extended Functional Groups",
    long_description=read("README.rst"),

    packages=find_packages(exclude=('tests',)),

    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
)
