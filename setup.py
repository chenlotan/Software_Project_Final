from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    packages=find_packages(),
    headers=['spkmeans.h', 'kmeans.h'],
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['spkmeansmodule.c', 'kmeans.c', 'spkmeans.c']
        )
    ]
)