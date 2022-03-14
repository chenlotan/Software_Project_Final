from setuptools import setup, find_packages, Extension

setup(
    name='spkmeansmodule',
    packages=find_packages(),
    ext_modules=[
        Extension(
            'spkmeansmodule',
            ['spkmeansmodule.c', 'kmeans.c', 'spkmeans.c']
        )
    ]
)