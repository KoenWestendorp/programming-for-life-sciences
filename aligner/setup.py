from setuptools import setup

setup(
    name='aligner',
    version='0.1.5',
    description='A sequence aligner that can be used to resolve and visualize homologies in protein sequences.',
    author='Koen Westendorp',
    author_email='koensswestendorp@gmail.com',
    packages=['aligner'],
    install_requires=['numpy', 'matplotlib']
)
