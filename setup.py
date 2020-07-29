from setuptools import setup, find_packages

with open("README.md", 'r') as fh:
    long_description = fh.read()

setup(
    url = "https://github.com/mesquita-daniel/StabBB",
    author = 'Daniel Mesquita',
    name="stabbb",
    version="0.0.0.1",
    description="Implementation of stabilized Barzilai-Borwein Method for gradient descent problems",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    keywords = ['Optimization', 'Gardient Descent', 'Line Search'],
    packages=find_packages(),
    classifiers=[
        "Environment :: Console",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Cython",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)


