from setuptools import setup
import io, re, os

setup(name='pycovis-postprocess',
      url='https://github.com/COVIS-Sonar/postprocessing',
      author='Aaron Marburg',
      author_email='amarburg@apl.washington.edu',
      # version=version,
      # description='Module for retrieving CamHD data through a LazyCache server',
      # long_description='README.rst',
      license='MIT',
      python_requires='>=3',
      packages=['pycovis.postprocess'],
      install_requires=["python-decouple","pyunpack","patool"]
)
