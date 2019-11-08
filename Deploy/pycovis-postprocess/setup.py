from setuptools import setup
import io, re, os

# def read(*names, **kwargs):
#     with io.open(
#         os.path.join(os.path.dirname(__file__), *names),
#         encoding=kwargs.get("encoding", "utf8")
#     ) as fp:
#         return fp.read()
#
# def find_version(*file_paths):
#     version_file = read(*file_paths)
#     version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
#                               version_file, re.M)
#     if version_match:
#         return version_match.group(1)
#     raise RuntimeError("Unable to find version string.")

#version = find_version('pycamhd', 'lazycache', '__init__.py')

setup(name='pycovis-postprocess',
      # version=version,
      # description='Module for retrieving CamHD data through a LazyCache server',
      # long_description='README.rst',
      # url='https://github.com/CamHD-Analysis/pycamhd-lazycache',
      author='Aaron Marburg',
      author_email='amarburg@apl.washington.edu',
      license='MIT',
      python_requires='>=3',
      packages=['pycovis.postprocess'],
      install_requires=["minio", "python-decouple"]
      # install_requires=['requests','pillow'],
      # setup_requires=['pytest-runner'],
      # tests_require=['pytest', 'numpy', 'pillow']
)
