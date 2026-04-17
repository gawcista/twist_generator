import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="TwistGenerator",
  version="1.0.0",
  author="GAO Yifan",
  author_email="gawcista@gmail.com",
  description="Python3 packages for twisted system generation",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/gawcista/twist_generator",
  packages=setuptools.find_packages(),
  install_requires=[
  "numpy",
  "spglib",
  ],
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
)
